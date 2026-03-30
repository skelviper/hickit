#include <cuda_runtime.h>
#include <cstdio>
#include <cmath>
#include <cstdint>
#include <cstddef>
#include <new>
#include <limits>
#include <float.h>
#include <algorithm>
#include <cstring>

#include "fdg_gpu.h"
#include "hkpriv.h"


extern int hk_verbose;

struct hk_gpu_vec3_t {
	float x, y, z;
};

static_assert(sizeof(hk_gpu_vec3_t) == sizeof(fvec3_t), "GPU vector size mismatch");

#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ >= 350
#define FDG_LDG_FLOAT(ptr) __ldg(ptr)
#define FDG_LDG_INT(ptr) __ldg(ptr)
#define FDG_LDG_UINT64(ptr) __ldg(ptr)
#else
#define FDG_LDG_FLOAT(ptr) (*(ptr))
#define FDG_LDG_INT(ptr) (*(ptr))
#define FDG_LDG_UINT64(ptr) (*(ptr))
#endif

__device__ __forceinline__ hk_gpu_vec3_t fdg_load_vec3(const hk_gpu_vec3_t *__restrict__ pos, int idx)
{
	hk_gpu_vec3_t v;
	v.x = FDG_LDG_FLOAT(&pos[idx].x);
	v.y = FDG_LDG_FLOAT(&pos[idx].y);
	v.z = FDG_LDG_FLOAT(&pos[idx].z);
	return v;
}

__device__ __forceinline__ void fdg_store_vec3(hk_gpu_vec3_t *__restrict__ arr, int idx, float x, float y, float z)
{
	arr[idx].x = x;
	arr[idx].y = y;
	arr[idx].z = z;
}

struct hk_fdg_gpu_ctx {
	int32_t n_beads;
	size_t bead_capacity;
	size_t capacity_pairs;
	size_t active_pairs;
	size_t block_table_cap;
	size_t block_keys_cap;
	size_t cell_capacity;
	size_t cell_hash_cap;
	hk_gpu_vec3_t *d_pos;
	hk_gpu_vec3_t *d_prev_pos;
	hk_gpu_vec3_t *d_force;
	struct hk_fdg_gpu_pair *d_pairs;
	struct hk_fdg_gpu_stats *d_stats;
	uint64_t *d_block_table;
	uint64_t *d_block_keys;
	uint64_t *d_cell_hash_keys;
	int *d_cell_hash_vals;
	int *d_bead_next;
	float4 *d_grid_params;
	float *d_bounds;
	float *h_bounds_init_pinned; /* pinned buffer for bounds initialization */
	int n_cells;
	cudaStream_t stream;
	uint32_t pair_counts[HK_FDG_PAIR_TYPE_COUNT];
	int pairs_ready;
	struct hk_fdg_gpu_stats *h_stats_pinned;
};

static int hk_fdg_gpu_cuda_check(cudaError_t err, const char *msg);

template <typename T>
static void hk_fdg_gpu_free_device(T **ptr)
{
	if (ptr && *ptr) {
		cudaFree(*ptr);
		*ptr = nullptr;
	}
}

template <typename T>
static int hk_fdg_gpu_resize_device(T **ptr, size_t count, const char *tag)
{
	if (count == 0) {
		hk_fdg_gpu_free_device(ptr);
		return 0;
	}
	T *tmp = nullptr;
	cudaError_t err = cudaMalloc(&tmp, count * sizeof(T));
	if (err != cudaSuccess) {
		if (hk_verbose >= 1)
			std::fprintf(stderr, "[E::fdg-gpu] failed to allocate %s (%zu elements): %s\n",
						 tag, count, cudaGetErrorString(err));
		return -1;
	}
	hk_fdg_gpu_free_device(ptr);
	*ptr = tmp;
	return 0;
}

static inline int hk_fdg_gpu_check_last_error(const char *msg)
{
	return hk_fdg_gpu_cuda_check(cudaGetLastError(), msg);
}

static int hk_fdg_gpu_cuda_check(cudaError_t err, const char *msg)
{
	if (err != cudaSuccess) {
		if (hk_verbose >= 1)
			std::fprintf(stderr, "[E::fdg-gpu] %s: %s\n", msg, cudaGetErrorString(err));
		return -1;
	}
	return 0;
}

static inline size_t next_pow2_size(size_t v)
{
	if (v == 0) return 1;
	--v;
	v |= v >> 1;
	v |= v >> 2;
	v |= v >> 4;
	v |= v >> 8;
	v |= v >> 16;
	if (sizeof(size_t) == 8)
		v |= v >> 32;
	return v + 1;
}

__host__ __device__ static inline uint64_t fdg_hash64(uint64_t key)
{
	key = ~key + (key << 21);
	key = key ^ key >> 24;
	key = (key + (key << 3)) + (key << 8);
	key = key ^ key >> 14;
	key = (key + (key << 2)) + (key << 4);
	key = key ^ key >> 28;
	key = key + (key << 31);
	return key;
}

__device__ static inline uint64_t pack_cell(int cx, int cy, int cz)
{
	return (uint64_t(cx) << 42) | (uint64_t(cy) << 21) | uint64_t(cz);
}

static __device__ __forceinline__ int fdg_cell_find_slot(const uint64_t *hash_keys, uint32_t mask, uint64_t key)
{
	const uint64_t EMPTY_KEY = 0xffffffffffffffffULL;
	uint32_t slot = uint32_t(fdg_hash64(key)) & mask;
	for (uint32_t iter = 0; iter <= mask; ++iter) {
		uint64_t stored = FDG_LDG_UINT64(&hash_keys[slot]);
		if (stored == EMPTY_KEY) return -1;
		if (stored == key) return slot;
		slot = (slot + 1u) & mask;
	}
	return -1;
}

__device__ static inline void atomicMinFloat(float *addr, float value)
{
	int *addr_as_int = reinterpret_cast<int*>(addr);
	int old = *addr_as_int;
	while (__int_as_float(old) > value) {
		int assumed = old;
		old = atomicCAS(addr_as_int, assumed, __float_as_int(value));
		if (assumed == old) break;
	}
}

__device__ static inline void atomicMaxFloat(float *addr, float value)
{
	int *addr_as_int = reinterpret_cast<int*>(addr);
	int old = *addr_as_int;
	while (__int_as_float(old) < value) {
		int assumed = old;
		old = atomicCAS(addr_as_int, assumed, __float_as_int(value));
		if (assumed == old) break;
	}
}

__device__ static inline double atomicAdd_double(double *addr, double val)
{
#if __CUDA_ARCH__ >= 600
	return atomicAdd(addr, val);
#else
	unsigned long long int *address_as_ull = reinterpret_cast<unsigned long long int*>(addr);
	unsigned long long int old = *address_as_ull, assumed;
	do {
		assumed = old;
		double new_val = __longlong_as_double(assumed) + val;
		old = atomicCAS(address_as_ull, assumed, __double_as_longlong(new_val));
	} while (assumed != old);
	return __longlong_as_double(old);
#endif
}

__device__ static inline bool block_contains(const uint64_t *table, uint32_t mask, uint64_t key)
{
	uint32_t slot = uint32_t(fdg_hash64(key)) & mask;
	for (uint32_t iter = 0; iter <= mask; ++iter) {
		uint64_t stored = FDG_LDG_UINT64(&table[slot]);
		if (stored == 0ULL) return false;
		if (stored == key) return true;
		slot = (slot + 1u) & mask;
	}
	return false;
}

__device__ __forceinline__ void warp_sum_match3(unsigned active_mask,
												 unsigned match_mask,
												 float vx,
												 float vy,
												 float vz,
												 float &sum_x,
												 float &sum_y,
												 float &sum_z)
{
	sum_x = 0.0f;
	sum_y = 0.0f;
	sum_z = 0.0f;
	unsigned lanes = match_mask;
	while (lanes) {
		int src_lane = __ffs(lanes) - 1;
		sum_x += __shfl_sync(active_mask, vx, src_lane);
		sum_y += __shfl_sync(active_mask, vy, src_lane);
		sum_z += __shfl_sync(active_mask, vz, src_lane);
		lanes &= (lanes - 1u);
	}
}

/* Optimized bounds kernel: warp+block reduction, then one atomic per block */
__global__ static void bounds_kernel(const hk_gpu_vec3_t *__restrict__ pos, int32_t n, float *__restrict__ bounds)
{
	float lmin_x = FLT_MAX, lmin_y = FLT_MAX, lmin_z = FLT_MAX;
	float lmax_x = -FLT_MAX, lmax_y = -FLT_MAX, lmax_z = -FLT_MAX;

	for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < n; i += blockDim.x * gridDim.x) {
		hk_gpu_vec3_t p = pos[i];
		lmin_x = fminf(lmin_x, p.x); lmax_x = fmaxf(lmax_x, p.x);
		lmin_y = fminf(lmin_y, p.y); lmax_y = fmaxf(lmax_y, p.y);
		lmin_z = fminf(lmin_z, p.z); lmax_z = fmaxf(lmax_z, p.z);
	}

	/* warp reduction */
	for (int offset = 16; offset > 0; offset >>= 1) {
		lmin_x = fminf(lmin_x, __shfl_down_sync(0xffffffff, lmin_x, offset));
		lmin_y = fminf(lmin_y, __shfl_down_sync(0xffffffff, lmin_y, offset));
		lmin_z = fminf(lmin_z, __shfl_down_sync(0xffffffff, lmin_z, offset));
		lmax_x = fmaxf(lmax_x, __shfl_down_sync(0xffffffff, lmax_x, offset));
		lmax_y = fmaxf(lmax_y, __shfl_down_sync(0xffffffff, lmax_y, offset));
		lmax_z = fmaxf(lmax_z, __shfl_down_sync(0xffffffff, lmax_z, offset));
	}

	__shared__ float s_min[3][8];
	__shared__ float s_max[3][8];
	int lane = threadIdx.x & 31;
	int warp_id = threadIdx.x >> 5;
	int n_warps = blockDim.x >> 5;

	if (lane == 0) {
		s_min[0][warp_id] = lmin_x; s_min[1][warp_id] = lmin_y; s_min[2][warp_id] = lmin_z;
		s_max[0][warp_id] = lmax_x; s_max[1][warp_id] = lmax_y; s_max[2][warp_id] = lmax_z;
	}
	__syncthreads();

	if (threadIdx.x == 0) {
		for (int w = 1; w < n_warps; ++w) {
			lmin_x = fminf(lmin_x, s_min[0][w]);
			lmin_y = fminf(lmin_y, s_min[1][w]);
			lmin_z = fminf(lmin_z, s_min[2][w]);
			lmax_x = fmaxf(lmax_x, s_max[0][w]);
			lmax_y = fmaxf(lmax_y, s_max[1][w]);
			lmax_z = fmaxf(lmax_z, s_max[2][w]);
		}
		atomicMinFloat(&bounds[0], lmin_x);
		atomicMinFloat(&bounds[1], lmin_y);
		atomicMinFloat(&bounds[2], lmin_z);
		atomicMaxFloat(&bounds[3], lmax_x);
		atomicMaxFloat(&bounds[4], lmax_y);
		atomicMaxFloat(&bounds[5], lmax_z);
	}
}

__global__ static void insert_block_keys_kernel(uint64_t *table, uint32_t mask, const uint64_t *keys, size_t n_keys)
{
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	for (size_t i = idx; i < n_keys; i += blockDim.x * (size_t)gridDim.x) {
		uint64_t key = keys[i];
		if (key == 0ULL) continue;
		uint32_t slot = uint32_t(fdg_hash64(key)) & mask;
		while (true) {
			uint64_t prev = atomicCAS(reinterpret_cast<unsigned long long*>(&table[slot]), 0ULL, key);
			if (prev == 0ULL || prev == key) break;
			slot = (slot + 1u) & mask;
		}
	}
}


__global__ static void fdg_insert_beads_kernel(const hk_gpu_vec3_t *__restrict__ pos,
							   int32_t n_beads,
							   const float *__restrict__ bounds,
							   float unit,
							   float rep_radius,
							   float4 *__restrict__ grid_params,
							   uint64_t *__restrict__ cell_hash_keys,
							   int *__restrict__ cell_heads,
							   int *__restrict__ bead_next,
							   uint32_t cell_hash_mask)
{
	const uint64_t EMPTY_KEY = 0xffffffffffffffffULL;

	/* Thread 0 computes grid params (replaces prepare_grid<<<1,1>>>) */
	__shared__ float4 s_grid;
	if (threadIdx.x == 0) {
		float cell_size = unit * rep_radius;
		float inv_cell = (cell_size > 0.0f) ? 1.0f / cell_size : 0.0f;
		float min_x = bounds[0] == FLT_MAX ? 0.0f : bounds[0];
		float min_y = bounds[1] == FLT_MAX ? 0.0f : bounds[1];
		float min_z = bounds[2] == FLT_MAX ? 0.0f : bounds[2];
		s_grid = make_float4(min_x, min_y, min_z, inv_cell);
		grid_params[0] = s_grid;
	}
	__syncthreads();

	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx >= n_beads) return;

	float inv_cell = s_grid.w;
	if (inv_cell <= 0.0f) {
		bead_next[idx] = -1;
		return;
	}
	float3 origin = make_float3(s_grid.x, s_grid.y, s_grid.z);
	hk_gpu_vec3_t p = pos[idx];
	int cx = int(floorf((p.x - origin.x) * inv_cell));
	int cy = int(floorf((p.y - origin.y) * inv_cell));
	int cz = int(floorf((p.z - origin.z) * inv_cell));
	if (cx < 0) cx = 0;
	if (cy < 0) cy = 0;
	if (cz < 0) cz = 0;
	uint64_t key = pack_cell(cx, cy, cz);
	uint32_t slot = uint32_t(fdg_hash64(key)) & cell_hash_mask;
	while (true) {
		unsigned long long *addr = reinterpret_cast<unsigned long long*>(&cell_hash_keys[slot]);
		uint64_t prev = atomicCAS(addr, EMPTY_KEY, key);
		if (prev == EMPTY_KEY || prev == key) {
			int old_head = atomicExch(&cell_heads[slot], idx);
			bead_next[idx] = old_head;
			break;
		}
		slot = (slot + 1u) & cell_hash_mask;
	}
}

__global__ static void fdg_force_kernel(const hk_gpu_vec3_t *__restrict__ pos,
										hk_gpu_vec3_t *__restrict__ force,
										const struct hk_fdg_gpu_pair *__restrict__ pairs,
										size_t n_pairs,
										struct hk_fdg_conf opt,
										float unit,
										struct hk_fdg_gpu_stats *__restrict__ stats,
										int32_t n_beads)
{
	size_t idx = blockIdx.x * (size_t)blockDim.x + threadIdx.x;
	if (idx >= n_pairs) return;

	struct hk_fdg_gpu_pair pair = pairs[idx];
	if (pair.type >= HK_FDG_PAIR_TYPE_COUNT || pair.type == HK_FDG_PAIR_TYPE_REPEL) return;
	if (pair.d_scale <= 0.0f) return;

	int i = pair.i;
	int j = pair.j;
	if (i < 0 || j < 0 || i >= n_beads || j >= n_beads) return;

	hk_gpu_vec3_t pi = pos[i];
	hk_gpu_vec3_t pj = pos[j];
	float dx = pi.x - pj.x;
	float dy = pi.y - pj.y;
	float dz = pi.z - pj.z;
	float dist_sq = dx * dx + dy * dy + dz * dz;
	if (dist_sq == 0.0f) return;

	float dist = sqrtf(dist_sq);
	float inv_dist = 1.0f / dist;
	float dirx = dx * inv_dist;
	float diry = dy * inv_dist;
	float dirz = dz * inv_dist;
	float dist_norm = dist / unit;
	float r = dist_norm / pair.d_scale;
	float energy = 0.0f;
	float force_mag = 0.0f;

	if (pair.type == HK_FDG_PAIR_TYPE_BACKBONE) {
		if (r < opt.d_b1) {
			float t = opt.d_b1 - r;
			energy = pair.k * t * t;
			force_mag = 2.0f * pair.k * t;
		} else if (r > opt.d_b2) {
			float t = r - opt.d_b2;
			energy = pair.k * t * t;
			force_mag = -2.0f * pair.k * t;
		}
	} else {
		if (r < opt.d_c1) {
			float t = opt.d_c1 - r;
			energy = pair.k * t * t;
			force_mag = 2.0f * pair.k * t;
		} else if (r > opt.d_c2 && r <= opt.d_c3) {
			float t = r - opt.d_c2;
			energy = pair.k * t * t;
			force_mag = -2.0f * pair.k * t;
		} else if (r > opt.d_c3) {
			float t = r - opt.d_c2;
			float term = opt.c_c1 * (r - opt.d_c3) + opt.c_c2 / t;
			energy = pair.k * term;
			force_mag = -pair.k * (opt.c_c1 - opt.c_c2 / (t * t));
		}
	}

	if (stats != nullptr) {
		if (energy > 0.0f)
			atomicAdd(&stats->energy[pair.type], energy);
		if (force_mag != 0.0f)
			atomicAdd(&stats->dist_sum[pair.type], dist_norm / pair.d_scale);
	}

	if (force_mag != 0.0f) {
		float fx = dirx * force_mag;
		float fy = diry * force_mag;
		float fz = dirz * force_mag;
		atomicAdd(&force[i].x, fx);
		atomicAdd(&force[i].y, fy);
		atomicAdd(&force[i].z, fz);
		atomicAdd(&force[j].x, -fx);
		atomicAdd(&force[j].y, -fy);
		atomicAdd(&force[j].z, -fz);
		if (stats != nullptr)
			atomicAdd(&stats->active[pair.type], 1u);
	}
}

__global__ static void fdg_repulsion_kernel(const hk_gpu_vec3_t *__restrict__ pos,
							 hk_gpu_vec3_t *__restrict__ force,
							 int32_t n_beads,
							 float unit,
							 float rep_radius,
							 float k_rep,
							 const uint64_t *__restrict__ cell_hash_keys,
							 const int *__restrict__ cell_hash_vals,
							 const int *__restrict__ bead_next,
							 uint32_t cell_hash_mask,
							 const uint64_t *__restrict__ block_table,
							 uint32_t block_mask,
							 struct hk_fdg_gpu_stats *__restrict__ stats,
							 const float4 *__restrict__ grid_params)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx >= n_beads) return;

	__shared__ float4 grid_shared;
	if (threadIdx.x == 0)
		grid_shared = grid_params[0];
	__syncthreads();

	float4 grid = grid_shared;
	float inv_cell = grid.w;
	if (inv_cell <= 0.0f)
		return;
	float3 origin = make_float3(grid.x, grid.y, grid.z);

	hk_gpu_vec3_t pi = fdg_load_vec3(pos, idx);
	float fx_acc = 0.0f;
	float fy_acc = 0.0f;
	float fz_acc = 0.0f;
	float energy_acc = 0.0f;
	float dist_acc = 0.0f;
	uint32_t active_cnt = 0;

	const float inv_unit = 1.0f / unit;
	const float rep_limit = rep_radius * unit;
	const float rep_limit_sq = rep_limit * rep_limit;
	const float two_k_rep = 2.0f * k_rep;

	int cx = __float2int_rd((pi.x - origin.x) * inv_cell);
	int cy = __float2int_rd((pi.y - origin.y) * inv_cell);
	int cz = __float2int_rd((pi.z - origin.z) * inv_cell);
	if (cx < 0) cx = 0;
	if (cy < 0) cy = 0;
	if (cz < 0) cz = 0;

	for (int dz = -1; dz <= 1; ++dz) {
		for (int dy = -1; dy <= 1; ++dy) {
			for (int dx = -1; dx <= 1; ++dx) {
				int ncx = cx + dx;
				int ncy = cy + dy;
				int ncz = cz + dz;
				if (ncx < 0 || ncy < 0 || ncz < 0) continue;
				uint64_t nkey = pack_cell(ncx, ncy, ncz);
				int slot = fdg_cell_find_slot(cell_hash_keys, cell_hash_mask, nkey);
				if (slot < 0) continue;
				int head = FDG_LDG_INT(&cell_hash_vals[slot]);
				while (head >= 0) {
					int j = head;
					head = FDG_LDG_INT(&bead_next[j]);
					if (j <= idx) continue;
					hk_gpu_vec3_t pj = fdg_load_vec3(pos, j);
					float dx2 = pi.x - pj.x;
					float dy2 = pi.y - pj.y;
					float dz2 = pi.z - pj.z;
					float dist_sq = dx2 * dx2 + dy2 * dy2 + dz2 * dz2;
					if (dist_sq == 0.0f || dist_sq >= rep_limit_sq) continue;
					float inv_dist = rsqrtf(dist_sq);
					float dist = dist_sq * inv_dist;
					float dist_unit = dist * inv_unit;
					uint64_t pair_key = (uint64_t(idx) << 32) | uint32_t(j);
					if (block_table && block_contains(block_table, block_mask, pair_key)) continue;
					float t = rep_radius - dist_unit;
					float energy = k_rep * t * t;
					float force_mag = two_k_rep * t;
					float scale = force_mag * inv_dist;
					float fx = dx2 * scale;
					float fy = dy2 * scale;
					float fz = dz2 * scale;
					fx_acc += fx;
					fy_acc += fy;
					fz_acc += fz;

#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ >= 700
					float fx_other = -fx;
					float fy_other = -fy;
					float fz_other = -fz;
					unsigned mask = __activemask();
					unsigned match = __match_any_sync(mask, j);
					float fx_total, fy_total, fz_total;
					warp_sum_match3(mask, match, fx_other, fy_other, fz_other, fx_total, fy_total, fz_total);
					int lane = threadIdx.x & 31;
					int leader = __ffs(match) - 1;
					if (lane == leader) {
						atomicAdd(&force[j].x, fx_total);
						atomicAdd(&force[j].y, fy_total);
						atomicAdd(&force[j].z, fz_total);
					}
#else
					atomicAdd(&force[j].x, -fx);
					atomicAdd(&force[j].y, -fy);
					atomicAdd(&force[j].z, -fz);
#endif
					energy_acc += energy;
					dist_acc += dist_unit;
					++active_cnt;
				}
			}
		}
	}

	if (fx_acc != 0.0f || fy_acc != 0.0f || fz_acc != 0.0f) {
		atomicAdd(&force[idx].x, fx_acc);
		atomicAdd(&force[idx].y, fy_acc);
		atomicAdd(&force[idx].z, fz_acc);
	}

	if (stats != nullptr) {
		if (energy_acc != 0.0f)
			atomicAdd(&stats->energy[HK_FDG_PAIR_TYPE_REPEL], energy_acc);
		if (dist_acc != 0.0f)
			atomicAdd(&stats->dist_sum[HK_FDG_PAIR_TYPE_REPEL], dist_acc);
		if (active_cnt)
			atomicAdd(&stats->active[HK_FDG_PAIR_TYPE_REPEL], active_cnt);
	}
}

__global__ static void fdg_update_kernel(hk_gpu_vec3_t *__restrict__ pos,
										 hk_gpu_vec3_t *__restrict__ prev_pos,
										 hk_gpu_vec3_t *__restrict__ force,
										 int32_t n_beads,
										 float step,
										 float coef_moment,
										 float max_f,
										 struct hk_fdg_gpu_stats *__restrict__ stats)
{
	extern __shared__ double s_force_sq[];
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int lane = threadIdx.x;
	double force_sq_acc = 0.0;

	if (idx < n_beads) {
		hk_gpu_vec3_t curr = fdg_load_vec3(pos, idx);
		hk_gpu_vec3_t prev = fdg_load_vec3(prev_pos, idx);
		hk_gpu_vec3_t f = fdg_load_vec3(force, idx);

		float fx = f.x;
		float fy = f.y;
		float fz = f.z;
		float force_sq = fx * fx + fy * fy + fz * fz;

		if (force_sq > 0.0f) {
			if (max_f > 0.0f) {
				float max_f_sq = max_f * max_f;
				if (force_sq > max_f_sq) {
					float inv_force_mag = rsqrtf(force_sq);
					float scale = max_f * inv_force_mag;
					fx *= scale;
					fy *= scale;
					fz *= scale;
					force_sq = max_f_sq;
				}
			}
			force_sq_acc = (double)force_sq;
		}

		float dx_prev = curr.x - prev.x;
		float dy_prev = curr.y - prev.y;
		float dz_prev = curr.z - prev.z;

		float next_x = fmaf(step, fx, fmaf(coef_moment, dx_prev, curr.x));
		float next_y = fmaf(step, fy, fmaf(coef_moment, dy_prev, curr.y));
		float next_z = fmaf(step, fz, fmaf(coef_moment, dz_prev, curr.z));

		fdg_store_vec3(prev_pos, idx, curr.x, curr.y, curr.z);
		fdg_store_vec3(pos, idx, next_x, next_y, next_z);
		fdg_store_vec3(force, idx, 0.0f, 0.0f, 0.0f);
	}

	s_force_sq[lane] = force_sq_acc;
	__syncthreads();

	for (int offset = blockDim.x >> 1; offset > 0; offset >>= 1) {
		if (lane < offset)
			s_force_sq[lane] += s_force_sq[lane + offset];
		__syncthreads();
	}

	if (lane == 0 && stats != nullptr) {
		double block_sum = s_force_sq[0];
		if (block_sum > 0.0)
			atomicAdd_double(&stats->force_sq_sum, block_sum);
	}
}

int hk_fdg_gpu_is_available(void)
{
	int count = 0;
	cudaError_t err = cudaGetDeviceCount(&count);
	if (err != cudaSuccess || count <= 0)
		return 0;
	return 1;
}

struct hk_fdg_gpu_ctx *hk_fdg_gpu_create(int32_t n_beads)
{
	(void)n_beads;
	struct hk_fdg_gpu_ctx *ctx = new (std::nothrow) hk_fdg_gpu_ctx;
	if (ctx == nullptr)
		return nullptr;
	std::memset(ctx, 0, sizeof(*ctx));
	ctx->stream = 0;
	return ctx;
}

void hk_fdg_gpu_destroy(struct hk_fdg_gpu_ctx *ctx)
{
	if (!ctx) return;
	if (ctx->stream) cudaStreamDestroy(ctx->stream);
	hk_fdg_gpu_free_device(&ctx->d_pos);
	hk_fdg_gpu_free_device(&ctx->d_prev_pos);
	hk_fdg_gpu_free_device(&ctx->d_force);
	hk_fdg_gpu_free_device(&ctx->d_pairs);
	hk_fdg_gpu_free_device(&ctx->d_stats);
	hk_fdg_gpu_free_device(&ctx->d_block_table);
	hk_fdg_gpu_free_device(&ctx->d_block_keys);
	hk_fdg_gpu_free_device(&ctx->d_cell_hash_keys);
	hk_fdg_gpu_free_device(&ctx->d_cell_hash_vals);
	hk_fdg_gpu_free_device(&ctx->d_bead_next);
	hk_fdg_gpu_free_device(&ctx->d_grid_params);
	hk_fdg_gpu_free_device(&ctx->d_bounds);
	if (ctx->h_stats_pinned)
		cudaFreeHost(ctx->h_stats_pinned);
	if (ctx->h_bounds_init_pinned)
		cudaFreeHost(ctx->h_bounds_init_pinned);
	delete ctx;
}

static int hk_fdg_gpu_ensure_bead_capacity(struct hk_fdg_gpu_ctx *ctx, size_t needed)
{
	if (needed <= ctx->bead_capacity)
		return 0;
	size_t new_cap = next_pow2_size(std::max<size_t>(needed, 1));
	if (hk_fdg_gpu_resize_device(&ctx->d_pos, new_cap, "position buffer") != 0)
		return -1;
	if (hk_fdg_gpu_resize_device(&ctx->d_prev_pos, new_cap, "previous position buffer") != 0)
		return -1;
	if (hk_fdg_gpu_resize_device(&ctx->d_force, new_cap, "force buffer") != 0)
		return -1;
	if (hk_fdg_gpu_resize_device(&ctx->d_bead_next, new_cap, "cell linked list buffer") != 0)
		return -1;
	ctx->bead_capacity = new_cap;
	ctx->cell_capacity = new_cap;
	return 0;
}

int hk_fdg_gpu_prepare(struct hk_fdg_gpu_ctx *ctx, int32_t n_beads, size_t n_block_keys)
{
	if (!ctx) return -1;
	ctx->n_beads = n_beads;
	if (hk_fdg_gpu_ensure_bead_capacity(ctx, (size_t)n_beads) != 0)
		return -1;
	if (ctx->d_stats == nullptr) {
		if (hk_fdg_gpu_resize_device(&ctx->d_stats, 1, "stats buffer") != 0)
			return -1;
	}
	if (ctx->d_grid_params == nullptr) {
		if (hk_fdg_gpu_resize_device(&ctx->d_grid_params, 1, "grid parameter buffer") != 0)
			return -1;
	}
	if (ctx->d_bounds == nullptr) {
		if (hk_fdg_gpu_resize_device(&ctx->d_bounds, 6, "bounds buffer") != 0)
			return -1;
	}
	if (n_block_keys > ctx->block_keys_cap) {
		if (hk_fdg_gpu_resize_device(&ctx->d_block_keys, n_block_keys, "blocklist key buffer") != 0)
			return -1;
		ctx->block_keys_cap = n_block_keys;
	}
	if (ctx->h_stats_pinned == nullptr) {
		if (cudaHostAlloc(&ctx->h_stats_pinned, sizeof(*ctx->h_stats_pinned), cudaHostAllocPortable) != cudaSuccess) {
			if (hk_verbose >= 1)
				std::fprintf(stderr, "[E::fdg-gpu] failed to allocate pinned stats buffer\n");
			return -1;
		}
	}
	if (ctx->h_bounds_init_pinned == nullptr) {
		if (cudaHostAlloc(&ctx->h_bounds_init_pinned, 6 * sizeof(float), cudaHostAllocPortable) != cudaSuccess) {
			if (hk_verbose >= 1)
				std::fprintf(stderr, "[E::fdg-gpu] failed to allocate pinned bounds init buffer\n");
			return -1;
		}
		ctx->h_bounds_init_pinned[0] = FLT_MAX;
		ctx->h_bounds_init_pinned[1] = FLT_MAX;
		ctx->h_bounds_init_pinned[2] = FLT_MAX;
		ctx->h_bounds_init_pinned[3] = -FLT_MAX;
		ctx->h_bounds_init_pinned[4] = -FLT_MAX;
		ctx->h_bounds_init_pinned[5] = -FLT_MAX;
	}
	ctx->pairs_ready = 0;
	ctx->active_pairs = 0;
	for (int t = 0; t < HK_FDG_PAIR_TYPE_COUNT; ++t)
		ctx->pair_counts[t] = 0;
	return 0;
}

int hk_fdg_gpu_set_blocklist(struct hk_fdg_gpu_ctx *ctx, const uint64_t *keys, size_t n_keys)
{
	if (!ctx) return -1;
	if (n_keys == 0) {
		ctx->block_table_cap = 0;
		ctx->block_keys_cap = 0;
		hk_fdg_gpu_free_device(&ctx->d_block_table);
		hk_fdg_gpu_free_device(&ctx->d_block_keys);
		return 0;
	}
	if (n_keys > ctx->block_keys_cap) {
		if (hk_fdg_gpu_resize_device(&ctx->d_block_keys, n_keys, "blocklist key buffer") != 0)
			return -1;
		ctx->block_keys_cap = n_keys;
	}
	cudaError_t err = cudaMemcpyAsync(ctx->d_block_keys, keys, n_keys * sizeof(uint64_t), cudaMemcpyHostToDevice, ctx->stream);
	if (hk_fdg_gpu_cuda_check(err, "copy block keys") != 0)
		return -1;

	size_t table_cap = next_pow2_size(std::max<size_t>(n_keys * 2, 4));
	if (table_cap != ctx->block_table_cap) {
		if (hk_fdg_gpu_resize_device(&ctx->d_block_table, table_cap, "blocklist hash table") != 0)
			return -1;
		ctx->block_table_cap = table_cap;
	}
	err = cudaMemsetAsync(ctx->d_block_table, 0, ctx->block_table_cap * sizeof(uint64_t), ctx->stream);
	if (hk_fdg_gpu_cuda_check(err, "clear block table") != 0)
		return -1;
	int threads = 256;
	int blocks = (int)((n_keys + threads - 1) / threads);
	if (blocks > 0) {
		insert_block_keys_kernel<<<blocks, threads, 0, ctx->stream>>>(ctx->d_block_table, (uint32_t)(ctx->block_table_cap - 1),
																	  ctx->d_block_keys, n_keys);
		if (hk_fdg_gpu_check_last_error("insert block keys kernel") != 0)
			return -1;
	}
	err = cudaStreamSynchronize(ctx->stream);
	return hk_fdg_gpu_cuda_check(err, "sync blocklist setup");
}

int hk_fdg_gpu_upload_positions(struct hk_fdg_gpu_ctx *ctx, const fvec3_t *pos_host, int32_t n_beads)
{
	if (!ctx || !pos_host) return -1;
	if (n_beads > ctx->n_beads) return -1;
	size_t bytes = (size_t)n_beads * sizeof(hk_gpu_vec3_t);
	cudaError_t err = cudaMemcpyAsync(ctx->d_pos, pos_host, bytes, cudaMemcpyHostToDevice, ctx->stream);
	if (hk_fdg_gpu_cuda_check(err, "upload positions") != 0)
		return -1;
	err = cudaMemcpyAsync(ctx->d_prev_pos, pos_host, bytes, cudaMemcpyHostToDevice, ctx->stream);
	if (hk_fdg_gpu_cuda_check(err, "upload previous positions") != 0)
		return -1;
	err = cudaMemsetAsync(ctx->d_force, 0, ctx->bead_capacity * sizeof(hk_gpu_vec3_t), ctx->stream);
	if (hk_fdg_gpu_cuda_check(err, "zero force buffer") != 0)
		return -1;
	return 0;
}

int hk_fdg_gpu_download_positions(const struct hk_fdg_gpu_ctx *ctx, fvec3_t *pos_host, int32_t n_beads)
{
	if (!ctx || !pos_host) return -1;
	if (n_beads > ctx->n_beads) return -1;
	size_t bytes = (size_t)n_beads * sizeof(hk_gpu_vec3_t);
	cudaError_t err = cudaMemcpy(pos_host, ctx->d_pos, bytes, cudaMemcpyDeviceToHost);
	return hk_fdg_gpu_cuda_check(err, "download positions");
}

int hk_fdg_gpu_ensure_capacity(struct hk_fdg_gpu_ctx *ctx, size_t n_pairs)
{
	if (!ctx) return -1;
	if (n_pairs <= ctx->capacity_pairs)
		return 0;
	size_t new_cap = next_pow2_size(std::max<size_t>(n_pairs, 1));
	if (hk_fdg_gpu_resize_device(&ctx->d_pairs, new_cap, "pair buffer") != 0)
		return -1;
	ctx->capacity_pairs = new_cap;
	return 0;
}

static int hk_fdg_gpu_prepare_cell_hash(struct hk_fdg_gpu_ctx *ctx, int n_cells)
{
	if (n_cells <= 0) {
		ctx->n_cells = 0;
		return 0;
	}
	size_t hash_cap = next_pow2_size(std::max<size_t>(n_cells * 2, 4));
	if (hash_cap != ctx->cell_hash_cap) {
		if (hk_fdg_gpu_resize_device(&ctx->d_cell_hash_keys, hash_cap, "cell hash keys") != 0)
			return -1;
		if (hk_fdg_gpu_resize_device(&ctx->d_cell_hash_vals, hash_cap, "cell hash values") != 0)
			return -1;
		ctx->cell_hash_cap = hash_cap;
	}
	cudaError_t err = cudaMemsetAsync(ctx->d_cell_hash_keys, 0xff, ctx->cell_hash_cap * sizeof(uint64_t), ctx->stream);
	if (hk_fdg_gpu_cuda_check(err, "clear cell hash keys") != 0)
		return -1;
	err = cudaMemsetAsync(ctx->d_cell_hash_vals, 0xff, ctx->cell_hash_cap * sizeof(int), ctx->stream);
	if (hk_fdg_gpu_cuda_check(err, "clear cell hash values") != 0)
		return -1;
	err = cudaMemsetAsync(ctx->d_bead_next, 0xff, ctx->bead_capacity * sizeof(int), ctx->stream);
	if (hk_fdg_gpu_cuda_check(err, "clear bead adjacency") != 0)
		return -1;
	ctx->n_cells = n_cells;
	return 0;
}

static inline int div_up_int(int a, int b)
{
	return (a + b - 1) / b;
}

int hk_fdg_gpu_compute(struct hk_fdg_gpu_ctx *ctx,
					   const struct hk_fdg_conf *opt,
					   const struct hk_fdg_gpu_pair *pairs,
					   size_t n_pairs,
					   float unit,
					   float rel_rep_k,
					   float rep_radius,
					   struct hk_fdg_gpu_stats *stats,
					   double *rms_force,
					   int need_sync)
{
	if (!ctx || !opt || !stats || !rms_force) return -1;
	int n_beads = ctx->n_beads;
	if (n_beads <= 0) {
		std::memset(stats, 0, sizeof(*stats));
		*rms_force = 0.0;
		return 0;
	}

	bool upload_pairs = (pairs != nullptr);
	if (upload_pairs) {
		if (hk_fdg_gpu_ensure_capacity(ctx, n_pairs) != 0)
			return -1;
		if (n_pairs > 0) {
			cudaError_t err = cudaMemcpyAsync(ctx->d_pairs, pairs, n_pairs * sizeof(*pairs), cudaMemcpyHostToDevice, ctx->stream);
			if (hk_fdg_gpu_cuda_check(err, "copy pairs") != 0)
				return -1;
		}
		ctx->active_pairs = n_pairs;
		ctx->pairs_ready = 1;
	} else {
		n_pairs = ctx->active_pairs;
		if (n_pairs > ctx->capacity_pairs)
			return -1;
	}

	struct hk_fdg_gpu_stats *stats_dev = need_sync ? ctx->d_stats : nullptr;
	cudaError_t err;
	if (need_sync) {
		err = cudaMemsetAsync(ctx->d_stats, 0, sizeof(*ctx->d_stats), ctx->stream);
		if (hk_fdg_gpu_cuda_check(err, "zero stats") != 0)
			return -1;
	}

	const int threads = 256;
	if (n_pairs > 0) {
		int blocks = div_up_int((int)n_pairs, threads);
		fdg_force_kernel<<<blocks, threads, 0, ctx->stream>>>(ctx->d_pos,
															   ctx->d_force,
															   ctx->d_pairs,
															   n_pairs,
															   *opt,
															   unit,
															   stats_dev,
															   n_beads);
		if (hk_fdg_gpu_check_last_error("force kernel") != 0)
			return -1;
	}

	float rep_k = opt->k_rel_rep * rel_rep_k;
	bool do_repulsion = (rep_radius > 0.0f && rep_k > 0.0f && n_beads > 1);
	auto pick_thread_count = [](int count) -> int {
		if (count >= 4096) return 256;
		if (count >= 1024) return 128;
		if (count >= 256) return 64;
		return 32;
	};
	if (do_repulsion) {
		/* Use pinned memcpy instead of <<<1,1>>> kernel for bounds init */
		err = cudaMemcpyAsync(ctx->d_bounds, ctx->h_bounds_init_pinned, 6 * sizeof(float), cudaMemcpyHostToDevice, ctx->stream);
		if (hk_fdg_gpu_cuda_check(err, "init bounds") != 0)
			return -1;
		int bound_blocks = div_up_int(n_beads, threads);
		bounds_kernel<<<bound_blocks, threads, 0, ctx->stream>>>(ctx->d_pos, n_beads, ctx->d_bounds);
		if (hk_fdg_gpu_check_last_error("bounds kernel") != 0)
			return -1;
		if (hk_fdg_gpu_prepare_cell_hash(ctx, n_beads) != 0)
			return -1;
		int bead_blocks = div_up_int(n_beads, threads);
		fdg_insert_beads_kernel<<<bead_blocks, threads, 0, ctx->stream>>>(ctx->d_pos,
							     n_beads,
							     ctx->d_bounds,
							     unit,
							     rep_radius,
							     ctx->d_grid_params,
							     ctx->d_cell_hash_keys,
							     ctx->d_cell_hash_vals,
							     ctx->d_bead_next,
							     (uint32_t)(ctx->cell_hash_cap - 1));
		if (hk_fdg_gpu_check_last_error("insert beads kernel") != 0)
			return -1;
		uint32_t block_mask = ctx->block_table_cap ? (uint32_t)(ctx->block_table_cap - 1) : 0u;
		float rep_radius_norm = rep_radius;
		int rep_threads = pick_thread_count(n_beads);
		int rep_blocks = div_up_int(n_beads, rep_threads);
		fdg_repulsion_kernel<<<rep_blocks, rep_threads, 0, ctx->stream>>>(ctx->d_pos,
							      ctx->d_force,
							      n_beads,
							      unit,
							      rep_radius_norm,
							      rep_k,
							      ctx->d_cell_hash_keys,
							      ctx->d_cell_hash_vals,
							      ctx->d_bead_next,
							      (uint32_t)(ctx->cell_hash_cap - 1),
							      ctx->d_block_table,
							      block_mask,
							      stats_dev,
							      ctx->d_grid_params);
		if (hk_fdg_gpu_check_last_error("repulsion kernel") != 0)
			return -1;
	}

	float step = opt->step * unit;
	int update_threads = pick_thread_count(n_beads);
	int update_blocks = div_up_int(n_beads, update_threads);
	fdg_update_kernel<<<update_blocks, update_threads, static_cast<size_t>(update_threads) * sizeof(double), ctx->stream>>>(ctx->d_pos,
																  ctx->d_prev_pos,
																  ctx->d_force,
																  n_beads,
																  step,
																  opt->coef_moment,
																  opt->max_f,
																  stats_dev);
	if (hk_fdg_gpu_check_last_error("update kernel") != 0)
		return -1;

	if (need_sync) {
		err = cudaMemcpyAsync(ctx->h_stats_pinned, ctx->d_stats, sizeof(*ctx->h_stats_pinned), cudaMemcpyDeviceToHost, ctx->stream);
		if (hk_fdg_gpu_cuda_check(err, "copy stats to host") != 0)
			return -1;
		err = cudaStreamSynchronize(ctx->stream);
		if (hk_fdg_gpu_cuda_check(err, "sync hk_fdg_gpu_compute") != 0)
			return -1;
		std::memcpy(stats, ctx->h_stats_pinned, sizeof(*stats));

		if (n_beads > 0 && stats->force_sq_sum > 0.0)
			*rms_force = sqrt(stats->force_sq_sum / (double)n_beads);
		else
			*rms_force = 0.0;
	} else {
		/* No sync: caller doesn't need stats this iteration */
		std::memset(stats, 0, sizeof(*stats));
		*rms_force = 0.0;
	}

	return 0;
}

int hk_fdg_gpu_pairs_ready(const struct hk_fdg_gpu_ctx *ctx)
{
	return ctx ? ctx->pairs_ready : 0;
}

size_t hk_fdg_gpu_get_active_pairs(const struct hk_fdg_gpu_ctx *ctx)
{
	return ctx ? ctx->active_pairs : 0;
}

void hk_fdg_gpu_set_pair_totals(struct hk_fdg_gpu_ctx *ctx, const uint32_t totals[HK_FDG_PAIR_TYPE_COUNT])
{
	if (!ctx) return;
	if (totals) {
		for (int i = 0; i < HK_FDG_PAIR_TYPE_COUNT; ++i)
			ctx->pair_counts[i] = totals[i];
	} else {
		for (int i = 0; i < HK_FDG_PAIR_TYPE_COUNT; ++i)
			ctx->pair_counts[i] = 0;
	}
}

void hk_fdg_gpu_get_pair_totals(const struct hk_fdg_gpu_ctx *ctx, uint32_t totals[HK_FDG_PAIR_TYPE_COUNT])
{
	if (!totals) return;
	if (!ctx) {
		for (int i = 0; i < HK_FDG_PAIR_TYPE_COUNT; ++i)
			totals[i] = 0;
		return;
	}
	for (int i = 0; i < HK_FDG_PAIR_TYPE_COUNT; ++i)
		totals[i] = ctx->pair_counts[i];
}
