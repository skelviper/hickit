# add on node03
CUDA_HOME ?= /usr/local/cuda
CUDA_LIB  ?= $(CUDA_HOME)/targets/x86_64-linux/lib

CXX ?= g++
LDFLAGS  += -L$(CUDA_LIB) -Wl,-rpath,$(CUDA_LIB)
LDLIBS   += -lcudart -lm -

#
CFLAGS=		-g -Wall -O2 -Wc++-compat -ffast-math
CPPFLAGS=
INCLUDES=
OBJS=		sdict.o io.o pair.o count.o phase.o bin.o blind.o blind_p9016_cli.o fdg.o image.o view3d.o fdg_gpu_stub.o
PROG=		hickit
LIBS=		-lm -lz
LIBS_GL=
ASAN_FLAG=
CUDA_OBJS=
CXX?=		g++
LINK?=		$(CC)
NVCC?=		$(CUDA_HOME)/bin/nvcc
CUDAFLAGS?=	-O3 -std=c++17 --use_fast_math -lineinfo -gencode arch=compute_89,code=sm_89
TEST_COMMON_SRCS=sdict.c io.c pair.c count.c phase.c bin.c blind.c fdg.c image.c view3d.c fdg_gpu_stub.c
PYTHON_ANALYSIS=bash -lc 'source "$$(conda info --base)/etc/profile.d/conda.sh" && conda activate analysis && python "$$@"' --

ifneq ($(asan),)
	ASAN_FLAG = -fsanitize=address
endif

ifneq ($(gl),)
	CPPFLAGS += -DHAVE_GL
	ifeq ($(shell uname),Darwin)
		CFLAGS += -Wno-deprecated-declarations
		LIBS_GL = -framework OpenGL -framework GLUT
	else
		LIBS_GL = -Wl,-Bstatic -lglut -Wl,-Bdynamic -lGLU -lGL -lXi
	endif
endif

ifeq ($(gpu),1)
	OBJS := $(filter-out fdg_gpu_stub.o,$(OBJS))
	CUDA_OBJS = fdg_gpu.o
	LIBS += -lcudart -lstdc++
	LINK = $(CXX)
endif

.PHONY:all clean depend test test_blind_eval_py test_blind_p9016_cli_smoke test_blind_p9016_prefix_diag test_blind_p9016_prefix_dry_relax test_blind_p9016_prefix_iter_loop test_blind_p9016_prefix_scheduled_loop test_blind_p9016_prefix_output test_blind_p9016_prefix_scaling_output test_blind_p9016_prefix_raw_output_audit test_blind_output_auditor_negative run_blind_p9016_full_cpu hickit-blind-p9016 run_blind_p9016_full_cpu_matrix run_blind_p9016_candidate_rerun run_blind_p9016_prior_rho_2x2 run_blind_p9016_scaffold_stage_diag run_blind_p9016_dscale_ablation run_blind_p9016_mr4mb_scan run_blind_p9016_mr4mb_inter_gate run_blind_p9016_diagnostic_candidates run_blind_p9016_diagnostic_ablation_grid run_blind_p9016_full_cpu_sweep run_blind_p9016_repulsion_strength_sweep run_blind_p9016_repulsion_tradeoff_sweep run_blind_p9016_repulsion_depth_sweep run_blind_p9016_candidate_outputs run_blind_p9016_candidate_curves run_blind_p9016_repul_iter_relax_scan run_blind_p9016_large_repel_iter_relax_scan audit_blind_p9016_full_cpu_output eval_blind_p9016_rep1_fine_grid test_blind_multiresolution test_blind_multiresolution_pipeline
.SUFFIXES:.c .o .cu

.c.o:
		$(CC) -c $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) $< -o $@

.cu.o:
		$(NVCC) $(CUDAFLAGS) $(CPPFLAGS) -I. $(INCLUDES) -c $< -o $@

all:$(PROG)

hickit:$(OBJS) $(CUDA_OBJS) main.o
		$(LINK) $(LDFLAGS) -o $@ $^ $(ASAN_FLAG) $(LIBS_GL) $(LIBS)

test: test_blind_eval_py test_diploid_index test_blind_input test_blind_binning test_blind_multiresolution test_blind_multiresolution_pipeline test_fdg_energy test_blind_posterior test_blind_coord_posterior test_blind_set_update test_blind_smoke test_blind_wedge test_blind_wedge_list test_blind_wedge_force test_blind_prior test_blind_backbone test_blind_repulsion test_blind_relax_step test_blind_relax_loop test_blind_dry_relax test_blind_init_coords test_blind_init_smoke test_blind_homolog_sep test_blind_homolog_batch test_blind_init_homolog_smoke test_blind_gauge test_blind_gauge_bmap test_blind_gauge_smoke test_blind_single_iter_diag test_blind_single_iter_cpu test_blind_iter_loop_cpu test_blind_diag_fixture test_blind_output test_blind_grid_smoke test_blind_p9016_cli_smoke
		./test_diploid_index
		if ./test_diploid_index --invalid-copy >/dev/null 2>&1; then \
			echo "expected invalid copy assertion to fail"; \
			exit 1; \
		fi
		./test_blind_input
		./test_blind_binning
		./test_blind_multiresolution
		./test_blind_multiresolution_pipeline
		./test_fdg_energy
		./test_blind_posterior
		./test_blind_coord_posterior
		./test_blind_set_update
		./test_blind_smoke
		./test_blind_wedge
		./test_blind_wedge_list
		./test_blind_wedge_force
		./test_blind_prior
		./test_blind_backbone
		./test_blind_repulsion
		./test_blind_relax_step
		./test_blind_relax_loop
		./test_blind_dry_relax
		./test_blind_init_coords
		./test_blind_init_smoke
		./test_blind_homolog_sep
		./test_blind_homolog_batch
		./test_blind_init_homolog_smoke
		./test_blind_gauge
		./test_blind_gauge_bmap
		./test_blind_gauge_smoke
		./test_blind_single_iter_diag
		./test_blind_single_iter_cpu
		./test_blind_iter_loop_cpu
		if ./test_blind_iter_loop_cpu --invalid-schedule-temperature >/dev/null 2>&1; then \
			echo "expected invalid schedule temperature assertion to fail"; \
			exit 1; \
		fi
		if ./test_blind_iter_loop_cpu --invalid-schedule-rho >/dev/null 2>&1; then \
			echo "expected invalid schedule rho assertion to fail"; \
			exit 1; \
		fi
		./test_blind_diag_fixture
		./test_blind_output
		./test_blind_grid_smoke

test_blind_eval_py:
		$(PYTHON_ANALYSIS) test_blind_eval.py

test_blind_p9016_cli_smoke: hickit audit_blind_p9016_full_cpu_output.bin testdata/p9016_blind_smoke.pairs
		rm -rf /tmp/hk_blind_p9016_cli_smoke_test
		./hickit blind-p9016 -i testdata/p9016_blind_smoke.pairs -o /tmp/hk_blind_p9016_cli_smoke_test --bd-iter 1 --bd-relax-steps 1 --bd-heldout-frac 0.5 --bd-base-k-mode uniform
		./audit_blind_p9016_full_cpu_output.bin /tmp/hk_blind_p9016_cli_smoke_test

test_diploid_index: test_diploid_index.c hickit.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) $< -o $@

test_blind_input: test_blind_input.c testdata/p9016_blind_fixture.pairs $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_input.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_binning: test_blind_binning.c testdata/p9016_blind_fixture.pairs $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_binning.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_multiresolution: test_blind_multiresolution.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_multiresolution.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_multiresolution_pipeline: test_blind_multiresolution_pipeline.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_multiresolution_pipeline.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_fdg_energy: test_fdg_energy.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_fdg_energy.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_posterior: test_blind_posterior.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_posterior.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_coord_posterior: test_blind_coord_posterior.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_coord_posterior.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_set_update: test_blind_set_update.c testdata/p9016_blind_fixture.pairs $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_set_update.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_smoke: test_blind_smoke.c testdata/p9016_blind_smoke.pairs $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_smoke.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_wedge: test_blind_wedge.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_wedge.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_wedge_list: test_blind_wedge_list.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_wedge_list.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_wedge_force: test_blind_wedge_force.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_wedge_force.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_prior: test_blind_prior.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_prior.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_backbone: test_blind_backbone.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_backbone.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_repulsion: test_blind_repulsion.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_repulsion.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_relax_step: test_blind_relax_step.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_relax_step.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_relax_loop: test_blind_relax_loop.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_relax_loop.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_dry_relax: test_blind_dry_relax.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_dry_relax.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_init_coords: test_blind_init_coords.c testdata/p9016_blind_fixture.pairs $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_init_coords.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_init_smoke: test_blind_init_smoke.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_init_smoke.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_homolog_sep: test_blind_homolog_sep.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_homolog_sep.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_homolog_batch: test_blind_homolog_batch.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_homolog_batch.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_init_homolog_smoke: test_blind_init_homolog_smoke.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_init_homolog_smoke.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_gauge: test_blind_gauge.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_gauge.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_gauge_bmap: test_blind_gauge_bmap.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_gauge_bmap.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_gauge_smoke: test_blind_gauge_smoke.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_gauge_smoke.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_single_iter_diag: test_blind_single_iter_diag.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_single_iter_diag.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_single_iter_cpu: test_blind_single_iter_cpu.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_single_iter_cpu.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_iter_loop_cpu: test_blind_iter_loop_cpu.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_iter_loop_cpu.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_diag_fixture: test_blind_diag_fixture.c testdata/p9016_blind_diag.pairs $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_diag_fixture.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_output: test_blind_output.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_output.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_grid_smoke: test_blind_grid_smoke.c testdata/p9016_blind_fixture.pairs audit_blind_p9016_full_cpu_output.bin $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_grid_smoke.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_p9016_prefix_diag: test_blind_p9016_prefix_diag.bin
		./test_blind_p9016_prefix_diag.bin

test_blind_p9016_prefix_diag.bin: test_blind_p9016_prefix_diag.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_p9016_prefix_diag.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_p9016_prefix_dry_relax: test_blind_p9016_prefix_dry_relax.bin
		./test_blind_p9016_prefix_dry_relax.bin

test_blind_p9016_prefix_dry_relax.bin: test_blind_p9016_prefix_dry_relax.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_p9016_prefix_dry_relax.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_p9016_prefix_iter_loop: test_blind_p9016_prefix_iter_loop.bin
		./test_blind_p9016_prefix_iter_loop.bin

test_blind_p9016_prefix_iter_loop.bin: test_blind_p9016_prefix_iter_loop.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_p9016_prefix_iter_loop.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_p9016_prefix_scheduled_loop: test_blind_p9016_prefix_scheduled_loop.bin
		./test_blind_p9016_prefix_scheduled_loop.bin

test_blind_p9016_prefix_scheduled_loop.bin: test_blind_p9016_prefix_scheduled_loop.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_p9016_prefix_scheduled_loop.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_p9016_prefix_output: test_blind_p9016_prefix_output.bin
		./test_blind_p9016_prefix_output.bin

test_blind_p9016_prefix_output.bin: test_blind_p9016_prefix_output.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_p9016_prefix_output.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_p9016_prefix_scaling_output: test_blind_p9016_prefix_scaling_output.bin
		./test_blind_p9016_prefix_scaling_output.bin

test_blind_p9016_prefix_scaling_output.bin: test_blind_p9016_prefix_scaling_output.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_p9016_prefix_scaling_output.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_p9016_prefix_raw_output_audit: test_blind_p9016_prefix_raw_output_audit.bin audit_blind_p9016_full_cpu_output.bin
		./test_blind_p9016_prefix_raw_output_audit.bin

test_blind_p9016_prefix_raw_output_audit.bin: test_blind_p9016_prefix_raw_output_audit.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_p9016_prefix_raw_output_audit.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_output_auditor_negative: test_blind_output_auditor_negative.bin audit_blind_p9016_full_cpu_output.bin
		./test_blind_output_auditor_negative.bin

test_blind_output_auditor_negative.bin: test_blind_output_auditor_negative.c hickit.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_output_auditor_negative.c -o $@ $(LIBS)

run_blind_p9016_full_cpu: run_blind_p9016_full_cpu.bin
		./run_blind_p9016_full_cpu.bin

run_blind_p9016_full_cpu.bin: run_blind_p9016_full_cpu.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) run_blind_p9016_full_cpu.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

hickit-blind-p9016: run_blind_p9016_full_cpu.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) run_blind_p9016_full_cpu.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

run_blind_p9016_full_cpu_matrix: run_blind_p9016_full_cpu_matrix.bin audit_blind_p9016_full_cpu_output.bin
		./run_blind_p9016_full_cpu_matrix.bin

run_blind_p9016_full_cpu_matrix.bin: run_blind_p9016_full_cpu_matrix.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) run_blind_p9016_full_cpu_matrix.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

run_blind_p9016_candidate_rerun: run_blind_p9016_candidate_rerun.bin audit_blind_p9016_full_cpu_output.bin
		./run_blind_p9016_candidate_rerun.bin

run_blind_p9016_candidate_rerun.bin: run_blind_p9016_full_cpu_matrix.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) -DHK_BLIND_P9016_FULL_OUTPUT_CANDIDATE_RERUN run_blind_p9016_full_cpu_matrix.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

run_blind_p9016_prior_rho_2x2: run_blind_p9016_prior_rho_2x2.bin audit_blind_p9016_full_cpu_output.bin
		./run_blind_p9016_prior_rho_2x2.bin

run_blind_p9016_prior_rho_2x2.bin: run_blind_p9016_full_cpu_matrix.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) -DHK_BLIND_P9016_PRIOR_RHO_2X2 run_blind_p9016_full_cpu_matrix.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

run_blind_p9016_scaffold_stage_diag: run_blind_p9016_scaffold_stage_diag.bin audit_blind_p9016_full_cpu_output.bin
		./run_blind_p9016_scaffold_stage_diag.bin

run_blind_p9016_scaffold_stage_diag.bin: run_blind_p9016_full_cpu_matrix.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) -DHK_BLIND_P9016_SCAFFOLD_STAGE_DIAG run_blind_p9016_full_cpu_matrix.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

run_blind_p9016_dscale_ablation: run_blind_p9016_dscale_ablation.bin audit_blind_p9016_full_cpu_output.bin
		./run_blind_p9016_dscale_ablation.bin

run_blind_p9016_dscale_ablation.bin: run_blind_p9016_full_cpu_matrix.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) -DHK_BLIND_P9016_DSCALE_ABLATION run_blind_p9016_full_cpu_matrix.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

run_blind_p9016_mr4mb_scan: run_blind_p9016_mr4mb_scan.bin audit_blind_p9016_full_cpu_output.bin
		HK_BLIND_P9016_BIN_SIZE_BP=4000000 ./run_blind_p9016_mr4mb_scan.bin

run_blind_p9016_mr4mb_scan.bin: run_blind_p9016_full_cpu_matrix.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) -DHK_BLIND_P9016_MR4MB_SCAN run_blind_p9016_full_cpu_matrix.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

run_blind_p9016_mr4mb_inter_tune: run_blind_p9016_mr4mb_inter_tune.bin audit_blind_p9016_full_cpu_output.bin
		HK_BLIND_P9016_BIN_SIZE_BP=4000000 ./run_blind_p9016_mr4mb_inter_tune.bin

run_blind_p9016_mr4mb_inter_tune.bin: run_blind_p9016_full_cpu_matrix.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) -DHK_BLIND_P9016_MR4MB_INTER_TUNE run_blind_p9016_full_cpu_matrix.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

run_blind_p9016_mr4mb_inter_gate: run_blind_p9016_mr4mb_inter_gate.bin audit_blind_p9016_full_cpu_output.bin
		HK_BLIND_P9016_BIN_SIZE_BP=4000000 ./run_blind_p9016_mr4mb_inter_gate.bin

run_blind_p9016_mr4mb_inter_gate.bin: run_blind_p9016_full_cpu_matrix.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) -DHK_BLIND_P9016_MR4MB_INTER_GATE run_blind_p9016_full_cpu_matrix.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

run_blind_p9016_mr4mb_force_balance: run_blind_p9016_mr4mb_force_balance.bin audit_blind_p9016_full_cpu_output.bin
		HK_BLIND_P9016_BIN_SIZE_BP=4000000 ./run_blind_p9016_mr4mb_force_balance.bin

run_blind_p9016_mr4mb_force_balance.bin: run_blind_p9016_full_cpu_matrix.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) -DHK_BLIND_P9016_MR4MB_FORCE_BALANCE run_blind_p9016_full_cpu_matrix.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

run_blind_p9016_diagnostic_candidates: run_blind_p9016_diagnostic_candidates.bin audit_blind_p9016_full_cpu_output.bin
		./run_blind_p9016_diagnostic_candidates.bin

run_blind_p9016_diagnostic_candidates.bin: run_blind_p9016_full_cpu_matrix.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) -DHK_BLIND_P9016_DIAGNOSTIC_CANDIDATES run_blind_p9016_full_cpu_matrix.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

run_blind_p9016_diagnostic_ablation_grid: run_blind_p9016_diagnostic_ablation_grid.bin audit_blind_p9016_full_cpu_output.bin
		./run_blind_p9016_diagnostic_ablation_grid.bin

run_blind_p9016_diagnostic_ablation_grid.bin: run_blind_p9016_full_cpu_matrix.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) -DHK_BLIND_P9016_DIAGNOSTIC_ABLATION_GRID run_blind_p9016_full_cpu_matrix.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

run_blind_p9016_full_cpu_sweep: run_blind_p9016_full_cpu_sweep.bin
		./run_blind_p9016_full_cpu_sweep.bin

run_blind_p9016_full_cpu_sweep.bin: run_blind_p9016_full_cpu_sweep.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) run_blind_p9016_full_cpu_sweep.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

run_blind_p9016_repulsion_strength_sweep: run_blind_p9016_repulsion_strength_sweep.bin
		./run_blind_p9016_repulsion_strength_sweep.bin

run_blind_p9016_repulsion_strength_sweep.bin: run_blind_p9016_repulsion_strength_sweep.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) run_blind_p9016_repulsion_strength_sweep.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

run_blind_p9016_repulsion_tradeoff_sweep: run_blind_p9016_repulsion_tradeoff_sweep.bin
		./run_blind_p9016_repulsion_tradeoff_sweep.bin

run_blind_p9016_repulsion_tradeoff_sweep.bin: run_blind_p9016_repulsion_strength_sweep.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) -DHK_BLIND_REPULSION_TRADEOFF_SWEEP run_blind_p9016_repulsion_strength_sweep.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

run_blind_p9016_repulsion_depth_sweep: run_blind_p9016_repulsion_depth_sweep.bin
		./run_blind_p9016_repulsion_depth_sweep.bin

run_blind_p9016_repulsion_depth_sweep.bin: run_blind_p9016_repulsion_strength_sweep.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) -DHK_BLIND_REPULSION_DEPTH_SWEEP run_blind_p9016_repulsion_strength_sweep.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

run_blind_p9016_candidate_outputs: run_blind_p9016_candidate_outputs.bin audit_blind_p9016_full_cpu_output.bin
		./run_blind_p9016_candidate_outputs.bin

run_blind_p9016_candidate_outputs.bin: run_blind_p9016_candidate_outputs.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) run_blind_p9016_candidate_outputs.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

run_blind_p9016_candidate_curves: run_blind_p9016_candidate_curves.bin
		./run_blind_p9016_candidate_curves.bin

run_blind_p9016_candidate_curves.bin: run_blind_p9016_candidate_curves.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) run_blind_p9016_candidate_curves.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

run_blind_p9016_repul_iter_relax_scan: run_blind_p9016_repul_iter_relax_scan.bin
		./run_blind_p9016_repul_iter_relax_scan.bin

run_blind_p9016_repul_iter_relax_scan.bin: run_blind_p9016_repul_iter_relax_scan.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) run_blind_p9016_repul_iter_relax_scan.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

run_blind_p9016_large_repel_iter_relax_scan: run_blind_p9016_large_repel_iter_relax_scan.bin
		./run_blind_p9016_large_repel_iter_relax_scan.bin

run_blind_p9016_large_repel_iter_relax_scan.bin: run_blind_p9016_large_repel_iter_relax_scan.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
		$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) run_blind_p9016_large_repel_iter_relax_scan.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

audit_blind_p9016_full_cpu_output: audit_blind_p9016_full_cpu_output.bin
		@if [ -z "$(OUTDIR)" ]; then \
			echo "usage: make audit_blind_p9016_full_cpu_output OUTDIR=/tmp/hk_blind_p9016_full_cpu_..." >&2; \
			exit 1; \
		fi
		./audit_blind_p9016_full_cpu_output.bin "$(OUTDIR)"

audit_blind_p9016_full_cpu_output.bin: audit_blind_p9016_full_cpu_output.c hickit.h krng.h
		$(CC) $(filter-out -ffast-math,$(CFLAGS)) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) audit_blind_p9016_full_cpu_output.c -o $@ $(LIBS)

eval_blind_p9016_rep1_fine_grid:
		$(PYTHON_ANALYSIS) evaluate_blind_p9016_grid.py \
			--summary /tmp/hk_blind_p9016_rep1_fine_n_rs_grid_135846_1777436403_0/scan_summary.tsv \
			--tdg /shared/zliu/CHARM/CHARM_mesc/data/tdg/P9016.1m.3dg.gz \
			--out-root /tmp/hk_blind_p9016_rep1_fine_n_rs_grid_135846_1777436403_0/eval \
			--trans-sample-per-chrpair 2048 \
			--sample-seed 17

clean:
		rm -fr gmon.out *.o a.out $(PROG) hickit-blind-p9016 test_diploid_index test_blind_input test_blind_binning test_blind_multiresolution test_blind_multiresolution_pipeline test_fdg_energy test_blind_posterior test_blind_coord_posterior test_blind_set_update test_blind_smoke test_blind_wedge test_blind_wedge_list test_blind_wedge_force test_blind_backbone test_blind_repulsion test_blind_relax_step test_blind_relax_loop test_blind_dry_relax test_blind_init_coords test_blind_init_smoke test_blind_homolog_sep test_blind_homolog_batch test_blind_init_homolog_smoke test_blind_gauge test_blind_gauge_bmap test_blind_gauge_smoke test_blind_single_iter_diag test_blind_single_iter_cpu test_blind_iter_loop_cpu test_blind_diag_fixture test_blind_output test_blind_p9016_prefix_diag.bin test_blind_p9016_prefix_dry_relax.bin test_blind_p9016_prefix_iter_loop.bin test_blind_p9016_prefix_scheduled_loop.bin test_blind_p9016_prefix_output.bin test_blind_p9016_prefix_scaling_output.bin test_blind_p9016_prefix_raw_output_audit.bin test_blind_output_auditor_negative.bin run_blind_p9016_full_cpu.bin run_blind_p9016_full_cpu_matrix.bin run_blind_p9016_candidate_rerun.bin run_blind_p9016_prior_rho_2x2.bin run_blind_p9016_scaffold_stage_diag.bin run_blind_p9016_dscale_ablation.bin run_blind_p9016_mr4mb_scan.bin run_blind_p9016_mr4mb_inter_gate.bin run_blind_p9016_mr4mb_force_balance.bin run_blind_p9016_diagnostic_candidates.bin run_blind_p9016_diagnostic_ablation_grid.bin run_blind_p9016_full_cpu_sweep.bin run_blind_p9016_repulsion_strength_sweep.bin run_blind_p9016_repulsion_tradeoff_sweep.bin run_blind_p9016_repulsion_depth_sweep.bin run_blind_p9016_candidate_outputs.bin run_blind_p9016_candidate_curves.bin run_blind_p9016_repul_iter_relax_scan.bin run_blind_p9016_large_repel_iter_relax_scan.bin audit_blind_p9016_full_cpu_output.bin *.a *.dSYM hickit.aux hickit.log hickit.pdf

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(CPPFLAGS) -- *.c)

# DO NOT DELETE

bin.o: hkpriv.h hickit.h krng.h khash.h ksort.h
blind.o: hkpriv.h hickit.h krng.h
count.o: hkpriv.h hickit.h krng.h kavl.h klist.h ksort.h
fdg.o: hkpriv.h hickit.h krng.h ksort.h kavl.h khash.h
image.o: hkpriv.h hickit.h krng.h ksort.h stb_image_write.h
io.o: hickit.h krng.h hkpriv.h kseq.h
main.o: hickit.h krng.h
pair.o: hkpriv.h hickit.h krng.h ksort.h
phase.o: hkpriv.h hickit.h krng.h ksort.h
sdict.o: hkpriv.h hickit.h krng.h khash.h
view3d.o: hkpriv.h hickit.h krng.h
fdg_gpu_stub.o: fdg_gpu.h hickit.h
fdg_gpu.o: fdg_gpu.h hickit.h
