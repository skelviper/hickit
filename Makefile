CUDA_HOME ?= /usr/local/cuda
CUDA_LIB  ?= $(CUDA_HOME)/targets/x86_64-linux/lib
NVCC      ?= $(CUDA_HOME)/bin/nvcc

CFLAGS   ?= -g -Wall -O2 -Wc++-compat -fno-fast-math -fno-unsafe-math-optimizations
NVCCFLAGS ?= -O3
CPPFLAGS ?=
INCLUDES ?=
LDFLAGS  += -L$(CUDA_LIB) -Wl,-rpath,$(CUDA_LIB)
LIBS     := -lm -lz
AUDIT_LIBS := -lm -lz
LIBS_GL  :=
ASAN_FLAG :=
TEST_RES_ROOT ?= $(if $(HK_BLIND_TEST_RES_ROOT),$(HK_BLIND_TEST_RES_ROOT),$(CURDIR)/test_res)
SMOKE_RES_ROOT ?= $(TEST_RES_ROOT)/smoke
P9016_MINIMAL_SMOKE_ROOT := $(SMOKE_RES_ROOT)/hk_blind_p9016_minimal_audit_smoke
P9016_SEP_ON_SMOKE_ROOT := $(SMOKE_RES_ROOT)/hk_blind_p9016_sep_on_smoke
P9016_EXPECTED_SEP_SMOKE_ROOT := $(SMOKE_RES_ROOT)/hk_blind_p9016_expected_sep_smoke
P9016_RESOLUTION_CHAIN_SMOKE_ROOT := $(SMOKE_RES_ROOT)/hk_blind_p9016_resolution_chain_smoke

ifeq ($(gpu),1)
	FDG_GPU_OBJ := fdg_gpu.o
	FDG_GPU_STAMP := .fdg_gpu_backend.gpu
	LIBS += -lcudart -lstdc++
else
	FDG_GPU_OBJ := fdg_gpu_stub.o
	FDG_GPU_STAMP := .fdg_gpu_backend.cpu
endif

OBJS := sdict.o io.o pair.o count.o phase.o bin.o fdg.o image.o view3d.o $(FDG_GPU_OBJ)
PROG := hickit
BLIND_COMMON_SRCS := sdict.c io.c pair.c count.c phase.c bin.c blind.c fdg.c image.c view3d.c
BLIND_COMMON_OBJS := sdict.o io.o pair.o count.o phase.o bin.o blind.o fdg.o image.o view3d.o $(FDG_GPU_OBJ)
ifneq ($(asan),)
	ASAN_FLAG = -fsanitize=address,undefined -fno-omit-frame-pointer
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

.PHONY: all clean depend test smoke_blind_p9016_minimal smoke_blind_p9016_sep_on smoke_blind_p9016_expected_sep smoke_blind_p9016_resolution_chain audit_blind_p9016_full_cpu_output
.SUFFIXES: .c .o

.c.o:
	$(CC) -c $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) $< -o $@

fdg_gpu.o: fdg_gpu.cu fdg_gpu.h hickit.h
	$(NVCC) -c $(NVCCFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

$(FDG_GPU_STAMP):
	@rm -f .fdg_gpu_backend.*
	@touch $@

all: $(PROG)

hickit: $(OBJS) main.o
	$(CC) $(LDFLAGS) -o $@ $^ $(ASAN_FLAG) $(LIBS_GL) $(LIBS)

test: smoke_blind_p9016_minimal

run_blind_p9016_minimal.bin: run_blind_p9016_minimal.c $(BLIND_COMMON_OBJS) $(FDG_GPU_STAMP) hickit.h hkpriv.h krng.h
	$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) run_blind_p9016_minimal.c $(BLIND_COMMON_OBJS) -o $@ $(LDFLAGS) $(LIBS)

run_blind_p9016_fixed_posterior_relax.bin: run_blind_p9016_fixed_posterior_relax.c $(BLIND_COMMON_OBJS) $(FDG_GPU_STAMP) hickit.h hkpriv.h krng.h
	$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) run_blind_p9016_fixed_posterior_relax.c $(BLIND_COMMON_OBJS) -o $@ $(LDFLAGS) $(LIBS)

smoke_blind_p9016_minimal: run_blind_p9016_minimal.bin audit_blind_p9016_full_cpu_output.bin testdata/p9016_blind_smoke.pairs
	rm -rf $(P9016_MINIMAL_SMOKE_ROOT)
	mkdir -p $(SMOKE_RES_ROOT)
	env -i PATH="$$PATH" LD_LIBRARY_PATH="$$LD_LIBRARY_PATH" \
	HK_BLIND_P9016_PAIRS=testdata/p9016_blind_smoke.pairs \
	HK_BLIND_P9016_ALLOW_CUSTOM_PAIRS=1 \
	HK_BLIND_P9016_OUTPUT_ROOT=$(P9016_MINIMAL_SMOKE_ROOT) \
	HK_BLIND_P9016_CONFIG_NAME=minimal_soft_sep_off \
	HK_BLIND_P9016_BIN_SIZE_BP=1000000 \
	HK_BLIND_P9016_MINIMAL_N_ITER=1 \
	HK_BLIND_P9016_MINIMAL_RELAX_STEPS=1 \
	HK_BLIND_P9016_INIT_EPS=0.5 \
	HK_BLIND_P9016_INIT_NOISE_SCALE=0.0 \
	HK_BLIND_P9016_MIN_SEP_UNIT=0.0 \
	HK_BLIND_P9016_LAMBDA_SEP=0.0 \
	HK_BLIND_P9016_INIT_SEED=17 \
	./run_blind_p9016_minimal.bin
	HK_BLIND_AUDIT_ALLOW_CUSTOM_RAW_PAIRS=1 ./audit_blind_p9016_full_cpu_output.bin $(P9016_MINIMAL_SMOKE_ROOT)/minimal_soft_sep_off

smoke_blind_p9016_sep_on: run_blind_p9016_minimal.bin audit_blind_p9016_full_cpu_output.bin testdata/p9016_blind_smoke.pairs
	rm -rf $(P9016_SEP_ON_SMOKE_ROOT)
	mkdir -p $(SMOKE_RES_ROOT)
	env -i PATH="$$PATH" LD_LIBRARY_PATH="$$LD_LIBRARY_PATH" \
	HK_BLIND_P9016_PAIRS=testdata/p9016_blind_smoke.pairs \
	HK_BLIND_P9016_ALLOW_CUSTOM_PAIRS=1 \
	HK_BLIND_P9016_OUTPUT_ROOT=$(P9016_SEP_ON_SMOKE_ROOT) \
	HK_BLIND_P9016_CONFIG_NAME=minimal_soft_sep_on_smoke \
	HK_BLIND_P9016_BIN_SIZE_BP=1000000 \
	HK_BLIND_P9016_MINIMAL_N_ITER=1 \
	HK_BLIND_P9016_MINIMAL_RELAX_STEPS=1 \
	HK_BLIND_P9016_INIT_EPS=1.0 \
	HK_BLIND_P9016_INIT_NOISE_SCALE=0.05 \
	HK_BLIND_P9016_MIN_SEP_UNIT=2.0 \
	HK_BLIND_P9016_LAMBDA_SEP=0.05 \
	./run_blind_p9016_minimal.bin
	HK_BLIND_AUDIT_ALLOW_CUSTOM_RAW_PAIRS=1 ./audit_blind_p9016_full_cpu_output.bin $(P9016_SEP_ON_SMOKE_ROOT)/minimal_soft_sep_on_smoke
	awk -F'\t' 'function abs(x){return x<0?-x:x} $$1=="min_sep_unit" && abs(($$2+0)-2.0)<1e-6 { found=1 } END { exit found?0:1 }' $(P9016_SEP_ON_SMOKE_ROOT)/minimal_soft_sep_on_smoke/p9016_full.manifest.tsv
	awk -F'\t' 'function abs(x){return x<0?-x:x} $$1=="lambda_sep" && abs(($$2+0)-0.05)<1e-6 { found=1 } END { exit found?0:1 }' $(P9016_SEP_ON_SMOKE_ROOT)/minimal_soft_sep_on_smoke/p9016_full.manifest.tsv

smoke_blind_p9016_expected_sep: run_blind_p9016_minimal.bin audit_blind_p9016_full_cpu_output.bin testdata/p9016_blind_smoke.pairs
	rm -rf $(P9016_EXPECTED_SEP_SMOKE_ROOT)
	mkdir -p $(SMOKE_RES_ROOT)
	env -i PATH="$$PATH" LD_LIBRARY_PATH="$$LD_LIBRARY_PATH" \
	HK_BLIND_P9016_PAIRS=testdata/p9016_blind_smoke.pairs \
	HK_BLIND_P9016_ALLOW_CUSTOM_PAIRS=1 \
	HK_BLIND_P9016_OUTPUT_ROOT=$(P9016_EXPECTED_SEP_SMOKE_ROOT) \
	HK_BLIND_P9016_CONFIG_NAME=scaffold_expected_msep1p5_lsep0p5_smoke \
	HK_BLIND_P9016_INIT_MODE=unphased_scaffold_split \
	HK_BLIND_P9016_MINIMAL_N_ITER=1 \
	HK_BLIND_P9016_MINIMAL_RELAX_STEPS=1 \
	HK_BLIND_P9016_INIT_EPS=1.0 \
	HK_BLIND_P9016_INIT_NOISE_SCALE=0.05 \
	HK_BLIND_P9016_MIN_SEP_UNIT=1.5 \
	HK_BLIND_P9016_LAMBDA_SEP=0.5 \
	HK_BLIND_P9016_D_SCALE_MODE=expected_count \
	HK_BLIND_P9016_D_SCALE_EPS_COUNT=0.001 \
	./run_blind_p9016_minimal.bin
	HK_BLIND_AUDIT_ALLOW_CUSTOM_RAW_PAIRS=1 ./audit_blind_p9016_full_cpu_output.bin $(P9016_EXPECTED_SEP_SMOKE_ROOT)/scaffold_expected_msep1p5_lsep0p5_smoke
	awk -F'\t' '$$1=="d_scale_mode" && $$2=="posterior_count" { found=1 } END { exit found?0:1 }' $(P9016_EXPECTED_SEP_SMOKE_ROOT)/scaffold_expected_msep1p5_lsep0p5_smoke/p9016_full.manifest.tsv
	awk -F'\t' '$$1=="d_scale_mode_input_string" && $$2=="expected_count" { found=1 } END { exit found?0:1 }' $(P9016_EXPECTED_SEP_SMOKE_ROOT)/scaffold_expected_msep1p5_lsep0p5_smoke/p9016_full.manifest.tsv
	awk -F'\t' '$$1=="legacy_expected_count_alias_used" && $$2=="1" { found=1 } END { exit found?0:1 }' $(P9016_EXPECTED_SEP_SMOKE_ROOT)/scaffold_expected_msep1p5_lsep0p5_smoke/p9016_full.manifest.tsv
	awk -F'\t' 'function abs(x){return x<0?-x:x} $$1=="d_scale_eps_count" && abs(($$2+0)-0.001)<1e-8 { found=1 } END { exit found?0:1 }' $(P9016_EXPECTED_SEP_SMOKE_ROOT)/scaffold_expected_msep1p5_lsep0p5_smoke/p9016_full.manifest.tsv
	awk -F'\t' 'function abs(x){return x<0?-x:x} $$1=="min_sep_unit" && abs(($$2+0)-1.5)<1e-6 { found=1 } END { exit found?0:1 }' $(P9016_EXPECTED_SEP_SMOKE_ROOT)/scaffold_expected_msep1p5_lsep0p5_smoke/p9016_full.manifest.tsv
	awk -F'\t' 'function abs(x){return x<0?-x:x} $$1=="lambda_sep" && abs(($$2+0)-0.5)<1e-6 { found=1 } END { exit found?0:1 }' $(P9016_EXPECTED_SEP_SMOKE_ROOT)/scaffold_expected_msep1p5_lsep0p5_smoke/p9016_full.manifest.tsv
	awk -F'\t' '$$1=="init_mode" && $$2=="unphased_scaffold_split" { found=1 } END { exit found?0:1 }' $(P9016_EXPECTED_SEP_SMOKE_ROOT)/scaffold_expected_msep1p5_lsep0p5_smoke/p9016_full.manifest.tsv
	awk -F'\t' 'function abs(x){return x<0?-x:x} $$1=="init_eps_effective" && ($$2+0)>0.0 { found=1 } END { exit found?0:1 }' $(P9016_EXPECTED_SEP_SMOKE_ROOT)/scaffold_expected_msep1p5_lsep0p5_smoke/p9016_full.manifest.tsv
	awk -F'\t' 'function abs(x){return x<0?-x:x} $$1=="init_noise_scale_effective" && abs(($$2+0)-0.05)<1e-6 { found=1 } END { exit found?0:1 }' $(P9016_EXPECTED_SEP_SMOKE_ROOT)/scaffold_expected_msep1p5_lsep0p5_smoke/p9016_full.manifest.tsv
	test -s $(P9016_EXPECTED_SEP_SMOKE_ROOT)/scaffold_expected_msep1p5_lsep0p5_smoke/p9016_full.sep_diag.tsv

smoke_blind_p9016_resolution_chain: run_blind_p9016_minimal.bin audit_blind_p9016_full_cpu_output.bin testdata/p9016_blind_smoke.pairs
	rm -rf $(P9016_RESOLUTION_CHAIN_SMOKE_ROOT)
	mkdir -p $(SMOKE_RES_ROOT)
	HK_BLIND_P9016_PAIRS=testdata/p9016_blind_smoke.pairs \
	HK_BLIND_P9016_ALLOW_CUSTOM_PAIRS=1 \
	HK_BLIND_P9016_OUTPUT_ROOT=$(P9016_RESOLUTION_CHAIN_SMOKE_ROOT) \
	HK_BLIND_P9016_CHAIN=1 \
	HK_BLIND_P9016_RESOLUTION_CHAIN=4000000,1000000,200000 \
	HK_BLIND_P9016_MINIMAL_N_ITER=1 \
	HK_BLIND_P9016_MINIMAL_RELAX_STEPS=1 \
	./run_blind_p9016_minimal.bin
	HK_BLIND_AUDIT_ALLOW_CUSTOM_RAW_PAIRS=1 ./audit_blind_p9016_full_cpu_output.bin $(P9016_RESOLUTION_CHAIN_SMOKE_ROOT)/4m/minimal_soft_sep_off
	HK_BLIND_AUDIT_ALLOW_CUSTOM_RAW_PAIRS=1 ./audit_blind_p9016_full_cpu_output.bin $(P9016_RESOLUTION_CHAIN_SMOKE_ROOT)/1m/minimal_soft_sep_off
	HK_BLIND_AUDIT_ALLOW_CUSTOM_RAW_PAIRS=1 ./audit_blind_p9016_full_cpu_output.bin $(P9016_RESOLUTION_CHAIN_SMOKE_ROOT)/200k/minimal_soft_sep_off

audit_blind_p9016_full_cpu_output: audit_blind_p9016_full_cpu_output.bin

audit_blind_p9016_full_cpu_output.bin: audit_blind_p9016_full_cpu_output.c hickit.h krng.h
	$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) audit_blind_p9016_full_cpu_output.c -o $@ $(AUDIT_LIBS)

clean:
	rm -f $(PROG) *.o *.bin .fdg_gpu_backend.*

depend:
	( LC_ALL=C ; export LC_ALL; makedepend -Y -- $(CFLAGS) $(CPPFLAGS) -- *.c ) 2>/dev/null
