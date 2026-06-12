CUDA_HOME ?= /usr/local/cuda
CUDA_LIB  ?= $(CUDA_HOME)/targets/x86_64-linux/lib

CFLAGS   ?= -g -Wall -O2 -Wc++-compat -ffast-math
CPPFLAGS ?=
INCLUDES ?=
LDFLAGS  += -L$(CUDA_LIB) -Wl,-rpath,$(CUDA_LIB)
LIBS     := -lm -lz
LIBS_GL  :=
ASAN_FLAG :=

OBJS := sdict.o io.o pair.o count.o phase.o bin.o blind.o fdg.o image.o view3d.o fdg_gpu_stub.o
PROG := hickit
TEST_COMMON_SRCS := sdict.c io.c pair.c count.c phase.c bin.c blind.c fdg.c image.c view3d.c fdg_gpu_stub.c
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

.PHONY: all clean depend test smoke_blind_p9016_minimal audit_blind_p9016_full_cpu_output
.SUFFIXES: .c .o

.c.o:
	$(CC) -c $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) $< -o $@

all: $(PROG)

hickit: $(OBJS) main.o
	$(CC) $(LDFLAGS) -o $@ $^ $(ASAN_FLAG) $(LIBS_GL) $(LIBS)

test: test_diploid_index test_blind_input test_blind_binning test_fdg_energy test_blind_posterior test_blind_smoke test_blind_relax_loop test_blind_output test_blind_output_auditor_negative
	./test_diploid_index
	if ./test_diploid_index --invalid-copy >/dev/null 2>&1; then \
		echo "expected invalid copy assertion to fail"; \
		exit 1; \
	fi
	./test_blind_input
	./test_blind_binning
	./test_fdg_energy
	./test_blind_posterior
	./test_blind_smoke
	./test_blind_relax_loop
	./test_blind_output
	./test_blind_output_auditor_negative

test_blind_input: test_blind_input.c testdata/p9016_blind_fixture.pairs $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
	$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_input.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_binning: test_blind_binning.c testdata/p9016_blind_fixture.pairs $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
	$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_binning.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_diploid_index: test_diploid_index.c hickit.h krng.h
	$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) $< -o $@

test_fdg_energy: test_fdg_energy.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
	$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_fdg_energy.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_posterior: test_blind_posterior.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
	$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_posterior.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_smoke: test_blind_smoke.c testdata/p9016_blind_smoke.pairs $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
	$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_smoke.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_relax_loop: test_blind_relax_loop.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
	$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_relax_loop.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

test_blind_output: test_blind_output.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
	$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_output.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

run_blind_p9016_minimal.bin: run_blind_p9016_minimal.c $(TEST_COMMON_SRCS) hickit.h hkpriv.h krng.h
	$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) run_blind_p9016_minimal.c $(TEST_COMMON_SRCS) -o $@ $(LIBS)

smoke_blind_p9016_minimal: run_blind_p9016_minimal.bin audit_blind_p9016_full_cpu_output.bin testdata/p9016_blind_smoke.pairs
	rm -rf /tmp/hk_blind_p9016_minimal_audit_smoke
	HK_BLIND_P9016_PAIRS=testdata/p9016_blind_smoke.pairs \
	HK_BLIND_P9016_OUTPUT_ROOT=/tmp/hk_blind_p9016_minimal_audit_smoke \
	HK_BLIND_P9016_BIN_SIZE_BP=1000000 \
	HK_BLIND_P9016_MINIMAL_N_ITER=1 \
	HK_BLIND_P9016_MINIMAL_RELAX_STEPS=1 \
	HK_BLIND_P9016_MSTEP_GRAPH_MODE=raw_expected_soft_all \
	HK_BLIND_P9016_HARD_PC_SIZES=0 \
	./run_blind_p9016_minimal.bin
	./audit_blind_p9016_full_cpu_output.bin /tmp/hk_blind_p9016_minimal_audit_smoke/minimal_soft_sep_off

audit_blind_p9016_full_cpu_output: audit_blind_p9016_full_cpu_output.bin

audit_blind_p9016_full_cpu_output.bin: audit_blind_p9016_full_cpu_output.c hickit.h krng.h
	$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) audit_blind_p9016_full_cpu_output.c -o $@ $(LIBS)

test_blind_output_auditor_negative: test_blind_output_auditor_negative.c audit_blind_p9016_full_cpu_output.bin hickit.h krng.h
	$(CC) $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) test_blind_output_auditor_negative.c -o $@ $(LIBS)

clean:
	rm -f $(PROG) *.o *.bin test_diploid_index test_fdg_energy test_blind_posterior \
		test_blind_input test_blind_binning test_blind_smoke \
		test_blind_relax_loop test_blind_output \
		test_blind_output_auditor_negative

depend:
	( LC_ALL=C ; export LC_ALL; makedepend -Y -- $(CFLAGS) $(CPPFLAGS) -- *.c ) 2>/dev/null
