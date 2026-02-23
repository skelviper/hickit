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
OBJS=		sdict.o io.o pair.o count.o phase.o bin.o fdg.o image.o view3d.o fdg_gpu_stub.o
PROG=		hickit
LIBS=		-lm -lz
LIBS_GL=
ASAN_FLAG=
CUDA_OBJS=
CXX?=		g++
LINK?=		$(CC)
NVCC?=		$(CUDA_HOME)/bin/nvcc
CUDAFLAGS?=	-O3 -std=c++17

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

.PHONY:all clean depend
.SUFFIXES:.c .o .cu

.c.o:
		$(CC) -c $(CFLAGS) $(ASAN_FLAG) $(CPPFLAGS) $(INCLUDES) $< -o $@

.cu.o:
		$(NVCC) $(CUDAFLAGS) $(CPPFLAGS) -I. $(INCLUDES) -c $< -o $@

all:$(PROG)

hickit:$(OBJS) $(CUDA_OBJS) main.o
		$(LINK) $(LDFLAGS) -o $@ $^ $(ASAN_FLAG) $(LIBS_GL) $(LIBS)

clean:
		rm -fr gmon.out *.o a.out $(PROG) *.a *.dSYM hickit.aux hickit.log hickit.pdf

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(CPPFLAGS) -- *.c)

# DO NOT DELETE

bin.o: hkpriv.h hickit.h krng.h khash.h ksort.h
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
