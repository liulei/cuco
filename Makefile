# If we want to use cuda gramma, compiler should be switched to nvcc

CC			=	gcc
NVCC		=	nvcc

#OPTIMIZE	=	-O2 -W -Wall
OPTIMIZE	=	-O2

GSL_INCL	=	-I/usr/local/include
GSL_LIBS	=	-L/usr/local/lib

CUFFT_INCL	=	-I/usr/local/cuda/include
CUFFT_LIBS	=	-L/usr/local/cuda/lib

CUTIL_INCL	=	-I/home/liulei/NVIDIA_CUDA_SDK/common/inc
CUTIL_LIBS	=	-L/home/liulei/NVIDIA_CUDA_SDK/lib

CUDPP_INCL	=	-I/home/liulei/NVIDIA_CUDA_SDK/common/inc
CUDPP_LIBS	=	-L/home/liulei/NVIDIA_CUDA_SDK/common/lib/linux

OPTIONS		=	$(OPTIMIZE)

EXEC		=	cuco

C_SOURCES	=	allvars.c begrun.c init.c read_ic.c \
				main.c driftfac.c longrange.c pm_periodic.c \
				run.c predict.c accel.c io.c timestep.c

C_OBJS		=	$(patsubst %.c, %.c.o, $(C_SOURCES))

CPP_SOURCES	=	gravlist.cpp

CPP_OBJS	=	$(patsubst %.cpp, %.cpp.o, $(CPP_SOURCES))

CU_SOURCES	=	gravlist.cu radixsort.cu

CU_OBJS		=	$(patsubst %.cu, %.cu.o, $(CU_SOURCES))

OBJS		=	$(C_OBJS) $(CPP_OBJS) $(CU_OBJS)

INCL		=	allvars.h proto.h proto.cuh radixsort.h Makefile

CFLAGS		=	$(OPTIONS) $(GSL_INCL) $(CUFFT_INCL) \
				$(CUTIL_INCL) $(CUDPP_INCL)

LIBS		=	$(GSL_LIBS) -lgsl -lgslcblas -lm \
				$(CUFFT_LIBS) -lcudart -lcufft \
				$(CUTIL_LIBS) -lcutil \
				$(CUDPP_LIBS) -lcudpp

$(EXEC)	:	$(OBJS) $(INCL)
	$(CC) $(OBJS) $(LIBS) -o $(EXEC)

%.c.o	:	%.c
	$(CC) -c $^ $(CFLAGS) $(LIBS) -o $@

%.cpp.o	:	%.cpp
	$(CC) -c $^ $(CFLAGS) $(LIBS) -o $@

%.cu.o	:	%.cu
	$(NVCC) -c $^ $(CFLAGS) $(LIBS) -o $@

clean:
	rm -f $(C_OBJS) $(CU_OBJS) $(CPP_OBJS) $(EXEC) *.linkinfo
