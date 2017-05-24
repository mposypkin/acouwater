CPP = g++
NVCC = nvcc
INC1 = ../alglib/
INCDIRS = -I${INC1}
CPPFLAGS =  -O3 -Wall -std=c++11 ${INCDIRS} -g -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS 
NVCCFLAGS = -gencode arch=compute_61,code=sm_61 -gencode arch=compute_50,code=sm_50 -gencode arch=compute_20,code=sm_21 -O3 -D_FORCE_INLINES -g -std=c++11 --ptxas-options="-v -dlcm=cg" -I/usr/local/cuda/samples/common/inc ${INCDIRS}

ils_gpu.o: ils_gpu.cu
	${NVCC} ${NVCCFLAGS} ils_gpu.cu -c

pd-sat-gpu: main.o alglibinternal.o alglibmisc.o ap.o linalg.o specialfunctions.o sspemdd_sequential.o sspemdd_parallel.o ils_gpu.o
	${CPP} ${CPPFLAGS} main.o alglibinternal.o alglibmisc.o ap.o linalg.o specialfunctions.o sspemdd_sequential.o sspemdd_parallel.o ils_gpu.o -o SSPEMDD_gpu -pthread -L/usr/local/cuda/lib64  -lcudart
	
clean:
	rm -rf *.o
	rm SSPEMDD_parallel
	rm SSPEMDD
	rm SSPEMDD_gpu
	clear
