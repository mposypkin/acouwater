BACKEND = 1
ROOT = ..
OBJS = sspemdd_sequential.o 
LIBS = sspemdd.a
ADDIT_INC = -I$(ROOT)/alglib 
include $(ROOT)/all.inc
-include deps.inc

cpu: BACKEND=2
cpu: all
cpu: 
	ar r sspemdd.a ../gpu/ils_cpu.o

gpu: BACKEND=3
gpu: all
gpu: 
	ar r sspemdd.a ../gpu/ils_gpu.o
