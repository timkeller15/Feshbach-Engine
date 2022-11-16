CUDA_HOME = /opt/cuda/
GPU_ARCH	= sm_52
OS:=	$(shell uname)
ifeq ($(OS),Darwin)
CUDA_LIB	= $(CUDA_HOME)/lib
CUDA_HEADER	= $(CUDA_HOME)/include
CC		= $(CUDA_HOME)/bin/nvcc -ccbin /usr/bin/clang --ptxas-options=-v#-save-temps
CFLAGS		= -g -std=c++11 -Wno-deprecated-gpu-targets
else
CUDA_LIB	= $(CUDA_HOME)/lib64
CUDA_HEADER	= $(CUDA_HOME)/include
CC		= $(CUDA_HOME)/bin/nvcc --ptxas-options=-v --compiler-options -Wall #-save-temps
CHOSTFLAGS	= #-fopenmp
CFLAGS		= -g -O3 -std=c++11 -Xcompiler '-std=c++11' -Xcompiler '-fopenmp' #-L$(CUTT_DIR) -l:libcutt.a
endif

CUDA_FLAGS 	= -lcufft

CLINKER		= $(CC)
RM		= /bin/rm
INCFLAGS	= -I$(CUDA_HEADER)
LDFLAGS		= -L$(CUDA_LIB)
EXECS		= gpue_ramp # BINARY NAME HERE

DEPS = ./include/ds.h ./include/evolution.h ./include/fileIO.h ./include/init.h ./include/kernels.h ./include/operators.h ./include/parser.h ./include/split_op.h

OBJ = fileIO.o kernels.o split_op.o ds.o parser.o evolution.o init.o operators.o

%.o: ./src/%.cc $(DEPS)
	$(CC) -c -o $@ $(INCFLAGS) $(CFLAGS) $(LDFLAGS) -Xcompiler "-fopenmp" -arch=$(GPU_ARCH) $<

%.o: ./src/%.cu $(DEPS)
	$(CC) -c -o $@ $(INCFLAGS) $(CFLAGS) $(LDFLAGS) $(CUDA_FLAGS) -Xcompiler "-fopenmp" -arch=$(GPU_ARCH) $< -dc

gpue_ramp: $(OBJ)
	$(CC) -o $@ $(INCFLAGS) $(CFLAGS) $(LDFLAGS) $(CUDA_FLAGS) -Xcompiler "-fopenmp" -arch=$(GPU_ARCH) $^

clean:
	@-$(RM) -f r_0 Phi_0 E* px_* py_0* xPy* xpy* ypx* x_* y_* yPx* p0* p1* p2* EKp* EVr* gpot wfc* Tpot 0* V_* K_* Vi_* Ki_* 0i* k s_* si_* *.o *~ PI* $(EXECS) $(OTHER_EXECS) *.dat *.eps *.ii *.i *cudafe* *fatbin* *hash* *module* *ptx test* vort* v_opt*;
