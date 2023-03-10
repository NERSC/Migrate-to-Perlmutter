NVCCFLAGS = -arch=sm_80
NVCFLAGS = -cuda -gpu
# nvcc can be used to build cuda code (*.cu) for any compiler:
vec_add: vecAdd.cu
	nvcc $(NVCCFLAGS) vecAdd.cu -o $@
# the nvidia compilers recognize cuda code within a C or C++ source file,
# and can compile it directly (either via the Cray compiler wrapper, or with
# nvc or nvc++):
vec_add_nv: vecAdd.cpp
	CC $(NVCFLAGS) vecAdd.cpp -o vec_add
clean:
	rm -f vec_add vec_add_nv

