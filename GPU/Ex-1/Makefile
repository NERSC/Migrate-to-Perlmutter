NVCCFLAGS = -arch=sm_80
NVCFLAGS = -cuda -gpu=cc80
vec_add: vecAdd.cu
	nvcc $(NVCCFLAGS) vecAdd.cu -o $@
vec_add_nv: vecAdd.cu
	nvc++ $(NVCFLAGS) vecAdd.cpp -o vec_add
clean:
	rm vec_add vec_add_nv

