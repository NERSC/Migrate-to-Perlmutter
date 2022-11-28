NVCCFLAGS = -arch=sm_80
NVCFLAGS = -cuda -gpu=cc80

vec_add: vecAdd.cu
	CC $(NVCFLAGS) vecAdd.cu -o $@
clean:
	rm vec_add
