NVCCFLAGS = -arch=sm_80
NVCFLAGS = -cuda -gpu=cc80

vec_add: vecAdd.cu
	CC $(NVCFLAGS) vecAdd.cu -o $@
vec_add_gnu_cu: vecAdd.cu
	CC vecAdd.cu -o $@
vec_add_gnu_cpp: vecAdd.cpp
	CC vecAdd.cpp -o $@
clean:
	rm -f vec_add
