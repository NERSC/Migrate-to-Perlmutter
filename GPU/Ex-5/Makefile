NVCCFLAGS = -arch=sm_80
NVCFLAGS = -gpu=cc80

vec_add: kernels.o kernels.h vecAdd.cpp
	CC -o $@ -lhwloc vecAdd.cpp kernels.o
kernels.o: kernels.cu kernels.h
	nvcc $(NVCCFLAGS) -c kernels.cu -o $@
clean:
	rm -f *.o vec_add
