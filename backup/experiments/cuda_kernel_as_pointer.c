

    #include <cuda.h>
    #include <stdio.h>


    __global__ void foo(int* a) {
       *a=1;
    }


    void Call(void(*func)(int*), int Nob, int BlkSz, int* a, cudaStream_t stream) {
        (*func)<<<Nob, BlkSz, 0, stream>>>(a);
    }



    int main() {
        int *a, b;
        cudaMalloc(&a, sizeof(int));

	Call(foo, 1, 1, a, 0);

	cudaMemcpy(&b, a, sizeof(int), cudaMemcpyDeviceToHost);
        printf("result is %d\n", b);
        cudaFree(a);
        return 0;
    }

    $ nvcc -o foo foo.cu
    $ ./foo 
    result is 1