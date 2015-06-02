#ifndef __CUDA_KERNEL__
#define __CUDA_KERNEL__

extern float call_kernel( size_t nloops, size_t nblocks_x, size_t nblocks_y, size_t nthreads );

#endif /*__CUDA_KERNEL__*/