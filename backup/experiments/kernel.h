#ifndef __CUDA_KERNEL__
#define __CUDA_KERNEL__

#include "config.h"

#if defined( __CUDA__ )
	#define _cu_index() threadIdx.x+blockIdx.x*blockDim.x+blockIdx.y*blockDim.x*gridDim.x
#endif /*__CUDA__*/

extern float call_kernel( size_t nloops, size_t nblocks_x, size_t nblocks_y, size_t nthreads );

#endif /*__CUDA_KERNEL__*/