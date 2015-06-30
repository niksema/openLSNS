#ifndef __CUDA_ENGINE__
#define __CUDA_ENGINE__

#include "config.h"

#if defined( __CUDA__ )
	#define _cu_index() threadIdx.x+blockIdx.x*blockDim.x+blockIdx.y*blockDim.x*gridDim.x
#endif /*__CUDA__*/

#endif /*__CUDA_ENGINE__*/