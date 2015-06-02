#include <cuda.h>
#include <helper_cuda.h>

#include "gateproc.h"
#include "kernel.h"
#include "cuda_timer.h" 

typedef struct __lsns_align(16) __dev_chan_data{/*SoA data allocated in device memory*/
	float4 __lsns_align( 16 ) *devG;	// Gmax(x), g(y), I(z), GE(w)
	float4 __lsns_align( 16 ) *devMH;	// M(x), H(y), reserved(z), reserved(w)
	int4 __lsns_align( 16 ) *devTypeMH;	// type_m (x), par_m(y), type_h(z), par_h(w)
	int4 __lsns_align( 16 ) *devSharedVar;// indices of shared variables: V(x), E(y), reserved(z), reserved(w)
} dev_icdata;	/*data related to ions channel*/

#define NUM_GPARS 128

__constant__ __lsns_align(16) gate_par GatePar[NUM_GPARS];

/*
 * 192 threads
 * Threads Per Block			192
 * Registers Per Thread			20
 * Shared Memory Per Block (bytes)	6144 (8 float/int variables per thread)
 * =============================================
 * Occupancy of each Multiprocessor	100%
 */
// <<<blocks, threads>>>; blocks = 16X1, threads = 192X1
//	Coalesced access to device's memory for [1D grid X 1D block] metric:
//	shmem[threadIdx.x] = gmem[threadIdx.x + blockIdx.x*blockDim.x ];

//	Coalesced, 2D grid, 1D block:
//	shmem[threadIdx.x] = gmem[threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*blockDim.x*gridDim.x ]; 

//	Coalesced, 2D grid and block:
//	int x = threadIdx.x + blockIdx.x*blockDim.x;
//	int y = threadIdx.y + blockIdx.y*blockDim.y;
//	int elementPitch = blockDim.x*gridDim.x;
//	shmem[threadIdx.x+threadIdx.y*blockDim.x] = gmem[x + y*elementPitch];
__global__ void chan_kernel( float4 *dev_g, float4 *dev_mh, int4 *dev_typemh, int4 *dev_sharedvar )
{
	int index = threadIdx.x+blockIdx.x*blockDim.x+blockIdx.y*blockDim.x*gridDim.x;

	float4 mh = dev_mh[index];
	float4 g = dev_g[index];
	
	int4 typemh = dev_typemh[index];
	int4 sharedvar = dev_sharedvar[index];

	typemh.x = 4; typemh.y = 1; //activation
	typemh.z = 5; typemh.w = 2; //inactivation
	sharedvar.x = 0; sharedvar.y = 1; sharedvar.z = 2; sharedvar.w = 3;

//	float m = mh.x, h = mh.y;

	lsns_gate(typemh.x, GatePar[typemh.y], 0, 0.1, mh.x );
	lsns_gate(typemh.z, GatePar[typemh.w], 0, 0.1, mh.y );
/*
 	float4 __lsns_align( 16 ) *devG;	// Gmax(x), g(y), I(z), GE(w)
	float4 __lsns_align( 16 ) *devMH;	// M(x), H(y), reserved(z), reserved(w)
	int4 __lsns_align( 16 ) *devTypeMH;	// type_m (x), par_m(y), type_h(z), par_h(w)
	int4 __lsns_align( 16 ) *devSharedVar;// indices of shared variables: V(x), E(y), reserved(z), reserved(w)
 */
//	mh.x = m; mh.y = h;
	g.x = 1.; g.y = g.x*mh.x*mh.y; g.z = 0; g.w = 0;
	
	dev_g[index] = g;
	dev_mh[index] = mh;
}

bool allocate_devdata( dev_icdata *data, size_t size )
{
	// allocate data
	if( cudaSuccess != cudaMalloc(( void **)&( data->devG ), size*sizeof( float4 )))
		return false;
	if( cudaSuccess != cudaMalloc(( void **)&( data->devMH ), size*sizeof( float4 )))
		return false;
	if( cudaSuccess != cudaMalloc(( void **)&( data->devTypeMH ), size*sizeof( int4 )))
		return false;
	if( cudaSuccess != cudaMalloc(( void **)&( data->devSharedVar ), size*sizeof( int4 )))
		return false;
	return true;
}

void free_devdata( dev_icdata *data )
{
	if( data != NULL ){
		if( data->devG ){
			cudaFree( data->devG );
		}
		if( data->devMH ){
			cudaFree( data->devMH );
		}
		if( data->devTypeMH ){
			cudaFree( data->devTypeMH );
		}
		if( data->devSharedVar ){
			cudaFree( data->devSharedVar );
		}
	}
}

float call_kernel( size_t nloops, size_t nblocks_x, size_t nblocks_y, size_t nthreads )
{
	float time = -1.;
	int devID = gpuGetMaxGflopsDeviceId();
	cudaError_t error_id = ( cudaSetDevice( devID ));
	if( error_id == cudaSuccess ){
/*
todo:	allocate memory on host
	specify parameters to be processed      
	allocate memory on device
	copy data from host to device
*/
//		non-fragmented data related to channels/synapses
//		float *G;		// conductance
//		float *Gmax;	// maximal conductance
//		float *M;		// activation
//		float *H;		// inactivation/saturation
//		float *I;		// current
//		float *GE;		// the term which is necessary to implement algorithm of exponential Euler 
//		fragmented data related to channels/synapses (must be cached)
//		float *V;		// membrane potential
//		float *E;		// reversal potential
		size_t data_size = nblocks_x*nblocks_y*nthreads;
		dev_icdata devdata; devdata.devG = NULL; devdata.devMH = NULL; devdata.devTypeMH = NULL; devdata.devSharedVar = NULL;
		if( !allocate_devdata( &devdata, data_size ) ){
			printf("Memory allocation error\n" );
			exit( 0 );
		}
		printf("devG = %p\n", devdata.devG );
		time = 0.;
		cuda_timer timer;
		dim3 blocks( nblocks_x, nblocks_y );
		dim3 threads( nthreads, 1 );
		for( size_t i = 0; i < nloops; ++i ){
			timer.start();
			chan_kernel<<<blocks, threads>>>( devdata.devG, devdata.devMH, devdata.devTypeMH, devdata.devSharedVar );
			timer.stop();
			time += timer.elapsed();
		}
		free_devdata( &devdata );
	}
	else{
		printf("call_kernel returned %d\n-> %s\n", error_id, cudaGetErrorString( error_id ));
		exit(EXIT_FAILURE);
	}
	cudaDeviceReset();
	return time;
}

/*
 Q. Let's say that i have array of 32bit int in global memory and i want to copy it 
    to shared memory with coalesced access. 

 A. What you want ultimately depends on whether your input data is a 1D or 2D array, 
    and whether your grid and blocks are 1D or 2D. The simplest case is both 1D:

	shmem[threadIdx.x] = gmem[blockDim.x * blockIdx.x + threadIdx.x];

    This is coalesced. The rule of thumb I use is that the most rapidly varying coordinate 
    (the threadIdx) is added on as offset to the block offset (blockDim * blockIdx). The end 
    result is that the indexing stride between threads in the block is 1. If the stride gets 
    larger, then you lose coalescing.

    The simple rule (on Fermi and later GPUs) is that if the addresses for all threads in a 
    warp fall into the same aligned 128-byte range, then a single memory transaction will
    result (assuming caching is enabled for the load, which is the default). If they fall
    into two aligned 128-byte ranges, then two memory transactions result, etc.

    On GT2xx and earlier GPUs, it gets more complicated. But you can find the details of that
    in the programming guide.

    Additional examples:

    Not coalesced:

	shmem[threadIdx.x] = gmem[blockDim.x + blockIdx.x * threadIdx.x];

    Not coalesced, but not too bad on GT200 and later:

	stride = 2;
	shmem[threadIdx.x] = gmem[blockDim.x*blockIdx.x+stride*threadIdx.x];

    Not coalesced at all:

	stride = 32;
	shmem[threadIdx.x] = gmem[blockDim.x * blockIdx.x + stride * threadIdx.x];

    Coalesced, 2D grid, 1D block:

	int elementPitch = blockDim.x * gridDim.x;
	shmem[threadIdx.x] = gmem[blockIdx.y * elementPitch + blockIdx.x * blockDim.x + threadIdx.x]; 

    Coalesced, 2D grid and block:

	int x = blockIdx.x * blockDim.x + threadIdx.x;
	int y = blockIdx.y * blockDim.y + threadIdx.y;
	int elementPitch = blockDim.x * gridDim.x;
	shmem[threadIdx.y * blockDim.x + threadIdx.x] = gmem[y * elementPitch + x];
 */

/*
You seem to be a bit confused about the thread hierachy that CUDA has; in a nutshell, for a kernel there will be 1 grid, 
(which I always visualize as a 3-dimensional cube). Each of its elements is a block, such that a grid declared as dim3 
grid(10, 10, 2); would have 10*10*2 total blocks. In turn, each block is a 3-dimensional cube of threads.
With that said, it's common to only use the x-dimension of the blocks and grids, which is what it looks like the code 
in your question is doing. This is especially revlevant if you're working with 1D arrays. In that case, your 
tid+=blockDim.x * gridDim.x line would in effect be the unique index of each thread within your grid. This is because 
your blockDim.x would be the size of each block, and your gridDim.x would be the total number of blocks.

So if you launch a kernel with parameters

dim3 block_dim(128,1,1);
dim3 grid_dim(10,1,1);
kernel<<<grid_dim,block_dim>>>(...);

then in your kernel had threadIdx.x + blockIdx.x*blockDim.x you would effectively have:

threadIdx.x range from [0 ~ 128)
blockIdx.x range from [0 ~ 10)
blockDim.x equal to 128
gridDim.x equal to 10

Hence in calculating threadIdx.x + blockIdx.x*blockDim.x, you would have values within the range defined by: 
[0, 128) + 128 * [1, 10), which would mean your tid values would range from {0, 1, 2, ..., 1279}. This is useful 
for when you want to map threads to tasks, as this provides a unique identifier for all of your threads in your kernel.

However, if you have

int tid = threadIdx.x + blockIdx.x * blockDim.x;
tid += blockDim.x * gridDim.x;

then you'll essentially have: tid = [0, 128) + 128 * [1, 10) + (128 * 10), and your tid values would range from 
{1280, 1281, ..., 2559} I'm not sure where that would be relevant, but it all depends on your application and 
how you map your threads to your data. This mapping is pretty central to any kernel launch, and you're the one 
who determines how it should be done. When you launch your kernel you specify the grid and block dimensions, and 
you're the one who has to enforce the mapping to your data inside your kernel. As long as you don't exceed your 
hardware limits (for modern cards, you can have a maximum of 2^10 threads per block and 2^16 - 1 blocks per thread)
*/
