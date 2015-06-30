#include <cuda.h>
#include <helper_cuda.h>

#include "gateproc.h"
#include "pumpproc.h"
#include "kernel.h"
#include "cuda_timer.h" 

// dev_mh: M(x), H(y), Power M(z), Power H(w)
#define _gate_m( mh ) ( mh ).x
#define _gate_h( mh ) ( mh ).y
#define _gate_powm( mh ) ( mh ).z
#define _gate_powh( mh ) ( mh ).w
// dev_g: Gmax(x), g(y), I(z), GE(w)
#define _gate_gmax( g ) ( g ).x
#define _gate_g( g ) ( g ).y
#define _gate_i( g ) ( g ).z
#define _gate_ge( g ) ( g ).w
// dev_type: type_m (x), par_m(y), type_h(z), par_h(w)
#define _gate_typem( t ) ( t ).x
#define _gate_parm( t ) ( t ).y
#define _gate_typeh( t ) ( t ).z
#define _gate_parh( t ) ( t ).w
// dev_ve (indices of shared variables): V(x), E(y), ions inside the cell for activation(z), ions inside the cell for inactivation(w)
#define _gate_v( sh ) ( sh ).x
#define _gate_e( sh ) ( sh ).y
#define _gate_in_m( sh ) ( sh ).z
#define _gate_in_h( sh ) ( sh ).w
// dev_e: E(x), In(y), Out(z), RTFz(w)
#define _ions_eds( e ) ( e ).x
#define _ions_in( e ) ( e ).y
#define _ions_out( e ) ( e ).z
#define _ions_rtfz( e ) ( e ).w
// dev_i: Ipump(x), Ichan(y), reserved(z), reserved(w)
#define _ions_ipump( i ) ( i ).x
#define _ions_ichan( i ) ( i ).y
// dev_v: V(x), reserver(y), reserved(z), reserved(w)
#define _cell_v( v ) ( v ).x

typedef struct __lsns_align(16) __ichan_data{/*SoA data allocated in device memory*/
	float4 __lsns_align( 16 ) *G;		// Gmax(x), g(y), I(z), GE(w)
	float4 __lsns_align( 16 ) *MH;	// M(x), H(y), Power M(z), Power H(w)
	int4 __lsns_align( 16 ) *TypeMH;	// type_m (x), par_m(y), type_h(z), par_h(w)
	int4 __lsns_align( 16 ) *SharedVar;	// indices of shared variables: V(x), E(y), reserved(z), reserved(w)
	int4 __lsns_align( 16 ) *iG;		// size(x), index1(y), index2(z), index3(w), index4 (+1.x) ... index11(+2.w)
	float4 __lsns_align( 16 ) *V;		// V(x), reserved(y), reserved(z), reserved(w)
	float4 __lsns_align( 16 ) *E;		// E(x), In(y), Out(z), RTFz(w)
	float4 __lsns_align( 16 ) *I;		// E(x), In(y), Out(z), RTFz(w)
} icdata;	/*data related to ions channel*/

#define NUM_GPARS 128

__constant__ __lsns_align(16) gate_par Gates[NUM_GPARS];
__constant__ __lsns_align(16) float SimStep = 0.1;

__global__ void ions_kernel(  /*int4 *dev_type*/ float4 *dev_e, float4 *dev_i, float4 *dev_g, int4 *dev_ig )
{
	int index = _cu_index();				// 
	float step = SimStep;					// integration step
	float4 e = dev_e[index];				//
	float4 i = dev_i[index];				//
	__syncthreads();
	
	// calculate ions dynamics, ions current and pump current
	float ipump, ichan;					// pump current and ions current
	float eds = _ions_eds( e );				// resting potential 
	float in = _ions_in( e );				// concentation of ions inside the cell
	float out = _ions_out( e );				// concentation of ions outside cell
	float v = 0; 						//-->> get membrane potential
	float rtfz = 1;						//-->> get rtfz from parameters
	float T = 1;						//-->> get time constant for ions dynamics
	float gchan = 0; 						// sum of conductances of all ion currents which are involved to ions dynamics
	lsns_fastsum( _gate_g, dev_g, dev_ig+index*3, gchan );

	// 2-order Runge-Kutta integration method
	// sub-step 1
//	lsns_ipump( type, par, in, ipump )
//	lsns_ichan( gchan, v, in, out, rtfz, eds, ichan )
	// sub-step 2
	float in1 = in-step*( ichan+ipump )/( 2*T );
//	lsns_ipump( type, par, in1, ipump )
//	lsns_ichan( gchan, v, in1, out, rtfz, eds, ichan )
	// the final calculations for in, eds, ichan
	in = in-step*( ichan+ipump )/T;
//	lsns_ipump( type, par, in, ipump )
//	lsns_ichan( gchan, v, in, out, rtfz, eds, ichan )

	// store the results
	_ions_ipump( i ) = ipump;
	_ions_ichan( i ) = ichan;
	_ions_in( e ) = in;
	_ions_eds( e ) = eds;
	dev_e[index] = e;
	dev_i[index] = i;
}

// calculate ions (aka channels) and synaptic currents. 
__global__ void ichan_kernel( float4 *dev_g, float4 *dev_mh, int4 *dev_type, int4 *dev_ve, float4 *dev_v, float4 *dev_e )
{
	int index = _cu_index();
	float step = SimStep;					// integration step
	// load data
	int4 tp = dev_type[index];				// type of ion channel (generic, a-b, etc) and its parameters (half-voltage, slope, etc)
	int4 ve = dev_ve[index];				// external parameters (membrane potential, rest potential, etc)
	float4 v = dev_v[_gate_v( ve )];			// get the membrane potential
	float4 e = dev_e[_gate_e( ve )];			// get the resting potential
	float4 g = dev_g[index];				// properties of ions channel (conductance, current, etc)
//>>>> remove later, for testing purpose only
_gate_gmax( g ) = 1.;
_gate_typem( tp ) = 0;
_gate_parm( tp ) = 1;
_gate_typeh( tp ) = 0;
_gate_parh( tp ) = 2;
_gate_v( ve ) = threadIdx.x; 
_gate_e( ve ) = threadIdx.x;
//<<<< remove later, for testing purpose only
	// load properties of gate variables (activation, inactivation, etc) load if needed
	float4 mh = ( _gate_typem( tp )+_gate_typeh( tp ) != LSNS_NOGATE)? dev_mh[index]: float4();
	// load Mg- or Ca- concentration inside the cell for NMDA synapse or Z-channels if needed
	float in_m = ( _gate_typem( tp ) >= LSNS_ZGENERIC_INSTANT )? _ions_in( dev_e[_gate_in_m( ve )] ):0;
	// load Mg- or Ca- concentration inside the cell for NMDA synapse or Z-channels if needed
	float in_h = ( _gate_typeh( tp ) >= LSNS_ZGENERIC_INSTANT )? _ions_in( dev_e[_gate_in_h( ve )] ):0;
//>>>> remove later, for testing purpose only
_gate_m( mh ) = 0.3;
_gate_h( mh ) = 0.1;
_gate_powm( mh ) = 3;
_gate_powh( mh ) = 2; 
//<<<< remove later, for testing purpose only

	// calculating
	float eds = _ions_eds( e );				// get resting potential
	float vm = _cell_v( v );				// get membrane potential
//>>>> remove later, for testing purpose only
eds = 0; vm = 0;
//<<<< remove later, for testing purpose only
	float mp, hp; 
	proc_gate( _gate_typem( tp ), Gates[_gate_parm( tp )], in_m, vm, step, _gate_powm( mh ), _gate_m( mh ), mp );
	proc_gate( _gate_typeh( tp ), Gates[_gate_parh( tp )], in_h, vm, step, _gate_powh( mh ), _gate_h( mh ), hp );
	_gate_g( g ) = _gate_gmax( g )*mp*hp;		// g
	_gate_i( g ) = _gate_g( g )*( vm-eds );		// I
	_gate_ge( g ) = _gate_g( g )*eds;			// ge
	// store results
	dev_g[index] = g;
	if( _gate_typem( tp )+_gate_typeh( tp ) != LSNS_NOGATE){
		dev_mh[index] = mh;
	}
}

// dev_v: V(x), G(y), GE(z), I(w)
__global__ void neurons_kernel( float4 *dev_v )
{
	int index = _cu_index(); //threadIdx.x+blockIdx.x*blockDim.x+blockIdx.y*blockDim.x*gridDim.x;
	float4 v = dev_v[index];
	__syncthreads();
	// g = sum( all chan g );
	// ge = sum( all chan ge );
	// i = sum( all chan i );
	v.x = -60; 
	
	dev_v[index] = v;
}

bool allocate_devdata( size_t size, icdata &data )
{
	data.G = NULL; data.MH = NULL; data.TypeMH = NULL; data.SharedVar = NULL; data.iG = NULL; data.V = NULL; data.E = NULL;
	// allocate data
	if( cudaSuccess != cudaMalloc(( void **)&( data.G ), size*4*sizeof( float4 )))
		return false;
	if( cudaSuccess != cudaMalloc(( void **)&( data.MH ), size*4*sizeof( float4 )))
		return false;
	if( cudaSuccess != cudaMalloc(( void **)&( data.TypeMH ), size*4*sizeof( int4 )))
		return false;
	if( cudaSuccess != cudaMalloc(( void **)&( data.SharedVar ), size*4*sizeof( int4 )))
		return false;
	if( cudaSuccess != cudaMalloc(( void **)&( data.iG ), size*4*sizeof( int4 )))
		return false;
	if( cudaSuccess != cudaMalloc(( void **)&( data.V ), size*sizeof( float4 )))
		return false;
	if( cudaSuccess != cudaMalloc(( void **)&( data.E ), size*sizeof( float4 )))
		return false;
	if( cudaSuccess != cudaMalloc(( void **)&( data.I ), size*sizeof( float4 )))
		return false;
	
	cudaMemset( data.G, 0, size*4*sizeof( float4 ));
	cudaMemset( data.MH, 0, size*4*sizeof( float4 ));
	cudaMemset( data.TypeMH, 0, size*4*sizeof( int4 ));
	cudaMemset( data.SharedVar, 0, size*4*sizeof( int4 ));
	cudaMemset( data.iG, 0, size*4*sizeof( int4 ));
	cudaMemset( data.V, 0, size*sizeof( float4 ));
	cudaMemset( data.E, 0, size*sizeof( float4 ));
	cudaMemset( data.I, 0, size*sizeof( float4 ));
	return true;
}

void free_devdata( icdata &data )
{
	if( data.G ){
		cudaFree( data.G );
	}
	if( data.MH ){
		cudaFree( data.MH );
	}
	if( data.TypeMH ){
		cudaFree( data.TypeMH );
	}
	if( data.SharedVar ){
		cudaFree( data.SharedVar );
	}
	if( data.iG ){
		cudaFree( data.iG );
	}
	if( data.V ){
		cudaFree( data.V );
	}
	if( data.E ){
		cudaFree( data.E );
	}
	if( data.I ){
		cudaFree( data.E );
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
		icdata devdata;
		printf( "Data size %ld\n ", 60*data_size*4 );
		if( !allocate_devdata( data_size, devdata )){
			printf("Memory allocation error\n" );
			exit( 0 );
		}
		time = 0.;
		cuda_timer timer;
		dim3	blocks( nblocks_x, nblocks_y*4 ),
			blocks1( nblocks_x, nblocks_y ),
			threads( nthreads, 1 );
		for( size_t i = 0; i < nloops; ++i ){
			timer.start();
/*			
			ions_kernel<<<blocks1, threads>>>( devdata.E, devdata.I, devdata.G, devdata.iG );
			cudaDeviceSynchronize();
*/			
			ichan_kernel<<<blocks, threads>>>( devdata.G, devdata.MH, devdata.TypeMH, devdata.SharedVar, devdata.V, devdata.E );
/*			
			cudaDeviceSynchronize();
			neurons_kernel<<<blocks1, threads>>>( devdata.V );
*/			
			timer.stop();
			time += timer.elapsed();
		}
		free_devdata( devdata );
	}
	else{
		printf("call_kernel returned %d\n-> %s\n", error_id, cudaGetErrorString( error_id ));
		exit(EXIT_FAILURE);
	}
	cudaDeviceReset();
	return time;
}

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
