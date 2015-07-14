#ifndef __LSNS_CONFIG_H
#define __LSNS_CONFIG_H

// #define __CUDA__
// #define __LINUX__
// #define __GUI__
#define __CONSOLE__
#define __WINDOWS__
#define __LSNS_DEBUG__

#if defined( __LINUX__ )
	#undef __WINDOWS__
#elif defined( __WINDOWS__ )
	#undef __LINUX__
#else
	#define __LINUX__
	#undef __WINDOWS__
#endif

#if defined( __GUI__ )
	#undef __CONSOLE__
#elif defined( __CONSOLE__ )
	#undef __GUI__
#else
	#define __CONSOLE__
	#undef __GUI__
#endif

#if defined( __CUDACC__ ) // NVCC
	#define __lsns_inline __forceinline__ __device__
	#define __lsns_align( n ) __align__( n )
#else
	#define __lsns_inline inline
	#if defined( __GNUC__ ) // GCC
		#define __lsns_align( n ) __attribute__(( aligned( n )))
	#elif defined(_MSC_VER ) // MSVC
		#define __lsns_align( n ) __declspec( align( n ))
	#else
		#error "Please provide a definition for __lsns_align macro for your host compiler!"
	#endif

	typedef struct __int2{
		int x;
		int y;
	} int2;
	typedef struct __int4{
		int x;
		int y;
		int z;
		int w;
	} int4;
	typedef struct __float2{
		float x;
		float y;
	} float2;
	typedef struct __float4{
		float x;
		float y;
		float z;
		float w;
	} float4;
#endif

#endif /*__LSNS_CONFIG_H*/

