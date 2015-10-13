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
	#define __lsns_cached( p ) __ldg( p )
#else
	#define __lsns_inline inline
	#define __lsns_cached( p ) p
	#if defined( __GNUC__ ) // GCC
		#include <stdint.h>
		#define __lsns_align( n ) __attribute__(( aligned( n )))
		typedef int32_t __lsns_int32;
		typedef uint32_t __lsns_uint32;
	#elif defined(_MSC_VER ) // MSVC
		#define __lsns_align( n ) __declspec( align( n ))
		typedef __int32 __lsns_int32;
		typedef unsigned __int32 __lsns_uint32; 
	#else
		#define __lsns_align( n )
		typedef int __lsns_int32;
		typedef unsigned int __lsns_uint32; 
		#error "Please provide a definition for __lsns_align macro for your host compiler!"
		#error "Please provide a definition for __lsns_int32 macro for your host compiler!"
		#error "Please provide a definition for __lsns_uint32 macro for your host compiler!"
	#endif
	typedef struct __int2{
		__lsns_int32 x;
		__lsns_int32 y;
	} int2;
	typedef struct __uint2{
		__lsns_uint32 x;
		__lsns_uint32 y;
	} uint2;
	typedef struct __int4{
		__lsns_int32 x;
		__lsns_int32 y;
		__lsns_int32 z;
		__lsns_int32 w;
	} int4;
	typedef struct __uint4{
		__lsns_uint32 x;
		__lsns_uint32 y;
		__lsns_uint32 z;
		__lsns_uint32 w;
	} uint4;
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

