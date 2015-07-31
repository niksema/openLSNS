#ifndef __LSNS_MATH_H
#define __LSNS_MATH_H

#include "config.h"
#include "debug.h"

///////////////////////////////////////////////////////////////////////////////
// some useful math macroses
#if defined( __CUDA__ )
	#include <cuda.h>
	
	#define lsns_div( x, y ) __fdividef(( x ),( y ))
	#define lsns_pow( x, y ) __powf(( x ),( y ))
	#define lsns_exp( x ) __expf( x )
	#define lsns_cosh( x ) __coshf( x )
	#define lsns_log( x ) __logf( x )
#else
	#include <math.h>
	
	#define lsns_div( x, y ) ( x )/( y )
	#define lsns_pow( x, y ) pow(( x ),( y ))
	#define lsns_exp( x ) exp( x )
	#define lsns_cosh( x ) cosh( x )
	#define lsns_log( x ) log( x )
#endif
///////////////////////////////////////////////////////////////////////////////
// +'lsns_fastsum' unrolls the calculation of total sum of series that consists 
// of up-to-24 elements and stores the results into 'res'. The series is located 
// in 'data' array. 'fn' is the method to extract the actual number from 'data[i]'. 
// The 'lut' keeps the actual size (<24) of the series and reference to actual 
// position (offset) of the particular element in the array 'data'.
#define lsns_fastsum( fn, data, lut, res, N )\
{\
	int4 lut1 = ( lut )[0], lut2, lut3, lut4, lut5, lut6; \
	int n = lut1.x; \
	switch( n ){ \
		case 23: case 22: case 21: case 20: \
			lut6 = ( lut )[5]; \
		case 19: case 18: case 17: case 16: \
			lut5 = ( lut )[4]; \
		case 15: case 14: case 13: case 12: \
			lut4 = ( lut )[3]; \
		case 11: case 10: case 9: case 8: \
			lut3 = ( lut )[2]; \
		case 7: case 6: case 5: case 4: \
			lut2 = ( lut )[1]; \
	} \
	( res  ) = 0; \
	int i = 0; \
	switch( n ){ \
		case 23: \
			i = lut6.w; \
			__lsns_assert( i >= 0 && i < N );\
			res += fn(( data )[i]); \
		case 22: \
			i = lut6.z; \
			__lsns_assert( i >= 0 && i < N );\
			res += fn(( data )[i]); \
		case 21: \
			i = lut6.y; \
			__lsns_assert( i >= 0 && i < N );\
			res += fn(( data )[i]); \
		case 20: \
			i = lut6.x; \
			__lsns_assert( i >= 0 && i < N );\
			res += fn(( data )[i]); \
		case 19: \
			i = lut5.w; \
			__lsns_assert( i >= 0 && i < N );\
			res += fn(( data )[i]); \
		case 18: \
			i = lut5.z; \
			__lsns_assert( i >= 0 && i < N );\
			res += fn(( data )[i]); \
		case 17: \
			i = lut5.y; \
			__lsns_assert( i >= 0 && i < N );\
			res += fn(( data )[i]); \
		case 16: \
			i = lut5.x; \
			__lsns_assert( i >= 0 && i < N );\
			res += fn(( data )[i]); \
		case 15: \
			i = lut4.w; \
			__lsns_assert( i >= 0 && i < N );\
			res += fn(( data )[i]); \
		case 14: \
			i = lut4.z;\
			__lsns_assert( i >= 0 && i < N );\
			res += fn(( data )[i] ); \
		case 13: \
			i = lut4.y; \
			__lsns_assert( i >= 0 && i < N );\
			res += fn(( data )[i] ); \
		case 12: \
			i = lut4.x; \
			__lsns_assert( i >= 0 && i < N );\
			res += fn(( data )[i] ); \
		case 11: \
			i = lut3.w; \
			__lsns_assert( i >= 0 && i < N );\
			res += fn(( data )[i]); \
		case 10: \
			i = lut3.z;\
			__lsns_assert( i >= 0 && i < N );\
			res += fn(( data )[i] ); \
		case 9: \
			i = lut3.y; \
			__lsns_assert( i >= 0 && i < N );\
			res += fn(( data )[i] ); \
		case 8: \
			i = lut3.x; \
			__lsns_assert( i >= 0 && i < N );\
			res += fn(( data )[i] ); \
		case 7: \
			i = lut2.w; \
			__lsns_assert( i >= 0 && i < N );\
			res += fn(( data )[i] ); \
		case 6: \
			i = lut2.z; \
			__lsns_assert( i >= 0 && i < N );\
			res += fn(( data )[i] ); \
		case 5: \
			i = lut2.y; \
			__lsns_assert( i >= 0 && i < N );\
			res += fn(( data )[i] ); \
		case 4: \
			i = lut2.x; \
			__lsns_assert( i >= 0 && i < N );\
			res += fn(( data )[i] ); \
		case 3: \
			i = lut1.w; \
			__lsns_assert( i >= 0 && i < N );\
			res += fn(( data )[i] ); \
		case 2: \
			i = lut1.z; \
			__lsns_assert( i >= 0 && i < N );\
			res += fn(( data )[i] ); \
		case 1: \
			i = lut1.y; \
			__lsns_assert( i >= 0 && i < N );\
			res += fn(( data )[i] ); \
		default:; \
	} \
}
///////////////////////////////////////////////////////////////////////////////
// +'lsns_fastsum2' unrolls the calculation of total sum of series that consists 
// of up-to-24 elements and stores the results into 'res1' and 'res2'. The series 
// is located in 'data' array. 'fn1' and 'fn2' are the methods to extract the 
// actual numbers from 'data[i]'. The 'lut' keeps the actual size (<24) of the 
// series and reference to actual position (offset) of the particular element in 
// the array 'data'.
#define lsns_fastsum2( fn1, fn2, data, lut, res1, res2, N )\
{\
	int4 lut1 = ( lut )[0], lut2, lut3, lut4, lut5, lut6; \
	int n = lut1.x; \
	switch( n ){ \
		case 23: case 22: case 21: case 20: \
			lut6 = ( lut )[5]; \
		case 19: case 18: case 17: case 16: \
			lut5 = ( lut )[4]; \
		case 15: case 14: case 13: case 12: \
			lut4 = ( lut )[3]; \
		case 11: case 10: case 9: case 8: \
			lut3 = ( lut )[2]; \
		case 7: case 6: case 5: case 4: \
			lut2 = ( lut )[1]; \
	} \
	( res1  ) = 0; \
	( res2  ) = 0; \
	int i = 0; \
	switch( n ){ \
		case 23: \
			i = lut6.w; \
			__lsns_assert( i >= 0 && i < N );\
			res1 += fn1(( data )[i]); \
			res2 += fn2(( data )[i]); \
		case 22: \
			i = lut6.z; \
			__lsns_assert( i >= 0 && i < N );\
			res1 += fn1(( data )[i]); \
			res2 += fn2(( data )[i]); \
		case 21: \
			i = lut6.y; \
			__lsns_assert( i >= 0 && i < N );\
			res1 += fn1(( data )[i]); \
			res2 += fn2(( data )[i]); \
		case 20: \
			i = lut6.x; \
			__lsns_assert( i >= 0 && i < N );\
			res1 += fn1(( data )[i]); \
			res2 += fn2(( data )[i]); \
		case 19: \
			i = lut5.w; \
			__lsns_assert( i >= 0 && i < N );\
			res1 += fn1(( data )[i]); \
			res2 += fn2(( data )[i]); \
		case 18: \
			i = lut5.z; \
			__lsns_assert( i >= 0 && i < N );\
			res1 += fn1(( data )[i]); \
			res2 += fn2(( data )[i]); \
		case 17: \
			i = lut5.y; \
			__lsns_assert( i >= 0 && i < N );\
			res1 += fn1(( data )[i]); \
			res2 += fn2(( data )[i]); \
		case 16: \
			i = lut5.x; \
			__lsns_assert( i >= 0 && i < N );\
			res1 += fn1(( data )[i]); \
			res2 += fn2(( data )[i]); \
		case 15: \
			i = lut4.w; \
			__lsns_assert( i >= 0 && i < N );\
			res1 += fn1(( data )[i]); \
			res2 += fn2(( data )[i]); \
		case 14: \
			i = lut4.z;\
			__lsns_assert( i >= 0 && i < N );\
			res1 += fn1(( data )[i]); \
			res2 += fn2(( data )[i]); \
		case 13: \
			i = lut4.y; \
			__lsns_assert( i >= 0 && i < N );\
			res1 += fn1(( data )[i]); \
			res2 += fn2(( data )[i]); \
		case 12: \
			i = lut4.x; \
			__lsns_assert( i >= 0 && i < N );\
			res1 += fn1(( data )[i]); \
			res2 += fn2(( data )[i]); \
		case 11: \
			i = lut3.w; \
			__lsns_assert( i >= 0 && i < N );\
			res1 += fn1(( data )[i]); \
			res2 += fn2(( data )[i]); \
		case 10: \
			i = lut3.z;\
			__lsns_assert( i >= 0 && i < N );\
			res1 += fn1(( data )[i]); \
			res2 += fn2(( data )[i]); \
		case 9: \
			i = lut3.y; \
			__lsns_assert( i >= 0 && i < N );\
			res1 += fn1(( data )[i]); \
			res2 += fn2(( data )[i]); \
		case 8: \
			i = lut3.x; \
			__lsns_assert( i >= 0 && i < N );\
			res1 += fn1(( data )[i]); \
			res2 += fn2(( data )[i]); \
		case 7: \
			i = lut2.w; \
			__lsns_assert( i >= 0 && i < N );\
			res1 += fn1(( data )[i]); \
			res2 += fn2(( data )[i]); \
		case 6: \
			i = lut2.z; \
			__lsns_assert( i >= 0 && i < N );\
			res1 += fn1(( data )[i]); \
			res2 += fn2(( data )[i]); \
		case 5: \
			i = lut2.y; \
			__lsns_assert( i >= 0 && i < N );\
			res1 += fn1(( data )[i]); \
			res2 += fn2(( data )[i]); \
		case 4: \
			i = lut2.x; \
			__lsns_assert( i >= 0 && i < N );\
			res1 += fn1(( data )[i]); \
			res2 += fn2(( data )[i]); \
		case 3: \
			i = lut1.w; \
			__lsns_assert( i >= 0 && i < N );\
			res1 += fn1(( data )[i]); \
			res2 += fn2(( data )[i]); \
		case 2: \
			i = lut1.z; \
			__lsns_assert( i >= 0 && i < N );\
			res1 += fn1(( data )[i]); \
			res2 += fn2(( data )[i]); \
		case 1: \
			i = lut1.y; \
			__lsns_assert( i >= 0 && i < N );\
			res1 += fn1(( data )[i]); \
			res2 += fn2(( data )[i]); \
		default:; \
	} \
}
///////////////////////////////////////////////////////////////////////////////
// +'lsns_euler' 1-step Euler method: y = y+step*f
#define lsns_euler( y, f, step ) \
	( y )+( step )*( f )
///////////////////////////////////////////////////////////////////////////////
// +'lsns_exp_euler' 1-step exponential Euler method: y = exp(-step/t )*( y-y0 )+y0
#define lsns_exp_euler( y, y0, step, t ) \
	( __lsns_assert( t > 0.0 ), lsns_exp(-lsns_div( step, t ))*(( y )-( y0 ))+( y0 ))

#endif /*__LSNS_MATH_H*/

