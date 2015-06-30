#ifndef __LSNS_MATH_H
#define __LSNS_MATH_H

#include "config.h"

#if defined( __CUDA__ )
	#include <cuda.h>
	
	#define lsns_div( x, y ) fdividef(( x ),( y ))
	#define lsns_pow( x, y ) powf(( x ),( y ))
	#define lsns_exp( x ) expf( x )
	#define lsns_cosh( x ) coshf( x )
	#define lsns_log( x ) __logf( x )
#else
	#include <math.h>
	
	#define lsns_div( x, y ) ( x )/( y )
	#define lsns_pow( x, y ) pow(( x ),( y ))
	#define lsns_exp( x ) exp( x )
	#define lsns_cosh( x ) cosh( x )
	#define lsns_log( x ) log( x )
#endif

// Unroll the calculation of total sum of series that consists of up-to-16 elements and
// stores the results into 'res'. The series is located in 'data' array. 'fn' is the method 
// to extract the acual number from 'data[i]'. The 'lut' keeps the actual size of the series
// and reference to actual position of the particular element in the array 'data'.
#define lsns_fastsum( fn, data, lut, res )\
{\
	int4 lut1 = ( lut )[0], lut2, lut3, lut4; \
	int n = lut1.x; \
	switch( n ){ \
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
		case 15: \
			i = lut4.w; \
			res += fn(( data )[i]); \
		case 14: \
			i = lut4.z;\
			res += fn(( data )[i] ); \
		case 13: \
			i = lut4.y; \
			res += fn(( data )[i] ); \
		case 12: \
			i = lut4.x; \
			res += fn(( data )[i] ); \
		case 11: \
			i = lut3.w; \
			res += fn(( data )[i]); \
		case 10: \
			i = lut3.z;\
			res += fn(( data )[i] ); \
		case 9: \
			i = lut3.y; \
			res += fn(( data )[i] ); \
		case 8: \
			i = lut3.x; \
			res += fn(( data )[i] ); \
		case 7: \
			i = lut2.w; \
			res += fn(( data )[i] ); \
		case 6: \
			i = lut2.z; \
			res += fn(( data )[i] ); \
		case 5: \
			i = lut2.y; \
			res += fn(( data )[i] ); \
		case 4: \
			i = lut2.x; \
			res += fn(( data )[i] ); \
		case 3: \
			i = lut1.w; \
			res += fn(( data )[i] ); \
		case 2: \
			i = lut1.z; \
			res += fn(( data )[i] ); \
		case 1: \
			i = lut1.y; \
			res += fn(( data )[i] ); \
		default:; \
	} \
}

// 1-step Euler method: y = y+step*f
#define lsns_euler( y, f, step ) \
	( y )+( step )*( f )
// 1-step exponential Euler method : y = exp(-step/t )*( y-y0 )+y0
#define lsns_exp_euler( y, y0, step, t ) \
	lsns_exp(-lsns_div( step, t ))*(( y )-( y0 ))+( y0 )

#endif /*__LSNS_MATH_H*/
