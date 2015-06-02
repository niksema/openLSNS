#ifndef __LSNS_MATH_H
#define __LSNS_MATH_H

#include "config.h"

#if defined( __CUDA__ )
	#include <cuda.h>
	
	#define lsns_div( x, y ) fdividef(( x ),( y ))
	#define lsns_pow( x, y ) powf(( x ),( y ))
	#define lsns_exp( x ) expf( x )
	#define lsns_cosh( x ) coshf( x )
#else
	#include <math.h>
	
	#define lsns_div( x, y ) ( x )/( y )
	#define lsns_pow( x, y ) pow(( x ),( y ))
	#define lsns_exp( x ) exp( x )
	#define lsns_cosh( x ) cosh( x )
#endif

// 1-step Euler method: y = y+step*f
#define lsns_euler( y, f, step ) \
	( y )+( step )*( f )
// 1-step exponential Euler method : y = exp(-step/t )*( y-y0 )+y0
#define lsns_exp_euler( y, y0, step, t ) \
	lsns_exp(-lsns_div( step, t ))*(( y )-( y0 ))+( y0 )

#endif /*__LSNS_MATH_H*/