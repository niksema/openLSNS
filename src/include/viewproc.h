#ifndef __VIEWPROC_H
#define __VIEWPROC_H

#include "lsnsmath.h"
#include "views.h"

///////////////////////////////////////////////////////////////////////////////
// +'store_data': saves the data from array 'src' to 'res' (both are variables  
// in float4 format) according to look-up-table defined in 'lut'.
// LUT format is: bits 31..30 are coding the field in each float4 structure
// (00 - x, 01 - y, 10 - z, 11 - w) in array 'src'; 
// bits 29..0 are coding the offset in an array 'src'.
__lsns_inline float4 get_data( float4 *src, int4 lut, int imax )
{
	int4vw *_lut = ( int4vw *)&lut;
	float4 d = {0};
	float4 res = {0};
	float4vw *_res = ( float4vw *)&res;
	float4vw *_src = ( float4vw *)&d;
	for( int i = 0, pos = -1, index = 0, _index = -1; i < 4; _index = index, ++i ){
		pos = ( _lut->pDat[i] ) >> 30;						// get offset (00 - x, 01 - y, 10 - z, 11 - w)
		index = ( _lut->pDat[i] ) & 0x3FFFFFFF;					// get index
		__lsns_assert( pos >= 0 && pos < 4 );					// DEBUG: check the range for 'pos' variable
		__lsns_assert( index >= 0 && index < imax );				// DEBUG: check the range for 'index' variable
		d = ( _index != index )? src[index]: d;					// load data from src array to variable d if needed
		_res->pDat[i] = _src->pDat[pos];					// store loaded data
	}
	return res;
}

#endif /*__VIEWPROC_H*/
