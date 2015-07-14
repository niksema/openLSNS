#ifndef __VIEWPROC_H
#define __VIEWPROC_H

#include "lsnsmath.h"
#include "views.h"

///////////////////////////////////////////////////////////////////////////////
// 
__lsns_inline void store_data( float4 *dat, int4 lut, float4 &res, int imax )
{
	int4vw *_lut = ( int4vw *)&lut;
	float4vw *_res = ( float4vw *)&res;
	float4vw *_dat = 0L;
	float4 d = {0};
	for( int i = 0, _index = -1; i < 4; ++i ){
		int pos = ( _lut->pDat[i] ) >> 30;
		int index = ( _lut->pDat[i] ) & 0x3FFFFFFF;
		__lsns_assert( index >= 0 && index < imax );
		d = ( _index != index )? dat[index]: d;
		_dat = ( float4vw *)&d;
		_res->pDat[i] = _dat->pDat[pos];
		_index = index;
	}
}

#endif /*__VIEWPROC_H*/
