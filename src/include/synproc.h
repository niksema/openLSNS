#ifndef __SYNPROC_H
#define __SYNPROC_H

#include "lsnsmath.h"
#include "synapses.h"

__lsns_inline float proc_synsum1( float4 &v, float4 &w  )
{
	return v.x*w.x+v.y*w.y+v.z*w.z+v.w*w.w;
}

__lsns_inline float proc_synsum2( float4 &v, float4 &w  )
{
	/*lsns_sigma( v );*/
	return v.x*w.x+v.y*w.y+v.z*w.z+v.w*w.w;
}

///////////////////////////////////////////////////////////////////////////////
// 'lsns_synsum' calculates the total sum of weighted inputs from 4 presynaptic
// neurons.
//
#define lsns_synsum4( fn_sum, fn_v, v_data, lut_data, w_data, res )\
{\
	/* load synaptic weights for 4 presynaptic neurons */\
	float4 w = w_data;\
	/*load look-up-table for 4 presynaptic units (neurons/drives/outputs/feedbacks etc)*/\
	int4 vlut = lut_data;\
	/*load raw data from 4 presynaptic units*/\
	float4 v_raw[4] = { v_data[vlut.x], v_data[vlut.y], v_data[vlut.z], v_data[vlut.w] };\
	/*extract input signal from raw data*/\
	float4 v = { fn_v( v_raw[0] ), fn_v( v_raw[1] ), fn_v( v_raw[2] ), fn_v( v_raw[3] ) };\
	/*calculate the summation*/\
	res += fn_sum( v, w );\
}

#endif /*__SYNPROC_H*/
