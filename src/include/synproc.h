#ifndef __SYNPROC_H
#define __SYNPROC_H

#include "lsnsmath.h"
#include "synapses.h"

__lsns_inline float proc_synsum( float4 &v, float4 &w  )
{
	return v.x*w.x+v.y*w.y+v.z*w.z+v.w*w.w;
}

///////////////////////////////////////////////////////////////////////////////
// 'lsns_sig_syn' calculates the neurotransmitter release for model of sigma 
// synapse. Parameters v12 and slope are hidden
#define lsns_sig_syn( v )
	lsns_msigm( v, v12, slope )
///////////////////////////////////////////////////////////////////////////////
// 'lsns_synsum4' calculates the total sum of weighted inputs from 4 presynaptic
// neurons.
//
// pCellV - v_data
// pCellLUT[j] - lut_data
// pWall[j] w_data
// LSNS_WSUM_SYN -  fn_vproc=none fn_v=_cell_v
// LSNS_SIGMA_SYN - fn_vproc=lsns_msigm fn_v=_cell_v
// LSNS_PULSE_SYN - fn_vproc=none fn_v=_cell_spike
#define lsns_synsum4( fn_vproc, fn_v, v_data, lut_data, w_data, res )\
	/* load synaptic weights for 4 presynaptic neurons */\
	float4 w = w_data;\
	/*load look-up-table for 4 presynaptic units (neurons/drives/outputs/feedbacks etc)*/\
	int4 vlut = lut_data;\
	/*load raw data from 4 presynaptic units*/\
	float4 v_raw[4] = { v_data[vlut.x], v_data[vlut.y], v_data[vlut.z], v_data[vlut.w] };\
	/*extract input signal from raw data*/\
	float4 v = { fn_v( v_raw[0] ), fn_v( v_raw[1] ), fn_v( v_raw[2] ), fn_v( v_raw[3] ) };\
	/*calculate the summation*/\
	res += fn_vproc( v.x )*w.x+fn_vprocv( v.y )*w.y+fn_vproc( v.z )*w.z+fn_vproc( v.w )*w.w;

#endif /*__SYNPROC_H*/
