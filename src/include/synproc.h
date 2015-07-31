#ifndef __SYNPROC_H
#define __SYNPROC_H

#include "lsnsmath.h"
#include "synapses.h"

__lsns_inline float proc_synsum( float4 &v, float4 &w  )
{
	return v.x*w.x+v.y*w.y+v.z*w.z+v.w*w.w;
}

#endif /*__SYNPROC_H*/
