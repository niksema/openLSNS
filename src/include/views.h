#ifndef __VIEWS_H
#define __VIEWS_H

#include "config.h"

///////////////////////////////////////////////////////////////////////////////
// +'float4vw' is dirty hack that provides easy access to the fields of int4 
// structure
typedef union __lsns_align( 16 ) __int4_vw{
	int4 Dat;
	int pDat[4];
} int4vw;
///////////////////////////////////////////////////////////////////////////////
// +'float4vw' is dirty hack that provides easy access to the fields of float4 
// structure
typedef union __lsns_align( 16 ) __float4_vw{
	float4 Dat;
	float pDat[4];
} float4vw;

#endif /*__VIEWS_H*/
