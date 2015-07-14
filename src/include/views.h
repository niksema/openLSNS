#ifndef __VIEWS_H
#define __VIEWS_H

#include "config.h"

typedef union __lsns_align( 16 ) __int4_vw{
	int4 Dat;
	int pDat[4];
} int4vw;

typedef union __lsns_align( 16 ) __float4_vw{
	float4 Dat;
	float pDat[4];
} float4vw;

#endif /*__VIEWS_H*/
