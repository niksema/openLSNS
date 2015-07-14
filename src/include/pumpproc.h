#ifndef __PUMPPROC_H
#define __PUMPPROC_H

#include "lsnsmath.h"
#include "pumps.h"

///////////////////////////////////////////////////////////////////////////////
// Performes no calculation. Always returns 0.
//=============================================================================
__lsns_inline float calc_nopump( pumppar &par, float in )
{
	return 0;
}

///////////////////////////////////////////////////////////////////////////////
// Calculate current of Ca-pump (description+reference)
//=============================================================================
__lsns_inline float calc_capump( pumppar &par, float in )
{
	return _pump_rpump( par )*( in-_pump_in0( par ));
}

///////////////////////////////////////////////////////////////////////////////
// Calculate current of Na-pump (description+reference)
//=============================================================================
__lsns_inline float calc_napump( pumppar &par, float in )
{
	float kp = lsns_pow( _napump_kp( par ), 3 );
	in = lsns_pow( in, 3 );
	return _pump_rpump( par )*( in/( in+kp )-_napump_in03( par ));
}

///////////////////////////////////////////////////////////////////////////////
// Calculate the current for any type of ions pump 
//----------------------------------------------------------------------------
// The general description of ions dynamics is calculated as follow:
// T*dIn/dt = -Rchan*Ichan-Rpump*Ipump
// The method 'lsns_ipump' declared here calulates the second term of the equation,
// namely: Rpump*Ipump for either Na- or Ca- pumps.
//=============================================================================
#define lsns_ipump( type, par, in, apump, ipump ) \
	switch( type ){ \
		case LSNS_CA_PUMP: \
			apump = 1; /*todo: calculate the apump*/ \
			ipump = calc_capump( par, in ); \
			break; \
		case LSNS_NA_PUMP: \
			apump = 1; /*todo: calculate the apump*/ \
			ipump = calc_napump( par, in ); \
			break; \
		default: \
			apump = 0; \
			ipump = 0;\
	} 

#endif /*__PUMPPROC_H*/
