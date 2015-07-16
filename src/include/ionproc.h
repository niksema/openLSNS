#ifndef __IONPROC_H
#define __IONPROC_H

#include "lsnsmath.h"
#include "ions.h"

#define KCA 1.E-5f/( 2.f*9.648f )
///////////////////////////////////////////////////////////////////////////////
// +'lsns_eds' calculates Nernst potential E = (R*T)/(z*F)*ln([Out]/[In])
#define lsns_eds( rtfz, in, out ) \
	( rtfz )*lsns_log(( out )/( in ))
///////////////////////////////////////////////////////////////////////////////
// +'lsns_ichan' calculates ions current that is involved into ions dynamics
#define lsns_ichan( apump, v, eds, gchan, ichan ) \
	ichan = apump*(( v )-( eds ))*( gchan );
//=================== pump calculation macros =================================
///////////////////////////////////////////////////////////////////////////////
// +performes no calculation. Always returns 0.
__lsns_inline float proc_nopump( ionspar &par, float in )
{
	return 0;
}
///////////////////////////////////////////////////////////////////////////////
// +calculates current of Ca-pump
__lsns_inline float proc_capump( ionspar &par, float in )
{
	return _pump_rpump( par )*( in-_pump_in0( par ));
}
///////////////////////////////////////////////////////////////////////////////
// +calculates current of Na-pump
__lsns_inline float proc_napump( ionspar &par, float in )
{
	in = lsns_pow( in, 3 );
	return _pump_rpump( par )*( in/( in+_napump_kp3( par ))-_napump_in03( par ));
}
///////////////////////////////////////////////////////////////////////////////
// +calculates the ipump and apump for any type of ions dynamics
//----------------------------------------------------------------------------
// The general description of ions dynamics is calculated as follow:
// T*dIn/dt = -aPump*iChan-rPump*iPump
// The macros 'lsns_ipump' declared here calulates the second term of the equation,
// namely: Rpump*Ipump for either Na- or Ca- pumps.
#define lsns_ipump( type, par, in, apump, ipump, kca ) \
	switch( type ){ \
		case LSNS_NA_PUMP: \
			apump = _pump_apump( par ); \
			ipump = proc_napump( par, in ); \
			break; \
		case LSNS_CA_PUMP_STD: \
			apump = _pump_apump( par ); \
			ipump = proc_capump( par, in ); \
			break; \
		case LSNS_CA_PUMP_SHELL: \
			apump =  _pump_tau( par )*(in+_pump_k( par ))/(in+_pump_k( par )+_pump_b( par ))*( kca ); \
			ipump = proc_capump( par, in ); \
			break; \
		default: \
			apump = 0; \
			ipump = 0;\
	} 

#endif /*__IONPROC_H*/
