#ifndef __IONPROC_H
#define __IONPROC_H

#include "lsnsmath.h"
#include "ions.h"

// calculate Nernst potential E = (R*T)/(z*F)*ln([Out]/[In])
#define lsns_eds( rtfz, in, out ) \
	( rtfz )*lsns_log(( out )/( in ))

// calculate ions current that is involved into ions dynamics
#define lsns_ichan( a, v, eds, gchan, ichan ) \
	ichan = a*(( v )-( eds ))*( gchan );

#endif /*__IONPROC_H*/
