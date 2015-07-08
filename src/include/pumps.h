#ifndef __PUMPS_H
#define __PUMPS_H

#include "config.h"

///////////////////////////////////////////////////////////////////////////////
// Pumps ids and descriptions
//=============================================================================
enum __lsns_pump_types{
	LSNS_NOPUMP				= 0,
	LSNS_CA_PUMP			= 1,
	LSNS_NA_PUMP			= 2,
	LSNS_MAX_PUMP,
};

///////////////////////////////////////////////////////////////////////////////
// Extract parameters from the structure 'pump_par'
//=============================================================================
// get parameter In0 for all-type-pumps
#define _pump_in0( par ) ( par ).Par1.x
// get parameter Rpump for all-type-pumps
#define _pump_rpump( par ) ( par ).Par1.y
// get parameter Kp for Na-pump
#define _napump_kp( par ) ( par ).Par1.z
// get parameter pre-calculated parameter In03 = In0^3/(In0^3-Kp^3) for Na-pump 
#define _napump_in03( par ) ( par ).Par1.w

///////////////////////////////////////////////////////////////////////////////
// Structure 'gate_par' to store all constant parameters for gata variables
// of all types
//=============================================================================
typedef struct __lsns_align( 16 ) __pump_par{
	float4 Par1;
	float4 Par2;
} pump_par;

extern const char *lsns_pump_types[LSNS_MAX_PUMP];

#endif /*__PUMPS_H*/
