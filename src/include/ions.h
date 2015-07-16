#ifndef __IONS_H
#define __IONS_H

#include "config.h"

#define LSNS_MAX_IONPARS		16
///////////////////////////////////////////////////////////////////////////////
// +ions ids and descriptions
enum __lsns_eds_types{
	LSNS_NOSPEC_EDS			= 0,						// non-specific ions
	LSNS_NA_IONS			= 1,						// na ions
	LSNS_K_IONS			= 2,						// k ions
	LSNS_CA_IONS			= 3,						// ca ions
	LSNS_CL_IONS			= 4,						// cl ions
	LSNS_MG_IONS			= 5,						// mg ions
	LSNS_MAX_IONS,
};
extern const char *lsns_ions_types[LSNS_MAX_IONS];
///////////////////////////////////////////////////////////////////////////////
// +pumps ids and descriptions
enum __lsns_pump_types{
	LSNS_NO_DYN			= 0,						// none
	LSNS_NA_PUMP			= 1,						// apump = 1
	LSNS_CA_PUMP_STD		= 2,						// apump = 0.0045
	LSNS_CA_PUMP_SHELL		= 3,						// apump = T*(Ca+K)/(Ca+K+B )*KCa
	LSNS_MAX_PUMP,
};
extern const char *lsns_pump_types[LSNS_MAX_PUMP];
///////////////////////////////////////////////////////////////////////////////
// Extract parameters from the structure 'ionspar'
//=============================================================================
// get parameter In0 for all-type-pumps
#define _pump_in0( par ) ( par ).Par1.x
// get parameter Rpump for all-type-pumps
#define _pump_rpump( par ) ( par ).Par1.y
// get parameter Kp^3 for Na-pump
#define _napump_kp3( par ) ( par ).Par1.z
// get parameter pre-calculated parameter In03 = In0^3/(In0^3-Kp^3) for Na-pump 
#define _napump_in03( par ) ( par ).Par1.w
// get time constant of ions dynamics
#define _pump_tau( par ) ( par ).Par2.x
// get apump coefficient of ions dynamics
#define _pump_apump( par ) ( par ).Par2.y
// get k coefficient of ions dynamics
#define _pump_k( par ) ( par ).Par2.z
// get b coefficient of ions dynamics
#define _pump_b( par ) ( par ).Par2.w

///////////////////////////////////////////////////////////////////////////////
// Structure 'gate_par' to store all constant parameters for gata variables
// of all types
typedef struct __lsns_align( 16 ) __ions_par{
	float4 Par1;
	float4 Par2;
} ionspar;


#endif /*__IONS_H*/
