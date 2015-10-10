#ifndef __GATES_H
#define __GATES_H

#include "config.h"

#define LSNS_MAX_GATEPARS 64	/*maximal number of parameters for different types of gate variables*/
///////////////////////////////////////////////////////////////////////////////
// +gates ids and brief descriptions
enum __lsns_gate_types{
	LSNS_NOGATE			= 0,						// none (gate is always 1)
	LSNS_BYPASSGATE			= 1,						// bypass (bypass calulation)
	LSNS_GENERIC_INSTANT		= 2,						// generic instant
	LSNS_GENERIC_T			= 3,						// generic standard
	LSNS_GENERIC_TMOD		= 4,						// generic modified
	LSNS_GENERIC_TAMOD		= 5,						// generic A-current
	LSNS_ALPHABETA_INSTANT		= 6,						// alpha/beta instant
	LSNS_ALPHABETA_T		= 7,						// alpha/beta standard
	LSNS_ZGENERIC_INSTANT		= 8,						// z-current instant
	LSNS_ZGENERIC_T			= 9,						// z-current standard
	LSNS_ZAPHABETA_INSTANT		= 10,						// z-current (alpha/beta) instant
	LSNS_ZAPHABETA_T		= 11,						// z-current (alpha/beta) standard
	LSNS_PS_NMDA			= 12,						// nmda post-synapse
	LSNS_SYNAPSE			= 13,						// synaptic current
	LSNS_MAX_GATES,
};
extern const char *lsns_gate_types[LSNS_MAX_GATES];
///////////////////////////////////////////////////////////////////////////////
// +extract parameters for all generic descriptions from the structure 'gatepar'
//=============================================================================
// get half-voltage for gate variable 
#define _ggatev12( par ) ( par ).Par1.x
// get slope for gate variable 
#define _ggateslp( par ) ( par ).Par1.y
// get voltage threshold for time constant
#define _ggatevtr( par ) ( par ).Par3.x
// get t0 for time constant
#define _ggatet0( par ) ( par ).Par1.z
// get tmax for time constant
#define _ggatetmax( par ) ( par ).Par1.w
// get tup for time constant
#define _ggatetup( par ) ( par ).Par3.y
// get half-voltage for time constant
#define _ggatev12t( par ) ( par ).Par2.x
// get slope for time constant
#define _ggateslpt( par ) ( par ).Par2.y
// get half-voltage-2 for time constant
#define _ggatev12t2( par ) ( par ).Par2.z
// get slope-2 for time constant
#define _ggateslpt2( par ) ( par ).Par2.w
///////////////////////////////////////////////////////////////////////////////
// +extract parameters for alpha/beta descriptions from the structure 'gatepar'
//=============================================================================
//	Alpha/Beta = (A*v+B)/( exp(( v+C )/D )+E)
// get parameter T0 for alpha component of gate variable 
#define _abgatet0( par ) ( par ).Par1.x
// get parameter Tmax for beta component of gate variable 
#define _abgatetmax( par ) ( par ).Par1.y
// get A for alpha component of gate variable 
#define _abgateAa( par ) ( par ).Par1.z
// get B for alpha component of gate variable 
#define _abgateBa( par ) ( par ).Par1.w
// get C for alpha component of gate variable 
#define _abgateCa( par ) ( par ).Par2.x
// get D for alpha component of gate variable 
#define _abgateDa( par ) ( par ).Par2.y
// get E for alpha component of gate variable 
#define _abgateEa( par ) ( par ).Par2.z
// get A for beta component of gate variable 
#define _abgateAb( par ) ( par ).Par2.w
// get B for beta component of gate variable 
#define _abgateBb( par ) ( par ).Par3.x
// get C for beta component of gate variable 
#define _abgateCb( par ) ( par ).Par3.y
// get D for beta component of gate variable 
#define _abgateDb( par ) ( par ).Par3.z
// get E for beta component of gate variable 
#define _abgateEb( par ) ( par ).Par3.w
///////////////////////////////////////////////////////////////////////////////
// +extract parameters for all z-descriptions from the structure 'gatepar'
//=============================================================================
// get parameter Alpha z-gate variable
#define _zgateA( par ) ( par ).Par1.z
// get parameter Beta z-gate variable
#define _zgateB( par ) ( par ).Par1.w
// get parameter Lymbda z-gate variable
#define _zgateL( par ) ( par ).Par1.x
// get halh-valtage z-gate variable
#define _zgatev12( par ) ( par ).Par1.x
// get parameter Gamma z-gate variable
#define _zgateG( par ) ( par ).Par1.y
// get slope z-gate variable
#define _zgateslp( par ) ( par ).Par1.y
// get parameter T0 z-gate variable
#define _zgatet0( par ) ( par ).Par2.x
// get parameter Tmax z-gate variable
#define _zgatetmax( par ) ( par ).Par2.y
///////////////////////////////////////////////////////////////////////////////
// +structure 'gatepar' to store all constant parameters for gata variables
// of all types
typedef struct __lsns_align( 16 ) __gate_par{
	float4 Par1;
	float4 Par2;
	float4 Par3;
} gatepar;

#endif /*__GATES_H*/
