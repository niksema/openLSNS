#ifndef __GATES_H
#define __GATES_H

#include "config.h"

///////////////////////////////////////////////////////////////////////////////
// Gates ids and descriptions
//=============================================================================
enum __lsns_gate_types{
	LSNS_NOGATE				= 0,
	LSNS_BYPASSGATE			= 1,
	LSNS_GENERIC_INSTANT		= 2,
	LSNS_GENERIC_T			= 3,
	LSNS_GENERIC_TMOD			= 4,
	LSNS_GENERIC_TAMOD		= 5,
	LSNS_ALPHABETA_INSTANT		= 6,
	LSNS_ALPHABETA_T			= 7,
	LSNS_ZGENERIC_INSTANT		= 8,
	LSNS_ZGENERIC_T			= 9,
	LSNS_ZAPHABETA_INSTANT		= 10,
	LSNS_ZAPHABETA_T			= 11,
	LSNS_PS_NMDA			= 12,
	LSNS_MAX_GATE,
};

///////////////////////////////////////////////////////////////////////////////
// Extract parameters for all generic descriptions from the structure 'gate_par'
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
// Extract parameters for alpha/beta descriptions from the structure 'gate_par'
//=============================================================================
// get half-voltage for alpha component of gate variable 
#define _abgatev12a( par ) ( par ).Par1.x
// get half-voltage for beta component of gate variable 
#define _abgatev12b( par ) ( par ).Par2.x
// get slope for alpha component of gate variable 
#define _abgateslpa( par ) ( par ).Par1.y
// get slope for beta component of gate variable 
#define _abgateslpb( par ) ( par ).Par2.y
// get parameter A for alpha component of gate variable 
#define _abgateAa( par ) ( par ).Par1.z
// get parameter A for beta component of gate variable 
#define _abgateAb( par ) ( par ).Par2.z
// get parameter B for alpha component of gate variable 
#define _abgateBa( par ) ( par ).Par1.w
// get parameter B for beta component of gate variable 
#define _abgateBb( par ) ( par ).Par2.w
// get parameter C for alpha component of gate variable 
#define _abgateCa( par ) ( par ).Par3.x
// get parameter C for beta component of gate variable 
#define _abgateCb( par ) ( par ).Par3.y
// get parameter T0 for alpha component of gate variable 
#define _abgatet0( par ) ( par ).Par3.z
// get parameter Tmax for beta component of gate variable 
#define _abgatetmax( par ) ( par ).Par3.w

///////////////////////////////////////////////////////////////////////////////
// Extract parameters for all z-descriptions from the structure 'gate_par'
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
// Structure 'gate_par' to store all constant parameters for gata variables
// of all types
//=============================================================================
typedef struct __gate_par{
	float4 Par1;
	float4 Par2;
	float4 Par3;
} gate_par;

extern const char *lsns_gate_types[LSNS_MAX_GATE];

#endif /*__GATES_H*/
