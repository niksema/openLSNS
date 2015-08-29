#ifndef __SYNAPSES_H
#define __SYNAPSES_H

#include "config.h"

#define LSNS_MAX_SYNPARS 64
///////////////////////////////////////////////////////////////////////////////
// +gates ids and brief descriptions
enum __lsns_synapse_types{
	LSNS_NOSYN			= 0,						// none 
	LSNS_BYPASS_SYN			= 1,						// bypass model of synapse (for network unit like drives, outputs etc): dt = 1; edt = 0;
	LSNS_PULSE_SYN			= 2,						// pulse model of synapse: dt = 1; edt = exp( -step/tmax )
	LSNS_MAX_SYNS,
};

///////////////////////////////////////////////////////////////////////////////
// +extract parameters for all synaptic current from the structure 'synpar'
//=============================================================================
// get rate of transmitter release
#define _synA( par ) ( par ).Par1.x
// get parameter Edt: exp( -step/T ) for fast synapse or 1-step/T for other types of synapses
#define _synEdt( par ) ( par ).Par1.y
// get parameter Dt: 1 for fast synapse or step/T for other types of synapses
#define _synDt( par ) ( par ).Par1.z
// get parameter Dt: 1 for fast synapse or step/T for other types of synapses
#define _synTmax( par ) ( par ).Par1.w
///////////////////////////////////////////////////////////////////////////////
// structure 'synpar' to store all constant parameters for synapses
// of all types
typedef struct __lsns_align( 16 ) __synapres_par{
	float4 Par1;
} synpar;

#endif /*__SYNAPSES_H*/
