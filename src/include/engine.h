#ifndef __LSNS_ENGINE_H
#define __LSNS_ENGINE_H

#include "config.h"
#include "gateproc.h"
#include "ionproc.h"
#include "pumpproc.h"

#define MAX_CHAN_PER_CELL	15				/*<24*/
#define MAX_CHAN_PER_PUMP	7				/*<24*/
#define MAX_IPUMP_PER_CELL	3				/*<24*/

///////////////////////////////////////////////////////////////////////////////
// gate types of ion channel
//	x - type of activation;
//	y - parameters of activation;
//	z - type of inactivation;
//	w - parameters of inactivation;
#define _gate_typem( t ) ( t ).x
#define _gate_parm( t ) ( t ).y
#define _gate_typeh( t ) ( t ).z
#define _gate_parh( t ) ( t ).w

///////////////////////////////////////////////////////////////////////////////
// parameters of gate variables of ion channel
//	x - activation;
//	y - power of activation;
//	z - inactivation;
//	w - power of inactivation.
#define _gate_m( mh ) ( mh ).x
#define _gate_h( mh ) ( mh ).y
#define _gate_powm( mh ) ( mh ).z
#define _gate_powh( mh ) ( mh ).w

///////////////////////////////////////////////////////////////////////////////
// parameters of the ion channel
//	x - maximal conductance;
//	y - conductance;
//	z - current;
//	w - G*Eds production
#define _chan_gmax( g ) ( g ).x
#define _chan_g( g ) ( g ).y
#define _chan_ge( g ) ( g ).z
#define _chan_i( g ) ( g ).w

///////////////////////////////////////////////////////////////////////////////
// shared parameters (indices of variables into correspondent arrays) for ion channels:
//	x - membrane potential;
//	y - resting potential;
//	z - concentration of ions inside the cell for activation;
//	w - concentration of ions inside the cell for inactivation
#define _chan_lut_v( sh ) ( sh ).x
#define _chan_lut_e( sh ) ( sh ).y
#define _chan_lut_inm( sh ) ( sh ).z
#define _chan_lut_inh( sh ) ( sh ).w

///////////////////////////////////////////////////////////////////////////////
// type of ion dynamics: (channels)
//	x - type of ions pump (Na-pump, Ca-pump, Ca-pump_1 etc);
//	y - type of reversal potential (non-specific, specific etc);
//	z - reserved;
//	w - reserved. 
#define _ions_typepump( tp ) ( tp ).x 
#define _ions_typeeds( tp ) ( tp ).y

///////////////////////////////////////////////////////////////////////////////
// parameters of ions dynamics:
//	x - reversal potential (Eds);
//	y - concentration of ions inside the cell;
//	z - concentration of ions outside the cell;
//	w - RT/Fz constant for specific ions
#define _ions_eds( e ) ( e ).x
#define _ions_in( e ) ( e ).y
#define _ions_out( e ) ( e ).z
#define _ions_rtfz( e ) ( e ).w

///////////////////////////////////////////////////////////////////////////////
// parameters of ion current:
//	x - pump current;
//	y - channels current;
//	z - time constant of ions dynamics;
//	w - reserved.
#define _ions_ipump( i ) ( i ).x
#define _ions_ichan( i ) ( i ).y
#define _ions_tau( i ) ( i ).z

///////////////////////////////////////////////////////////////////////////////
// shared parameters (indices of variables into correspondent arrays) for ion dynamics:
//	x - membrane potential;
//	y - reserved;
//	z - reserved;
//	w - reserved
#define _ions_lut_v( sh ) ( sh ).x

///////////////////////////////////////////////////////////////////////////////
// cell_v: V(x), C (y), spike onset(z), reserved(w)
#define _cell_v( v ) ( v ).x
#define _cell_c( v ) ( v ).y
#define _cell_spike( v ) ( v ).z
#define _cell_iadd( v ) ( v ).w

///////////////////////////////////////////////////////////////////////////////
// 
typedef struct __ions_proc{
	// local variables (read-only)
	int4 *ions_type;					// type of ions: x - pump type, y - eds type
	int4 *ions_shared;					// indices of shared variables: x - cell_v array
	int4 *gchan_lut;					// look-up-table of channel current: x - counter, the rest are actual indices of chan_g array
	// local variables (read/write)
	float4 *ions_e;						// reversal potential
	float4 *ions_i;						// ion currents
	// shared variables
	float4 *chan_g;						// conductance
	float4 *cell_v;						// membrane potential
} ionproc;

///////////////////////////////////////////////////////////////////////////////
// 
typedef struct __channel_proc{
	// local variables (read-only)
	int4 *chan_type;					// type of channel
	int4 *chan_shared;					// indices of shared variables: x - cell_v, y - ions_e/eds, z - ions_e/in for m, w - ions_e/in for h
	// local variables (read-write)
	float4 *chan_g;						// conductance
	float4 *chan_mh;					// gate variables
	// shared variables
	float4 *cell_v;						// membrane potential
	float4 *ions_e;						// reversal potential
} chanproc;

///////////////////////////////////////////////////////////////////////////////
// 
typedef struct __cell_proc{
	// local variables (read-only)
	int4 *gchan_lut;					// look-up-table of channel current: x - counter, the rest are actual indices of chan_g array
	int4 *ipump_lut;					// look-up-table of pump current: x - counter, the rest are actual indices of ions_i array
	// local variables (read-write)
	float4 *cell_v;						// membrane potential
	// shared variables
	float4 *chan_g;						// channel currents
	float4 *ions_i;						// pump currents
} cellproc;


///////////////////////////////////////////////////////////////////////////////
// 
typedef struct __network_proc{
	float Step
	int MaxIons;
	int MaxChan;
	int MaxCells;
	gatepar Gates[????];
	ionproc Ions;
	chanproc Channels;
	cellproc Cells;
} netproc;

#if !defined (__CUDA__)

///////////////////////////////////////////////////////////////////////////////
// calculate the properties of ions dynamics such as pump current, concentration 
// of ions inside the cell, etc.
extern void ions_kernel( float step, int index, ionproc *data );

///////////////////////////////////////////////////////////////////////////////
// calculate the properties of channel and synaptic currents such as conductance, 
// current etc.
extern void chan_kernel( float step, int index, chanproc *data );

///////////////////////////////////////////////////////////////////////////////
// calculate the properties of Hodgkin-Huxley cell such as membrane potential,
// onset of the spike etc.
extern void cell_kernel( float step, int index, cellproc *data );

#endif

#endif /*__LSNS_ENGINE_H*/

