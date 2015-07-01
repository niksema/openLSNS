#ifndef __LSNS_ENGINE_H
#define __LSNS_ENGINE_H

#include "config.h"
#include "gateproc.h"
#include "ionproc.h"
#include "pumpproc.h"

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
// shared parameters (indices of variables into correspondent arrays) for ion channels:
//	x - membrane potential;
//	y - resting potential;
//	z - concentration of ions inside the cell for activation;
//	w - concentration of ions inside the cell for inactivation
#define _gate_v( sh ) ( sh ).x
#define _gate_e( sh ) ( sh ).y
#define _gate_inm( sh ) ( sh ).z
#define _gate_inh( sh ) ( sh ).w

///////////////////////////////////////////////////////////////////////////////
// parameters of the ion channel
//	x - maximal conductance;
//	y - conductance;
//	z - current;
//	w - G*Eds production
#define _gate_gmax( g ) ( g ).x
#define _gate_g( g ) ( g ).y
#define _gate_ge( g ) ( g ).z
#define _gate_i( g ) ( g ).w

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
//	z - reserved;
//	w - reserved. 
#define _ions_ipump( i ) ( i ).x
#define _ions_ichan( i ) ( i ).y

///////////////////////////////////////////////////////////////////////////////
// type of ion dynamics: (channels)
//	x - type of ions pump (Na-pump, Ca-pump, Ca-pump_1 etc);
//	y - type of reversal potential (non-specific, specific etc);
//	z - reserved;
//	w - reserved. 
#define _ions_typepump( tp ) ( tp ).x 
#define _ions_typeeds( tp ) ( tp ).y

///////////////////////////////////////////////////////////////////////////////
// cell_v: V(x), C (y), spike onset(z), reserved(w)
#define _cell_v( v ) ( v ).x
#define _cell_c( v ) ( v ).y
#define _cell_spike( v ) ( v ).z
#define _cell_iadd( v ) ( v ).w

#endif /*__LSNS_ENGINE_H*/

