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
// CellV: V(x), C (y), spike onset(z), reserved(w)
#define _cell_v( v ) ( v ).x
#define _cell_c( v ) ( v ).y
#define _cell_spike( v ) ( v ).z
#define _cell_iadd( v ) ( v ).w

///////////////////////////////////////////////////////////////////////////////
// 
typedef struct __lsns_align( 16 ) __ions_data{
	// local variables (read-only)
	int4 __lsns_align( 16 ) *IonsType;			// type of ions: x - pump type, y - eds type
	int4 __lsns_align( 16 ) *IonsShared;			// indices of shared variables: x - CellV array
	int4 __lsns_align( 16 ) *GchanLUT;			// look-up-table of channel current: x - counter, the rest are actual indices of ChanG array
	// local variables (read/write)
	float4 __lsns_align( 16 ) *IonsE;			// reversal potential
	float4 __lsns_align( 16 ) *IonsI;			// ion currents
	// shared variables
	float4 __lsns_align( 16 ) *ChanG;			// conductance
	float4 __lsns_align( 16 ) *CellV;			// membrane potential
} iondat;

///////////////////////////////////////////////////////////////////////////////
// 
typedef struct __lsns_align( 16 ) __channel_data{
	// local variables (read-only)
	int4 __lsns_align( 16 ) *ChanType;			// type of channel
	int4 __lsns_align( 16 ) *ChanShared;			// indices of shared variables: x - CellV, y - IonsE/eds, z - IonsE/in for m, w - IonsE/in for h
	// local variables (read-write)
	float4 __lsns_align( 16 ) *ChanG;			// conductance
	float4 __lsns_align( 16 ) *ChanMH;			// gate variables
	// shared variables
	float4 __lsns_align( 16 ) *CellV;			// membrane potential
	float4 __lsns_align( 16 ) *IonsE;			// reversal potential
} chandat;

///////////////////////////////////////////////////////////////////////////////
// 
typedef struct __lsns_align( 16 ) __cell_data{
	// local variables (read-only)
	int4 __lsns_align( 16 ) *GchanLUT;			// look-up-table of channel current: x - counter, the rest are actual indices of ChanG array
	int4 __lsns_align( 16 ) *IpumpLUT;			// look-up-table of pump current: x - counter, the rest are actual indices of IonsI array;
	// local variables (read-write)
	float4 __lsns_align( 16 ) *CellV;			// membrane potential
	// shared variables
	float4 __lsns_align( 16 ) *ChanG;			// channel currents
	float4 __lsns_align( 16 ) *IonsI;			// ion currents
} celldat;

///////////////////////////////////////////////////////////////////////////////
// 
typedef struct __lsns_align( 16 ) __iobuf{
	// local variables (read-write)
	float4 __lsns_align( 16 ) *IOData;			// data to display
	// local variables (read-only). LUT format: from
	int4 __lsns_align( 16 ) *IonsILUT;			// look-up-table for IonsI
	int4 __lsns_align( 16 ) *IonsELUT;			// look-up-table for IonsE
	int4 __lsns_align( 16 ) *ChanGLUT;			// look-up-table for ChanG
	int4 __lsns_align( 16 ) *ChanMHLUT;			// look-up-table for ChanMH
	int4 __lsns_align( 16 ) *CellVLUT;			// look-up-table for CellV
	// shared variables
	float4 __lsns_align( 16 ) *IonsI;			// ion and pump currents
	float4 __lsns_align( 16 ) *IonsE;			// reversal potential, concentration of ions inside and ouside cell 
	float4 __lsns_align( 16 ) *ChanG;			// channel conguctance, currents etc
	float4 __lsns_align( 16 ) *ChanMH;			// gate variables
	float4 __lsns_align( 16 ) *CellV;			// membrane potential, spike onset etc

} iodat;


///////////////////////////////////////////////////////////////////////////////
// 
typedef struct __lsns_align( 16 ) __network_data{
	// global memory on host
	iondat Ions;						// ions
	chandat Channels;					// channels
	celldat Cells;						// neurons (compartments)
	iodat IOData;						// data buffer to display the results of simulation
	// global memory on device
	struct __lsns_align( 16 ) __network_data *DevMap;	// 
} netdat;

///////////////////////////////////////////////////////////////////////////////
// 
typedef struct __lsns_align( 16 ) __network_parameters{
	// constant memory
	float Step
	int MaxIons;
	int MaxChan;
	int MaxCells;
	gatepar Gates[LSNS_MAX_GPARS];
} netpar; 


///////////////////////////////////////////////////////////////////////////////
// interface 
extern netdat *lsns_alloc( netpar &par );			// allocate memory on host
extern void lsns_free( netdat *data );				// free memory both on host and device
extern bool lsns_map2dev( netdat *data, netpar &parpar );	// allocate memory on device and copy the network configuration from host to device.
extern bool lsns_run( floar step = -1. );			// run simulation

#endif /*__LSNS_ENGINE_H*/

