#ifndef __LSNS_ENGINE_H
#define __LSNS_ENGINE_H

#include "config.h"
#include "gateproc.h"
#include "ionproc.h"
#include "viewproc.h"

#define MAX_CHAN_PER_PUMP	7				/*<24*/
#define MAX_CHAN_PER_CELL	15				/*<24*/
#define MAX_IPUMP_PER_CELL	3				/*<24*/
#define MAX_STORED_STEPS	16
//=================== syndat macroses =========================================
///////////////////////////////////////////////////////////////////////////////
// type of synapes (decodes SynType):
//	x - type of synapse;
//	y - parameters
//	z - size
//	w - initial index in Wall&CellLUT arrays
#define _syn_type( tp ) ( tp ).x
#define _syn_par( tp ) ( tp ).y
#define _syn_size( tp ) ( tp ).z
#define _syn_lut( tp ) ( tp ).w
///////////////////////////////////////////////////////////////////////////////
// parameters of the synaptic weight (decodes Wsyn):
//	x - total weight of all connections;
//	y - production: ( rate of transmitter release )*( plasticity )
//	z - 1 for pulse model or step/T other models;
//	w -  exp( step/T ) for pulse model or 1-step/T for othe models
#define _wsyn_total( wsyn ) ( wsyn ).x
#define _wsyn_alpha( wsyn ) ( wsyn ).y
#define _wsyn_dt( wsyn ) ( wsyn ).z
#define _wsyn_edt( wsyn ) ( wsyn ).w
///////////////////////////////////////////////////////////////////////////////
// syndat maps data which are related to synaptic properties onto global memory
typedef struct __lsns_align( 16 ) __synapses_data{
	// look-up-tables for shared variables (read-only)
	int4 __lsns_align( 16 ) *SynType;			//  x - type of synapse, y - parameters, z - size, w - initial index in Wall&CellLUT arrays
	int4 __lsns_align( 16 ) *CellLUT;			// look-up-table for all neurons which are converged onto particular synapse
	float4 __lsns_align( 16 ) *Wall;			// all weights for particular synapse
	// local variables (read/write)
	float4 __lsns_align( 16 ) *Wsyn;			// x - total sum, y - ( rate of transmitter release )*( plasticity ), z - 1 for pulse model or step/T other models, w - exp( step/T ) for pulse model or 1-step/T for othe models
	// shared variables
	float4 __lsns_align( 16 ) *CellV;			// cell properties: x - membrane potential, y - membrane capacitance, z - spike onset, w - injected current
} syndat;
//=================== iondat macroses =========================================
///////////////////////////////////////////////////////////////////////////////
// type of ion dynamics (decodes IonsType):
//	x - type of ions pump (Na-pump, Ca-pump, Ca-pump_1 etc);
//	y - type of reversal potential (non-specific, specific etc);
//	z - parameters of ions pump;
//	w - reserved. 
#define _ions_typepump( tp ) ( tp ).x 
#define _ions_typeeds( tp ) ( tp ).y
#define _ions_parpump( tp ) ( tp ).z
///////////////////////////////////////////////////////////////////////////////
// shared parameters for ion dynamics (decodes IonsLUT):
//	x - membrane potential;
//	y - reserved;
//	z - reserved;
//	w - reserved
#define _ions_lut_v( sh ) ( sh ).x
///////////////////////////////////////////////////////////////////////////////
// parameters of ions dynamics (decodes IonsE):
//	x - reversal potential (Eds);
//	y - concentration of ions inside the cell;
//	z - concentration of ions outside the cell;
//	w - RT/Fz constant for specific ions
#define _ions_eds( e ) ( e ).x
#define _ions_in( e ) ( e ).y
#define _ions_out( e ) ( e ).z
#define _ions_rtfz( e ) ( e ).w
///////////////////////////////////////////////////////////////////////////////
// parameters of ion current (decodes IonsI):
//	x - pump current;
//	y - channels current;
//	z - reserved;
//	w - reserved.
#define _ions_ipump( i ) ( i ).x
#define _ions_ichan( i ) ( i ).y
///////////////////////////////////////////////////////////////////////////////
// iondat maps data which are related to ions dynamics onto global memory
typedef struct __lsns_align( 16 ) __ions_data{
	// look-up-tables for shared variables (read-only)
	int4 __lsns_align( 16 ) *IonsType;			// type of ions: x - pump type, y - eds type, z - pump parameters, w - reserved
	int4 __lsns_align( 16 ) *IonsLUT;			// indices of shared variables: x - cell properties, y, z, w - reserved
	int4 __lsns_align( 16 ) *ChanGLUT;			// look-up-table of channel current: x - counter, the rest are actual indices of ChanG array
	// local variables (read/write)
	float4 __lsns_align( 16 ) *IonsE;			// ions properties: x - reversal potential (Eds), y - concentration of ions inside the cell, z - concentration of ions outside the cell, w - RT/Fz constant for specific ions
	float4 __lsns_align( 16 ) *IonsI;			// ion currents: x - pump current, y - channels current, z - time constant of ions dynamics, w - reserved.
	// shared variables
	float4 __lsns_align( 16 ) *ChanG;			// channel current: x - maximal conductance, y - conductance, z - current, w - G*Eds production
	float4 __lsns_align( 16 ) *CellV;			// cell properties: x - membrane potential, y - membrane capacitance, z - spike onset, w - injected current
} iondat;

//=================== chandat macroses ========================================
///////////////////////////////////////////////////////////////////////////////
// gate types of ion channel (decodes ChanType):
//	x - type of activation;
//	y - parameters of activation;
//	z - type of inactivation;
//	w - parameters of inactivation;
#define _gate_typem( t ) ( t ).x
#define _gate_parm( t ) ( t ).y
#define _gate_typeh( t ) ( t ).z
#define _gate_parh( t ) ( t ).w
///////////////////////////////////////////////////////////////////////////////
// indices of shared parameters for ion channels (decodes ChanLUT):
//	x - membrane potential;
//	y - resting potential;
//	z - concentration of ions inside the cell for activation for z-channels/or total weight of connections for the synapse (for m component);
//	w - concentration of ions inside the cell for inactivation for z-channels/or total weight of connections for the synapse (for h component)
#define _chan_lutv( lut ) ( lut ).x
#define _chan_lute( lut ) ( lut ).y
#define _chan_lutm( lut ) ( lut ).z
#define _chan_luth( lut ) ( lut ).w
///////////////////////////////////////////////////////////////////////////////
// parameters of gate variables of ion channel (decodes ChanMH):
//	x - activation;
//	y - inactivation;
//	z - power of activation for channels;
//	w - power of inactivation for channels
#define _gate_m( mh ) ( mh ).x
#define _gate_h( mh ) ( mh ).y
#define _gate_powm( mh ) ( mh ).z
#define _gate_powh( mh ) ( mh ).w
///////////////////////////////////////////////////////////////////////////////
// parameters of the ion channel (decodes ChanG):
//	x - maximal conductance;
//	y - conductance;
//	z - current;
//	w - G*Eds production
#define _chan_gmax( g ) ( g ).x
#define _chan_g( g ) ( g ).y
#define _chan_ge( g ) ( g ).z
#define _chan_i( g ) ( g ).w
///////////////////////////////////////////////////////////////////////////////
// chandat maps data which are related to ion channels onto global memory
typedef struct __lsns_align( 16 ) __channel_data{
	// local variables (read-only)
	int4 __lsns_align( 16 ) *ChanType;			// type of channel: x - type of activation, y - parameters of activation, z - type of inactivation, w - parameters of inactivation.
	int4 __lsns_align( 16 ) *ChanLUT;			// indices of shared variables: x - CellV, y - IonsE/eds, z - IonsE/in for M (z-channel)/or W total for synapse, w - IonsE/in for H (z-channel)/or W total for synapse
	// local variables (read-write)
	float4 __lsns_align( 16 ) *ChanMH;			// gate variables: x - activation, y - inactivation, z - power of activation, w - power of inactivation
	float4 __lsns_align( 16 ) *ChanG;			// channel current: x - maximal conductance, y - conductance, z - current, w - G*Eds production
	// shared variables
	float4 __lsns_align( 16 ) *CellV;			// cell properties: x - membrane potential, y - membrane capacitance, z - spike onset, w - injected current
	float4 __lsns_align( 16 ) *IonsE;			// ions properties: x - reversal potential (Eds), y - concentration of ions inside the cell, z - concentration of ions outside the cell, w - RT/Fz constant for specific ions
	float4 __lsns_align( 16 ) *Wsyn;			// synaptic weights: x - total sum, y - ( rate of transmitter release )*( plasticity ), z - 1 for pulse model or step/T other models, w - exp( step/T ) for pulse model or 1-step/T for othe models
} chandat;

//=================== celldat macroses ========================================
///////////////////////////////////////////////////////////////////////////////
// parameters of cell (decodes CellV):
//	x - membrane potential; 
//	y - membrane capacitance;
//	z - spike onset;
//	w - injected current
#define _cell_v( v ) ( v ).x
#define _cell_c( v ) ( v ).y
#define _cell_spike( v ) ( v ).z
#define _cell_iadd( v ) ( v ).w
///////////////////////////////////////////////////////////////////////////////
// celldat maps data which are related to neurons onto global memory
typedef struct __lsns_align( 16 ) __cell_data{
	// local variables (read-only)
	int4 __lsns_align( 16 ) *ChanGLUT;			// look-up-table of channel current: x - counter, the rest are actual indices of ChanG array
	int4 __lsns_align( 16 ) *IonsILUT;			// look-up-table of pump current: x - counter, the rest are actual indices of IonsI array;
	// local variables (read-write)
	float4 __lsns_align( 16 ) *CellV;			// cell properties: x - membrane potential, y - membrane capacitance, z - spike onset, w - injected current
	// shared variables
	float4 __lsns_align( 16 ) *ChanG;			// channel current: x - maximal conductance, y - conductance, z - current, w - G*Eds production
	float4 __lsns_align( 16 ) *IonsI;			// ion currents: x - pump current, y - channels current, z - time constant of ions dynamics, w - reserved.
} celldat;

//=============================================================================
///////////////////////////////////////////////////////////////////////////////
// iodat contains the pointers to both host and device memory which are related 
// to data to be displayed
typedef struct __lsns_align( 16 ) __iobuf{
	// local variables (read-only). 
	// LUT format: bits 31..30 are coding the offset in each float4 variable (00 - x, 01 - y, 10 - z, 11 - w); 
	// bits 29..0 are coding the offset in the global array 'GlobalData'
	int4 __lsns_align( 16 ) *ViewLUT;			// look-up-table for data needed to be stored (located both in device and host memory)
	// local variables (read-write)
	float4 __lsns_align( 16 ) *ViewData[MAX_STORED_STEPS];	// data to display (located in device memory)
	float4 __lsns_align( 16 ) *GlobalViewData;		// data to display (pinned memory located in host memory)
	// shared variables
	float4 __lsns_align( 16 ) *GlobalData;			// array of global data (located in device memory)
} iodat;

//=============================================================================
///////////////////////////////////////////////////////////////////////////////
// netdat contains the pointers to global memory which are related to network
// parameters and structure
typedef struct __lsns_align( 16 ) __network_data{
	// global memory
	int4 __lsns_align( 16 ) *GlobalLUT;			// look-up-tables (located in both device and host memory); size = MaxGlobalLUT
	float4 __lsns_align( 16 ) *GlobalData;			// array of global data (located in both device and host memory); size = MaxGlobalData
	int4 __lsns_align( 16 ) *GlobalViewLUT;			// look-up-table for data needed to be stored (located both in device and host memory)
	float4 __lsns_align( 16 ) *GlobalViewData;		// array of data to be displayed (pinned memory located in host)
	// network elements (channels, ions, cells, etc) mapped on global memory
	iondat Ions;						// ions data
	chandat Channels;					// channels data
	celldat Cells;						// cells (compartments) data
	iodat IOData;						// data buffer to display the results of simulations
	// global memory on device
	struct __lsns_align( 16 ) __network_data *DevMap;	// pointer to device-specific memory (for non-cuda version this is a pointer to the same data structure)
} netdat;

//=============================================================================
///////////////////////////////////////////////////////////////////////////////
// netpar structure
typedef struct __lsns_align( 16 ) __network_parameters{
	// data which have to be initialized before memory allocation by 'lsns_alloc' method
	float Threshold;					// threshold for spike discrimination
	float Step;						// integration step
	int MaxIons;						// total number of ions parameters
	int MaxChan;						// total number of ion channels (includeing synapses)
	int MaxCells;						// total numbed of parameters to be displayed
	int MaxViewPars;					// total number of network parameters
	gatepar Gates[LSNS_MAX_GATEPARS];			// properties of ion channels
	ionspar Ions[LSNS_MAX_IONPARS];				// properties of ion dynamics
	// data which will initialized after memory allocation by 'lsns_alloc' method
	int MaxGlobalData;					// = 2*MaxIons+2*MaxChan+MaxCells
	int MaxGlobalLUT;					// = MaxIons*( 3+MAX_CHAN_PER_PUMP/4 )+2*MaxChan+MaxCells*(2+MAX_CHAN_PER_CELL/4+MAX_IPUMP_PER_CELL/4)
	netdat *Data;						// pointer to network data
} netpar; 

//=============================================================================
///////////////////////////////////////////////////////////////////////////////
// interface methods
extern netdat *lsns_alloc( netpar &par );			// allocate memory on host
extern void lsns_free( netdat *data );				// free memory both on host and device
extern bool lsns_map2dev( netdat *data, netpar &parpar );	// allocate memory on device and copy the network configuration from host to device.
extern bool lsns_run( float step = -1. );			// run simulation

#endif /*__LSNS_ENGINE_H*/

