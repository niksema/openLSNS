#include "precompile.h"

#include "engine.h"
#include "debug.h"

#if !defined ( __CUDA__ )

#include <stdlib.h>
#include <string.h>

///////////////////////////////////////////////////////////////////////////////
// macroses which define the dimension the network data
#define MAX_GATES LSNS_MAX_GATE
#define MAX_GPARS LSNS_MAX_GPARS
#define MAX_CHAN MaxChan
#define MAX_IONS MaxIons
#define MAX_CELLS MaxCells
#define MAX_PARS MaxGlobalData
#define MAX_VIEW_PARS MaxViewPars

///////////////////////////////////////////////////////////////////////////////
// macroses to get access to networks data
#define pIonsType ( data->IonsType )
#define pIonsE ( data->IonsE )
#define pIonsI ( data->IonsI )
#define pIonsLUT ( data->IonsLUT )
#define pChanType ( data->ChanType )
#define pChanG ( data->ChanG )
#define pChanMH ( data->ChanMH )
#define pChanLUT ( data->ChanLUT )
#define pCellV ( data->CellV )
#define pChanGLUT ( data->ChanGLUT )
#define pIonsILUT ( data->IonsILUT )
#define pDevDat( index ) ( data->DevData[index] )
#define pHostDat ( data->HostData )
#define pGlobalDat ( data->GlobalData )
#define pGlobalLUT ( data->GlobalLUT )


///////////////////////////////////////////////////////////////////////////////
// common data for all networks elements (read-only)
static float Step = -1.f;
static float Threshold = -10.f;
static int Counter = 0;
static int MaxChan = 0;
static int MaxIons = 0;
static int MaxCells = 0;
static int MaxViewPars = 0;
static int MaxGlobalData = 0;
static gatepar Gates[LSNS_MAX_GPARS] = {0};
static pumppar Pumps[LSNS_MAX_PARPUMPS] = {0};

///////////////////////////////////////////////////////////////////////////////
// control_kernel: 
// Input parameters:
// Output parameters:
void control_kernel( int index /*, chandat *data*/ )
{
	if( index == 0 ){
		Counter++;
	}
}

///////////////////////////////////////////////////////////////////////////////
// ions_kernel: the engine kernel to calculate the properties of ions dynamics.
// such as pump current, concentration of ions inside the cell, etc.
// Input parameters:
// Output parameters:
void ions_kernel( int index, iondat *data )
{
	__lsns_assert( index>= 0 && index < MAX_IONS );			// DEBUG: check the range for 'index' variable
	// load parameters for simulation
	float step = Step;
	// load type of ions (Na, K, Ca etc, and type of ion dynamics)
	int4 tp = pIonsType[index];
	int type_dyn = _ions_typepump( tp );
	int type_eds = _ions_typeeds( tp );
	if( type_dyn != LSNS_NO_DYN ){
		// load ions parameters, resting potential, ions concentrations if needed
		pumppar par = Pumps[_ions_parpump( tp )];	// load pump properties
		int4 lut = pIonsLUT[index];			// load references to external parameters (channel conductances etc)
		float4 e = pIonsE[index];			// load ions parameters
		float4 i = pIonsI[index];			// load ions parameters
//--->>> possible CUDA optimization (try to use shared variables)
		float v = _cell_v( pCellV[_ions_lut_v( lut )] );// get membrane potential
//---<<< possible CUDA optimization
		float tau = _ions_tau( i );			// get time constant for ions dynamics
		float out = _ions_out( e );			// concentation of ions outside cell
		float rtfz = _ions_rtfz( e );			// rtfz for specific neurons
		// start calculations
		float gchan = 0; 				// total conductances of all ion currents which are involved to ions dynamics
		lsns_fastsum( _chan_g, pChanG, pChanGLUT+index*( MAX_CHAN_PER_PUMP/4+1 ), gchan, MAX_CHAN ); // sum all gchan
		// 2-order Runge-Kutta integration method
		float ipump = 0, apump = 0, ichan = 0;
		// sub-step 1
		float in = _ions_in( e );			// concentation of ions inside the cell
		float eds = _ions_eds( e );			// resting potential
		lsns_ipump( type_dyn, par, in, apump, ipump );	// calculate pump current
		lsns_ichan( apump, v, eds, gchan, ichan );	// calculate channels current
		// sub-step 2
		float in1 = in-step*( ichan+ipump )/( 2*tau );
		eds = ( type_eds != LSNS_NOSPEC_EDS )? lsns_eds( rtfz, in1, out ): eds; // calculate resting potential if needed
		lsns_ipump( type_dyn, par, in1, apump, ipump );	// calculate pump current
		lsns_ichan( apump, v, eds, gchan, ichan );	// calculate channels current
		// the final calculations for in, eds, ichan
		in = in-step*( ichan+ipump )/tau;
		eds = ( type_eds != LSNS_NOSPEC_EDS )? lsns_eds( rtfz, in, out ): eds; // calculate resting potential if needed
		lsns_ipump( type_dyn, par, in, apump, ipump );	// calculate pump current
		lsns_ichan( apump, v, eds, gchan, ichan );	// calculate channels current
		
		// store the results
		_ions_in( e ) = in;
		_ions_out( e ) = out;
		_ions_eds( e ) = eds;
		_ions_ipump( i ) = ipump;
		_ions_ichan( i ) = ichan;
		pIonsE[index] = e;
		pIonsI[index] = i;
	}
	else if( type_eds != LSNS_NOSPEC_EDS ){			// calculate reversal potential only (no ions dynamics)
		// load ions parameters
		float4 e = pIonsE[index];			// load ions parameters
		float in = _ions_in( e );			// concentation of ions inside the cell
		float out = _ions_out( e );			// concentation of ions outside cell
		float rtfz = _ions_rtfz( e );			// rtfz for specific ions
		// calculate resting potential
		_ions_eds( e ) = lsns_eds( rtfz, in, out );
		// store the results
		pIonsE[index] = e;
	}
}

///////////////////////////////////////////////////////////////////////////////
// connect_kernel: 
// Input parameters:
// Output parameters:
void connect_kernel( int index /*, chandat *data*/ )
{
}

///////////////////////////////////////////////////////////////////////////////
// chan_kernel: the engine kernel to calculate the properties of channel and 
// synaptic currents such as conductance, current etc.
// Input parameters:
// Output parameters:
void chan_kernel( int index, chandat *data )
{
	__lsns_assert( index>= 0 && index < MAX_CHAN );		// DEBUG: check the range for 'index' variable 
	float step = Step;
	// load type of ion channel (generic, a-b, etc) and its parameters (half-voltage, slope, etc)
	int4 tp = pChanType[index];
	// load references to external parameters (membrane potential, rest potential, etc)
	int4 sh = pChanLUT[index];
	__lsns_assert( _gate_typem( tp ) < MAX_GATES );		// DEBUG: check the range for _gate_typem( tp )
	__lsns_assert( _gate_typeh( tp ) < MAX_GATES );		// DEBUG: check the range for  _gate_typeh( tp ) 
	__lsns_assert( _gate_parm( tp ) < MAX_GPARS );		// DEBUG: check the range for  _gate_parm( tp ) 
	__lsns_assert( _gate_parh( tp ) < MAX_GPARS );		// DEBUG: check the range for  _gate_parh( tp ) 
	__lsns_assert( _chan_lut_v( sh ) < MAX_CELLS );		// DEBUG: check the range for _chan_lut_v( sh )
	__lsns_assert( _chan_lut_e( sh ) < MAX_IONS );		// DEBUG: check the range for _chan_lut_e( sh )
	__lsns_assert( _chan_lut_inm( sh ) < MAX_IONS );	// DEBUG: check the range for _chan_lut_inm( sh )
	__lsns_assert( _chan_lut_inh( sh ) < MAX_IONS );	// DEBUG: check the range for _chan_lut_inh( sh )
	// load properties of ions channel (conductance, current, etc)
	float4 g = pChanG[index];
	// load properties of gate variables (activation, inactivation, etc) if needed
	float4 mh = ( _gate_typem( tp )+_gate_typeh( tp ) != LSNS_NOGATE)? pChanMH[index]: float4();
//--->>> possible CUDA optimization (try to use shared variables)
	// load shared variables (resting potential, membrane potential, etc)
	float eds = _ions_eds( pIonsE[_chan_lut_e( sh )] );// extract resting potential from 'IonsE'
	float vm = _cell_v( pCellV[_chan_lut_v( sh )] );	// extract membrane potential from 'CellV'
	// load Mg- or Ca- concentration inside the cell for NMDA synapse or Z-channels from 'ions_e' if needed
	float in_m = ( _gate_typem( tp ) >= LSNS_ZGENERIC_INSTANT )? _ions_in( pIonsE[_chan_lut_inm( sh )] ):0;
	// load Mg- or Ca- concentration inside the cell for NMDA synapse or Z-channels from 'ions_e' if needed
	float in_h = ( _gate_typeh( tp ) >= LSNS_ZGENERIC_INSTANT )? _ions_in( pIonsE[_chan_lut_inh( sh )] ):0;
//---<<< possible CUDA optimization
	// perform calculations
	float mp, hp;
	proc_gate( _gate_typem( tp ), Gates[_gate_parm( tp )], in_m, vm, step, _gate_powm( mh ), _gate_m( mh ), mp );
	proc_gate( _gate_typeh( tp ), Gates[_gate_parh( tp )], in_h, vm, step, _gate_powh( mh ), _gate_h( mh ), hp );
	_chan_g( g ) = _chan_gmax( g )*mp*hp;			// g
	_chan_ge( g ) = _chan_g( g )*eds;			// ge
	_chan_i( g ) = _chan_g( g )*( vm-eds );			// I
	// save properties of ions channel (conductance, current, etc)
	pChanG[index] = g;
	if( _gate_typem( tp )+_gate_typeh( tp ) != LSNS_NOGATE ){
		// save properties of gate variables (activation, inactivation, etc) if needed
		pChanMH[index] = mh;
	}
}

///////////////////////////////////////////////////////////////////////////////
// cell_kernel: the engine kernel to calculate the properties of Hodgkin-Huxley
// cell such as membrane potential, onset of the spike etc.
// Input parameters:
// Output parameters:
void cell_kernel( int index, celldat *data )
{
	__lsns_assert( index>= 0 && index < MAX_CELLS );	// DEBUG: check the range for 'index' variable 
	// load properties of the cell (membrane potential, etc)
	float step = Step;
	float threshold = Threshold;
	float4 v = pCellV[index];
	// preprocessing
	float g = 1, ge = 0, ipump = 0;
	lsns_fastsum( _ions_ipump, pIonsI, pIonsILUT+index*( MAX_IPUMP_PER_CELL/4+1 ), ipump, MAX_IONS ); // sum all ipump
	lsns_fastsum2( _chan_g, _chan_ge, pChanG, pChanGLUT+index*( MAX_CHAN_PER_CELL/4+1 ), g, ge, MAX_CHAN ); // sum all g and ge
	__lsns_assert( g != 0.0 );				// DEBUG: check if conductance is not zero
	// perform calculations
	float time = _cell_c( v )/g;
	float v_inf = ( ge+ipump+_cell_iadd( v ))/g;
	float vm = lsns_exp_euler( _cell_v( v ), v_inf, step, time );
	float spike = ( vm > threshold && _cell_v( v )  <= threshold )? 1.f: 0.f;
	// store the results of simulation
	_cell_v( v ) = vm;
	_cell_spike( v ) = spike;
	pCellV[index] = v;
}

///////////////////////////////////////////////////////////////////////////////
// store_kernel: 
// Input parameters:
// Output parameters:
///////////////////////////////////////////////////////////////////////////////
// store2dev_kernel: Stores the preliminary saved results to device memory.
// The results of the simulation store to device memory into 
// arrays float4 *DevData[MAX_STORED_STEPS]. These arrays organized in 
// such a way that guaranteed continues memory allocation to rich the maximal 
// performance while I/O operations: 
//	DevData[0] = alloc(size*MAX_STORED_STEPS);
//	for( i = 1; i < MAX_STORED_STEPS; ++i ){
//		DevData[i] = DevData[0]+i*size;
//	}
// Input parameters:
// Output parameters:
// MAX_VIEW_PARS = ( MAX_IIONS_VIEWS+MAX_EIONS_VIEWS+MAX_GCHAN_VIEWS+MAX_MHCHAN_VIEWS+MAX_CELLV_VIEWS )/4;
void store2dev_kernel( int index, iodat *data )
{
	__lsns_assert( index>= 0 && index < MAX_VIEW_PARS );	// DEBUG: check the range for 'index' variable
	float4 res = {0};
	float4 *src = pGlobalDat;
	int4 lut = pGlobalLUT[index];
	store_data( src, lut, res, MAX_PARS );
	pDevDat(Counter%MAX_STORED_STEPS)[index] = res;
}

///////////////////////////////////////////////////////////////////////////////
// store2host_kernel: Stores the preliminary saved results to host memory.
// Initially, the results of the simulation store to device memory into 
// arrays float4 *DevData[MAX_STORED_STEPS]. These arrays organized in 
// such a way that guaranteed continues memory allocation to rich the maximal 
// performance while I/O operations: 
//	DevData[0] = alloc(size*MAX_STORED_STEPS);
//	for( i = 1; i < MAX_STORED_STEPS; ++i ){
//		DevData[i] = DevData[0]+i*size;
//	}
// Input parameters:
// Output parameters:
void store2host_kernel( int index , iodat *data )
{
	__lsns_assert( index>= 0 && index < MAX_VIEW_PARS*MAX_STORED_STEPS );// DEBUG: check the range for 'index' variable
	pHostDat[index] = pDevDat(0)[index];			// copy data from device memory to host memory
}

///////////////////////////////////////////////////////////////////////////////
// interface 
 
// allocate memory on host
netdat *lsns_alloc( netpar &par )
{
	// setup the constant parameters of the network
	netdat *data = ( netdat *)malloc( sizeof( netdat )); __lsns_assert( data != NULL );
	data->Ions.MaxIons = par.MaxIons; 
	data->Channels.MaxChan = par.MaxChan;
	data->Cells.MaxCells = par.MaxCells;
	Step = par.Step; Threshold = par.Threshold; MaxIons = par.MaxIons; MaxChan = par.MaxChan; MaxCells = par.MaxCells;
	MaxGlobalData = par.MaxGlobalData = 2*MaxIons+2*MaxChan+MaxCells;
	memcpy ( Gates, par.Gates, sizeof( netpar )*LSNS_MAX_GPARS );
	memcpy ( Pumps, par.Pumps, sizeof( pumppar )*LSNS_MAX_PARPUMPS );
	// allocate memory for actual parameters
	data->IOData.GlobalData = data->GlobalData = ( float4 *)malloc( sizeof( float4 )*( MaxGlobalData ));
	// initialize the arrays of parameters
	data->Ions.IonsI = data->Cells.IonsI = data->GlobalData;
	data->Ions.IonsE = data->Channels.IonsE = data->GlobalData+MaxIons;
	data->Channels.ChanG = data->Cells.ChanG = data->Ions.ChanG = data->GlobalData+2*MaxIons;
	data->Channels.ChanMH = data->GlobalData+2*MaxIons+MaxChan;
	data->Cells.CellV = data->Channels.CellV = data->Ions.CellV = data->GlobalData+2*MaxIons+2*MaxChan;
	// allocate memory for look-up-tables
	data->Ions.IonsType = ( int4 * )malloc( sizeof( int4 )*MaxIons ); 
	data->Ions.IonsLUT = ( int4 * )malloc( sizeof( int4 )*MaxIons );
	data->Ions.ChanGLUT = ( int4 * )malloc( sizeof( int4 )*MaxIons*( MAX_CHAN_PER_PUMP/4+1 ) );

	data->Channels.ChanType = ( int4 * )malloc( sizeof( int4 )*MaxChan );
	data->Channels.ChanLUT = ( int4 * )malloc( sizeof( int4 )*MaxChan );

	data->Cells.ChanGLUT = ( int4 *)malloc( sizeof( int4 )*MaxCells*( MAX_CHAN_PER_CELL/4+1 ) ); 
	data->Cells.IonsILUT = ( int4 *)malloc( sizeof( int4 )*MaxCells*( MAX_IPUMP_PER_CELL/4+1 ) );
	// device-specific memory (for non-cuda version this is a pointer to the same data structure)
	data->DevMap = data;
	return data;
}

// free memory both on host and device
void lsns_free( netdat *data )
{
	// free memory for actual parameters
	free( data->GlobalData ); 
	data->GlobalData = data->IOData.GlobalData = NULL;
	data->Ions.IonsI = data->Cells.IonsI = NULL;
	data->Ions.IonsE = data->Channels.IonsE = NULL;
	data->Channels.ChanG = data->Cells.ChanG = data->Ions.ChanG = NULL;
	data->Channels.ChanMH = NULL;
	data->Cells.CellV = data->Channels.CellV = data->Ions.CellV = NULL;
	// free memory for look-up-tables
	free( data->Ions.IonsType );
	free( data->Ions.IonsLUT );
	free( data->Ions.ChanGLUT );
	free( data->Channels.ChanType );
	free( data->Channels.ChanLUT );
	free( data->Cells.ChanGLUT );
	free( data->Cells.IonsILUT );
	// device-specific memory (for non-cuda version this is a pointer to the same data structure)
	data->DevMap = NULL;
	free( data );
	// clean up the constant parameters of the network
	Step = -1; MaxIons = 0; MaxChan = 0; MaxCells = 0;
	memset( Gates, 0, sizeof( netpar )*LSNS_MAX_GPARS );
	memset( Pumps, 0, sizeof( pumppar )*LSNS_MAX_PARPUMPS );
}

// allocate memory on device and copy the network configuration from host to device.
bool lsns_map2dev( netdat *data, netpar &parpar )
{
	// device-specific memory (for non-cuda version this is a pointer to the same data structure)
	data->DevMap = data;
	return true;
}

// run simulation
bool lsns_run( float step )
{
	return true;
}

#endif

