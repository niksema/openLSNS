#include "precompile.h"

#include "engine.h"
#include "debug.h"

#if !defined ( __CUDA__ )

#include <stdlib.h>
#include <string.h>

///////////////////////////////////////////////////////////////////////////////
// macroses which define the dimension the network data
#define MAX_WSYNS 1	/*todo: implement synaptic summation*/
#define MAX_GATES LSNS_MAX_GATES
#define MAX_GPARS LSNS_MAX_GATEPARS
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
#define pWsyn ( data->Wsyn )
#define pChanGLUT ( data->ChanGLUT )
#define pIonsILUT ( data->IonsILUT )
#define pDevDat( index ) ( data->ViewData[index] )
#define pHostDat ( data->GlobalViewData )
#define pGlobalDat ( data->GlobalData )
#define pGlobalViewLUT ( data->ViewLUT )
///////////////////////////////////////////////////////////////////////////////
// common data for all networks elements (read-only)
static float Threshold = -10.f;
static float Step = -1.f;
static int StepCounter = 0;
static int MaxIons = 0;
static int MaxChan = 0;
static int MaxCells = 0;
static int MaxViewPars = 0;
static int MaxGlobalData = 0;
static gatepar Gates[LSNS_MAX_GATEPARS] = {0};
static ionspar Ions[LSNS_MAX_IONPARS] = {0};

///////////////////////////////////////////////////////////////////////////////
// control_kernel: 
// Input parameters:
// Output parameters:
void control_kernel( int index /*, ctrldat *data*/ )
{
	if( index == 0 ){
		++StepCounter;
	}
}

///////////////////////////////////////////////////////////////////////////////
// connect_kernel: 
// Input parameters:
// Output parameters:
void connect_kernel( int index /*, connectdat *data*/ )
{
	__lsns_assert( index >= 0 && index < MAX_WSYNS );	// DEBUG: check the range for 'index' variable
	/*
	wtab - synapse type, parameters, index for wsyn&wlut, index for placticity
	wsyn is float4 array (x - total sum, rest)
	wlut is look-up-table (x - size, rest indecis for cellv array)
	*/
}

///////////////////////////////////////////////////////////////////////////////
// +'ions_kernel': the engine kernel to calculate the properties of ions dynamics.
// such as pump current, concentration of ions inside the cell, etc.
// Input parameters:
// Output parameters:
void ions_kernel( int index, iondat *data )
{
	__lsns_assert( index >= 0 && index < MAX_IONS );	// DEBUG: check the range for 'index' variable
	// load parameters for simulation
	float step = Step;
	// load type of ions (Na, K, Ca etc, and type of ion dynamics)
	int4 tp = pIonsType[index];
	int type_dyn = _ions_typepump( tp );
	int type_eds = _ions_typeeds( tp );
	if( type_dyn != LSNS_NO_DYN ){
		// load ions parameters, resting potential, ions concentrations if needed
		ionspar par = Ions[_ions_parpump( tp )];	// load pump properties
		int4 lut = pIonsLUT[index];			// load references to external parameters (channel conductances etc)
		float4 e = pIonsE[index];			// load ions parameters
		float4 i = pIonsI[index];			// load ions parameters
		float4 v = pCellV[_ions_lut_v( lut )];		// load cell parameters
//--->>> possible CUDA optimization (try to use shared variables)
		float vm = _cell_v( v );			// get membrane potential
//---<<< possible CUDA optimization
		float tau = _pump_tau( par );			// get time constant for ions dynamics
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
		lsns_ipump( type_dyn, par, in, apump, ipump, KCA );// calculate pump current
		lsns_ichan( apump, vm, eds, gchan, ichan );	// calculate channels current
		// sub-step 2
		float in1 = in-step*( ichan+ipump )/( 2*tau );
		eds = ( type_eds != LSNS_NOSPEC_EDS )? lsns_eds( rtfz, in1, out ): eds; // calculate resting potential if needed
		lsns_ipump( type_dyn, par, in1, apump, ipump, KCA );// calculate pump current
		lsns_ichan( apump, vm, eds, gchan, ichan );	// calculate channels current
		// the final calculations for in, eds, ichan
		in = in-step*( ichan+ipump )/tau;
		eds = ( type_eds != LSNS_NOSPEC_EDS )? lsns_eds( rtfz, in, out ): eds; // calculate resting potential if needed
		lsns_ipump( type_dyn, par, in, apump, ipump, KCA );// calculate pump current
		lsns_ichan( apump, vm, eds, gchan, ichan );	// calculate channels current
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
// +'chan_kernel': the engine kernel to calculate the properties of channel and 
// synaptic currents such as conductance, current etc.
void chan_kernel( int index, chandat *data )
{
	__lsns_assert( index >= 0 && index < MAX_CHAN );				// DEBUG: check the range for 'index' variable 
	float step = Step;
	// load type of ion channel (generic, a-b, etc) and its parameters (half-voltage, slope, etc)
	int4 tp = pChanType[index];
	// load references to external parameters (membrane potential, rest potential, etc)
	int4 sh = pChanLUT[index];
	__lsns_assert( _gate_typem( tp ) >= 0 && _gate_typem( tp ) < MAX_GATES );	// DEBUG: check the range for _gate_typem( tp )
	__lsns_assert( _gate_typeh( tp ) >= 0 && _gate_typeh( tp ) < MAX_GATES );	// DEBUG: check the range for  _gate_typeh( tp ) 
	__lsns_assert( _gate_parm( tp ) >= 0 && _gate_parm( tp ) < MAX_GPARS );		// DEBUG: check the range for  _gate_parm( tp ) 
	__lsns_assert( _gate_parh( tp ) >= 0 && _gate_parh( tp ) < MAX_GPARS );		// DEBUG: check the range for  _gate_parh( tp ) 
	__lsns_assert( _chan_lutv( sh ) >= 0 &&_chan_lutv( sh ) < MAX_CELLS );		// DEBUG: check the range for _chan_lut_v( sh )
	__lsns_assert( _chan_lute( sh ) >= 0 && _chan_lute( sh ) < MAX_IONS );		// DEBUG: check the range for _chan_lut_e( sh )
	__lsns_assert( _gate_typem( tp ) <= LSNS_PS_NMDA && _chan_lutm( sh ) >= 0 && _chan_lutm( sh ) < MAX_IONS );	// DEBUG: check the range for _chan_lut_m( sh )
	__lsns_assert( _gate_typeh( tp ) <= LSNS_PS_NMDA && _chan_luth( sh ) >= 0 && _chan_luth( sh ) < MAX_IONS );	// DEBUG: check the range for _chan_lut_h( sh )
	__lsns_assert( _gate_typem( tp ) > LSNS_PS_NMDA && _chan_lutm( sh ) >= 0 && _chan_lutm( sh ) < MAX_WSYNS );	// DEBUG: check the range for _chan_lut_m( sh )
	__lsns_assert( _gate_typeh( tp ) > LSNS_PS_NMDA && _chan_luth( sh ) >= 0 && _chan_luth( sh ) < MAX_WSYNS );	// DEBUG: check the range for _chan_lut_h( sh )
	// load properties of ions channel (conductance, current, etc)
	float4 g = pChanG[index];
	// load properties of gate variables (activation, inactivation, etc) if needed
	float4 mh = ( _gate_typem( tp )+_gate_typeh( tp ) != LSNS_NOGATE)? pChanMH[index]: float4();
//todo: { possible CUDA optimization (try to use shared variables)
	// load shared variables (resting potential, membrane potential, etc) to shared memory
	float eds = _ions_eds( pIonsE[_chan_lute( sh )] );				// extract resting potential from 'IonsE'
	float vm = _cell_v( pCellV[_chan_lutv( sh )] );					// extract membrane potential from 'CellV'
//todo: } possible CUDA optimization
	// perform calculations
	float mp, hp;
	proc_gate( _gate_typem( tp ), Gates[_gate_parm( tp )], step, vm, _chan_lutm( sh ), pIonsE, pWsyn, _gate_modm( mh ), _gate_m( mh ), mp );
	proc_gate( _gate_typeh( tp ), Gates[_gate_parh( tp )], step, vm, _chan_luth( sh ), pIonsE, pWsyn, _gate_modh( mh ), _gate_h( mh ), hp );
	_chan_g( g ) = _chan_gmax( g )*mp*hp;						// g
	_chan_ge( g ) = _chan_g( g )*eds;						// ge
	_chan_i( g ) = _chan_g( g )*( vm-eds );						// I
	// save properties of ions channel (conductance, current, etc)
	pChanG[index] = g;
	// save properties of gate variables (activation, inactivation, etc) if needed
	if( _gate_typem( tp )+_gate_typeh( tp ) != LSNS_NOGATE ){
		pChanMH[index] = mh;
	}
}

///////////////////////////////////////////////////////////////////////////////
// +'cell_kernel': the engine kernel to calculate the properties of Hodgkin-Huxley
// cell such as membrane potential, onset of the spike etc.
void cell_kernel( int index, celldat *data )
{
	__lsns_assert( index >= 0 && index < MAX_CELLS );				// DEBUG: check the range for 'index' variable 
	// load properties of the cell (membrane potential, etc)
	float step = Step;
	float threshold = Threshold;
	float4 v = pCellV[index];
	// preprocessing
	float g_sum = 0.f, ge_sum = 0.f, ipump = 0.f;
	lsns_fastsum( _ions_ipump, pIonsI, pIonsILUT+index*( MAX_IPUMP_PER_CELL/4+1 ), ipump, MAX_IONS );		// sum all ipump
	lsns_fastsum2( _chan_g, _chan_ge, pChanG, pChanGLUT+index*( MAX_CHAN_PER_CELL/4+1 ), g_sum, ge_sum, MAX_CHAN );	// sum all g and ge
	__lsns_assert( g_sum != 0.0 );							// DEBUG: check if total conductance is not zero
	// perform calculations
	float time = _cell_c( v )/g_sum;
	float v_inf = ( ge_sum+ipump+_cell_iadd( v ))/g_sum;
	float vm = lsns_exp_euler( _cell_v( v ), v_inf, step, time );
	float spike = ( vm > threshold && _cell_v( v )  <= threshold )? 1.f: 0.f;
	// store the results of simulation
	_cell_v( v ) = vm;
	_cell_spike( v ) = spike;
	pCellV[index] = v;
}

///////////////////////////////////////////////////////////////////////////////
// +'store2dev_kernel': stores the simulation results to device memory.
// The results of the simulation store to device memory into 
// arrays 'float4 *ViewData[MAX_STORED_STEPS]'. 
void store2dev_kernel( int index, iodat *data, int counter )
{
	__lsns_assert( index >= 0 && index < MAX_VIEW_PARS/4+1 );			// DEBUG: check the range for 'index' variable
	__lsns_assert( counter >= 0 && counter < MAX_STORED_STEPS );			// DEBUG: check the range for 'counter' variable
	int4 lut = pGlobalViewLUT[index];
	float4 *src = pGlobalDat;
	float4 res = get_data( src, lut, MAX_PARS );
	pDevDat( counter )[index] = res;
}

///////////////////////////////////////////////////////////////////////////////
// +'store2host_kernel': copies the preliminary saved results to host memory.
// Initially, the simulation results for 'MAX_STORED_STEPS' steps are stored 
// into device memory 
void store2host_kernel( int index , iodat *data )
{
	__lsns_assert( index >= 0 && index < ( MAX_VIEW_PARS/4+1 )*MAX_STORED_STEPS );	// DEBUG: check the range for 'index' variable
	pHostDat[index] = pDevDat(0)[index];						// copy data from device memory to host (pinned) memory
}

///////////////////////////////////////////////////////////////////////////////
// interface 
//=================== global methods ==========================================
///////////////////////////////////////////////////////////////////////////////
//+'lsns_alloc' allocates memory on host 
netdat *lsns_alloc( netpar &par )
{
	// allocate memory for 'netdat' structure
	par.Data = ( netdat *)malloc( sizeof( netdat )); __lsns_assert( par.Data != NULL );
	// calculate the dimension of arrays of actual parameters and look-up-tables
	par.MaxGlobalData = 2*par.MaxIons+2*par.MaxChan+par.MaxCells;
	par.MaxGlobalLUT = ( 3+MAX_CHAN_PER_PUMP/4 )*par.MaxIons+2*par.MaxChan+(2+MAX_CHAN_PER_CELL/4+MAX_IPUMP_PER_CELL/4)*par.MaxCells;
	// allocate memory for actual parameters and look-up-tables
	par.Data->GlobalLUT = ( int4 * )malloc( sizeof( int4 )*par.MaxGlobalLUT ); __lsns_assert( par.Data->GlobalLUT != NULL );
	par.Data->GlobalData = ( float4 *)malloc( sizeof( float4 )*par.MaxGlobalData ); __lsns_assert( par.Data->GlobalData != NULL );
	par.Data->GlobalViewLUT = ( int4 * )malloc( sizeof( int4 )*( par.MaxViewPars/4+1 )); __lsns_assert( par.Data->GlobalViewLUT != NULL );
	par.Data->GlobalViewData = ( float4 * )malloc( sizeof( float4 )*( par.MaxViewPars/4+1 )*MAX_STORED_STEPS ); __lsns_assert( par.Data->GlobalViewData != NULL ); // pinned-memory for cuda version
	// map the arrays of specific parameters onto the global memory
	par.Data->Cells.IonsI			= 
	par.Data->Ions.IonsI			= par.Data->GlobalData;
	par.Data->Channels.IonsE		= 
	par.Data->Ions.IonsE			= par.Data->Ions.IonsI+par.MaxIons;
	par.Data->Cells.ChanG			=
	par.Data->Ions.ChanG			=
	par.Data->Channels.ChanG		= par.Data->Ions.IonsE+par.MaxIons;
	par.Data->Channels.ChanMH		= par.Data->Channels.ChanG+par.MaxChan;
	par.Data->Cells.CellV			=
	par.Data->Ions.CellV			=
	par.Data->Channels.CellV		= par.Data->Channels.ChanMH+par.MaxChan;
	par.Data->IOData.GlobalData		= par.Data->GlobalData;
	par.Data->IOData.GlobalViewData		= par.Data->GlobalViewData;	// pinned-memory for cuda version
	for( int i = 0; i < MAX_STORED_STEPS; ++i ){
		par.Data->IOData.ViewData[i]	= par.Data->IOData.GlobalViewData+i*(par.MaxViewPars/4+1);
	}
	// map the arrays of specific look-up-tables onto the global memory
	par.Data->Ions.IonsType			= par.Data->GlobalLUT;
	par.Data->Ions.IonsLUT			= par.Data->Ions.IonsType+par.MaxIons;
	par.Data->Ions.ChanGLUT			= par.Data->Ions.IonsLUT+par.MaxIons;
	par.Data->Channels.ChanType		= par.Data->Ions.ChanGLUT+par.MaxIons*( MAX_CHAN_PER_PUMP/4+1 );
	par.Data->Channels.ChanLUT		= par.Data->Channels.ChanType+par.MaxChan;
	par.Data->Cells.ChanGLUT		= par.Data->Channels.ChanLUT+par.MaxChan;
	par.Data->Cells.IonsILUT		= par.Data->Cells.ChanGLUT+par.MaxCells*( MAX_CHAN_PER_CELL/4+1 );
	par.Data->IOData.ViewLUT		= par.Data->GlobalViewLUT;
	// reset device-specific memory
	par.Data->DevMap = NULL;
	return par.Data;
}

///////////////////////////////////////////////////////////////////////////////
// +'lsns_free': frees both host and device memory
bool lsns_free( netpar &par )
{
	__lsns_assert( par.Data != NULL );
	// free device-specific memory for actual parameters and look-up-tables 
	// (for non-cuda version this is a pointer to the same data structure, so don't do anything)
	par.Data->DevMap = NULL;
	// free host-specific memory for both actual parameters and look-up-tables
	__lsns_assert( par.Data->GlobalData != NULL ); free( par.Data->GlobalData ); 
	__lsns_assert( par.Data->GlobalLUT != NULL ); free( par.Data->GlobalLUT );
	// free IO buffer memory for both look-up tables and actual data (for non-cuda version its the memory on host only)
	__lsns_assert( par.Data->GlobalViewData != NULL ); free( par.Data->GlobalViewData ); 
	__lsns_assert( par.Data->GlobalViewLUT != NULL ); free( par.Data->GlobalViewLUT );
	// free memory for netdat structure
	free( par.Data );
	// reset the constant parameters of the network
	Step = -1.f; Threshold = -10.f; StepCounter = 0;
	MaxIons = MaxChan = MaxCells = MaxViewPars = MaxGlobalData = 0;
	memset( Gates, 0, sizeof( netpar )*LSNS_MAX_GATEPARS );
	memset( Ions, 0, sizeof( ionspar )*LSNS_MAX_IONPARS );
	return true;
}

///////////////////////////////////////////////////////////////////////////////
// +'lsns_map2dev': allocates memory on device and copy the network configuration 
// from host to device setup the constant parameters of the network
bool lsns_map2dev( netpar &par )
{
	// setup the constant parameters of the network which are stored in constant memody for cuda version
	Step = par.Step; Threshold = par.Threshold; StepCounter = 0;
	MaxIons = par.MaxIons; MaxChan = par.MaxChan; MaxCells = par.MaxCells; MaxGlobalData = par.MaxGlobalData;
	MaxViewPars = par.MaxViewPars;
	memcpy ( Gates, par.Gates, sizeof( gatepar )*LSNS_MAX_GATEPARS );
	memcpy ( Ions, par.Ions, sizeof( ionspar )*LSNS_MAX_IONPARS );
	// allocate device-specific memory and copy initialized arrays from host memory to device memory
	// (for non-cuda version this is a pointer to the same data structure, so don't do anything)
	par.Data->DevMap = par.Data;
	///////////////////////////////////////////////////////////////////////////////
	// cuda version notes:
	// the array 'GlobalViewData' must be allocated ones as pinned-memory on host machine;
	// the arrays 'ViewData' must be allocated on device memory and should be organized in 
	// specific was that guarantees the maximal performance while coping the result of 
	// simulation (see 'store2host_kernel' method):
	//-------------------------------------------------------------------------------
	//	ViewData[0] = alloc_on_device(size*MAX_STORED_STEPS);
	//	for( i = 1; i < MAX_STORED_STEPS; ++i ){
	//		ViewData[i] = ViewData[0]+i*size;
	//	}
	return true;
}

///////////////////////////////////////////////////////////////////////////////
// +'lsns_run' runs simulation
bool lsns_run( netpar &par, int max_step )
{
	for( int step = 0; step < max_step; ++step ){
		control_kernel( 0 ); // don't do anything now
		connect_kernel( 0 ); // don't do anything now
		for( int ion = 0; ion < MaxIons; ++ion ){
			ions_kernel( ion, &( par.Data->DevMap->Ions ));
		}
		for( int chan = 0; chan < MaxChan; ++chan ){
			chan_kernel( chan, &( par.Data->DevMap->Channels ));
		}
		for( int cell = 0; cell < MaxCells; ++cell ){
			cell_kernel( cell, &( par.Data->DevMap->Cells ));
		}
		for( int view = 0; view < par.MaxViewPars/4+1; ++view ){
			store2dev_kernel( view, &( par.Data->DevMap->IOData ), step%MAX_STORED_STEPS );
		}
		if( step > 0 && step%MAX_STORED_STEPS == 0 ){
			for( int view = 0; view < ( par.MaxViewPars/4+1 )*MAX_STORED_STEPS; ++view ){
				store2host_kernel( view, &( par.Data->DevMap->IOData ));
			}
		}
	}
	return true;
}

#endif
