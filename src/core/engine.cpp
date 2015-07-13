#include "engine.h"
#include "debug.h"

#if !defined ( __CUDA__ )

#include <stdlib.h>
#include <string.h>

#define MAX_GATES LSNS_MAX_GATE
#define MAX_GPARS LSNS_MAX_GPARS
#define MAX_CHAN MaxChan
#define MAX_IONS MaxIons
#define MAX_CELLS MaxCells
#define MAX_VIEW_PARS MaxViews

#define lsns_ions_type ( data->IonsType )
#define lsns_ions_e ( data->IonsE )
#define lsns_ions_i ( data->IonsI )
#define lsns_ions_shared ( data->IonsShared )

#define lsns_chan_type ( data->ChanType )
#define lsns_chan_g ( data->ChanG )
#define lsns_chan_mh ( data->ChanMH )
#define lsns_chan_shared ( data->ChanShared )

#define lsns_cell_v ( data->CellV )
#define lsns_gchan_lut ( data->GchanLUT )
#define lsns_ipump_lut ( data->IpumpLUT )
/*
#define lsns_store_data( res, lut, dat ) \
	{
	
		int pos = ( lut.x >> 30 )
		int index = lut.x & 0x3FFFFFFF;
		res.x = 
		index = lut.y & 0x3FFFFFFF;
		res.y = 
		index = lut.z & 0x3FFFFFFF;
		res.z = 
		index = lut.w & 0x3FFFFFFF;
		res.w = 
	}
*/
#define lsns_dev_dat ( data->DevData[0] )
#define lsns_host_dat ( data->HostData )

static float Step = -1.f;
static float Threshold = -10.f;

static int Counter = 0;
static int MaxChan = 0;
static int MaxIons = 0;
static int MaxCells = 0;
static int MaxViews = 0;
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
	__lsns_assert( index < MAX_IONS );			// DEBUG: check the range for 'index' variable
	// load 
	float step = Step;
	// load type of ions (Na, K, Ca etc, and type of ion dynamics)
	int4 tp = lsns_ions_type[index];
	int type_dyn = _ions_typepump( tp );
	int type_eds = _ions_typeeds( tp );
	if( type_dyn != LSNS_NO_DYN ){
		// load ions parameters, resting potential, ions concentrations if needed
		pumppar par = Pumps[_ions_parpump( tp )];	// load pump properties
		int4 sh = lsns_ions_shared[index];		// load references to external parameters (channel conductances etc)
		float4 e = lsns_ions_e[index];			// load ions parameters
		float4 i = lsns_ions_i[index];			// load ions parameters
//--->>> possible CUDA optimization (try to use shared variables)
		float v = _cell_v( lsns_cell_v[_ions_lut_v( sh )] );// get membrane potential
//---<<< possible CUDA optimization
		float tau = _ions_tau( i );			// get time constant for ions dynamics
		float out = _ions_out( e );			// concentation of ions outside cell
		float rtfz = _ions_rtfz( e );			// rtfz for specific neurons
		// start calculations
		float gchan = 0; 				// total conductances of all ion currents which are involved to ions dynamics
		lsns_fastsum( _chan_g, lsns_chan_g, lsns_gchan_lut+index*( MAX_CHAN_PER_PUMP/4+1 ), gchan, MAX_CHAN ); // sum all gchan
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
		lsns_ions_e[index] = e;
		lsns_ions_i[index] = i;
	}
	else if( type_eds != LSNS_NOSPEC_EDS ){			// calculate reversal potential only (no ions dynamics)
		// load ions parameters
		float4 e = lsns_ions_e[index];			// load ions parameters
		float in = _ions_in( e );			// concentation of ions inside the cell
		float out = _ions_out( e );			// concentation of ions outside cell
		float rtfz = _ions_rtfz( e );			// rtfz for specific ions
		// calculate resting potential
		_ions_eds( e ) = lsns_eds( rtfz, in, out );
		// store the results
		lsns_ions_e[index] = e;
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
	__lsns_assert( index < MAX_CHAN );			// DEBUG: check the range for 'index' variable 
	float step = Step;
	// load type of ion channel (generic, a-b, etc) and its parameters (half-voltage, slope, etc)
	int4 tp = lsns_chan_type[index];
	// load references to external parameters (membrane potential, rest potential, etc)
	int4 sh = lsns_chan_shared[index];
	__lsns_assert( _gate_typem( tp ) < MAX_GATES );		// DEBUG: check the range for _gate_typem( tp )
	__lsns_assert( _gate_typeh( tp ) < MAX_GATES );		// DEBUG: check the range for  _gate_typeh( tp ) 
	__lsns_assert( _gate_parm( tp ) < MAX_GPARS );		// DEBUG: check the range for  _gate_parm( tp ) 
	__lsns_assert( _gate_parh( tp ) < MAX_GPARS );		// DEBUG: check the range for  _gate_parh( tp ) 
	__lsns_assert( _chan_lut_v( sh ) < MAX_CELLS );		// DEBUG: check the range for _chan_lut_v( sh )
	__lsns_assert( _chan_lut_e( sh ) < MAX_IONS );		// DEBUG: check the range for _chan_lut_e( sh )
	__lsns_assert( _chan_lut_inm( sh ) < MAX_IONS );	// DEBUG: check the range for _chan_lut_inm( sh )
	__lsns_assert( _chan_lut_inh( sh ) < MAX_IONS );	// DEBUG: check the range for _chan_lut_inh( sh )
	// load properties of ions channel (conductance, current, etc)
	float4 g = lsns_chan_g[index];
	// load properties of gate variables (activation, inactivation, etc) if needed
	float4 mh = ( _gate_typem( tp )+_gate_typeh( tp ) != LSNS_NOGATE)? lsns_chan_mh[index]: float4();
//--->>> possible CUDA optimization (try to use shared variables)
	// load shared variables (resting potential, membrane potential, etc)
	float eds = _ions_eds( lsns_ions_e[_chan_lut_e( sh )] );// extract resting potential from 'IonsE'
	float vm = _cell_v( lsns_cell_v[_chan_lut_v( sh )] );	// extract membrane potential from 'CellV'
	// load Mg- or Ca- concentration inside the cell for NMDA synapse or Z-channels from 'ions_e' if needed
	float in_m = ( _gate_typem( tp ) >= LSNS_ZGENERIC_INSTANT )? _ions_in( lsns_ions_e[_chan_lut_inm( sh )] ):0;
	// load Mg- or Ca- concentration inside the cell for NMDA synapse or Z-channels from 'ions_e' if needed
	float in_h = ( _gate_typeh( tp ) >= LSNS_ZGENERIC_INSTANT )? _ions_in( lsns_ions_e[_chan_lut_inh( sh )] ):0;
//---<<< possible CUDA optimization
	// perform calculations
	float mp, hp;
	proc_gate( _gate_typem( tp ), Gates[_gate_parm( tp )], in_m, vm, step, _gate_powm( mh ), _gate_m( mh ), mp );
	proc_gate( _gate_typeh( tp ), Gates[_gate_parh( tp )], in_h, vm, step, _gate_powh( mh ), _gate_h( mh ), hp );
	_chan_g( g ) = _chan_gmax( g )*mp*hp;			// g
	_chan_ge( g ) = _chan_g( g )*eds;			// ge
	_chan_i( g ) = _chan_g( g )*( vm-eds );			// I
	// save properties of ions channel (conductance, current, etc)
	lsns_chan_g[index] = g;
	if( _gate_typem( tp )+_gate_typeh( tp ) != LSNS_NOGATE ){
		// save properties of gate variables (activation, inactivation, etc) if needed
		lsns_chan_mh[index] = mh;
	}
}

///////////////////////////////////////////////////////////////////////////////
// cell_kernel: the engine kernel to calculate the properties of Hodgkin-Huxley
// cell such as membrane potential, onset of the spike etc.
// Input parameters:
// Output parameters:
void cell_kernel( int index, celldat *data )
{
	__lsns_assert( index < MAX_CELLS );			// DEBUG: check the range for 'index' variable 
	// load properties of the cell (membrane potential, etc)
	float step = Step;
	float threshold = Threshold;
	float4 v = lsns_cell_v[index];
	// preprocessing
	float g = 1, ge = 0, ipump = 0;
	lsns_fastsum( _ions_ipump, lsns_ions_i, lsns_ipump_lut+index*( MAX_IPUMP_PER_CELL/4+1 ), ipump, MAX_IONS ); // sum all ipump
	lsns_fastsum2( _chan_g, _chan_ge, lsns_chan_g, lsns_gchan_lut+index*( MAX_CHAN_PER_CELL/4+1 ), g, ge, MAX_CHAN ); // sum all g and ge
	__lsns_assert( g != 0.0 );				// DEBUG: check if conductance is not zero
	// perform calculations
	float time = _cell_c( v )/g;
	float v_inf = ( ge+ipump+_cell_iadd( v ))/g;
	float vm = lsns_exp_euler( _cell_v( v ), v_inf, step, time );
	float spike = ( vm > threshold && _cell_v( v )  <= threshold )? 1.f: 0.f;
	// store the results of simulation
	_cell_v( v ) = vm;
	_cell_spike( v ) = spike;
	lsns_cell_v[index] = v;
}

///////////////////////////////////////////////////////////////////////////////
// store_kernel: 
// Input parameters:
// Output parameters:
/*
	// local variables (read-write)
	float4 __lsns_align( 16 ) *DevData[MAX_STORED_STEPS];	// data to display (located in device memory)
	float4 __lsns_align( 16 ) *HostData;			// data to display (pinned memory located in the host)

	devdata = alloc(size*MAX_STORED_STEPS);
	for( i = 0; i < MAX_STORED_STEPS; ++i ){
		DevData[i] = devdata+i*size;
	}
	HostData = alloc(size*MAX_STORED_STEPS);

	// local variables (read-only). 
	// LUT format: bits 31..30 are coding the offset in each float4 variable (00 - x, 01 - y, 10 - z, 11 - w); 
	// bits 29..0 are coding the offset in an arrays of shared variables
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

*/

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
void store2dev_kernel( int index /*, iodat *data*/ )
{
	__lsns_assert( index < MAX_VIEW_PARS );			// DEBUG: check the range for 'index' variable
	/*
	int4 lut;
	float4 data;
	if( index < MAX_IIONS_VIEWS ){
		int i = index;
		int4 lut = IonsILUT[i];
		lsns_store_data( data, lut, IonsI ); 
	} 
	else if( index < MAX_IIONS_VIEWS+MAX_EIONS_VIEWS ){
		int i = index-MAX_IIONS_VIEWS;
		int4 lut = IonsELUT[i];
		lsns_store_data( data, lut, IonsI ); 
	}
	else if( index < MAX_IIONS_VIEWS+MAX_EIONS_VIEWS+MAX_GCHAN_VIEWS ){
		int i = index-MAX_IIONS_VIEWS-MAX_EIONS_VIEWS;
		int4 lut = ChanGLUT[i];
		lsns_store_data( data, lut, IonsI ); 
	}
	else if( index < MAX_IIONS_VIEWS+MAX_EIONS_VIEWS+MAX_GCHAN_VIEWS+MAX_MHCHAN_VIEWS ){
		int i = index-MAX_IIONS_VIEWS-MAX_EIONS_VIEWS-MAX_GCHAN_VIEWS;
		int4 lut = ChanMHLUT[i];
		lsns_store_data( data, lut, IonsI ); 
	}
	else if( index < MAX_IIONS_VIEWS+MAX_EIONS_VIEWS+MAX_GCHAN_VIEWS+MAX_MHCHAN_VIEWS+MAX_CELLV_VIEWS ){
		int i = index-MAX_IIONS_VIEWS-MAX_EIONS_VIEWS-MAX_GCHAN_VIEWS-MAX_MHCHAN_VIEWS;
		int4 lut = CellVLUT[i];
		lsns_store_data( data, lut, IonsI ); 
	}
	DevData[index+offset] = data;
	}
	*/
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
	__lsns_assert( index < MAX_VIEW_PARS );			// DEBUG: check the range for 'index' variable
	lsns_host_dat[index] = lsns_dev_dat[index];		// copy data from device memory to host memory
}

///////////////////////////////////////////////////////////////////////////////
// interface 
 
// allocate memory on host
netdat *lsns_alloc( netpar &par )
{
	netdat *data = ( netdat *)malloc( sizeof( netdat )); __lsns_assert( data != NULL );
	data->Ions.MaxIons = par.MaxIons; 
	data->Channels.MaxChan = par.MaxChan;
	data->Cells.MaxCells = par.MaxCells;
	// allocate memory for actual parameters
	data->Ions.IonsI = data->Cells.IonsI = data->IOData.IonsI = ( float4 *)malloc( sizeof( float4 )*MaxIons );
	data->Ions.IonsE = data->Channels.IonsE = data->IOData.IonsE = ( float4 *)malloc( sizeof( float4 )*MaxIons );
	data->Channels.ChanG = data->Cells.ChanG = data->Ions.ChanG = data->IOData.ChanG = ( float4 *)malloc( sizeof( float4 )*MaxChan );
	data->Channels.ChanMH = data->IOData.ChanMH = ( float4 *)malloc( sizeof( float4 )*MaxChan );
	data->Cells.CellV = data->Channels.CellV = data->Ions.CellV = data->IOData.CellV = ( float4 *)malloc( sizeof( float4 )*MaxCells );
	// allocate memory for look-up-tables
	data->Ions.IonsType = ( int4 * )malloc( sizeof( int4 )*MaxIons ); 
	data->Ions.IonsShared = ( int4 * )malloc( sizeof( int4 )*MaxIons );
	data->Ions.GchanLUT = ( int4 * )malloc( sizeof( int4 )*MaxIons*( MAX_CHAN_PER_PUMP/4+1 ) );
	data->Channels.ChanType = ( int4 * )malloc( sizeof( int4 )*MaxChan );
	data->Channels.ChanShared = ( int4 * )malloc( sizeof( int4 )*MaxChan );
	data->Cells.GchanLUT = ( int4 *)malloc( sizeof( int4 )*MaxCells*( MAX_CHAN_PER_CELL/4+1 ) ); 
	data->Cells.IpumpLUT = ( int4 *)malloc( sizeof( int4 )*MaxCells*( MAX_IPUMP_PER_CELL/4+1 ) );
	// device-specific memory (for non-cuda version this is a pointer to the same data structure)
	data->DevMap = data;
	// setup the constant parameters of the network
	Step = par.Step; Threshold = par.Threshold; MaxIons = par.MaxIons; MaxChan = par.MaxChan; MaxCells = par.MaxCells;
	memcpy ( Gates, par.Gates, sizeof( netpar )*LSNS_MAX_GPARS );
	memcpy ( Pumps, par.Pumps, sizeof( pumppar )*LSNS_MAX_PARPUMPS );
	return data;
}

// free memory both on host and device
void lsns_free( netdat *data )
{
	// free memory for actual parameters
	free( data->Ions.IonsI ); data->Ions.IonsI = data->Cells.IonsI = data->IOData.IonsI = NULL;
	free( data->Ions.IonsE ); data->Ions.IonsE = data->Channels.IonsE = data->IOData.IonsE = NULL;
	free( data->Channels.ChanG ); data->Channels.ChanG = data->Cells.ChanG = data->Ions.ChanG = data->IOData.ChanG = NULL;
	free( data->Channels.ChanMH ); data->Channels.ChanMH = data->IOData.ChanMH = NULL;
	free( data->Cells.CellV ); data->Cells.CellV = data->Channels.CellV = data->Ions.CellV = data->IOData.CellV = NULL;
	// free memory for look-up-tables
	free( data->Ions.IonsType );
	free( data->Ions.IonsShared );
	free( data->Ions.GchanLUT );
	free( data->Channels.ChanType );
	free( data->Channels.ChanShared );
	free( data->Cells.GchanLUT );
	free( data->Cells.IpumpLUT );
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

