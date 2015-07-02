#include "engine.h"
#include "debug.h"

#if !defined ( __CUDA__ )

#define MAX_GATES LSNS_MAX_GATE
#define MAX_GPARS LSNS_MAX_GPARS
#define MAX_CHAN MaxChan
#define MAX_IONS MaxIons
#define MAX_CELLS MaxCells

#define lsns_ions_type data->ions_type
#define lsns_ions_e data->ions_e
#define lsns_ions_i data->ions_i
#define lsns_ions_shared data->ions_shared

#define lsns_chan_type data->chan_type
#define lsns_chan_g data->chan_g
#define lsns_chan_mh data->chan_mh
#define lsns_chan_shared data->chan_shared

#define lsns_cell_v data->cell_v
#define lsns_gchan_lut data->gchan_lut
#define lsns_ipump_lut data->ipump_lut

static int MaxChan;
static int MaxIons;
static int MaxCells;

///////////////////////////////////////////////////////////////////////////////
// ions_kernel: the engine kernel to calculate the properties of ions dynamics.
// such as pump current, concentration of ions inside the cell, etc.
// Input parameters:
// Output parameters:
void ions_kernel( float step, int index, ionproc *data )
{
	__lsns_assert( index < MAX_IONS );			// DEBUG: check the range for 'index' variable
	// load 
	int4 tp = lsns_ions_type[index];			// load type of ions 
	int type_dyn = _ions_typepump( tp );
	int type_eds = _ions_typeeds( tp );
	if( type_dyn != LSNS_NO_DYN ){
		// load ions parameters, resting potential, ions concentrations if needed
		int4 sh = lsn_ions_shared[index];		// load references to external parameters (channel conductances etc)
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
		lsns_fastsum( _chan_g, lsns_chan_g, lsn_gchan_lut+index*( MAX_CHAN_PER_PUMP/4+1 ), gchan, MAX_CHAN ); // sum all gchan
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
// chan_kernel: the engine kernel to calculate the properties of channel and 
// synaptic currents such as conductance, current etc.
// Input parameters:
// Output parameters:
void chan_kernel( float step, int index, chanproc *data )
{
	__lsns_assert( index < MAX_CHAN );			// DEBUG: check the range for 'index' variable 
	int4 tp = lsns_chan_type[index];			// load type of ion channel (generic, a-b, etc) and its parameters (half-voltage, slope, etc)
	int4 sh = lsns_chan_shared[index];			// load references to external parameters (membrane potential, rest potential, etc)
	float4 g = lsns_chan_g[index];				// load properties of ions channel (conductance, current, etc)
	__lsns_assert( _gate_typem( tp ) < MAX_GATES );		// DEBUG: check the range for _gate_typem( tp )
	__lsns_assert( _gate_typeh( tp ) < MAX_GATES );		// DEBUG: check the range for  _gate_typeh( tp ) 
	__lsns_assert( _gate_parm( tp ) < MAX_GPARS );		// DEBUG: check the range for  _gate_parm( tp ) 
	__lsns_assert( _gate_parh( tp ) < MAX_GPARS );		// DEBUG: check the range for  _gate_parh( tp ) 
	__lsns_assert( _chan_lut_v( sh ) < MAX_CELLS );		// DEBUG: check the range for _chan_lut_v( sh )
	__lsns_assert( _chan_lut_e( sh ) < MAX_IONS );		// DEBUG: check the range for _chan_lut_e( sh )
	__lsns_assert( _chan_lut_inm( sh ) < MAX_IONS );	// DEBUG: check the range for _chan_lut_inm( sh )
	__lsns_assert( _chan_lut_inh( sh ) < MAX_IONS );	// DEBUG: check the range for _chan_lut_inh( sh )
	// load properties of gate variables (activation, inactivation, etc) if needed
	float4 mh = ( _gate_typem( tp )+_gate_typeh( tp ) != LSNS_NOGATE)? lsns_chan_mh[index]: float4();
//--->>> possible CUDA optimization (try to use shared variables)
	// load shared variables (resting potential, membrane potential, etc)
	float eds = _ions_eds( lsns_ions_e[_chan_lut_e( sh )] );// extract resting potential from 'ions_e'
	float vm = _cell_v( lsns_cell_v[_chan_lut_v( sh )] );	// extract membrane potential from 'cell_v'
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
	// store the results of simulation
	lsns_chan_g[index] = g;					// save properties of ions channel (conductance, current, etc)
	if( _gate_typem( tp )+_gate_typeh( tp ) != LSNS_NOGATE ){
		lsns_chan_mh[index] = mh;			// save properties of gate variables (activation, inactivation, etc) if needed
	}
}

///////////////////////////////////////////////////////////////////////////////
// cell_kernel: the engine kernel to calculate the properties of Hodgkin-Huxley
// cell such as membrane potential, onset of the spike etc.
// Input parameters:
// Output parameters:
void cell_kernel( float step, int index, cellproc *data )
{
	__lsns_assert( index < MAX_CELLS );			// DEBUG: check the range for 'index' variable 
	// load properties of the cell (membrane potential, etc)
	float4 v = lsns_cell_v[index];
	// preprocessing
	float g = 1, ge = 0, ipump = 0;
	lsns_fastsum( _ions_ipump, lsns_ions_i, lsns_ipump_lut+index*( MAX_IPUMP_PER_CELL/4+1 ), ipump, MAX_IONS ); // sum all ipump
	lsns_fastsum2( _chan_g, _chan_ge, lsns_chan_g, lsn_gchan_lut+index*( MAX_CHAN_PER_CELL/4+1 ), g, ge, MAX_CHAN ); // sum all g and ge
	__lsns_assert( g != 0.0 );				// DEBUG: check if conductance is not zero
	// perform calculations
	float time = _cell_c( v )/g;
	float v_inf = ( ge+ipump+_cell_iadd( v ))/g;
	float vm = lsns_exp_euler( _cell_v( v ), v_inf, step, time );
	float spike = ( vm > thresh && _cell_v( v )  <= thresh )? 1.: 0.;
	// store the results of simulation
	_cell_v( v ) = vm;
	_cell_spike( v ) = spike;
	lsns_cell_v[index] = v;
}

#endif

