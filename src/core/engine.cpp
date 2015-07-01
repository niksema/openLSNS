#if !defined ( __CUDA__ )

#include <assert.h>
#include "engine.h"

typedef struct __ions_proc{
	int4 *ions_type;
	int4 *ions_shared;
	float4 *ions_e;
	float4 *ions_i;
	float4 *chan_g;
	int4 *chan_ig;
} ionsproc;

///////////////////////////////////////////////////////////////////////////////
// ions_kernel: the engine kernel to calculate the properties of ions dynamics.
// such as pump current, concentration of ions inside the cell, etc.
// Input parameters:
// Output parameters:
void ions_kernel( float step, int index, ionsproc *data )
{
	assert( index < MaxIons );				// DEBUG: check the range for 'index' variable

	// load 
	int4 tp = ions_type[index];				// load type of ions 
	int type_dyn = ???;
	int type_eds = ???;
	// load data
	if( type_dyn != LSNS_NO_DYN ){
		// load ions parameters, resting potential, ions concentrations if needed
		int4 sh = ions_shared[index];			// load references to external parameters (channel conductances etc)
		
		float4 e = data->ions_e[index];			// load ions parameters
		float4 i;					// ions currents (no need to load)
		float v = 0; 					//-->> get membrane potential
		float tau = 1;					//-->> get time constant for ions dynamics

		float out = _ions_out( e );			// concentation of ions outside cell
		float rtfz = _ions_rtfz( e );			// rtfz for specific neurons
		// start calculations
		float gchan = 0; 				// total conductances of all ion currents which are involved to ions dynamics
		lsns_fastsum( _gate_g, chan_g, chan_ig+index*3, gchan );
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
		
		_ions_in( e ) = in;
		_ions_out( e ) = out;
		_ions_eds( e ) = eds;
		_ions_ipump( i ) = ipump;
		_ions_ichan( i ) = ichan;
		// store the results
		data->ions_e[index] = e;
		data->ions_i[index] = i;
	}
	else if( type_eds != LSNS_NOSPEC_EDS ){			// calculate reversal potential only (no ions dynamics)
		// load ions parameters
		float4 e = data->ions_e[index];			// load ions parameters
		float in = _ions_in( e );			// concentation of ions inside the cell
		float out = _ions_out( e );			// concentation of ions outside cell
		float rtfz = _ions_rtfz( e );			// rtfz for specific ions
		// calculate resting potential
		_ions_eds( e ) = lsns_eds( rtfz, in, out );
		// store the results
		data->ions_e[index] = e;
	}
}

typedef struct __chan_proc{
	int4 *chan_type;
	int4 *chan_shared;
	float4 *chan_g;
	float4 *chan_mh;
	float4 *cell_v;
	float4 *ions_e;
} chanproc;

typedef struct __cell_proc{
	int4 *cell_shared;
	float4 *cell_v;
	float4 *chan_g;
	float4 *ions_i;
} cellproc;

///////////////////////////////////////////////////////////////////////////////
// chan_kernel: the engine kernel to calculate the properties of channel and 
// synaptic currents such as conductance, current etc.
// Input parameters:
// Output parameters:
void chan_kernel( float step, int index, chanproc *data )
{
	assert( index < MaxChan );				// DEBUG: check the range for 'index' variable 
	int4 tp = data->chan_type[index];			// load type of ion channel (generic, a-b, etc) and its parameters (half-voltage, slope, etc)
	int4 sh = data->chan_shared[index];			// load references to external parameters (membrane potential, rest potential, etc)
	float4 g = data->chan_g[index];				// load properties of ions channel (conductance, current, etc)
	assert( _gate_typem( tp ) < LSNS_MAX_GATE );		// DEBUG: check the range for _gate_typem( tp )
	assert( _gate_typeh( tp ) < LSNS_MAX_GATE );		// DEBUG: check the range for  _gate_typeh( tp ) 
	assert( _gate_v( sh ) < MaxCells );			// DEBUG: check the range for _gate_v( sh )
	assert( _gate_e( sh ) < MaxIons );			// DEBUG: check the range for _gate_e( sh )
	assert( _gate_inm( sh ) < MaxIons );			// DEBUG: check the range for _gate_inm( sh )
	assert( _gate_inh( sh ) < MaxIons );			// DEBUG: check the range for _gate_inh( sh )
	// load properties of gate variables (activation, inactivation, etc) if needed
	float4 mh = ( _gate_typem( tp )+_gate_typeh( tp ) != LSNS_NOGATE)? data->chan_mh[index]: float4();
	// load shared variables (resting potential, membrane potential, etc)
	float eds = _ions_eds( data->ions_e[_gate_e( sh )] );	// extract resting potential from 'ions_e'
	float vm = _cell_v( data->cell_v[_gate_v( sh )] );	// extract membrane potential from 'cell_v'
	// load Mg- or Ca- concentration inside the cell for NMDA synapse or Z-channels from 'ions_e' if needed
	float in_m = ( _gate_typem( tp ) >= LSNS_ZGENERIC_INSTANT )? _ions_in( _ions_in( data->ions_e[_gate_inm( sh )] ):0;
	// load Mg- or Ca- concentration inside the cell for NMDA synapse or Z-channels from 'ions_e' if needed
	float in_h = ( _gate_typeh( tp ) >= LSNS_ZGENERIC_INSTANT )? _ions_in( _ions_in( data->ions_e[_gate_inh( sh )] ):0;
	// perform calculations
	float mp, hp;
	proc_gate( _gate_typem( tp ), Gates[_gate_parm( tp )], in_m, vm, step, _gate_powm( mh ), _gate_m( mh ), mp );
	proc_gate( _gate_typeh( tp ), Gates[_gate_parh( tp )], in_h, vm, step, _gate_powh( mh ), _gate_h( mh ), hp );
	_gate_g( g ) = _gate_gmax( g )*mp*hp;			// g
	_gate_ge( g ) = _gate_g( g )*eds;			// ge
	_gate_i( g ) = _gate_g( g )*( vm-eds );			// I
	// store the results of simulation
	data->chan_g[index] = g;				// save properties of ions channel (conductance, current, etc)
	if( _gate_typem( tp )+_gate_typeh( tp ) != LSNS_NOGATE ){
		data->chan_mh[index] = mh;			// save properties of gate variables (activation, inactivation, etc) if needed
	}
}

///////////////////////////////////////////////////////////////////////////////
// cell_kernel: the engine kernel to calculate the properties of Hodgkin-Huxley
// cell such as membrane potential, onset of the spike etc.
// Input parameters:
// Output parameters:
void cell_kernel( float step, int index, cellproc *data )
{
	assert( index < MaxCell );				// DEBUG: check the range for 'index' variable 
	float4 v = data->cell_v[index];				// load properties of the cell (membrane potential, etc)
	// start calculations
	float g = 1, ge = 0, ipump = 0;
	lsns_fastsum( _ions_ipump, data->ions_i, data->ions_ii+index, ipump ); // sum all ipump
	lsns_fastsum2( _gate_g, _gate_ge, data->chan_g, data->chan_ig+index*4, g, ge ); // sum all g and ge
	assert( g != 0.0 );					// DEBUG: check if conductance is not zero
	//--- 2. perform calculations
	float time = _cell_c( v )/g;
	float v_inf = ( ge+ipump+_cell_iadd( v ))/g;
	float vm = lsns_exp_euler( _cell_v( v ), v_inf, step, time );
	float spike = ( vm > thresh && _cell_v( v )  <= thresh )? 1.: 0.;
	//--- 3.  store the results of simulation
	_cell_v( v ) = vm;
	_cell_spike( v ) = spike;
	// store the results of simulation
	data->cell_v[index] = v;
}

#endif

