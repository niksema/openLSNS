#include <cuda.h>
#include <helper_cuda.h>

#include "engine.h"

#if defined( __CUDA__ )
	#define _cu_index() threadIdx.x+blockIdx.x*blockDim.x+blockIdx.y*blockDim.x*gridDim.x
#endif /*__CUDA__*/

__constant__ __lsns_align(16) float SimStep = 0.1;
__constant__ __lsns_align(16) gate_par Gates[LSNS_MAX_GPARS];

// calculate ions (aka channels) and synaptic currents. 
__global__ void ichan_kernel( int4 *chan_type, int4 *chan_shared, float4 *chan_g, float4 *chan_mh, float4 *cell_v, float4 *ions_e )
{
	//--- 1. load data
	int index = _cu_index();
	float step = SimStep;					// get integration step
	//--- 1.1. load local variables
// todo: check the range for 'index' variable in debug mode [0...MaxChan]
	int4 tp = chan_type[index];				// load type of ion channel (generic, a-b, etc) and its parameters (half-voltage, slope, etc)
	int4 sh = chan_shared[index];				// load references to external parameters (membrane potential, rest potential, etc)
	float4 g = chan_g[index];				// load properties of ions channel (conductance, current, etc)
// todo: check the range for _gate_typem( tp ) and _gate_typeh( tp ) in debug mode [LSNS_NOGATE...LSNS_MAX_GATE]
	// load properties of gate variables (activation, inactivation, etc) if needed
	float4 mh = ( _gate_typem( tp )+_gate_typeh( tp ) != LSNS_NOGATE)? chan_mh[index]: float4();
	//--- 1.2. load shared variables
// todo: check the range for _gate_v( sh ) and _gate_e( sh ) indices for cell_v [0...MaxCells] and ions_e [0...MaxIons] variables in debug mode 
// todo: the possible optimization is to use shared variables to load v and e parameters
	float eds = _ions_eds( ions_e[_gate_e( sh )] );		// extract resting potential from 'ions_e'
	float vm = _cell_v( cell_v[_gate_v( sh )] );		// extract membrane potential from 'cell_v'
// todo: check the range for _gate_in_m( sh ) and _gate_in_m( sh ) indices for ions_e variable in debug mode [0...MaxIons]
	// load Mg- or Ca- concentration inside the cell for NMDA synapse or Z-channels from 'ions_e' if needed
	float in_m = ( _gate_typem( tp ) >= LSNS_ZGENERIC_INSTANT )? _ions_in( e ):0;
	// load Mg- or Ca- concentration inside the cell for NMDA synapse or Z-channels from 'ions_e' if needed
	float in_h = ( _gate_typeh( tp ) >= LSNS_ZGENERIC_INSTANT )? _ions_in( e ):0;
	//--- 2. perform calculations
	float mp, hp;
	proc_gate( _gate_typem( tp ), Gates[_gate_parm( tp )], in_m, vm, step, _gate_powm( mh ), _gate_m( mh ), mp );
	proc_gate( _gate_typeh( tp ), Gates[_gate_parh( tp )], in_h, vm, step, _gate_powh( mh ), _gate_h( mh ), hp );
	_gate_g( g ) = _gate_gmax( g )*mp*hp;			// g
	_gate_i( g ) = _gate_g( g )*( vm-eds );			// I
	_gate_ge( g ) = _gate_g( g )*eds;			// ge
	//--- 3.  store the results of simulation
	chan_g[index] = g;					// save properties of ions channel (conductance, current, etc)
	if( _gate_typem( tp )+_gate_typeh( tp ) != LSNS_NOGATE ){
		chan_mh[index] = mh;				// save properties of gate variables (activation, inactivation, etc) if needed
	}
}

// calculate ions dynamics. 
__global__ void ions_kernel( int4 *ions_type, int4 *ions_shared, float4 *ions_e, float4 *ions_i, float4 *chan_g, int4 *chan_ig )
{
	int index = _cu_index();				// 
	float step = SimStep;					// integration step
// todo: check the range for 'index' variable in debug mode [0...MaxIons]
	// load 
	int4 tp = ions_type[index];				// load type of ions 
	int type_dyn = ???;
	int type_eds = ???;
	// load data
	if( type_dyn != LSNS_NO_DYN ){
		// load ions parameters, resting potential, ions concentrations if needed
		int4 sh = ions_shared[index];				// load references to external parameters (channel conductances etc)
		
		float4 e = ions_e[index];				// load ions parameters
		float v = 0; 						//-->> get membrane potential
		float tau = 1;						//-->> get time constant for ions dynamics
		float in = _ions_in( e );				// concentation of ions inside the cell
		float out = _ions_out( e );				// concentation of ions outside cell
		float rtfz = _ions_rtfz( e );				// rtfz for specific neurons
		float eds = _ions_eds( e );				// resting potential
		// start calculations
		float gchan = 0; 						// total conductances of all ion currents which are involved to ions dynamics
		lsns_fastsum( _gate_g, chan_g, chan_ig+index*3, gchan );
		// 2-order Runge-Kutta integration method
		float ipump = 0, apump = 0, ichan = 0;
		// sub-step 1
		esd = ( type_eds != LSNS_NOSPEC_EDS )? lsns_eds( rtfz, in, out ): eds; // calculate resting potential
		lsns_ipump( type_dyn, par, in, apump, ipump );	// calculate pump current
		lsns_ichan( apump, v, eds, gchan, ichan );	// calculate channels current
		// sub-step 2
		float in1 = in-step*( ichan+ipump )/( 2*tau );
		esd = ( type_eds != LSNS_NOSPEC_EDS )? lsns_eds( rtfz, in1, out ): eds; // calculate resting potential
		lsns_ipump( type_dyn, par, in1, apump, ipump );	// calculate pump current
		lsns_ichan( apump, v, eds, gchan, ichan );	// calculate channels current
		// the final calculations for in, eds, ichan
		in = in-step*( ichan+ipump )/tau;
		esd = ( type_eds != LSNS_NOSPEC_EDS )? lsns_eds( rtfz, in, out ): eds; // calculate resting potential
		lsns_ipump( type_dyn, par, in, apump, ipump );	// calculate pump current
		lsns_ichan( apump, v, eds, gchan, ichan );	// calculate channels current
		
		_ions_in( e ) = in;
		_ions_out( e ) = out;
		_ions_eds( e ) = eds;
		// store the results
		float4 i;
		_ions_ipump( i ) = ipump;
		_ions_ichan( i ) = ichan;
		ions_e[index] = e;
		ions_i[index] = i;
	}
	else if( type_eds != LSNS_NOSPEC_EDS ){
		// load ions parameters
		float4 e = dev_e[index];				// load ions parameters
		float in = _ions_in( e );				// concentation of ions inside the cell
		float out = _ions_out( e );				// concentation of ions outside cell
		float rtfz = _ions_rtfz( e );				// rtfz for specific neurons
		// calculate resting potential
		_ions_eds( e ) = lsns_eds( rtfz, in, out );
		// store the results
		ions_e[index] = e;
	}
}

// dev_v: V(x), G(y), GE(z), I(w)
__global__ void neurons_kernel( int4 *cell_shared, float4 *cell_v, float4 *chan_g )
{
	int index = _cu_index();
	
	//--- load local variables
	float4 v = cell_v[index];
	v.x = -60; 
	// start calculations
	float c_cell = 1; 					// membrane capacitance
	float g_cell = 1, ge_cell = 0, ipump_cell = 0, iadd_cell = 0;
// todo sum all g
// todo sum all ge
// todo sum all ipump
	//--- 2. perform calculations
	float time = c_cell/g_cell;
	float v_inf = ( ge_cell+ipump_cell+iadd_cell )/g_cell;
	float vm = lsns_exp_euler( v.x, v_inf, step, time );
	float spike = ( vm > thresh && v.x  <= thresh )? 1.: 0.;
	//--- 3.  store the results of simulation
	v.x = vm;
	v.y = spike;
	dev_v[index] = v;
}
