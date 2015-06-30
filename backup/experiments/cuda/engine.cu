#include <cuda.h>
#include <helper_cuda.h>

#include "engine.h"
#include "gateproc.h"

// type: type_m (x), par_m(y), type_h(z), par_h(w)
#define _gate_typem( t ) ( t ).x
#define _gate_parm( t ) ( t ).y
#define _gate_typeh( t ) ( t ).z
#define _gate_parh( t ) ( t ).w

// shared (indices of shared variables): V(x), E(y), ions inside the cell for activation(z), ions inside the cell for inactivation(w)
#define _gate_v( sh ) ( sh ).x
#define _gate_e( sh ) ( sh ).y
#define _gate_in_m( sh ) ( sh ).z
#define _gate_in_h( sh ) ( sh ).w

// chan_g: Gmax(x), g(y), I(z), GE(w)
#define _gate_gmax( g ) ( g ).x
#define _gate_g( g ) ( g ).y
#define _gate_i( g ) ( g ).z
#define _gate_ge( g ) ( g ).w

// chan_mh: M(x), H(y), Power M(z), Power H(w)
#define _gate_m( mh ) ( mh ).x
#define _gate_h( mh ) ( mh ).y
#define _gate_powm( mh ) ( mh ).z
#define _gate_powh( mh ) ( mh ).w

// cell_v: V(x), reserver(y), reserved(z), reserved(w)
#define _cell_v( v ) ( v ).x

// ions_e: E(x), In(y), Out(z), RTFz(w)
#define _ions_eds( e ) ( e ).x
#define _ions_in( e ) ( e ).y
#define _ions_out( e ) ( e ).z
#define _ions_rtfz( e ) ( e ).w

// ions_i: Ipump(x), Ichan(y), reserved(z), reserved(w)
#define _ions_ipump( i ) ( i ).x
#define _ions_ichan( i ) ( i ).y

// ions_type: pump type(x), eds type(y), reserved(z), reserved(w)
#define _ions_typepump( tp ) ( tp ).x 
#define _ions_typeeds( tp ) ( tp ).y


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
	float4 v = cell_v[_gate_v( sh )];			// load cells properties like: membrane potential, etc 
	float4 e = ions_e[_gate_e( sh )];			// load properties of ions like: resting potential, ions inside the cell, ions ouside of cell, etc
// todo: check the range for _gate_in_m( sh ) and _gate_in_m( sh ) indices for ions_e variable in debug mode [0...MaxIons]
	// load Mg- or Ca- concentration inside the cell for NMDA synapse or Z-channels if needed
	float in_m = ( _gate_typem( tp ) >= LSNS_ZGENERIC_INSTANT )? _ions_in( ions_e[_gate_in_m( sh )] ):0;
	// load Mg- or Ca- concentration inside the cell for NMDA synapse or Z-channels if needed
	float in_h = ( _gate_typeh( tp ) >= LSNS_ZGENERIC_INSTANT )? _ions_in( ions_e[_gate_in_h( sh )] ):0;
	//--- 2. perform calculations
	float eds = _ions_eds( e );				// extract resting potential from 'e'
	float vm = _cell_v( v );				// extract membrane potential from 'v'
	float mp, hp;
	proc_gate( _gate_typem( tp ), Gates[_gate_parm( tp )], in_m, vm, step, _gate_powm( mh ), _gate_m( mh ), mp );
	proc_gate( _gate_typeh( tp ), Gates[_gate_parh( tp )], in_h, vm, step, _gate_powh( mh ), _gate_h( mh ), hp );
	_gate_g( g ) = _gate_gmax( g )*mp*hp;		// g
	_gate_i( g ) = _gate_g( g )*( vm-eds );		// I
	_gate_ge( g ) = _gate_g( g )*eds;			// ge
	//--- 3.  store the results of simulation
	chan_g[index] = g;					// save properties of ions channel (conductance, current, etc)
	if( _gate_typem( tp )+_gate_typeh( tp ) != LSNS_NOGATE ){
		chan_mh[index] = mh;				// save properties of gate variables (activation, inactivation, etc) if needed
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
