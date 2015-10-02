#ifndef __GATEPROC_H
#define __GATEPROC_H

#include "lsnsmath.h"
#include "gates.h"

///////////////////////////////////////////////////////////////////////////////
// +implementation of some usefull macroses to calculate the gate variables
//=============================================================================
// generic gate variable: 1/( 1+exp(-(v-v12 )/slp ))
#define lsns_ggate( v, v12, slp ) \
	lsns_div( 1, 1+lsns_exp( -lsns_div(( v )-( v12 ), slp )))
// time constant for generic gate variable: t0+t/cosh( -( v-v12 )/slp ))
#define lsns_ggate_t( t0, t, v, v12, slp ) \
	( t0 )+lsns_div( t, lsns_cosh( -lsns_div(( v )-( v12 ), slp )))
// time constant for modified generic gate variable: t0+2*t/( exp(-( v-v12 )/slp )+exp(( v-v12_2 )/slp2 ))
#define lsns_ggate_tmod( t0, t, v, v12, slp, v12_2, slp_2 ) \
	( t0 )+lsns_div( 2*( t ), lsns_exp( -lsns_div(( v )-( v12 ), slp ))+lsns_exp( lsns_div(( v )-( v12_2 ), slp_2 )))
// alha/beta gate variable: a*( b*v-v12 )/( exp(( v-v12 )/slp )-c )
#define lsns_abgate( v, v12, slp, a, b, c ) \
	( a )*lsns_div(( b )*( v )-( v12 ), lsns_exp( lsns_div(( v )-( v12 ), slp ))-( c ))
///////////////////////////////////////////////////////////////////////////////
// +implementation of the gate variables of all types
//=============================================================================
///////////////////////////////////////////////////////////////////////////////
// +performes no calculation. Always returns 1.
__lsns_inline float proc_nogate( gatepar &par, float v, float step, float gate )
{
	return 1;
}
///////////////////////////////////////////////////////////////////////////////
// +bypass the gate calculation. Always returns 'gate'
__lsns_inline float proc_passgate( gatepar &par, float v, float step, float gate )
{
	return gate;
}
///////////////////////////////////////////////////////////////////////////////
// +implements generic description of gate variable [M/H] for voltage-activated 
// ion channel. T*d[M/H]/dt = [M/H]inf-[M/H]; 
//-----------------------------------------------------------------------------
// four subtypes of generic gate variable for voltage activated ion channel are:
// 	1) 'ggate1' instant (time constant = 0):
//		[M/H]inf = 1/(1+exp(-(V-V12)/Slope))
//		[M/H] = [M/H]inf;
//	2) 'ggate2' generic description for time constant !=0:
//		[M/H]inf = 1/(1+exp(-(V-V12)/Slope))
//		T = T0+Tmax/cosh((V-V12T)/SlopeT)
//		T*d[M/H]/dt = [M/H]inf-[M/H];
//	3) 'ggate3' modified generic description for time constant !=0 :
//		[M/H]inf = 1/(1+exp(-(V-V12)/Slope))
//		T = T0+2*Tmax/(exp((V-V12T)/SlopeT)+exp(-(V-V12T2)/SlopeT2))
//		T*d[M/H]/dt = [M/H]inf-[M/H];
//	4) 'ggate4' modified generic description for time constant !=0 (A-current):
//		[M/H]inf = 1/(1+exp(-(V-V12)/Slope))
//		T = Tup;
//		if( v < Vthreshold ){
//			T = T0+2*Tmax/(exp((V-V12T)/SlopeT)+exp(-(V-V12T2)/SlopeT2))
//		}
//		T*d[M/H]/dt = [M/H]inf-[M/H];
//=========================== ggate1 ==========================================
__lsns_inline float proc_ggate1( gatepar &par, float v, float step, float gate )
{
	return lsns_ggate( v, _ggatev12( par ), _ggateslp( par ));
}
//=========================== ggate2 ==========================================
__lsns_inline float proc_ggate2( gatepar &par, float v, float step, float gate )
{
	float gate_inf = lsns_ggate( v, _ggatev12( par ), _ggateslp( par ));
	float time = lsns_ggate_t( _ggatet0( par ), _ggatetmax( par ), v, _ggatev12t( par ), _ggateslpt( par ));
	return lsns_exp_euler( gate, gate_inf, step, time );
}
//=========================== ggate3 ==========================================
__lsns_inline float proc_ggate3( gatepar &par, float v, float step, float gate )
{
	float gate_inf = lsns_ggate( v, _ggatev12( par ), _ggateslp( par ));
	float time = lsns_ggate_tmod( _ggatet0( par ), _ggatetmax( par ), v, _ggatev12t( par ), _ggateslpt( par ), _ggatev12t2( par ), _ggateslpt2( par ));
	return lsns_exp_euler( gate, gate_inf, step, time );
}
//=========================== ggate4 ==========================================
__lsns_inline float proc_ggate4( gatepar &par, float v, float step, float gate )
{
	float gate_inf = lsns_ggate( v, _ggatev12( par ), _ggateslp( par ));
	float time = ( v < _ggatevtr( par ))? 
			lsns_ggate_tmod( _ggatet0( par ), _ggatetmax( par ), v, _ggatev12t( par ), _ggateslpt( par ), _ggatev12t2( par ), _ggateslpt2( par )):
			_ggatetup( par ); // need to be optimized;
	return lsns_exp_euler( gate, gate_inf, step, time );
}
///////////////////////////////////////////////////////////////////////////////
// +implements alpa/beta description of gate variable [M/H] for voltage-
// activated ion channel. T*d[M/H]/dt = Alpha*(1-[M/H])-Beta*[M/H];
//-----------------------------------------------------------------------------
// two subtypes of alpha/beta gate variable for voltage activated ion channel are:
// 	1) 'abgate1' instant (time constant is 0):
//		[M/H] = Alpha/(Alpha+Beta);
//	2) 'abgate2' alpha/beta description for time constant:
//		T = T0+Tmax/(Alpha+Beta);
//		[M/H]inf = Alpha/(Alpha+Beta);
//		T*d[M/H]/dt = [M/H]inf-[M/H];
//-----------------------------------------------------------------------------
//	Alpha/Beta = A*( B*V-V12 )/( C+exp( -( v-V12 )/Slp ));
//=========================== abgate1 =========================================
__lsns_inline float proc_abgate1( gatepar &par, float v, float step, float gate )
{
	float alpha = lsns_abgate( v, _abgatev12a( par ), _abgateslpa( par ), _abgateAa( par ), _abgateBa( par ), _abgateCa( par ));
	float beta = lsns_abgate( v, _abgatev12b( par ), _abgateslpb( par ), _abgateAb( par ), _abgateBb( par ), _abgateCb( par ));
	return lsns_div( alpha, alpha+beta );
}
//=========================== abgate2 =========================================
__lsns_inline float proc_abgate2( gatepar &par, float v, float step, float gate )
{
	float alpha = lsns_abgate( v, _abgatev12a( par ), _abgateslpa( par ), _abgateAa( par ), _abgateBa( par ), _abgateCa( par ));
	float beta = lsns_abgate( v, _abgatev12b( par ), _abgateslpb( par ), _abgateAb( par ), _abgateBb( par ), _abgateCb( par ));
	float time = lsns_div( 1.f, alpha+beta );
	float gate_inf = alpha*time;
	time = _abgatet0( par )+_abgatetmax( par )*time;
	return lsns_exp_euler( gate, gate_inf, step, time );
}
///////////////////////////////////////////////////////////////////////////////
// +implements z-description of gate variable [M/H] for calcium-
// activated ion channel. T*d[M/H]/dt = [M/H]inf-[M/H];
//-----------------------------------------------------------------------------
// four subtypes of z-gate variable for voltage activated ion channel are:
// 	1) 'zgate1' instant generic description(time constant is 0):
//		k = Alpha*(Beta*[Ca_in]^Lymbda];
//		[M/H] = k/(1+k);
//	2) 'zgate2' z-gate generic description for time constant:
//		k = Alpha*(Beta*[Ca_in]^Lymbda];
//		[M/H]inf = 1/(1+k);
//		T = T0+Tmax/(1+Gamma*k);
//		T*d[M/H]/dt = [M/H]inf-[M/H];
// 	3) 'zgate3' instant alpha/beta description (time constant is 0):
//		alpha = [Ca_in]*Alpha*exp( -( v-V12 )/Slp );
//		beta = Beta*exp(( v-V12 )/Slp );
//		[M/H] = alpha/(alpha+beta);
//	4) 'zgate4' z-gate alpha/beta description for time constant:
//		alpha = [Ca_in]*Alpha*exp( -( v-V12 )/Slp );
//		beta = Beta*exp(( v-V12 )/Slp );
//		[M/H]inf = alpha/(alpha+beta);
//		T = T0+Tmax/(alpha+beta);
//		T*d[M/H]/dt = [M/H]inf-[M/H];
//=========================== zgate1 ==========================================
__lsns_inline float proc_zgate1( gatepar &par, float in, float v, float step, float gate )
{
	float alpha = _zgateA( par )*lsns_pow( _zgateB( par )*in, _zgateL( par ));
	return lsns_div( alpha, 1.f+alpha );
}
//=========================== zgate2 ==========================================
__lsns_inline float proc_zgate2( gatepar &par, float in, float v, float step, float gate )
{
	float alpha = _zgateA( par )*lsns_pow( _zgateB( par )*in, _zgateL( par ));
	float gate_inf = lsns_div( alpha, 1.f+alpha );
	float time = _zgatet0( par )+lsns_div( _zgatetmax( par ), 1.f+_zgateG( par )*alpha );
	return lsns_exp_euler( gate, gate_inf, step, time );
}
//=========================== zgate3 ==========================================
__lsns_inline float proc_zgate3( gatepar &par, float in, float v, float step, float gate )
{
	float alpha = in*_zgateA( par )*lsns_exp( lsns_div( -( v-_zgatev12( par )), _zgateslp( par )));
	float beta = _zgateB( par )*lsns_exp( lsns_div( v-_zgatev12( par ), _zgateslp( par )));
	return lsns_div( alpha, alpha+beta );
}
//=========================== zgate4 ==========================================
__lsns_inline float proc_zgate4( gatepar &par, float in, float v, float step, float gate )
{
	float alpha = in*_zgateA( par )*lsns_exp( lsns_div( -( v-_zgatev12( par )), _zgateslp( par )));
	float beta = _zgateB( par )*lsns_exp( lsns_div( v-_zgatev12( par ), _zgateslp( par )));
	float time = lsns_div( 1.f, alpha+beta );
	float gate_inf = alpha*time;
	time = _zgatet0( par )+_zgatetmax( par )*time;
	return lsns_exp_euler( gate, gate_inf, step, time );
}
///////////////////////////////////////////////////////////////////////////////
// +implements synaptic plasticity on post-synaptic neuron as gate variable
// that modulates the synaptic conductance.
//-----------------------------------------------------------------------------
// one subtype of gate variable for synaptic plasticity on post-synaptic neuron:
// 	1) 'psgate1' nmda synapse:
//		[M/H] = 1/(1+exp(−0.062*V)*[Mg]/3.57)
//=========================== psgate1 =========================================
__lsns_inline float proc_psgate1( gatepar &par, float in, float v, float step, float gate )
{
	return lsns_div( 1.f, 1.f+lsns_exp( -0.062f*v )*lsns_div( in, 3.57f ));
}

///////////////////////////////////////////////////////////////////////////////
// implements synaptic summation of gate variable [M/H]  for synaptic current
// [M/H] = [M/H]*Edt+A*h*w*Dt, where:
//	A is rate of transmitter release,
//	h is synaptic plasticity,
//	w is total sum of all pre-synaptic neurons converged on the synapse
//-----------------------------------------------------------------------------
// 	1) 'proc_syngate1' pulse synapse:
//		Edt = exp( step/T ) where T is time constant, Dt = 1 
//=========================== psgate1 =========================================
__lsns_inline float proc_syngate1( float w, float alpha, float edt, float dt, float gate )
{
	return gate*edt+w*alpha*dt;
}

///////////////////////////////////////////////////////////////////////////////
// +macros to calculate gate variable of any type
// BTW. It's funny but the current imlementation is faster then the one that
// uses the callback function
//-----------------------------------------------------------------------------
// 'type' is channel's type;
// 'par' are channels' parameters;
// 'step' is step of integration
// 'v' is membrane potential;
// 'lut' is either index for 'ions' aray or index for 'wsyn' array;
// 'ions' is array of ions' parameters of Ca- or Mg- ions which are used for NMDA synapse or Z-channels correspondingly;
// 'mod' is power of gate variable for channels
// 'gate' is gate variable;
// 'G' is 'gate'^'mod'  for channels or 'gate' for synapses
#define proc_gate( type, par, step, v, lut, ions, wsyn, mod, gate, G ) \
	switch( type ){ \
		case LSNS_NOGATE: \
			( gate ) = ( G ) = 1; \
			break; \
		case LSNS_BYPASSGATE: \
			( G ) = lsns_pow( gate, mod ); \
			break; \
		case LSNS_GENERIC_INSTANT: \
			( gate ) = proc_ggate1( par, v, step, gate ); \
			( G ) = lsns_pow( gate, mod ); \
			break; \
		case LSNS_GENERIC_T: \
			( gate ) = proc_ggate2( par, v, step, gate ); \
			( G ) = lsns_pow( gate, mod ); \
			break; \
		case LSNS_GENERIC_TMOD: \
			( gate ) = proc_ggate3( par, v, step, gate ); \
			( G ) = lsns_pow( gate, mod ); \
			break; \
		case LSNS_GENERIC_TAMOD: \
			( gate ) = proc_ggate4( par, v, step, gate ); \
			( G ) = lsns_pow( gate, mod ); \
			break; \
		case LSNS_ALPHABETA_INSTANT: \
			( gate ) = proc_abgate1( par, v, step, gate ); \
			( G ) = lsns_pow( gate, mod ); \
			break; \
		case LSNS_ALPHABETA_T: \
			( gate ) = proc_abgate2( par, v, step, gate ); \
			( G ) = lsns_pow( gate, mod ); \
			break; \
		case LSNS_ZGENERIC_INSTANT: \
			/* _ions_in( ions[lut] loads Ca- concentration inside the cell for Z-channels */\
			( gate ) = proc_zgate1( par, _ions_in( ions[lut] ), v, step, gate ); \
			( G ) = lsns_pow( gate, mod ); \
			break; \
		case LSNS_ZGENERIC_T: \
			/* _ions_in( ions[lut] loads Ca- concentration inside the cell for Z-channels */\
			( gate ) = proc_zgate2( par, _ions_in( ions[lut] ), v, step, gate ); \
			( G ) = lsns_pow( gate, mod ); \
			break; \
		case LSNS_ZAPHABETA_INSTANT: \
			/* _ions_in( ions[lut] loads Ca- concentration inside the cell for Z-channels */\
			( gate ) = proc_zgate3( par, _ions_in( ions[lut] ), v, step, gate ); \
			( G ) = lsns_pow( gate, mod ); \
			break; \
		case LSNS_ZAPHABETA_T: \
			/* _ions_in( ions[lut] loads Ca- concentration inside the cell for Z-channels */\
			( gate ) = proc_zgate4( par, _ions_in( ions[lut] ), v, step, gate ); \
			( G ) = lsns_pow( gate, mod ); \
			break; \
		case LSNS_PS_NMDA: \
			/* _ions_in( ions[lut] loads Mg- concentration inside the cell for NMDA synapse*/\
			( gate ) = proc_psgate1( par, _ions_in( ions[lut] ), v, step, gate ); \
			( G ) = lsns_pow( gate, mod ); \
			break; \
		case LSNS_SYNAPSE: \
			{\
			/*loads total weight of connections converged on the synapse*/\
			float4 w = wsyn[lut]; \
			( G ) = ( gate ) = proc_syngate1( _wsyn_total( w ), _wsyn_ah( w ), _wsyn_edt( w ), _wsyn_dt( w ), gate ); \
			}\
			break; \
		default: \
			G = gate = 0;\
	} 
	
#endif /*__GATEPROC_H*/
