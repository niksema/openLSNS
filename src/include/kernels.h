#ifndef __LSNS_KERNELS_H
#define __LSNS_KERNELS_H

#include "lsnsmath.h"
#include "engine.h"

#if defined( __LSNS_DEBUG__ )

__lsns_inline void __cellv_kernel( unsigned _index, float _step, float _threshold, int4 *_pchanglut, int4 *_pionilut, float4 *_pchang, float4 *_pioni, float4 *_pcellv )
{
	float4 v = (_pcellv)[_index];
	/* preprocessing */
	float g_sum = 0.f, ge_sum = 0.f, ipump = 0.f;
	/* sum all ipump */
	lsns_fastsum( _ions_ipump, (_pioni), (_pionilut)+(_index)*( MAX_IPUMP_PER_CELL/4+1 ), ipump, MAX_IONS );
	/* sum all g and ge */
	lsns_fastsum2( _chan_g, _chan_ge, (_pchang), (_pchanglut)+(_index)*( MAX_CHAN_PER_CELL/4+1 ), g_sum, ge_sum, MAX_CHAN );
	__lsns_assert( g_sum != 0.0 );							/* DEBUG: check if total conductance is not zero */
	/* perform calculations */
	float time = _cell_c( v )/g_sum;
	float v_inf = ( ge_sum+ipump+_cell_iadd( v ))/g_sum;
	float vm = lsns_exp_euler( _cell_v( v ), v_inf, (_step), time );
	float spike = ( vm > (_threshold) && _cell_v( v )  <= (_threshold) )? 1.f: 0.f;
	/* store the results of simulation */
	_cell_v( v ) = vm;
	_cell_spike( v ) = spike;
	(_pcellv)[_index] = v;
}

#else

#define __cellv_kernel( _index, _step, _threshold, _pchanglut, _pionilut, _pchang, _pioni, _pcellv ) \
	float4 v = (_pcellv)[_index];\
	/* preprocessing */\
	float g_sum = 0.f, ge_sum = 0.f, ipump = 0.f;\
	/* sum all ipump */\
	lsns_fastsum( _ions_ipump, (_pioni), (_pionilut)+(_index)*( MAX_IPUMP_PER_CELL/4+1 ), ipump, MAX_IONS );\
	/* sum all g and ge */\
	lsns_fastsum2( _chan_g, _chan_ge, (_pchang), (_pchanglut)+(_index)*( MAX_CHAN_PER_CELL/4+1 ), g_sum, ge_sum, MAX_CHAN );\
	__lsns_assert( g_sum != 0.0 );							/* DEBUG: check if total conductance is not zero */\
	/* perform calculations */\
	float time = _cell_c( v )/g_sum;\
	float v_inf = ( ge_sum+ipump+_cell_iadd( v ))/g_sum;\
	float vm = lsns_exp_euler( _cell_v( v ), v_inf, (_step), time );\
	float spike = ( vm > (_threshold) && _cell_v( v )  <= (_threshold) )? 1.f: 0.f;\
	/* store the results of simulation */\
	_cell_v( v ) = vm;\
	_cell_spike( v ) = spike;\
	(_pcellv)[_index] = v

#define __chan_kernel( _index, _step, _ptype, _plut, _pcellv, _pione, _pwsyn, _pchanmh, _pchang ) \
	/* load type of ion channel (generic, a-b, etc) and its parameters (half-voltage, slope, etc)*/\
	int4 tp = (_ptype)[_index];\
	/* load references to external parameters (membrane potential, rest potential, etc) */\
	int4 lut = (_plut)[_index];\
	__lsns_assert( _gate_typem( tp ) >= 0 && _gate_typem( tp ) < MAX_GATES );	/* DEBUG: check the range for _gate_typem( tp ) */\
	__lsns_assert( _gate_typeh( tp ) >= 0 && _gate_typeh( tp ) < MAX_GATES );	/* DEBUG: check the range for  _gate_typeh( tp ) */\
	__lsns_assert( _gate_parm( tp ) >= 0 && _gate_parm( tp ) < MAX_GPARS );		/* DEBUG: check the range for  _gate_parm( tp ) */\
	__lsns_assert( _gate_parh( tp ) >= 0 && _gate_parh( tp ) < MAX_GPARS );		/* DEBUG: check the range for  _gate_parh( tp ) */\
	__lsns_assert( _chan_lutv( lut ) >= 0 &&_chan_lutv( lut ) < MAX_CELLS );	/* DEBUG: check the range for _chan_lut_v( lut ) */\
	__lsns_assert( _chan_lute( lut ) >= 0 && _chan_lute( lut ) < MAX_IONS );	/* DEBUG: check the range for _chan_lut_e( lut ) */\
	__lsns_assert( _gate_typem( tp ) <= LSNS_PS_NMDA && _chan_lutm( lut ) >= 0 && _chan_lutm( lut ) < MAX_IONS ); /* DEBUG: check the range for _chan_lut_m( lut ) */\
	__lsns_assert( _gate_typeh( tp ) <= LSNS_PS_NMDA && _chan_luth( lut ) >= 0 && _chan_luth( lut ) < MAX_IONS ); /* DEBUG: check the range for _chan_lut_h( lut ) */\
	__lsns_assert( _gate_typem( tp ) > LSNS_PS_NMDA && _chan_lutm( lut ) >= 0 && _chan_lutm( lut ) < MAX_WSYNS ); /* DEBUG: check the range for _chan_lut_m( lut ) */\
	__lsns_assert( _gate_typeh( tp ) > LSNS_PS_NMDA && _chan_luth( lut ) >= 0 && _chan_luth( lut ) < MAX_WSYNS ); /* DEBUG: check the range for _chan_lut_h( lut ) */\
	/* load properties of ions channel (conductance, current, etc) */\
	float4 g = (_pghang)[_index];\
	 /* load properties of gate variables (activation, inactivation, etc) if needed */\
	float4 mh = ( _gate_typem( tp )+_gate_typeh( tp ) != LSNS_NOGATE)? (_pchanmh)[_index]: float4();\
/*todo: { possible CUDA optimization (try to use shared variables)
  todo: use __lsns_cached for pIonsE, pCellV */\
	/* load shared variables (resting potential, membrane potential, etc) to shared memory */\
	float eds = _ions_eds(( _pionse )[_chan_lute( lut )] );\
	float vm = _cell_v(( _pcellv )[_chan_lutv( lut )] );\
/*todo: } possible CUDA optimization */\
	float mp, hp;\
/*todo: keep time constants instead power in 'mh' variable.  powers for activation and inactivateion should be stored in Gates structure */\
	proc_gate( _gate_typem( tp ), Gates[_gate_parm( tp )], (_step), vm, _chan_lutm( lut ), (_pionse), (_pwsyn), _gate_powm( mh ), _gate_m( mh ), mp );\
	proc_gate( _gate_typeh( tp ), Gates[_gate_parh( tp )], (_step), vm, _chan_luth( lut ), (_pionse), (_pwsyn), _gate_powh( mh ), _gate_h( mh ), hp );\
	_chan_g( g ) = _chan_gmax( g )*mp*hp;\
	_chan_ge( g ) = _chan_g( g )*eds;\
	_chan_i( g ) = _chan_g( g )*( vm-eds );\
	/* save properties of ions channel (conductance, current, etc) */\
	pChanG[_index] = g;\
	/* save properties of gate variables (activation, inactivation, etc) if needed */\
	if( _gate_typem( tp )+_gate_typeh( tp ) != LSNS_NOGATE ){\
		pChanMH[_index] = mh;\
	}

#endif

#endif /*__LSNS_KERNELS_H*/
