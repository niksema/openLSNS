#include "precompile.h"

#include "gateproc.h"

const char *lsns_gate_types[LSNS_MAX_GATE] = {
	"None",
	"Bypass",
	"Generic instant",
	"Generic standard",
	"Generic modified",
	"Generic A-current",
	"Alpha-Beta instant",
	"Alpha-Beta standard",
	"Z-current instant",
	"Z-current standard",
	"Z-current (AB) instant",
	"Z-current (AB) standard",
	"NMDA post-synapse",
};

/*
#include <assert.h> 
#include <string>
#include <vector>

using std::string;
using std::vector;

///////////////////////////////////////////////////////////////////////////////
// Interface to specifie the parameters of a gate variable
//=============================================================================
class iogate{
	public: // constructors/destructor
		iogate( void );
		iogate( int type );
		iogate( const string &type );
		iogate( int type, const string &description );
		iogate( const string &type, const string &description );
		iogate( const iogate &gate );
		~iogate( void );
	public: // operators
		iogate &operator = ( const iogate &gate );
		float &operator[]( const string &name ); // get/set x 
		bool operator >> ( string &description ); // export to ASCII text in format "name1=a\nname2=b\n..."
		bool operator << ( const string &description ); // import from ASCII text in format "name1=a\nname2=b\n..."
	public: // methods
		bool create( int type );
		bool create( const string &type );
		bool cleate( int type, const string &description );
		bool cleate( const string &type, const string &description );
		string get_typename( void ) const{
			assert( Type >= 0 && Type < LSNS_MAX_GATE );
			return string( lsns_gate_types[Type] );
		};
		int get_type( void ) const{
			return Type; 
		};
		bool get_x( const string &name, float &x ) const;
		bool set_x( const string &name, float x );
		bool get_names( vector<string> &Names );
	public: // data
		int Type;
		string Gdescr;
		gate_par Gpar;
};

const char *lsns_gate_description[LSNS_MAX_GATE] = {
	"Function \'calc_nogate\'.\n\
	Performes no calculation. Always returns 1",

	"Function \'calc_passgate\'.\n\
	Bypass the gate calculation.",

	"Function \'calc_ggate1\'.\n\
	Generic description of gate variable [M/H] (time constant = 0)\n\
	[M/H] = 1/(1+exp(-(V-V12)/Slope))",

	"Function \'calc_ggate2\'.\n\
	Generic description of gate variable [M/H] (time constant != 0)\n\
	[M/H]inf = 1/(1+exp(-(V-V12)/Slope))\
	T = T0+Tmax/cosh((V-V12T)/SlopeT)\
	T*d[M/H]/dt = [M/H]inf-[M/H]",

	"Function \'calc_ggate3\'.\n\
	Modified generic description of gate variable [M/H] (time constant != 0)\n\
	[M/H]inf = 1/(1+exp(-(V-V12)/Slope))\n\
	T = T0+2*Tmax/(exp((V-V12T)/SlopeT)+exp(-(V-V12T2)/SlopeT2))\n\
	T*d[M/H]/dt = [M/H]inf-[M/H]\n",

	"Function \'calc_ggate4\'.\n\
	Modified generic description of gate variable for A-current [M/H] (time constant != 0)\n\
	[M/H]inf = 1/(1+exp(-(V-V12)/Slope))\n\
	T = Tup; if( v < Vthreshold )\n\
	\tT = T0+2*Tmax/(exp((V-V12T)/SlopeT)+exp(-(V-V12T2)/SlopeT2))\n\
	T*d[M/H]/dt = [M/H]inf-[M/H]",

	"Function \'calc_abgate1\'.\n\
	Alpha/Beta description of gate variable [M/H] (time constant = 0)\n\
	[M/H] = Alpha/(Alpha+Beta)",

	"Function \'calc_abgate2\'.\n\
	Alpha/Beta description of gate variable [M/H] (time constant != 0)\n\
	[M/H]inf = Alpha/(Alpha+Beta)\
	T = T0+Tmax/(Alpha+Beta)\
	T*d[M/H]/dt = [M/H]inf-[M/H]",

	"Function \'calc_zgate1\'.\n\
	Generic description of [Ca] activated gate variable [M/H] (time constant = 0)\n\
	k = Alpha*(Beta*[Ca]^Lymbda]\n\
	[M/H] = k/(1+k)",

	"Function \'calc_zgate2\'.\n\
	Generic description of [Ca] activated gate variable [M/H] (time constant != 0)\n\
	k = Alpha*(Beta*[Ca]^Lymbda]\n\
	[M/H] = k/(1+k)\n\
	T = T0+Tmax/(1+Gamma*k)\n\
	T*d[M/H]/dt = [M/H]inf-[M/H]",

	"Function \'calc_zgate3\'.\n\
	Alpha-Beta description of [Ca] activated gate variable [M/H] (time constant = 0)\n\
	alpha = [Ca]*Alpha*exp( -( v-V12 )/Slp )\n\
	beta = Beta*exp(( v-V12 )/Slp )\n\
	[M/H] = alpha/(alpha+beta)",

	"Function \'calc_zgate4\'.\n\
	Alpha-Beta description of [Ca] activated gate variable [M/H] (time constant != 0)\n\
	alpha = [Ca_in]*Alpha*exp( -( v-V12 )/Slp )\n\
	beta = Beta*exp(( v-V12 )/Slp )\n\
	[M/H]inf = alpha/(alpha+beta)\n\
	T = T0+Tmax/(alpha+beta)\n\
	T*d[M/H]/dt = [M/H]inf-[M/H]",

	"Function \'calc_psgate1\'.\n\
	Implementation of synaptic plasticity on post-synaptic neuron (nmda synapse)\n\
	[M/H] = 1/(1+exp(−0.062*V)*[Mg]/3.57)",
};
*/
