#ifndef __LOAD_NSM_H
#define __LOAD_NSM_H

#include "utilities.h"
#include "netunit.h"

/////////////////////////////////////////////////////////////////////////////
// class simdata
class simdata{
	public:
		simdata( void ) : Promt(), Seed(0), SimStep(0.05), SimTime(1000), UpdatingTime(10), 
			        Freq(0.3), HistBin(100), HistNorm(1000), BeginView(0), EndView(1000), TimeFactor(20){};
		simdata( const simdata &data ) : Promt(data.Promt), Seed(data.Seed), SimStep(data.SimStep), SimTime(data.SimTime), UpdatingTime(data.UpdatingTime), 
			        Freq(data.Freq), HistBin(data.HistBin), HistNorm(data.HistNorm), BeginView(data.BeginView), EndView(data.EndView), TimeFactor(data.TimeFactor){};
		~simdata( void ){};
	public:
		simdata &operator = ( const simdata &data );
	public:
		bool loadpar( istream &file, const string &parname );
		bool validate( const string &type );
	public:
		long Seed;
		// parameters of a simulation proc.
		double SimStep;			//
		double SimTime;			// msec
		unsigned long UpdatingTime;	// msec
		// parameters of visualization proc.
		double Freq;			// kHz 
		double HistBin;			// msec
		double HistNorm;		// msec
		double BeginView;		// msec
		double EndView;			// msec
		double TimeFactor;		// 
		string Promt;			// 
};

/////////////////////////////////////////////////////////////////////////////
// class iondata
class iondata{
	public:
		iondata( void ) : RTF( 26.54f ), Z( 1 ), In( .1f ), In0( 0.f ), Out( .1f ), T( 1.f ), IType( 0 ), 
			PType( 0 ), Eds( 0.f ), B( 0.03f ), K( 0.001f ), Rpmp( 1.f ), Rch( 1.f ), Kp( 15.f ){};
		iondata( const iondata &ion ) : RTF( 26.54f ), Z( ion.Z ), In( ion.In ), In0( ion.In0 ), Out( ion.Out ), T( ion.T ), IType( ion.IType ),
			PType( ion.PType ), Eds( ion.Eds ), B( ion.B ), K( ion.K ), Rpmp( ion.Rpmp ), Rch( ion.Rch ), Kp( ion.Kp ){};
		~iondata( void ){};
/*
//				nns<<"Rpmp = "<<pops[i].paraset.ions.ion_para[j].Rpmp<<endl;
//				nns<<"Rch = "<<pops[i].paraset.ions.ion_para[j].Rch<<endl;
//				nns << "Kp = " << Network().Ions[j]().K << endl;

*/
	public:
		iondata &operator = ( const iondata &ion );
	public:
		bool loadpar( istream &file, const string &parname );
		bool validate( const string &type );
		bool save( ostream &file );
	public:
		int IType;
		int PType;
		int Z;
		float RTF;
		hhnpair<float> In;
		hhnpair<float> In0;
		hhnpair<float> Out;
		hhnpair<float> Eds;

		float Rpmp;
		float Rch;
		float Kp;
		float T;
		float B;		// Ca specific
		float K;		// Ca specific
};

/////////////////////////////////////////////////////////////////////////////
// class syndata
class syndata{
	public:
		syndata( void ) : Eds( "0.0" ), Gmax( 0.f ), Ts( 25.f ), Vs( 0.5f ) {};
		syndata( const syndata &syn ) : Eds( syn.Eds ), Gmax( syn.Gmax ), Ts( syn.Ts ), Vs( syn.Vs ){};
		~syndata( void ){};
	public:
		syndata &operator = ( const syndata &syn );
	public:
		bool loadpar( istream &file, const string &parname );
		bool validate( const string &type );
	public:
		float Gmax;		// maximal conductance
		float Ts;		// time constant
		float Vs;		// amplitude of sinaptic stimulation
		string Eds;		// reversal potential
};

/////////////////////////////////////////////////////////////////////////////
// class dsyndata
class dsyndata{
	public:
		dsyndata( void ) : Kdep( 1.f ), Tdep( 500.f ){};
		dsyndata( const dsyndata &dsyn ) : Kdep( dsyn.Kdep ), Tdep( dsyn.Tdep ){};
		~dsyndata( void ){};
	public:
		dsyndata &operator = ( const dsyndata &dsyn );
	public:
		bool loadpar( istream &file, const string &parname );
		bool validate( const string &type );
	public:
		float Kdep;		// amplitude 
		float Tdep;		// time constant
};

/////////////////////////////////////////////////////////////////////////////
// class ssyndata
class ssyndata{
	public:
		ssyndata( void ) : Hv( 1.f ), Slp( 500.f ){};
		ssyndata( const ssyndata &ssyn ) : Hv( ssyn.Hv ), Slp( ssyn.Slp ){};
		~ssyndata( void ){};
	public:
		ssyndata &operator = ( const ssyndata &ssyn );
	public:
		bool loadpar( istream &file, const string &parname );
		bool validate( const string &type );
	public:
		float Hv;		// half-voltage
		float Slp;		// slope
};

/////////////////////////////////////////////////////////////////////////////
// class chandata
class chandata{
	public:
		chandata( void ) : M( 1.f, 0.f ), H( 1.f, 0.f ), Gmax( 1.f, 0.f ), Eds( -1.f, -1.f ), Pcl( 0.45f, 0.f ), Pna( 0.05f, 0.f ), Vg( "false" ){}; //(auto) Eds( -1.f, -1.f ) specific reversal potential, related to ion dynamics 
		chandata( const chandata &chdat ) : M( chdat.M ), H( chdat.H ), Gmax( chdat.Gmax ), Eds( chdat.Eds ), Pcl( chdat.Pcl ), Pna( chdat.Pna ), Tm( chdat.Tm ), Th( chdat.Th ), 
			PowM( chdat.PowM ), PowH( chdat.PowH ), V12m( chdat.V12m ), V12h( chdat.V12h ), Slpm( chdat.Slpm ), Slph( chdat.Slph ),	V12tm( chdat.V12tm ), V12th( chdat.V12th ), 
			Slptm( chdat.Slptm ), Slpth( chdat.Slpth ), Vg( chdat.Vg ){};
		~chandata( void ){};
	public:
		chandata &operator = ( const chandata &chdat );
	public:
		bool loadpar( istream &file, const string &parname );
		bool validate( const string &type );
	public:
		hhnpair<float> M;
		hhnpair<float> H;
		hhnpair<float> Tm;
		hhnpair<float> Th;
		hhnpair<float> PowM;
		hhnpair<float> PowH;
		hhnpair<float> V12m;
		hhnpair<float> V12h;
		hhnpair<float> Slpm;
		hhnpair<float> Slph;

		hhnpair<float> V12tm;
		hhnpair<float> V12th;
		hhnpair<float> Slptm;
		hhnpair<float> Slpth;

		hhnpair<float> Gmax;
		hhnpair<float> Eds;

		hhnpair<float> Pcl;		// leak related
		hhnpair<float> Pna;		// leak related
		string Vg;			// goldman reversal potential 
};

/////////////////////////////////////////////////////////////////////////////
// class cmpdata
class cmpdata{
	public:
		cmpdata( void ): V( -60.f, 0.f ){};
		cmpdata( const cmpdata &cdat ) : S( cdat.S ), V( cdat.V ), Chans( cdat.Chans ){};
		~cmpdata( void ){};
	public:
		cmpdata &operator = ( const cmpdata &cdat );
	public:
		bool loadpar( istream &file, const string &parname );
		bool validate( const string &type );
	public:
		hhnpair<float> S;
		hhnpair<float> V;
		vector<netunit<chandata >> Chans;
};

/////////////////////////////////////////////////////////////////////////////
// class hhndata
class hhndata{
	public:
		hhndata( void ) : C( 1.f, 0.f ){};
		hhndata( const hhndata &hhn ) : C( hhn.C ), P( hhn.P ), Ex( hhn.Ex ), Ex2( hhn.Ex2 ), InhI( hhn.InhI ), InhII( hhn.InhII ), Cmps( hhn.Cmps ){};
		~hhndata( void ){};
	public:
		hhndata &operator = ( const hhndata &hhn );
	public:
		bool loadpar( istream &file, const string &parname );
		bool validate( const string &type );
	public:
		hhnpair<float> C;
		hhnpair<float> P;
		hhnpair<float> Ex;
		hhnpair<float> Ex2;
		hhnpair<float> InhI;
		hhnpair<float> InhII;
		vector<netunit<cmpdata >> Cmps;
};

/////////////////////////////////////////////////////////////////////////////
// class popdata
class popdata{
	public:
		popdata( void ) : Size( 0 ){};
		popdata( const popdata &pdat ) : Size( pdat.Size ), Hhn( pdat.Hhn ){};
		~popdata( void ){};
	public:
		popdata &operator = ( const popdata &pdat );
	public:
		bool loadpar( istream &file, const string &parname );
		bool validate( const string &type );
	public:
		size_t Size;
		netunit<hhndata> Hhn;
};

/////////////////////////////////////////////////////////////////////////////
// class drivedata
class drivedata{
	public:
		drivedata( void ){};
		drivedata( const drivedata &ddat ) : Size( ddat.Size ){};
		~drivedata( void ){};
	public:
		drivedata &operator = ( const drivedata &ddat );
	public:
		bool loadpar( istream &file, const string &parname );
		bool validate( const string &type );
	public:
		size_t Size;
};

/////////////////////////////////////////////////////////////////////////////
// class outdata
class outdata{
	public:
		outdata( void ){};
		outdata( const outdata &odat ) : Size( odat.Size ), Tup( odat.Tup ), Tdown( odat.Tdown ), Threshold( odat.Threshold ){};
		~outdata( void ){};
	public:
		outdata &operator = ( const outdata &odat );
	public:
		bool loadpar( istream &file, const string &parname );
		bool validate( const string &type );
	public:
		size_t Size;
		float Tup;
		float Tdown;
		float Threshold;
};

/////////////////////////////////////////////////////////////////////////////
// class unitdata
class unitdata{
	public:
		unitdata( void ){};
		unitdata( const unitdata &udat ) : PData( udat.PData ), DData( udat.DData ), OData( udat.OData ){};
		~unitdata( void ){};
	public:
		unitdata &operator = ( const unitdata &udat );
	public:
		bool loadpar( istream &file, const string &parname );
		bool validate( const string &type );
	public:
		vector<netunit<popdata >> PData;	// list of populations
		vector<netunit<drivedata >> DData;	// list of drives
		vector<netunit<outdata >> OData;	// list of outputs
};

/////////////////////////////////////////////////////////////////////////////
// class wdata
class wdata{
	public:
		wdata( void ){};
		wdata( const wdata &wdat ) : Type( wdat.Type ), Prob( wdat.Prob ), Ex( wdat.Ex ), Ex2( wdat.Ex2 ), InhI( wdat.InhI ), InhII( wdat.InhII ){};
		~wdata( void ){};
	public:
		wdata &operator = ( const wdata &wdat );
	public:
		bool loadpar( istream &file, const string &parname );
		bool validate( const string &type );
	public:
		int Type;
		float Prob;
		hhnpair<float> Ex;
		hhnpair<float> Ex2;
		hhnpair<float> InhI;
		hhnpair<float> InhII;
};

/////////////////////////////////////////////////////////////////////////////
// class trgdata
class trgdata{
	public:
		trgdata( void ){};
		trgdata( const trgdata &tdat ) : WData( tdat.WData ){};
		~trgdata( void ){};
	public:
		trgdata &operator = ( const trgdata &sdat );
	public:
		bool loadpar( istream &file, const string &parname );
		bool validate( const string &type );
	public:
		vector<netunit<wdata >> WData;		// weight of connection
};

/////////////////////////////////////////////////////////////////////////////
// class srcdata
class srcdata{
	public:
		srcdata( void ){};
		srcdata( const srcdata &sdat ) : TrgData( sdat.TrgData ){};
		~srcdata( void ){};
	public:
		srcdata &operator = ( const srcdata &sdat );
	public:
		bool loadpar( istream &file, const string &parname );
		bool validate( const string &type );
	public:
		vector<netunit<trgdata >> TrgData;	// list of targets
};

/////////////////////////////////////////////////////////////////////////////
// class condata
class condata{
	public:
		condata( void ){};
		condata( const condata &cndat ) : SrcData( cndat.SrcData ){};
		~condata( void ){};
	public:
		condata &operator = ( const condata &cndat );
	public:
		bool loadpar( istream &file, const string &parname );
		bool validate( const string &type );
	public:
		vector<netunit<srcdata >> SrcData;	// list of sources
};

/////////////////////////////////////////////////////////////////////////////
// class ctrdata
class ctrdata{
	public:
		ctrdata( void ){};
		ctrdata( const ctrdata &ctdat ){};
		~ctrdata( void ){};
	public:
		ctrdata &operator = ( const ctrdata &ctdat );
	public:
		bool loadpar( istream &file, const string &parname );
		bool validate( const string &type );
	public:
};

/////////////////////////////////////////////////////////////////////////////
// class netdata
class netdata{
	public:
		netdata( void ) : Threshold( 0. ){};
		netdata( const netdata &net ) : Threshold( net.Threshold ), Ions( net.Ions ), Syns( net.Syns ), DSyns( net.DSyns ), SSyns( net.SSyns), Units( net.Units ){};
		~netdata( void ){};
	public:
		netdata &operator = ( const netdata &net );
	public:
		bool loadpar( istream &file, const string &parname );
		bool validate( const string &type );
	private:
		void clear( void );
	public:
		float Threshold;
		vector<netunit<iondata >> Ions;
		vector<netunit<syndata >> Syns;
		vector<netunit<dsyndata >> DSyns;
		vector<netunit<ssyndata >> SSyns;
		netunit<unitdata> Units;

		vector<netunit<condata >> Connections;
		vector<netunit<ctrdata >> Controls;
};

/////////////////////////////////////////////////////////////////////////////
// class nsm_model
class nsm_model{
	public:
		nsm_model( void ) : SimData(){};
		nsm_model( const nsm_model &model ) : SimData( model.SimData ){};
		~nsm_model( void ){};
	public:
		bool load_model( const char *filename );
		bool save_model( const char *filename );
	public:
		netunit<simdata> SimData;	// common parameters of simulation
		netunit<netdata> Network;	// network structure and parameters
};
#endif /*__LOAD_NSM_H*/
