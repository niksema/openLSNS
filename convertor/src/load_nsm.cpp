#include "stdafx.h"

#include "load_nsm.h"

float eds( float rtfz, float out, float in )
{
	return rtfz*log( out/in );
}

hhnpair<float> hack1( hhnpair<float>& e, const string &name, vector<netunit<iondata >> &ions )
{
	if( name[0] == 'N' && name[1] == 'a' ){
		for( size_t i = 0; i < ions.size(); ++i ){
			if( ions[i].Name == "Na" ){
				e = ions[i]().Eds;
			}
		}
	}
	if( name[0] == 'K' ){
		for( size_t i = 0; i < ions.size(); ++i ){
			if( ions[i].Name == "K" ){
				e = ions[i]().Eds;
			}
		}
	}
	if( name[0] == 'C' && name[1] == 'a' ){
		for( size_t i = 0; i < ions.size(); ++i ){
			if( ions[i].Name == "Ca" ){
				e = ions[i]().Eds;
			}
		}
	}
	return e;
}

bool hack2( const string &name, vector<netunit<chandata >> &chans )
{
	if( name[0] == 'C' && name[1] == 'a' ){
		for( size_t i = 0; i < chans.size(); ++i ){
			if( chans[i].Name.find( "Ca" ) != string::npos ){
				return true;
			}
		}
	}
	return false;
}

/////////////////////////////////////////////////////////////////////////////
// class simdata
//--- constructors/destructor
//--- public methods
simdata &simdata::operator = ( const simdata &data )
{
	Seed = data.Seed;
	SimStep = data.SimStep;
	SimTime  = data.SimTime;
	UpdatingTime = data.UpdatingTime;
	Freq = data.Freq;
	HistBin = data.HistBin;
	HistNorm = data.HistNorm;
	BeginView = data.BeginView;
	EndView = data.EndView;
	TimeFactor = data.TimeFactor;
	Promt = data.Promt;
	return *this;
};

bool simdata::loadpar( istream &file, const string &parname )
{
	string str;
	if( parname == "Seed"){
		file >> str >> Seed;
		return true;
	} 
	else if( parname == "SimStep"){
		file >> str >> SimStep;
		return true;
	}
	else if( parname == "SimNNStep"){
		file >> str >> SimStep;
		return true;
	}
	else if( parname == "SimTime"){
		file >> str >> SimTime;
		return true;
	}
	else if( parname == "UpdatingTime"){
		file >> str >> UpdatingTime;
		return true;
	}
	else if( parname == "Freq"){
		file >> str >> Freq;
		return true;
	}
	else if( parname == "HistBin"){
		file >> str >> HistBin;
		return true;
	}
	else if( parname == "HistNorm"){
		file >> str >> HistNorm;
		return true;
	}
	else if( parname == "BeginView"){
		file >> str >> BeginView;
		return true;
	}
	else if( parname == "EndView"){
		file >> str >> EndView;
		return true;
	}
	else if( parname == "TimeFactor"){
		file >> str >> TimeFactor;
		return true;
	}
	return false;
}

bool simdata::validate( const string &type )
{
	return true;
}

/////////////////////////////////////////////////////////////////////////////
// class iondata
//--- constructors/destructor
//--- public methods
iondata &iondata::operator = ( const iondata &ion )
{
	In = ion.In;
	In0 = ion.In0;
	Out = ion.Out;
	T = ion.T;
	Eds = ion.Eds;
	Kp = ion.Kp;
	Rpmp = ion.Rpmp;
	Rch = ion.Rch;
	B = ion.B;
	K = ion.K;
	IType = ion.IType;
	PType = ion.PType;
	return *this;
};

bool iondata::loadpar( istream &file, const string &parname )
{
	string str;
	if( parname == "In"){
		file >> str >> In;
		return true;
	}
	else if( parname == "In0"){
		file >> str >> In0;
		return true;
	}
	else if( parname == "Rpmp"){
		file >> str >> Rpmp;
		return true;
	}
	else if( parname == "Rch"){
		file >> str >> Rch;
		return true;
	}
	else if( parname == "Out"){
		file >> str >> Out;
		return true;
	}
	else if( parname == "T"){
		file >> str >> T;
		return true;
	}
	else if( parname == "Eds"){
		file >> str >> Eds;
		return true;
	}
	else if( parname == "B"){ // Ca dynamics
		file >> str >> B;
		return true;
	}
	else if( parname == "K"){ // Ca dynamics
		file >> str >> K;
		return true;
	}
	return false;
}

bool iondata::validate( const string &type )
{
	if( type == "Na" ){
		IType = 1;
		PType = 0;
		Z = 1;
		Eds.X = eds( RTF/Z, Out.X, In.X );
	}
	else if( type == "K" ){
		IType = 2;
		PType = 0;
		Z = 1;
		Eds.X = eds( RTF/Z, Out.X, In.X );
	}
	else if( type == "Ca" ){
		IType = 3;
		PType = 3;
		Z = 2;
		Eds.X = eds( RTF/Z, Out.X, In.X );
	}
	else if( type == "Cl" ){
		IType = 4;
		PType = 0;
		Z = -1;
		Eds.X = eds( RTF/Z, Out.X, In.X );
	}
	else{
		IType = 0;
		PType = 0;
		Z = 1;
		string mess = "unknown ions type:" + type;
		message( mess.c_str(), "Warning" );
		return false;
	}
	return true;
}

bool iondata::save( ostream &file )
{
	file << "IonType = " << IType << endl;
	file << "PumpType = " << PType << endl;
	if( IType == 0 ){ // non-spec Eds
		file << "E = " << Eds << endl;
	}
	else{ // Na, K, Ca, Cl etc
		file << "E = " << Eds << endl;
		file << "Z = " << Z << endl;
		file << "Cin = " << In << endl;
		file << "Cout = " << Out << endl;
	}
	if( PType != 0){ // ion pump
		file << "T = " << T << endl;		 // time constant
		file << "InEq = " << In0 << endl;
		file << "Rpmp = " << Rpmp << endl;
		file << "Rch = " << Rch << endl;
		switch( PType ){
			case 1: // Na pump
				file << "Kp = " << Kp << endl;
				break;
			case 2: // Ca standard
				break;
			case 3: // Ca shell
				file << "K = " << K << endl;
				file << "B = " << B << endl;
				break;
			default:;
		}
	}
	return true;
}

/////////////////////////////////////////////////////////////////////////////
// class syndata
//--- constructors/destructor
//--- public methods
syndata &syndata::operator = ( const syndata &syn )
{
	Eds = syn.Eds;
	Gmax =syn.Gmax;
	Ts = syn.Ts;
	Vs = syn.Vs;
	return *this;
}

bool syndata::loadpar( istream &file, const string &parname )
{
	string str;
	if( parname == "K" || parname == "Gmax"){
		file >> str >> Gmax;
		return true;
	}
	else if( parname == "Eds"){
		file >> str >> ws;
		getline( file, Eds );
		return true;
	}
	else if( parname == "TimeCostant"){
		file >> str >> Ts;
		return true;
	}
	else if( parname == "Amplitude"){
		file >> str >> Vs;
		return true;
	}
	return false;
}

bool syndata::validate( const string &type )
{
	if( !( type == "Excitatory" || type == "Excitatory2" || type == "Inhibitory I" || type == "Inhibitory II" )){
		return false;
	}
	return true;
}

/////////////////////////////////////////////////////////////////////////////
// class dsyndata
//--- constructors/destructor
//--- public methods
dsyndata &dsyndata::operator = ( const dsyndata &dsyn )
{
	Kdep = dsyn.Kdep;
	Tdep = dsyn.Tdep;
	return *this;
}

bool dsyndata::loadpar( istream &file, const string &parname )
{
	string str;
	if( parname == "Amplitude"){
		file >> str >> Kdep;
		return true;
	}
	else if( parname == "TimeCostant"){
		file >> str >> Tdep;
		return true;
	}
	return false;
}

bool dsyndata::validate( const string &type )
{
	if( !( type == "Excitatory" || type == "Excitatory2" || type == "Inhibitory I" || type == "Inhibitory II" )){
		return false;
	}
	return true;
}

/////////////////////////////////////////////////////////////////////////////
// class dsyndata
//--- constructors/destructor
//--- public methods
ssyndata &ssyndata::operator = ( const ssyndata &ssyn )
{
	Hv = ssyn.Hv;
	Slp = ssyn.Slp;
	return *this;
}

bool ssyndata::loadpar( istream &file, const string &parname )
{
	string str;
	if( parname == "Hv"){
		file >> str >> Hv;
		return true;
	}
	else if( parname == "Slp"){
		file >> str >> Slp;
		return true;
	}
	return false;
}

bool ssyndata::validate( const string &type )
{
	if( !( type == "Excitatory" || type == "Excitatory2" || type == "Inhibitory I" || type == "Inhibitory II" )){
		return false;
	}
	return true;
}

/////////////////////////////////////////////////////////////////////////////
// class chandata
//--- constructors/destructor
static const char *ChannelNames[] = {
	"Na fast",
	"NaP channel (generic)",
	"K",
	"KCa",
	"CaL",
	"Leak",
};

//--- public methods
chandata &chandata::operator = ( const chandata &chdat )
{
	M = chdat.M;
	H = chdat.H;
	Gmax = chdat.Gmax;
	Eds = chdat.Eds;
	Pcl = chdat.Pcl;
	Pna = chdat.Pna;
	Tm = chdat.Tm;
	Th = chdat.Th;
	V12m = chdat.V12m;
	V12h = chdat.V12m;
	Slpm = chdat.Slpm;
	Slph = chdat.Slph;
	V12tm = chdat.V12tm;
	V12th = chdat.V12tm;
	Slptm = chdat.Slptm;
	Slpth = chdat.Slpth;
	PowM = chdat.PowM;
	PowH = chdat.PowH;
	Vg = chdat.Vg;
	Family = chdat.Family;
	return *this;
}

bool chandata::loadpar( istream &file, const string &parname )
{
	string str;
	if( parname == "M" || parname == "Mna" || parname == "Mnap" || parname == "Mk" || parname == "Mcal" || parname == "Mkca" ){
		file >> str >> M;
		return true;
	}
	else if( parname == "H" || parname == "Hna" || parname == "Hnap" || parname == "Hk" ){
		file >> str >> H;
		return true;
	}
	else if( parname == "E"  || parname == "Eleak" ){
		file >> str >> Eds;
		return true;
	}
	else if( parname == "Gmax" ){
		file >> str >> Gmax;
		return true;
	}
	else if( parname == "Tm" || parname ==  "Tkca" ){
		file >> str >> Tm;
		return true;
	}
	else if( parname == "Th" ){
		file >> str >> Th;
		return true;
	}
	else if( parname == "V12m" ){
		file >> str >> V12m;
		return true;
	}
	else if( parname == "V12h" ){
		file >> str >> V12h;
		return true;
	}
	else if( parname == "Slpm" || parname ==  "Km" ){
		file >> str >> Slpm;
		return true;
	}
	else if( parname == "Slph" || parname ==  "Kh" ){
		file >> str >> Slph;
		return true;
	}
	else if( parname == "V12tm" ){
		file >> str >> V12tm;
		return true;
	}
	else if( parname == "V12th" ){
		file >> str >> V12th;
		return true;
	}
	else if( parname == "Slptm" || parname ==  "Ktm" ){
		file >> str >> Slptm;
		return true;
	}
	else if( parname == "Slpth" || parname ==  "Kth" ){
		file >> str >> Slpth;
		return true;
	}
	else if( parname == "PowM" || parname ==  "M_power" ){
		file >> str >> PowM;
		return true;
	}
	else if( parname == "PowH" || parname ==  "H_power" ){
		file >> str >> PowH;
		return true;
	}
	else if( parname == "Pcl" ){
		file >> str >> Pcl;
		return true;
	}
	else if( parname == "Pna" ){
		file >> str >> Pna;
		return true;
	}
	else if( parname == "AdjustableLeak" ){
		file >> str >> Vg;
		return true;
	}
	return false;
}

bool chandata::validate( const string &type )
{
	return true;
}

bool chandata::save( ostream &file, const string &name, int id )
{
	file << "Name = " << name << endl;
	file << "ID = " << id << endl;
	file << "g = " << Gmax << endl;
	file << "E = " << Eds << endl;
	if( Family == "generic" ){
		switch( type( name ) ){
			case 1:
				file << "mpower = "<< PowM.X << endl;
				file << "hpower = "<< PowH.X << endl;
				file << "Vhfm = "<< V12m << endl;
				file << "km = "<< Slpm <<endl;
				if( Slptm.X > 0.f ){
					file << "tm = "<< Tm << endl;
					file << "ktm = "<< Slptm << endl;
				}
				else{
					file << "tm0 = "<< Tm << endl;
				}
				file << "Vhfh = "<< V12h << endl;
				file << "kh = "<< Slph <<endl;
				if( Slpth.X > 0.f ){
					file << "th = "<< Th << endl;
					file << "kth = "<< Slpth << endl;
				}
				else{
					file << "th0 = "<< Th << endl;
				}
				break;
			case 2:
				file << "mpower = "<< PowM.X << endl;
				file << "Vhfm = "<< V12m << endl;
				file << "km = "<< Slpm <<endl;
				if( Slptm.X > 0.f ){
					file << "tm = "<< Tm << endl;
					file << "ktm = "<< Slptm << endl;
				}
				else{
					file << "tm0 = "<< Tm << endl;
				}
				break;
			case 3:
				file << "hpower = "<< PowH.X << endl;
				file << "Vhfh = "<< V12h << endl;
				file << "kh = "<< Slph <<endl;
				if( Slpth.X > 0.f ){
					file << "th = "<< Th << endl;
					file << "kth = "<< Slpth << endl;
				}
				else{
					file << "th0 = "<< Th << endl;
				}
				break;
			case 4:
				break;
			default:
				file << "unknown type of the channel" << endl;
		}
	}
	else if( Family == "ab-type" ){
		switch( type( name )){
			case 2:
				file << "mpower = " << PowM.X << endl;
				file << "AmA = " << -0.01 << endl;
				file << "BmA = " << 1 << endl;
				file << "CmA = " << 1 << endl;
				file << "VhfmA = " << -44 << endl;
				file << "kmA = " << -5 << endl;

				file << "AmB = " << 0.17/49.f << endl;
				file << "BmB = " << 0 << endl;
				file << "CmB = " << 0 << endl;
				file << "VhfmB = " << -49 << endl;
				file << "kmB = " << 40 << endl;
				break;
			default:
				file << "unknown type of the channel" << endl;
		}
	}
	else if( Family == "KCA" ){
		switch( type( name ) ){
			case 2:
				file << "mpower = "<< PowM.X << endl;
				file << "A = "<< 50000000 << endl;
				file << "B = "<< 1.f << endl;
				file << "Lymbda = "<< 2.f << endl;
				file << "Gamma = "<< 1.f << endl;
				file << "tm = "<< Tm << endl;
				break;
			default:
				file << "unknown type of the channel" << endl;
		}
	}
	else{
		return false;
	}
	return true;
}

int chandata::type( const string &name )
{
	Family = "";
	int nChan = sizeof( ChannelNames )/sizeof( char * );
	for( int i = 0; i < nChan; ++i ){
		if( name == ChannelNames[i] ){
			switch( i ){
				case 0: // Na fast
					Family = "generic";
					PowM( 3, 0 );
					PowH( 1, 0 );
					V12m( 43.8f, 0 );
					Slpm( 6.f, 0 );
					V12h( 67.5f, 0 );
					Slph( 10.8f, 0 );
					Tm( 0.252f, 0 );
					Th( 8.456f, 0);
					Slptm( 14.f, 0 );
					Slpth( 12.8f, 0 );
					return 1; // m & h
				case 1: // NaP channel (generic)
					Family = "generic";
					return 1; // m & h
				case 2: // K 
					Family = "ab-type";
					PowM( 4, 0 );
					Tm( 7.f, 0 );
					return 2; // m only
				case 3: // KCa
					Family = "KCA";
					PowM( 2, 0 );
					Tm.X *= 1000.f;
					return 2; // m only
				case 4: // CaL
					Family = "generic";
					PowM( 1, 0 );
					PowH( 1, 0 );
					V12m( 27.4f, 0 );
					Slpm( 5.7f, 0 );
					V12h( 52.4f, 0 );
					Slph( 5.2f, 0 );
					Tm( 0.5f, 0 ); // t0
					Th( 18.f, 0);  // t0
					return 1; // m only
				case 5: // Leak
					return 4; // no m/h
					Family = "leak";
			}
		}
	}
	return -1;
}

/////////////////////////////////////////////////////////////////////////////
// class cmpdata
//--- constructors/destructor
//--- public methods
cmpdata &cmpdata::operator = ( const cmpdata &cdat )
{
	S = cdat.S;
	V = cdat.V;
	Chans = cdat.Chans;
	return *this;
}

bool cmpdata::loadpar( istream &file, const string &parname )
{
	string str;
	if( parname == "S"){
		file >> str >> S;
		return true;
	}
	else if( parname == "V"){
		file >> str >> V;
		return true;
	}
	else if( parname == "<Channel"){
		load_unit( file, Chans, "Channel", "anyname" );
		return true;
	}
	return false;
}

bool cmpdata::validate( const string &type )
{
	return true;
}

/////////////////////////////////////////////////////////////////////////////
// class hhndata
//--- constructors/destructor
//--- public methods
hhndata &hhndata::operator = ( const hhndata &hhn )
{
	C = hhn.C;
	P = hhn.P;
	Ex = hhn.Ex;
	Ex2 = hhn.Ex2;
	InhI = hhn.InhI;
	InhII = hhn.InhII;
	Cmps = hhn.Cmps;
	return *this;
}

bool hhndata::loadpar( istream &file, const string &parname )
{
	string str;
	if( parname == "C"){
		file >> str >> C;
		return true;
	}
	else if( parname == "P"){
		file >> str >> P;
		return true;
	}
	else if( parname == "Ex"){
		file >> str >> Ex;
		return true;
	}
	else if( parname == "Ex2"){
		file >> str >> Ex2;
		return true;
	}
	else if( parname == "InhI"){
		file >> str >> InhI;
		return true;
	}
	else if( parname == "InhII"){
		file >> str >>InhII;
		return true;
	}
	else if( parname == "<Compartment"){
		load_unit( file, Cmps, "Compartment", "Soma|Dendrite" );
		return true;
	}
	return false;
}

bool hhndata::validate( const string &type )
{
	return true;
}


/////////////////////////////////////////////////////////////////////////////
// class popdata
//--- constructors/destructor
//--- public methods
popdata &popdata::operator = ( const popdata &pdat )
{
	Size = pdat.Size;
	Hhn = pdat.Hhn;
	return *this;
}

bool popdata::loadpar( istream &file, const string &parname )
{
	string str;
	if( parname == "Size"){
		file >> str >> Size;
		return true;
	}
	else if( parname == "<Neuron>" ){
		load_unit( file, Hhn, "Neuron" );
		return true;
	}
	return false;
}

bool popdata::validate( const string &type )
{
	return true;
}

/////////////////////////////////////////////////////////////////////////////
// class drivedata
//--- constructors/destructor
//--- public methods
drivedata &drivedata::operator = ( const drivedata &ddat )
{
	Size  = ddat.Size;
	return *this;
}

bool drivedata::loadpar( istream &file, const string &parname )
{
	string str;
	if( parname == "Size"){
		file >> str >> Size;
		return true;
	}
	return false;
}

bool drivedata::validate( const string &type )
{
	Size = 1;	// it's always = 1
	return true;
}

/////////////////////////////////////////////////////////////////////////////
// class outdata
//--- constructors/destructor
//--- public methods
outdata &outdata::operator = ( const outdata &odat )
{
	Size  = odat.Size;
	Tup = odat.Tup;
	Tdown = odat.Tdown;
	Threshold = odat.Threshold;
	return *this;
}

bool outdata::loadpar( istream &file, const string &parname )
{
	string str;
	if( parname == "Size"){
		file >> str >> Size;
		return true;
	}
	else if( parname == "Tup"){
		file >> str >> Tup;
		return true;
	}
	else if( parname == "Tdown"){
		file >> str >> Tdown;
		return true;
	}
	else if( parname == "Threshold"){
		file >> str >> Threshold;
		return true;
	}
	return false;
}

bool outdata::validate( const string &type )
{
	Size = 1;	// it's always = 1
	return true;
}

/////////////////////////////////////////////////////////////////////////////
// class unitdata
//--- constructors/destructor
//--- public methods
unitdata &unitdata::operator = ( const unitdata &udat )
{
	PData = udat.PData;
	DData = udat.DData;
	OData = udat.OData;
	return *this;
}

bool unitdata::loadpar( istream &file, const string &parname )
{
	if( parname == "<Population"){
		load_unit( file, PData, "Population", "anyname" );
		return true;
	}
	else if( parname == "<Drive"){
		load_unit( file, DData, "Drive", "anyname" );
		return true;
	}
	else if( parname == "<Output"){
		load_unit( file, OData, "Output", "anyname" );
		return true;
	}
	return false;
}

bool unitdata::validate( const string &type )
{
	return true;
}

/////////////////////////////////////////////////////////////////////////////
// class trgdata
//--- constructors/destructor
//--- public methods
wdata &wdata::operator = ( const wdata &wdat )
{
	Type = wdat.Type;
	Prob = wdat.Prob;
	Ex = wdat.Ex;
	Ex2 = wdat.Ex2;
	InhI = wdat.InhI;
	InhII = wdat.InhII;
	return *this;
}

bool wdata::loadpar( istream &file, const string &parname )
{
	string str;
	if( parname == "Type"){
		file >> str >> Type;
		return true;
	}
	else if( parname == "Probability"){
		file >> str >> Prob;
		return true;
	}
	else if( parname == "Ex"){
		file >> str >> Ex;
		return true;
	}
	else if( parname == "Ex2"){
		file >> str >> Ex2;
		return true;
	}
	else if( parname == "InhI"){
		file >> str >> InhI;
		return true;
	}
	else if( parname == "InhII"){
		file >> str >> InhII;
		return true;
	}
	return false;
}

bool wdata::validate( const string &type )
{
	return true;
}

/////////////////////////////////////////////////////////////////////////////
// class trgdata
//--- constructors/destructor
//--- public methods
trgdata &trgdata::operator = ( const trgdata &tdat )
{
	WData = tdat.WData;
	return *this;
}

bool trgdata::loadpar( istream &file, const string &parname )
{
	if( parname == "<Connect>"){
		load_unit( file, WData, "Connect" );
		return true;
	}
	return false;
}

bool trgdata::validate( const string &type )
{
	return true;
}

/////////////////////////////////////////////////////////////////////////////
// class srcdata
//--- constructors/destructor
//--- public methods
srcdata &srcdata::operator = ( const srcdata &sdat )
{
	TrgData = sdat.TrgData;
	return *this;
}

bool srcdata::loadpar( istream &file, const string &parname )
{
	if( parname == "<Target"){
		load_unit( file, TrgData, "Target", "anyname" );
		return true;
	}
	return false;
}

bool srcdata::validate( const string &type )
{
	return true;
}

/////////////////////////////////////////////////////////////////////////////
// class condata
//--- constructors/destructor
//--- public methods
condata &condata::operator = ( const condata &cndat )
{
	SrcData = cndat.SrcData;
	return *this;
}

bool condata::loadpar( istream &file, const string &parname )
{
	if( parname == "<Source"){
		load_unit( file, SrcData, "Source", "anyname" );
		return true;
	}
	return false;
}

bool condata::validate( const string &type )
{
	return true;
}

/////////////////////////////////////////////////////////////////////////////
// class ctrdata
//--- public methods
ctrdata &ctrdata::operator = ( const ctrdata &ctdat )
{
	return *this;
}

bool ctrdata::loadpar( istream &file, const string &parname )
{
	return false;
}

bool ctrdata::validate( const string &type )
{
	return true;
}

/////////////////////////////////////////////////////////////////////////////
// class netdata
//--- constructors/destructor
//--- public methods
netdata &netdata::operator = ( const netdata &net )
{
	Threshold = net.Threshold;
	Ions = net.Ions;
	Syns = net.Syns;
	DSyns = net.DSyns;
	SSyns = net.SSyns;
	Units = net.Units;
	return *this;
}

bool netdata::loadpar( istream &file, const string &parname )
{
	string str;
	if( parname == "Threshold"){
		file >> str >> Threshold;
		return true;
	}
	else if( parname == "<Ion"){
		load_unit( file, Ions, "Ion", "Na|K|Ca|Cl" );
		return true;
	}
	else if( parname == "<Synapse" ){
		load_unit( file,Syns, "Synapse", "Excitatory|Excitatory2|Inhibitory I|Inhibitory II" );
		return true;
	}
	else if( parname == "<Depression" ){
		load_unit( file, DSyns, "Depression", "Excitatory|Excitatory2|Inhibitory I|Inhibitory II" );
		return true;
	}
	else if( parname == "<Sigma" ){
		load_unit( file, SSyns, "Sigma", "Excitatory|Excitatory2|Inhibitory I|Inhibitory II" );
		return true;
	}
	else if( parname == "<Units>"){
		load_unit( file, Units, "Units"  );
		return true;
	}
	else if( parname == "<Connections>"){
		load_unit( file, Connections, "Connections"  );
		return true;
	}
/*
	else if( parname == "<Control"){
		load_unit( file, Controls, "Control", "anyname"  );
		return true;
	}
*/
	return false;
}
bool netdata::validate( const string &type )
{
	return true;
}

//--- private methods
void netdata::clear( void )
{
}

/////////////////////////////////////////////////////////////////////////////
// class nsm_model
//--- constructors/destructor
//--- public methods
bool nsm_model::load_model( const char *filename )
{
#ifdef __LINUX__
	remove_char( filename, 0x0D );
#endif
	bool LoadOK = true;
	string str;
	ifstream file( filename );
	while( file >> str ){
		if( str == "Generated" ){
			getline( file, str );
		}
		else if(str == "<SimulateData>" ){
			LoadOK = load_unit( file, SimData, "SimulateData" );
		}
		else if( str == "<Network" ){
			load_unit( file, Network, "Network", "Network" );
		}
/*
		else if( str == "<Views>" ){
//			load_unit( file, Views, "Views" );
		}
		else if( str == "<Records>" ){
//			load_unit( file, Records, "Records" );
		}
		else if( str == "<Biomechanics" ){
//			load_unit( file, Biomech, "Biomechanics", "Cat" );
		}
*/
		else{
			::unknown_string( file, str );
		}
	}
	return LoadOK;
}

bool nsm_model::save_model( const char *filename )
{
	ofstream nns( filename );
	if (nns.is_open()){
		nns << "Neural Network Simulator (NNS), 2015(C)" << endl;
		nns << endl;
		nns << "<VERSION>" << endl;
		nns << "NNS V1.0" << endl;
		nns << "</VERSION>" << endl; // the tag must be closed
		nns << endl;
		nns << "<CTR PARAM>" << endl;
		nns << "Seed = " << SimData().Seed << endl;
		nns << "SimStep = "<< SimData().SimStep << endl;
		nns << "SimDuration = "<< SimData().SimTime << endl;
		nns << "</CTR PARAM>" << endl;
		nns << endl;
		size_t npops = Network().Units().PData.size();
		nns << "<POPULATIONS>" << endl;					// Populations definition
		nns << "N of Pops = " << npops << endl;
		for( size_t i = 0; i < npops; ++i ){
			size_t nICH = Network().Units().PData[i]().Hhn().Cmps[0]().Chans.size();
			size_t nIONS = Network().Ions.size();
			nns << "<POP "<< i << ">" << endl;
			nns << "Name = " << Network().Units().PData[i].Name << endl;
			nns << "ID = " << i << endl;
			nns << "Size = " << Network().Units().PData[i]().Size << endl;
			nns << "Threshold = " << Network().Threshold << endl;
			nns << endl;
			nns << "<NEURON PARAM>" << endl;
			nns << "Cm = " << Network().Units().PData[i]().Hhn().C << endl; // Cm must be parameter of compartment
			nns << "Vm = " << Network().Units().PData[i]().Hhn().Cmps[0]().V << endl;
			nns << "Type = " << 1 << endl;					// ????
			nns << "<ICHANNELS PARAM>" << endl;
			nns << "nICH = " << nICH << endl;
			nns << "nICHTypeList = {";
			for( size_t j = 0; j < nICH-1; ++j ){
				string name = Network().Units().PData[i]().Hhn().Cmps[0]().Chans[j].Name;
				nns << Network().Units().PData[i]().Hhn().Cmps[0]().Chans[j]().type( name ) << ", "; // must be type instead of name
			}
			if( nICH > 0 ){
				string name = Network().Units().PData[i]().Hhn().Cmps[0]().Chans[nICH-1].Name;
				nns << Network().Units().PData[i]().Hhn().Cmps[0]().Chans[nICH-1]().type( name ); // must be type instead of name
			}
			nns << "}" << endl;
			nns << "nICHNameList = {";
			for( size_t j = 0; j < nICH-1; ++j ){
				nns << Network().Units().PData[i]().Hhn().Cmps[0]().Chans[j].Name << ", ";
			}
			if( nICH > 0 ){
				nns << Network().Units().PData[i]().Hhn().Cmps[0]().Chans[nICH-1].Name;
			}
			nns << "}" << endl << endl;
			for( size_t j = 0; j < nICH; ++j ){
				string name = Network().Units().PData[i]().Hhn().Cmps[0]().Chans[j].Name;
				hhnpair<float> eds = Network().Units().PData[i]().Hhn().Cmps[0]().Chans[j]().Eds;
				eds = hack1( eds, name, Network().Ions );
				Network().Units().PData[i]().Hhn().Cmps[0]().Chans[j]().Eds = eds;
				nns << "<ICH " << j << ">" << endl;
				Network().Units().PData[i]().Hhn().Cmps[0]().Chans[j]().save( nns, name, j );
				nns << "</ICH " << j << ">" << endl;
				nns << endl;
			}
			nns << "</ICHANNELS PARAM>" << endl;
			nns << "</NEURON PARAM>" << endl;
			nns << endl;
			nns << "<IONS PARAM>" << endl;
			size_t n_ions = 0;
			for( size_t j = 0; j < nIONS; ++j ){
				if( hack2( Network().Ions[j].Name, Network().Units().PData[i]().Hhn().Cmps[0]().Chans )){
					n_ions=1;
				}
			}
			nns << "nIONS = " << n_ions << endl;
			for( size_t j = 0; j < nIONS; ++j ){
				if( hack2( Network().Ions[j].Name, Network().Units().PData[i]().Hhn().Cmps[0]().Chans )){
					/*check if Ca channels are present then print it*/
					nns << "<ION " << j << ">" << endl;
					nns << "Name = " << Network().Ions[j].Name << endl;
					nns << "ID = " << j << endl;
					Network().Ions[j]().save( nns );
					nns << "</ION " << j << ">" << endl;
					nns << endl;
				}
			};
			nns << "</IONS PARAM>" << endl;
			nns << endl;
			nns << "</POP "<< i << ">" << endl;
			nns<<endl;
		}
		nns << "</POPULATIONS>" << endl;
		nns << "<NETWORKS>" << endl;						// All Networks definition
		nns<<endl;
/*
		nns<<"N of Nets = "<<nnets<<endl;
		nns<<"N of MNets = "<<nmnets<<endl;
		nns<<endl;

		for (int i=0;i<nnets;i++) {

			nns<<"<NET "<<i<<">"<<endl;
			nns<<"Name = "<<nets[i].name<<endl;
			nns<<"OnOff = "<<nets[i].OnOff<<endl;
			nns<<"nList = "<<nets[i].nList<<endl;

			nns<<"pnameList = {";
			for (int j=0;j<nets[i].nList;j++) {
				if (j != nets[i].nList-1)
					nns<<nets[i].pnameList[j]<<", ";
				else
					nns<<nets[i].pnameList[j]<<"}"<<endl;
			};

			nns<<"pidList = {";
			for (int j=0;j<nets[i].nList;j++) {
				if (j != nets[i].nList-1)
					nns<<nets[i].pidList[j]<<", ";
				else
					nns<<nets[i].pidList[j]<<"}"<<endl;
			};

			for (int j=0;j<nets[i].pconn.size();j++) {
				nns<<"<PCONN "<<j<<">"<<endl;
				nns<<"Source Pop = "<<nets[i].pconn[j].spop<<endl;
				nns<<"Source Prob = "<<nets[i].pconn[j].sprob.val<<" ("<<nets[i].pconn[j].sprob.var<<")"<<endl;
				nns<<"Target Pop = "<<nets[i].pconn[j].tpop<<endl;
				nns<<"Target Prob = "<<nets[i].pconn[j].tprob.val<<" ("<<nets[i].pconn[j].tprob.var<<")"<<endl;
				nns<<"Weight = "<<nets[i].pconn[j].w.val<<" ("<<nets[i].pconn[j].w.var<<")"<<endl;
				nns<<"</PCONN "<<j<<">"<<endl;
			};
			nns<<"</NET "<<i<<">"<<endl;
			nns<<endl;
		};

		for (int i=0;i<nmnets;i++) {
			
			nns<<"<MNET "<<i<<">"<<endl;
			nns<<"Name = "<<mnets[i].name<<endl;
			nns<<"OnOff = "<<mnets[i].OnOff<<endl;
			nns<<"nList = "<<mnets[i].nList<<endl;

			nns<<"nnameList = {";
			for (int j=0;j<mnets[i].nList;j++) {
				if (j != mnets[i].nList-1)
					nns<<mnets[i].nnameList[j]<<", ";
				else
					nns<<mnets[i].nnameList[j]<<"}"<<endl;
			};

			nns<<"nidList = {";
			for (int j=0;j<mnets[i].nList;j++) {
				if (j != mnets[i].nList-1)
					nns<<mnets[i].nidList[j]<<", ";
				else
					nns<<mnets[i].nidList[j]<<"}"<<endl;
			};
			
			for (int j=0;j<mnets[i].nconn.size();j++) {
				nns<<"<NCONN "<<j<<">"<<endl;
				nns<<"Source Net = "<<mnets[i].nconn[j].snet<<endl;
				nns<<"Source Pop = "<<mnets[i].nconn[j].spop<<endl;
				nns<<"Source Prob = "<<mnets[i].nconn[j].sprob.val<<" ("<<mnets[i].nconn[j].sprob.var<<")"<<endl;
				nns<<"Target Net = "<<mnets[i].nconn[j].tnet<<endl;
				nns<<"Target Pop = "<<mnets[i].nconn[j].tpop<<endl;
				nns<<"Target Prob = "<<mnets[i].nconn[j].tprob.val<<" ("<<mnets[i].nconn[j].tprob.var<<")"<<endl;
				nns<<"Weight = "<<mnets[i].nconn[j].w.val<<" ("<<mnets[i].nconn[j].w.var<<")"<<endl;
				nns<<"</NCONN "<<j<<">"<<endl;
			};

			nns<<"</MNET "<<i<<">"<<endl;
			nns<<endl;

		};
		
*/
		nns << "</NETWORKS>" << endl;

		nns.close();
	}
	return false;
}
