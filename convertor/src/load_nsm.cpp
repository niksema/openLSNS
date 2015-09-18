#include "stdafx.h"

#include "load_nsm.h"

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
	Out = ion.Out;
	T = ion.T;
	Eds = ion.Eds;
	B = ion.B;
	K = ion.K;
	return *this;
};

bool iondata::loadpar( istream &file, const string &parname )
{
	string str;
	if( parname == "In"){
		file >> str >> In;
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
	if( type == "Na" || type == "K" ){
		Z = 1;
	}
	else if( type == "Ca" ){
		Z = 2;
	}
	else if( type == "Cl" ){
		Z = -1;
	}
	else{
		Z = 1;
		string mess = "unknown ions type:" + type;
		message( mess.c_str(), "Warning" );
		return false;
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
		file >> str >> Slpm;
		return true;
	}
	else if( parname == "Slpth" || parname ==  "Kth" ){
		file >> str >> Slph;
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
		nns << "<POPULATIONS>" << endl;
		nns << "N of Pops = " << npops << endl;
		nns << "</POPULATIONS>" << endl;
/*
		for( size_t i = 0; i < Network().Units.size(); ++i ){

		}
*/
/*
		Network.save( nns );
bool netdata::save( ostream &file )
{
	float Threshold;
	vector<netunit<iondata >> Ions;
	vector<netunit<syndata >> Syns;
	vector<netunit<dsyndata >> DSyns;
	vector<netunit<ssyndata >> SSyns;
	vector<netunit<unitdata >> Units;
	vector<netunit<condata >> Connections;
	vector<netunit<ctrdata >> Controls;

	for( size_t i = 0; i < Units.size(); ++i ){
		Units[i].save( file );
		bool unitdata::save( ostream &file )
		{
				vector<netunit<popdata >> PData;	// list of populations
				vector<netunit<drivedata >> DData;	// list of drives
				vector<netunit<outdata >> OData;	// list of outputs
			file << "<POPULATIONS>" << endl;
			file << "N of Pops = " << PData.size() << endl;

			for( size_t i = 0; i < PData.size(); ++i ){
			}
			file << "</POPULATIONS>" << endl;
			return true;
		}
	}
	return true;
}
*/
		nns.close();
	}
	return false;
}
