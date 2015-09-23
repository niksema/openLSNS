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
	}
	else if( type == "K" ){
		IType = 2;
		PType = 0;
		Z = 1;
	}
	else if( type == "Ca" ){
		IType = 3;
		PType = 3;
		Z = 2;
	}
	else if( type == "Cl" ){
		IType = 4;
		PType = 0;
		Z = -1;
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
		float eds = ( RTF/Z )*log( Out.X/In.X );
		file << "E = " << eds << endl;
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
		nns << "<POPULATIONS>" << endl;					// Populations definition
		nns << "N of Pops = " << npops << endl;
		for( size_t i = 0; i < npops; ++i ){
//			nichs = pops[i].neu_s.nICh;
//			nions = pops[i].paraset.ions.nIons;
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
			size_t nICH = Network().Units().PData[i]().Hhn().Cmps[0]().Chans.size();
			nns << "nICH = " << nICH << endl;
			nns << "nICHTypeList = {";
			for( size_t j = 0; j < nICH-1; ++j ){
				nns << Network().Units().PData[i]().Hhn().Cmps[0]().Chans[j].Name << ", "; // must be type instead of name
			}
			if( nICH > 0 ){
				nns << Network().Units().PData[i]().Hhn().Cmps[0]().Chans[nICH-1].Name; // must be type instead of name
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
				nns << "<ICH " << j << ">" << endl;
/*
				switch (pops[i].neu_s.IChTypeList[j]){
				case 1:
					nns<<"Name = "<<pops[i].paraset.ichs.g_para.at(tc1).name<<endl;
					nns<<"ID = "<<pops[i].paraset.ichs.g_para.at(tc1).ID<<endl;

					nns<<"g = "<<pops[i].paraset.ichs.g_para.at(tc1).g.val<<" ("<<pops[i].paraset.ichs.g_para.at(tc1).g.var<<")"<<endl;
					nns<<"E = "<<pops[i].paraset.ichs.g_para.at(tc1).E.val<<" ("<<pops[i].paraset.ichs.g_para.at(tc1).E.var<<")"<<endl;
					nns<<"mpower = "<<pops[i].paraset.ichs.g_para.at(tc1).mpower<<endl;
					nns<<"hpower = "<<pops[i].paraset.ichs.g_para.at(tc1).hpower<<endl;

					nns<<"Vhfm = "<<pops[i].paraset.ichs.g_para.at(tc1).Vhfm.val<<" ("<<pops[i].paraset.ichs.g_para.at(tc1).Vhfm.var<<")"<<endl;
					nns<<"km = "<<pops[i].paraset.ichs.g_para.at(tc1).km.val<<" ("<<pops[i].paraset.ichs.g_para.at(tc1).km.var<<")"<<endl;
					nns<<"tm = "<<pops[i].paraset.ichs.g_para.at(tc1).tm.val<<" ("<<pops[i].paraset.ichs.g_para.at(tc1).tm.var<<")"<<endl;
					nns<<"ktm = "<<pops[i].paraset.ichs.g_para.at(tc1).ktm.val<<" ("<<pops[i].paraset.ichs.g_para.at(tc1).ktm.var<<")"<<endl;

					nns<<"Vhfh = "<<pops[i].paraset.ichs.g_para.at(tc1).Vhfh.val<<" ("<<pops[i].paraset.ichs.g_para.at(tc1).Vhfh.var<<")"<<endl;
					nns<<"kh = "<<pops[i].paraset.ichs.g_para.at(tc1).kh.val<<" ("<<pops[i].paraset.ichs.g_para.at(tc1).kh.var<<")"<<endl;
					nns<<"th = "<<pops[i].paraset.ichs.g_para.at(tc1).th.val<<" ("<<pops[i].paraset.ichs.g_para.at(tc1).th.var<<")"<<endl;
					nns<<"kth = "<<pops[i].paraset.ichs.g_para.at(tc1).kth.val<<" ("<<pops[i].paraset.ichs.g_para.at(tc1).kth.var<<")"<<endl;

					tc1++;
					break;
				case 2:
					nns<<"Name = "<<pops[i].paraset.ichs.m_para.at(tc2).name<<endl;
					nns<<"ID = "<<pops[i].paraset.ichs.m_para.at(tc2).ID<<endl;

					nns<<"g = "<<pops[i].paraset.ichs.m_para.at(tc2).g.val<<" ("<<pops[i].paraset.ichs.m_para.at(tc2).g.var<<")"<<endl;
					nns<<"E = "<<pops[i].paraset.ichs.m_para.at(tc2).E.val<<" ("<<pops[i].paraset.ichs.m_para.at(tc2).E.var<<")"<<endl;
					nns<<"mpower = "<<pops[i].paraset.ichs.m_para.at(tc2).mpower<<endl;

					nns<<"Vhfm = "<<pops[i].paraset.ichs.m_para.at(tc2).Vhfm.val<<" ("<<pops[i].paraset.ichs.m_para.at(tc2).Vhfm.var<<")"<<endl;
					nns<<"km = "<<pops[i].paraset.ichs.m_para.at(tc2).km.val<<" ("<<pops[i].paraset.ichs.m_para.at(tc2).km.var<<")"<<endl;
					nns<<"tm = "<<pops[i].paraset.ichs.m_para.at(tc2).tm.val<<" ("<<pops[i].paraset.ichs.m_para.at(tc2).tm.var<<")"<<endl;
					nns<<"ktm = "<<pops[i].paraset.ichs.m_para.at(tc2).ktm.val<<" ("<<pops[i].paraset.ichs.m_para.at(tc2).ktm.var<<")"<<endl;

					tc2++;
					break;
				case 3:
					nns<<"Name = "<<pops[i].paraset.ichs.h_para.at(tc3).name<<endl;
					nns<<"ID = "<<pops[i].paraset.ichs.h_para.at(tc3).ID<<endl;

					nns<<"g = "<<pops[i].paraset.ichs.h_para.at(tc3).g.val<<" ("<<pops[i].paraset.ichs.h_para.at(tc3).g.var<<")"<<endl;
					nns<<"E = "<<pops[i].paraset.ichs.h_para.at(tc3).E.val<<" ("<<pops[i].paraset.ichs.h_para.at(tc3).E.var<<")"<<endl;
					nns<<"hpower = "<<pops[i].paraset.ichs.h_para.at(tc3).hpower<<endl;

					nns<<"Vhfh = "<<pops[i].paraset.ichs.h_para.at(tc3).Vhfh.val<<" ("<<pops[i].paraset.ichs.h_para.at(tc3).Vhfh.var<<")"<<endl;
					nns<<"kh = "<<pops[i].paraset.ichs.h_para.at(tc3).kh.val<<" ("<<pops[i].paraset.ichs.h_para.at(tc3).kh.var<<")"<<endl;
					nns<<"th = "<<pops[i].paraset.ichs.h_para.at(tc3).th.val<<" ("<<pops[i].paraset.ichs.h_para.at(tc3).th.var<<")"<<endl;
					nns<<"kth = "<<pops[i].paraset.ichs.h_para.at(tc3).kth.val<<" ("<<pops[i].paraset.ichs.h_para.at(tc3).kth.var<<")"<<endl;

					tc3++;
					break;
				case 4:
					nns<<"Name = "<<pops[i].paraset.ichs.ng_para.at(tc4).name<<endl;
					nns<<"ID = "<<pops[i].paraset.ichs.ng_para.at(tc4).ID<<endl;

					nns<<"g = "<<pops[i].paraset.ichs.ng_para.at(tc4).g.val<<" ("<<pops[i].paraset.ichs.ng_para.at(tc4).g.var<<")"<<endl;
					nns<<"E = "<<pops[i].paraset.ichs.ng_para.at(tc4).E.val<<" ("<<pops[i].paraset.ichs.ng_para.at(tc4).E.var<<")"<<endl;

					tc4++;
					break;
				
				};
*/
				nns << "</ICH " << j << ">" << endl;
				nns << endl;
			}
			nns << "</ICHANNELS PARAM>" << endl;
			nns << "</NEURON PARAM>" << endl;
			nns << endl;
			nns << "<IONS PARAM>" << endl;
			size_t nIONS = Network().Ions.size();
			nns << "nIONS = " << nIONS << endl;
			for( size_t j = 0; j < nIONS; ++j ){
				nns << "<ION " << j << ">" << endl;
				nns << "Name = " << Network().Ions[j].Name << endl;
				nns << "ID = " << j << endl;
				Network().Ions[j]().save( nns );
				nns << "</ION " << j << ">" << endl;
				nns << endl;
			};
			nns << "</IONS PARAM>" << endl;
			nns << endl;
			nns << "</POP "<< i << ">" << endl;
			nns<<endl;
		}
		nns << "</POPULATIONS>" << endl;
		nns << "<NETWORKS>" << endl;						// All Networks definition
/*
		nns<<endl;
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