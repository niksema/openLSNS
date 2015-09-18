#ifndef __NETUNIT_H_
#define __NETUNIT_H_

#include <string>
#include <vector>

using std::string;
using std::vector;

class base_unit{
	public:
		base_unit( void ){};
		base_unit( const base_unit &unit ){};
		~base_unit(){};
	public:
		base_unit &operator = ( const base_unit &unit ){ return *this; };

};
/////////////////////////////////////////////////////////////////////////////
// class netunit
template <class type_, class auto_ = base_unit>
class netunit : public auto_{
	public:
		netunit( const char *tag = "" ) : Name(""), Tag( tag ), Unit(), auto_(){};
		netunit( const netunit<type_, auto_> &unit ) : Name( unit.Name ), Tag( unit.Tag ), Unit( unit.Unit ), auto_( unit ){};
		~netunit( void ){};
	public:
		netunit<type_,auto_> &operator = ( const netunit<type_,auto_> &unit ){
			auto_::operator = ( unit );
			Name = unit.Name;
			Tag = unit.Tag;
			Unit = unit.Unit;
			return *this;
		};
		type_ &operator()( void ){
			return Unit;
		};
	public:
		bool load( istream &file, const char *name ){
			return load( file, string( name ));
		};
		bool load( istream &file, const string &name ){
			string endtag = ( !Tag.empty() )? string("</") + Tag+ string( ">" ) : "";
			string str;
			while( file >> str){
				if( Unit.loadpar( file, str )){
				}
				else if( !endtag.empty() && str == endtag ){
					if( !Unit.validate( name )){
						str = name+" is not valid in <"+Tag + ">";
						message( str.c_str(), "Warning" );
					}
					Name = name;
					return true;
				}
				else{
					::unknown_string( file, str);
				}
			}
			return false;
		}
	public:
		string Name;
		string Tag;
		type_ Unit;
};

template <class type_, class auto_> 
bool load_unit( istream &file, netunit<type_,auto_> &unit, const char *tag = "", const char *name = "noname" )
{
	string name_( name );
	if( !name_.empty() && name_ != "noname" ){
		file >> ws;
		getline( file, name_, '>' );
	}
	netunit<type_,auto_> unit_( tag );
	if( unit_.load( file, name_ )){
		unit = unit_;
		return true;
	}
	return false;
};

template <class type_, class auto_> 
bool load_unit( istream &file, vector<netunit<type_, auto_ >> &units, const char *tag = "", const char *name = "noname" )
{
	string name_( name );
	if( !name_.empty() && name_ != "noname" ){
		file >> ws;
		getline( file, name_, '>' );
	}
	netunit<type_, auto_> unit( tag );
	if( unit.load( file, name_ )){
		units.push_back( unit );
		return true;
	}
	return false;
};

template <class type_, class auto_> 
bool save_unit( ostream &file, netunit<type_,auto_> &unit )
{
/*
	string name_( name );
	if( !name_.empty() && name_ != "noname" ){
		file >> ws;
		getline( file, name_, '>' );
	}
	netunit<type_,auto_> unit_( tag );
	if( unit_.load( file, name_ )){
		unit = unit_;
		return true;
	}
*/
	return false;
};

#endif /*__NETUNIT_H_*/
