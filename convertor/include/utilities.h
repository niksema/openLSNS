#ifndef __UTILITIES_H_
#define __UTILITIES_H_

#include <math.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using std::cout;
using std::endl;
using std::ifstream;
using std::istream;
using std::ios_base;
using std::getline;
using std::ofstream;
using std::ostream;
using std::string;
using std::vector;
using std::ws;

template<class T>
class hhnpair{
	public:
		hhnpair( void ) : X( 0 ), Y( 0 ){};
		hhnpair( const T &x ) : X( x ), Y( 0 ){};
		hhnpair( const T &x, const T &y ) : X( x ), Y( y ){};
		hhnpair( const hhnpair &xy )  : X( xy.X ), Y( xy.Y ){};
		~hhnpair( void ){};
	public:
		hhnpair &operator = ( const hhnpair &xy ){
			X = xy.X;
			Y = xy.Y;
			return *this;
		};
	public:
		T X;
		T Y;
};

template<class T>
ostream &operator << ( ostream &s, const hhnpair<T> &xy )
{
	s << xy.X;
	if( xy.Y != T( 0 ))
		s << "(" << xy.Y << ")";
	return  s;
};

template<class T>
istream &operator >> ( istream &s, hhnpair<T> &xy )
{
	char ch = 0;
	bool OK = true;

	T x, y = T( 0 );
	s >> ws >> x >> ws;
	s.get( ch );
	if( ch == '(' ){
		s >> y >> ws;
		s.get( ch );
		OK = ( ch == ')' );
	}
	else if( ch == '\n' ){
		y = T( 0 );
		OK = true;
	}
	else{
		s.putback( ch );
	}
	if( s && OK ){
		xy.X = x;
		xy.Y = y;
	}
	return s;
}

extern void message( const char *mess, const char *type );
extern bool remove_char( const char *filename, char ch );
extern void unknown_string( istream &file, const string &str );

#endif /*__UTILITIES_H_*/
