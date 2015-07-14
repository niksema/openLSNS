#include "precompile.h"

#include "views.h"

int main( void )
{
	float4vw *a;
	float4 b;
	b.x = 1;
	b.y = 2;
	b.z = 3;
	b.w = 4;
	a = ( float4vw *)&b;
	int k1 = 0x3FFFFFFF;
	unsigned int k2 = 0xFFFFFFFF;

	unsigned int n1 = k1 >> 29;
	unsigned int n2 = k1 >> 30;
	unsigned int n3 = k2 >> 29;
	unsigned int n4 = k2 >> 30;
	unsigned int n5 = k2 >> 31;
	unsigned int n6 = k2 >> 32;

	return 0;
}