#ifndef __LSNS_ASSERT_H
#define __LSNS_ASSERT_H

#include "config.h"
#include <assert.h>

#if defined( __LSNS_DEBUG__ )
	#define __lsns_assert( c ) assert( c )
#else
	#define __lsns_assert( c )
#endif



#endif /*__LSNS_ASSERT_H*/

