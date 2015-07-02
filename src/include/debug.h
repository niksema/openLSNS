#ifndef __LSNS_ASSERT_H
#define __LSNS_ASSERT_H

#include <assert.h>
#include "config.h"

#if defined( __LSNS_DEBUG__ )
	__lsns_assert( cond ) assert( cond )
#else
	__lsns_assert( cond )
#endif



#endif /*__LSNS_ASSERT_H*/

