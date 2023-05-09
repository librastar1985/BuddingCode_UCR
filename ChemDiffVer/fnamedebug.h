/*
*			fnamedebug.h
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*/

#if !defined(_FNAMEDEBUG_H)
#define _FNAMEDEBUG_H
#include "cdecs.h"

#if defined(IGNORE_ERRORS)
#define	Check_return(func,fname)					\
	(void) (func);
#else /* defined(IGNORE_ERRORS) */
#define	Check_return(func,fname)					 \
	{								 \
	    if (! (func))						 \
	    {							         \
	    	(void) printf("ERROR in %s(), %s failed\n",#fname,#func);\
	    	clean_up(ERROR);				         \
	    }							         \
	}
#endif /* defined(IGNORE_ERRORS) */

#if defined(DEBUG_STRING)
#define DEBUG_ENTER(fname)						\
	debug_print(DEBUG_STRING,"Entered %s()\n",#fname);
#define DEBUG_LEAVE(fname)						\
	debug_print(DEBUG_STRING,"Left %s()\n",#fname);
#undef DEBUG
#define	DEBUG		debugging(DEBUG_STRING)
#else /* defined(DEBUG_STRING) */
#define DEBUG_ENTER(fname)
#define DEBUG_LEAVE(fname)
#undef DEBUG
#define	DEBUG		NO
#endif /* defined(DEBUG_STRING) */

#include "navdecs.h"

#endif /* !defined(_FNAMEDEBUG_H) */
