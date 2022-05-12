/* maxflow.cpp */


#include <stdio.h>
#include "graph.h"


/*
	special constants for node->parent
*/
#define TERMINAL ( (arc *) 1 )		/* to terminal */
#define ORPHAN   ( (a