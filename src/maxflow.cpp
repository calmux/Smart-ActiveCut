/* maxflow.cpp */


#include <stdio.h>
#include "graph.h"


/*
	special constants for node->parent
*/
#define TERMINAL ( (arc *) 1 )		/* to terminal */
#define ORPHAN   ( (arc *) 2 )		/* orphan */


#define INFINITE_D ((int)(((unsigned)-1)/2))		/* infinite distance to the terminal */

/***********************************************************************/

/*
	Functions for processing active list.
	i->next points to the next node in the list
	(or to i, if i is the last node in the list).
	If i->next is NULL iff i is not in the list.

	There are two queues. Active nodes a