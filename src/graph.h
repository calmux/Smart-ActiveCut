
/* graph.h */
/*
	This software library implements the maxflow algorithm
	described in

		"An Experimental Comparison of Min-Cut/Max-Flow Algorithms for Energy Minimization in Vision."
		Yuri Boykov and Vladimir Kolmogorov.
		In IEEE Transactions on Pattern Analysis and Machine Intelligence (PAMI), 
		September 2004

	This algorithm was developed by Yuri Boykov and Vladimir Kolmogorov
	at Siemens Corporate Research. To make it available for public use,
	it was later reimplemented by Vladimir Kolmogorov based on open publications.

	If you use this software for research purposes, you should cite
	the aforementioned paper in any resulting publication.

	----------------------------------------------------------------------

	REUSING TREES:

	Starting with version 3.0, there is a also an option of reusing search
	trees from one maxflow computation to the next, as described in

		"Efficiently Solving Dynamic Markov Random Fields Using Graph Cuts."
		Pushmeet Kohli and Philip H.S. Torr
		International Conference on Computer Vision (ICCV), 2005

	If you use this option, you should cite
	the aforementioned paper in any resulting publication.
*/
	


/*
	For description, license, example usage see README.TXT.
*/

#ifndef __GRAPH_H__
#define __GRAPH_H__

#include <string.h>
#include "block.h"

#include <assert.h>
// NOTE: in UNIX you need to use -DNDEBUG preprocessor option to supress assert's!!!



// captype: type of edge capacities (excluding t-links)
// tcaptype: type of t-links (edges between nodes and terminals)
// flowtype: type of total flow
//
// Current instantiations are in instances.inc
template <typename captype, typename tcaptype, typename flowtype> class Graph
{
public:
	typedef enum
	{
		SOURCE	= 0,
		SINK	= 1
	} termtype; // terminals 
	typedef int node_id;

	/////////////////////////////////////////////////////////////////////////
	//                     BASIC INTERFACE FUNCTIONS                       //
    //              (should be enough for most applications)               //
	/////////////////////////////////////////////////////////////////////////

	// Constructor. 
	// The first argument gives an estimate of the maximum number of nodes that can be added
	// to the graph, and the second argument is an estimate of the maximum number of edges.
	// The last (optional) argument is the pointer to the function which will be called 
	// if an error occurs; an error message is passed to this function. 
	// If this argument is omitted, exit(1) will be called.
	//
	// IMPORTANT: It is possible to add more nodes to the graph than node_num_max 
	// (and node_num_max can be zero). However, if the count is exceeded, then 
	// the internal memory is reallocated (increased by 50%) which is expensive. 
	// Also, temporarily the amount of allocated memory would be more than twice than needed.
	// Similarly for edges.
	// If you wish to avoid this overhead, you can download version 2.2, where nodes and edges are stored in blocks.
	Graph(int node_num_max, int edge_num_max, void (*err_function)(char *) = NULL);

	// Destructor
	~Graph();

	// Adds node(s) to the graph. By default, one node is added (num=1); then first call returns 0, second call returns 1, and so on. 
	// If num>1, then several nodes are added, and node_id of the first one is returned.
	// IMPORTANT: see note about the constructor 
	node_id add_node(int num = 1);

	// Adds a bidirectional edge between 'i' and 'j' with the weights 'cap' and 'rev_cap'.