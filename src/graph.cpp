/* graph.cpp */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "graph.h"


template <typename captype, typename tcaptype, typename flowtype> 
	Graph<captype, tcaptype, flowtype>::Graph(int node_num_max, int edge_num_max, void (*err_function)(char *))
	: node_num(0),
	  nodeptr_block(NULL),
	  error_function(err_function)
{
	if (node_num_max < 16) node_num_max = 16;
	if (edge_num_max < 16) edge_num_max = 16;

	nodes = (node*) malloc(node_num_max*sizeof(node));
	arcs = (arc*) malloc(2*edge_num_max*sizeof(arc));
	if (!nodes || !arcs) { if (error_function) (*error_function)("Not enough memory!"); exit(1); }

	node_last = nodes;
	node_max = nodes + node_num_ma