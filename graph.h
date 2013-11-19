/*
 * graph.h
 *
 *  Created on: Nov 14, 2013
 *      Author: olga
 */

#include <iostream>
#include <list>
#include <limits.h>

#ifndef GRAPH_H_
#define GRAPH_H_

using namespace std;

class Graph
{
    int V;    			// No. of vertices
    list<int> *adj;		// Pointer to an array containing adjacency lists
    bool isCyclicUtil(int v, bool visited[], bool *rs);  // used by isCyclic()

public:
    Graph(int V);   // Constructor
    void addEdge(int v, int w);   // to add an edge to graph
    bool isCyclic();    // returns true if there is a cycle in this graph
};

#endif
