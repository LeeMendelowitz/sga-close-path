// Author: Lee Mendelowitz
// Date:   9/28/2012
// Graph Search Functionality

#ifndef EDGE_GEN_H
#define EDGE_GEN_H
#include "SGUtil.h"

// Perform a bounded BFS to collect edges starting from pVertex in direction dir.
EdgePtrVec boundedBFS(Vertex * pVertex, EdgeDir dir, size_t maxNodes, int maxDistance);
StringGraph * makePathGraph(StringGraph * pGraph, Vertex * pX, EdgeDir dX, Vertex * pY, EdgeDir dY, int maxDistance);

#endif
