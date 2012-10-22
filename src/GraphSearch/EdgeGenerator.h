// Author: Lee Mendelowitz
// Date:   9/28/2012
// Graph Search Functionality

#ifndef EDGE_GEN_H
#define EDGE_GEN_H
#include "SGUtil.h"

// Perform a bounded BFS to collect edges starting from pVertex in direction dir, up to a maximum distance.
// NOTE: This search differs from the search functionality provided in SGSearch.
// Search distance is measured from the beginning of the start vertex.
// |----------> A
// 
// | 
// 0
//
//          <-----| B
// |<------>| distance from A to B.
//
// Distance 0 is at the start of vertex A
// A vertex reachable from A is measured as the distance between 0 and the start of that vertex.
EdgePtrVec boundedBFS(Vertex * pVertex, EdgeDir dir, int maxDistance);

// This version copies StringGraphAttributes from pGraph
// Make a subgraph of nodes that are on gauranteed to be on paths from Vertex pX to pY
// with distance less than maxDistance.
// See description for boundedBFS for how distance is measured.

// dX: direction of edge out of pX on walk to pY.                                                                                                                                        
// dY: direction of edge out of pY on walk to pX.
// Case 1: pX Forward, pY Reverse, then dX = ED_SENSE, dY = ED_SENSE |--->.......<----|
// Case 2: pX Forward, pY Forward, then dX = ED_SENSE, dY = ED_ANTISENSE     |--->.......|---->
// Case 3: pX Reverse, pX Forward, then dX = ED_ANTISENSE, dY = ED_ANTISENSE <---|......|----->
// Case 4: pX Reverse, pY Reverse, then dX = ED_ANTISENSE, dY = ED_SENSE  <----|......<----|
StringGraph * makePathGraph(StringGraph * pGraph, Vertex * pX, EdgeDir dX, Vertex * pY, EdgeDir dY, int maxDistanceX);

#endif
