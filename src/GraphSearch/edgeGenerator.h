// Author: Lee Mendelowitz
// Date:   9/28/2012
// Graph Search Functionality
#ifndef EDGE_GEN_H
#define EDGE_GEN_H

#include "SGUtil.h"
#include <set>
#include <map>

typedef std::set<Edge *> EdgePtrSet;
typedef std::pair<const Vertex *, EdgeDir> VDirPair;
typedef std::map<VDirPair, int> DistanceMap;

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
EdgePtrVec boundedBFS(const Vertex * pVertex, EdgeDir dir, int maxDistance);

// Perform a bounded BFS to collect edges starting from pVertex in direction dir, up to a maximum distance.
// Use only the edges in the allowableEdges set.
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
EdgePtrVec boundedBFS(const Vertex * pVertex, EdgeDir dir, int maxDistance, const EdgePtrSet& allowableEdges);
void boundedBFS(const Vertex * pVertex, EdgeDir dir, int halfDistance, int maxDistance,
                const Vertex * pGoal, EdgeDir goalDir,
                EdgePtrVec& halfDistEdges, EdgePtrVec& maxDistEdges);

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
StringGraph * makePathGraph(const StringGraph * pGraph, const Vertex * pX, EdgeDir dX, const Vertex * pY, EdgeDir dY, int maxDistanceX);

// Collect edges that are on gauranteed to be on paths from Vertex pX to pY
// with distance less than maxDistance.
// Return the edges in sorted order.
//
// See description for boundedBFS for how distance is measured.
// Return a pointer to the subgraph, or NULL if the subgraph is empty.

// dX: direction of edge out of pX on walk to pY.
// dY: direction of edge out of pY on walk to pX.
// Case 1: pX Forward, pY Reverse, then dX = ED_SENSE, dY = ED_SENSE |--->.......<----|
// Case 2: pX Forward, pY Forward, then dX = ED_SENSE, dY = ED_ANTISENSE     |--->.......|---->
// Case 3: pX Reverse, pX Forward, then dX = ED_ANTISENSE, dY = ED_ANTISENSE <---|......|----->
// Case 4: pX Reverse, pY Reverse, then dX = ED_ANTISENSE, dY = ED_SENSE  <----|......<----|
EdgePtrVec getPathEdges(const Vertex * pX, EdgeDir dX, const Vertex * pY, EdgeDir dY, int maxDistanceX);
EdgePtrVec getPathEdges2(const Vertex * pX, EdgeDir dX, const Vertex * pY, EdgeDir dY, int maxDistanceX);
EdgePtrVec getPathEdges3(const Vertex * pX, EdgeDir dX, const Vertex * pY, EdgeDir dY, int maxDistanceX);
EdgePtrVec dijkstra(const Vertex * pVertex, EdgeDir dir, int maxDistance, DistanceMap& distMap);
EdgePtrVec dijkstra(const Vertex * pVertex, EdgeDir dir, int maxDistance, DistanceMap& distMap, const EdgePtrSet& allowableEdges);
#endif
