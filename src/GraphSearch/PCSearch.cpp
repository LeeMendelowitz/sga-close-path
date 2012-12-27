#include <queue>

#include "PCSearch.h"
#include "edgeGenerator.h"
#include "dfs.h"

//#include <iostream>

#define PCSEARCH_DEBUG 0

// Declarations
Vertex * getOrigVertex(const StringGraph * pGraph, const Vertex * pNewVertex);
Edge * getOrigEdge(const StringGraph * pGraph, const Edge * pNewEdge);
SGWalk convertWalk(const StringGraph * pGraph, const SGWalk& walkIn);

// This algorithm will first create a subgraph that consists only of vertexes and edges that could be used
// on a path from pX to pY on a path with length less than the prescribed maxDistance.
// Then it searches for all possible valid paths in the subgraph
//
// Return true if the search completed, else return false
bool PCSearch::findWalks(StringGraph * pGraph, SGSearchParams params, bool exhaustive, SGWalkVector& outWalks)
{
    outWalks.clear();
    ///////////////////////////////////////////////////
    // DEBUG
    #if PCSEARCH_DEBUG > 0
    std::cout << "PCSearch with params:\n";
    params.print();
    #endif
    ///////////////////////////////////////////////////

    assert(params.maxDistance >= 0);
    assert(params.minDistance >= 0);
    assert(params.maxDistance >= params.minDistance);

    // Create a subgraph with nodes that are gauranteed to be on a path satisfying the gap constraints
    Vertex * pX = params.pStartVertex;
    VertexID pXid = pX->getID();
    Vertex * pY = params.pEndVertex;
    VertexID pYid = pY->getID();

    // TODO: Need to correctly handle the case if pX == pY. For now, return nothing.
    if( pX == pY)
    {
    //    return true;
    }

    // Note: params.maxDistance is the distance from the start of X to the start of Y:
    // |-----------> X               Y <----------------|
    // |<----------------------------->| maxDistance
    StringGraph * pSubgraph = makePathGraph(pGraph, pX, params.searchDir, pY, !params.goalDir, params.maxDistance);

    // makePathGraph returns NULL if there X & Y are not connected by a path satisfying path constraints
    if (pSubgraph == NULL)
    {
        // Do not add to outwalks
        return true;
    }

    assert(pSubgraph != NULL);
    #if PCSEARCH_DEBUG > 0
    std::cout << "Writing subgraph to file: temp.asqg.gz\n";
    pSubgraph->writeASQG("temp.asqg.gz");
    #endif

    // Modify the search tree parameters to operate on the subgraph
    SGSearchParams sgParams(params);
    sgParams.pStartVertex = pSubgraph->getVertex(pXid);
    sgParams.pEndVertex = pSubgraph->getVertex(pYid);
    assert(sgParams.pStartVertex != NULL);
    assert(sgParams.pEndVertex != NULL);

    // Convert the maxDistance and minDistance from to the distance expected by SGSearch,
    // which is the number of bases from the end of the start vertex to the end of the last vertex.
    size_t lY = pY->getSeqLen();
    size_t lX = pX->getSeqLen();
    sgParams.minDistance = params.minDistance -lX + lY;
    sgParams.maxDistance = params.maxDistance -lX + lY;
    sgParams.startDistance = 0;
    sgParams.allowGoalRepeat = true;
    sgParams.goalOriented = true;
    sgParams.minDistanceEnforced = true;
    sgParams.maxDistanceEnforced = true;
    bool searchComplete = SGSearch::findWalks(sgParams, exhaustive, outWalks);

    delete pSubgraph;
    return searchComplete;
}

// This algorithm will first create a subgraph that consists only of vertexes and edges that could be used
// on a path from pX to pY on a path with length less than the prescribed maxDistance.
// Then it searches for all possible valid paths in the subgraph
//
// Return true if the search completed, else return false
bool PCSearch::findWalks2(StringGraph * pGraph, SGSearchParams params, bool exhaustive, SGWalkVector& outWalks)
{
    outWalks.clear();
    ///////////////////////////////////////////////////
    // DEBUG
    #if PCSEARCH_DEBUG > 0
    std::cout << "PCSearch with params:\n";
    params.print();
    #endif
    ///////////////////////////////////////////////////

    assert(params.maxDistance >= 0);
    assert(params.minDistance >= 0);
    assert(params.maxDistance >= params.minDistance);

    // Create a subgraph with nodes that are gauranteed to be on a path satisfying the gap constraints
    Vertex * pX = params.pStartVertex;
    VertexID pXid = pX->getID();
    Vertex * pY = params.pEndVertex;
    VertexID pYid = pY->getID();

    // Note: params.maxDistance is the distance from the start of X to the start of Y:
    // |-----------> X               Y <----------------|
    // |<----------------------------->| maxDistance
    int dummy;
    EdgePtrVec allowedEdges = getPathEdges3(pX, params.searchDir, pY, !params.goalDir, params.maxDistance, dummy);
    //EdgePtrVec allowedEdges = getPathEdges2(pX, params.searchDir, pY, !params.goalDir, params.maxDistance);
   // EdgePtrVec allowedEdges = getPathEdges(pX, params.searchDir, pY, !params.goalDir, params.maxDistance);

    if (allowedEdges.size() == 0)
    {
        return true; // Search completed, found no paths
    }

    // Modify the search parameters for SGSearch:
    //  - Convert the maxDistance and minDistance from to the distance expected by SGSearch,
    //    which is the number of bases from the end of the start vertex to the end of the last vertex.
    SGSearchParams sgParams(params);
    size_t lY = pY->getSeqLen();
    size_t lX = pX->getSeqLen();
    sgParams.minDistance = params.minDistance -lX + lY;
    sgParams.maxDistance = params.maxDistance -lX + lY;
    sgParams.startDistance = 0;
    sgParams.allowGoalRepeat = true;
    sgParams.goalOriented = true;
    sgParams.minDistanceEnforced = true;
    sgParams.maxDistanceEnforced = true;
    sgParams.enforceAllowedEdges = true;
    sgParams.pAllowedEdges = &allowedEdges;
    bool searchComplete = SGSearch::findWalks(sgParams, exhaustive, outWalks);
    return searchComplete;
}


// Perform a one sided BFS from pX to pY to find all possible walks
// consistent with the minDistance and maxDistance.
//
// Return true if the search completed, else return false
bool PCSearch::findWalksOneSidedBFS(StringGraph * pGraph, SGSearchParams params, bool exhaustive, SGWalkVector& outWalks)
{
    outWalks.clear();

    ///////////////////////////////////////////////////
    // DEBUG
    #if PCSEARCH_DEBUG > 0
    std::cout << "PCSearch with params:\n";
    params.print();
    #endif
    ///////////////////////////////////////////////////

    assert(params.maxDistance >= 0);
    assert(params.minDistance >= 0);
    assert(params.maxDistance >= params.minDistance);

    // Create a subgraph with nodes that are gauranteed to be on a path satisfying the gap constraints
    Vertex * pX = params.pStartVertex;
    VertexID pXid = pX->getID();
    Vertex * pY = params.pEndVertex;
    VertexID pYid = pY->getID();

    // Note: params.maxDistance is the distance from the start of X to the start of Y:
    // |-----------> X               Y <----------------|
    // |<----------------------------->| maxDistance

    // Modify the search parameters for SGSearch:
    //  - Convert the maxDistance and minDistance from to the distance expected by SGSearch,
    //    which is the number of bases from the end of the start vertex to the end of the last vertex.
    SGSearchParams sgParams(params);
    size_t lY = pY->getSeqLen();
    size_t lX = pX->getSeqLen();
    sgParams.minDistance = params.minDistance -lX + lY;
    sgParams.maxDistance = params.maxDistance -lX + lY;
    sgParams.startDistance = 0;
    sgParams.allowGoalRepeat = true;
    sgParams.goalOriented = true;
    sgParams.minDistanceEnforced = true;
    sgParams.maxDistanceEnforced = true;
    bool searchComplete = SGSearch::findWalks(sgParams, exhaustive, outWalks);
    return searchComplete;
}

// NOTE: The minDistance and maxDistance in params is the GAP size between the start and end vertex.
// This is a different definition of distance than that which is used in SGSearch.
//
// This algorithm will first create a subgraph that consists only of vertexes and edges that could be used
// on a path from pX to pY on a path with length less than the prescribed maxDistance.
// Then it searches for all possible valid paths in the subgraph
//
// Return true if the search completed, else return false
bool PCSearch::findWalksDFS(StringGraph * pGraph, SGSearchParams params, bool exhaustive, SGWalkVector& outWalks, int& shortestDistance)
{

    outWalks.clear();

    ///////////////////////////////////////////////////
    // DEBUG
    #if PCSEARCH_DEBUG > 0
    std::cout << "PCSearch with params:\n";
    params.print();
    #endif
    ///////////////////////////////////////////////////

    assert(params.maxDistance >= 0);
    assert(params.minDistance >= 0);
    assert(params.maxDistance >= params.minDistance);

    // Create a subgraph with nodes that are gauranteed to be on a path satisfying the gap constraints
    Vertex * pX = params.pStartVertex;
    VertexID pXid = pX->getID();
    Vertex * pY = params.pEndVertex;
    VertexID pYid = pY->getID();

    // Note: params.maxDistance is the distance from the start of X to the start of Y:
    // |-----------> X               Y <----------------|
    // |<----------------------------->| maxDistance
    EdgePtrVec allowedEdges = getPathEdges3(pX, params.searchDir, pY, !params.goalDir, params.maxDistance, shortestDistance);

    if (allowedEdges.size() == 0)
    {
        return true; // Search completed, found no paths
    }

    size_t MAX_STEPS = 10000;

    DFS::DFSearch dfSearch(params.pStartVertex, params.searchDir, params.pEndVertex,
                      params.goalDir, params.maxDistance, params.minDistance, params.pDFSAllocator, allowedEdges);

    while(dfSearch.stepOnce() && dfSearch.getNumSteps() < MAX_STEPS)
    { }

    #if PCSEARCH_DEBUG > 0 
    if(!dfSearch.foundAll())
    {
        std::cout << pXid << "-" << pYid
                  << " numSteps=" << dfSearch.getNumSteps()
                  << " numWalks=" << dfSearch.getNumWalks() << "\n";
    }
    #endif

    bool foundAll = dfSearch.foundAll();
    dfSearch.swapWalks(outWalks);
    return foundAll;
}



// Helper functions to recover the vertex/edge from original graph
// using the vertex/edge in subgraph

Vertex * getOrigVertex(const StringGraph * pGraph, const Vertex * pNewVertex)
{
    VertexID id = pNewVertex->getID();
    return pGraph->getVertex(id);
}

Edge * getOrigEdge(const StringGraph * pGraph, const Edge * pNewEdge)
{
    Vertex * pStartNew = pNewEdge->getStart();
    Vertex * pEndNew = pNewEdge->getEnd();
    Vertex * pStartOrig = getOrigVertex(pGraph, pStartNew);
    Vertex * pEndOrig = getOrigVertex(pGraph, pEndNew);

    // Update the description of the edge by replacing the pointer to the ending Vertex
    EdgeDesc desc = pNewEdge->getDesc();
    desc.pVertex = pEndOrig;

    // Find the edge in the original starting Vertex
    Edge * newEdge = pStartOrig->getEdge(desc);
    return newEdge;
}

// Convert a walk on the subgraph to a walk in the original graph
SGWalk convertWalk(const StringGraph * pGraph, const SGWalk& walkIn)
{
    EdgePtrVec newEdges = walkIn.getEdges();
    EdgePtrVec origEdges(newEdges);
    for(size_t i = 0; i < newEdges.size(); i++)
    {
        origEdges[i] = getOrigEdge(pGraph, newEdges[i]);
    }
    SGWalk origWalk = SGWalk(origEdges, walkIn.isIndexed());
    return origWalk;
}
