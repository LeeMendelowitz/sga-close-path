#include "PCSearch.h"
#include "EdgeGenerator.h"
#include <queue>
//#include <iostream>


// TO DO: 
// We must convert any walks discovered to point to the nodes in the original string graph,
// since we delete the subGraph (which frees the nodes and edges).

// Return true if the search completed, else return false
// The minDistance and maxDistance in params is the gap size between the start and end vertex.
// This is a different definition of distance than that which is used in SGSearch.
// This algorithm will first create a subgraph that consists only of vertexes and edges that could be used
// on a path from pX to pY on a path with length less than the prescribed maxDistance.
// Then it searches for all possible valid paths in the subgraph
bool PCSearch::findWalks(StringGraph * pGraph, SGSearchParams params, bool exhaustive, SGWalkVector& outWalks)
{
    ///////////////////////////////////////////////////
    // DEBUG
    std::cout << "PCSearch with params:\n";
    params.print();
    ///////////////////////////////////////////////////

    // Create a subgraph with nodes that are gauranteed to be on a path satisfying the gap constraints
    Vertex * pX = params.pStartVertex;
    VertexID pXid = pX->getID();
    Vertex * pY = params.pEndVertex;
    VertexID pYid = pY->getID();

    // Note: params.maxDistance is the gap between pX and pY. A negative gap implies an overlap.
    StringGraph * pSubgraph = makePathGraph(pGraph, pX, params.searchDir, pY, !params.goalDir, params.maxDistance);

    // Modify the search tree parameters to operate on the subgraph
    SGSearchParams sgParams(params);
    sgParams.pStartVertex = pSubgraph->getVertex(pXid);
    sgParams.pEndVertex = pSubgraph->getVertex(pYid);

    // If the start or end vertex is not in the subgraph, then no valid path exists.
    if ( sgParams.pStartVertex == NULL ||
         sgParams.pEndVertex == NULL)
    {
        delete pSubgraph;
        return true;
    }

    // Convert the maxDistance and minDistance from gap sizes to the distance expected by SGSearch,
    // which is the number of bases from the end of the start vertex to the end of the last vertex.
    size_t lY = pY->getSeqLen();
    sgParams.minDistance = params.minDistance + lY;
    sgParams.maxDistance = params.maxDistance + lY;
    /*
    sgParams.allowGoalRepeat = true;
    sgParams.goalOriented = true;
    sgParams.minDistanceEnforced = true;
    sgParams.maxDistanceEnforced = true;
    */
    bool searchComplete = SGSearch::findWalks(sgParams, exhaustive, outWalks);

    delete pSubgraph;
    return searchComplete;
}
