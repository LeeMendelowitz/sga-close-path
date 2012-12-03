// PairClosure - Algorithms and data structures
// for searching a string graph for paths between two vertexes.


#ifndef PAIRCLOSURE_H
#define PAIRCLOSURE_H

#include "SGUtil.h"
#include "SGSearch.h"
#include "SGWalk.h"
#include "GraphSearchTree.h"
#include <deque>

namespace PCSearch
{
    // Find walks between two vertexes with constraints on the allowed gap
    // size between two vertexes.
    // The params should contain the maxDistance and minDistance for the gap size.
    bool findWalks(StringGraph * pGraph,
                   SGSearchParams params,
                   bool exhaustive,
                   SGWalkVector& outWalks);

    // Find walks without creating a subgraph. Instead, use the allowableEdgeSet in
    // params
    bool findWalks2(StringGraph * pGraph,
                    SGSearchParams params,
                    bool exhaustive,
                    SGWalkVector& outWalks);
};

#endif
