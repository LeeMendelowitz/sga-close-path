//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldSearch - Find walks through a scaffold graph
//
#include "ScaffoldSearch.h"
#include "ScaffoldVertex.h"
#include "SimpleAllocator.h"

typedef SimpleAllocator<ScaffoldSearchNode> ScaffoldSearchNodeAllocator;

//
ScaffoldWalkBuilder::ScaffoldWalkBuilder(ScaffoldWalkVector& outWalks) : m_outWalks(outWalks), m_pCurrWalk(NULL)
{
    
}

//
ScaffoldWalkBuilder::~ScaffoldWalkBuilder()
{
    // The pointer to the current walk should be NULL
    // or else finishCurrentWalk() was not called for the last
    // walk and the graph search algorithm has a bug
    assert(m_pCurrWalk == NULL);
}

//
void ScaffoldWalkBuilder::startNewWalk(ScaffoldVertex* pStartVertex)
{
    m_outWalks.push_back(ScaffoldWalk(pStartVertex));
    m_pCurrWalk = &m_outWalks.back();
}

//
void ScaffoldWalkBuilder::addEdge(ScaffoldEdge* pEdge)
{
    m_pCurrWalk->addEdge(pEdge);
}

//
void ScaffoldWalkBuilder::finishCurrentWalk()
{
    m_pCurrWalk = NULL;
}

//
void ScaffoldWalkBuilder::buildWalkFromEdges(ScaffoldEdgePtrVector& edgeVec)
{
    m_outWalks.push_back(ScaffoldWalk(NULL));
    m_outWalks.back().setEdges(edgeVec);
}

//
void ScaffoldSearch::findVariantWalks(ScaffoldVertex* pX, 
                                      EdgeDir initialDir, 
                                      int maxDistance,
                                      size_t maxWalks, 
                                      ScaffoldWalkVector& outWalks)
{
    (void)maxWalks;
    
    //ScaffoldSearchNodeAllocator * pAllocator = new ScaffoldSearchNodeAllocator;
    ScaffoldSearchParams params(pX, NULL, initialDir, maxDistance);
    params.nodeLimit = 10000;
    ScaffoldSearchTree searchTree(params);

    // Iteravively perform the BFS using the search tree. After each step
    // we check if the search has collapsed to a single vertex.
    bool done = false;
    while(!done)
    {
        done = !searchTree.stepOnce();
        if(searchTree.wasSearchAborted())
            break;

        ScaffoldVertex* pCollapsedVertex;
        bool isCollapsed = searchTree.hasSearchConverged(pCollapsedVertex);
        if(isCollapsed)
        {
            assert(pCollapsedVertex != NULL);
            ScaffoldWalkBuilder builder(outWalks);
            searchTree.buildWalksContainingVertex(pCollapsedVertex, &builder);
            return;           
        }
    }



    // no collapsed walk found, return empty set
    outWalks.clear();
}

int ScaffoldSearch::findCoveringWalk(const ScaffoldWalkVector& allWalks, 
                                      const ScaffoldVertexPtrVector& coverVector)
{
    for(size_t i = 0; i < allWalks.size(); ++i)
    {
        bool containsAll = true;
        for(size_t j = 0; j < coverVector.size(); ++j)
        {
            int idxInWalk = allWalks[i].findVertex(coverVector[j]);
            if(idxInWalk == -1)
            {
                containsAll = false;
            }
        }
        
        if(containsAll)
            return i;
    }
    return -1;
}


void ScaffoldSearch::findPrimaryWalks(ScaffoldVertex* pX,
                                      ScaffoldVertex* pY,
                                      EdgeDir initialDir,
                                      int maxDistance,
                                      size_t maxNodes, 
                                      ScaffoldWalkVector& outWalks)
{
    //ScaffoldSearchNodeAllocator * pAllocator = new ScaffoldSearchNodeAllocator;
    ScaffoldSearchParams params(pX, pY, initialDir, maxDistance);
    params.nodeLimit = maxNodes;
    //params.pNodeAllocator = pAllocator;
    ScaffoldSearchTree searchTree(params);

    // Iteravively perform the BFS using the search tree.
    while(searchTree.stepOnce()) { }

    // If the search was aborted, do not return any walks
    // because we do not know if there are more valid paths from pX
    // to pY that we could not find because the search space was too large
    if(!searchTree.wasSearchAborted())
    {
        // Extract the walks from the graph as a vector of edges
        ScaffoldWalkBuilder builder(outWalks);
        searchTree.buildWalksToGoal(&builder);
    }
   
   //delete pAllocator;
}

//
void ScaffoldSearch::printWalks(const ScaffoldWalkVector& walkVector)
{
    for(size_t i = 0; i < walkVector.size(); ++i)
    {
        std::cout << "walk " << i << " ";
        walkVector[i].print();
    }
}
