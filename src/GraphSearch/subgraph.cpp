#include <set>


#include "subgraph.h"
#include "SGUtil.h"
#include "SGAlgorithms.h"


// Code to handle creating a subgraph StringGraph from an existing StringGraph
StringGraph * Subgraph::copyGraph(const StringGraph * pGraph)
{
    StringGraph * pCopy = new StringGraph;

    // Set the graph parameters to match the main graph
    pCopy->setContainmentFlag(pGraph->hasContainment());
    pCopy->setTransitiveFlag(pGraph->hasTransitive());
    pCopy->setMinOverlap(pGraph->getMinOverlap());
    pCopy->setErrorRate(pGraph->getErrorRate());

    return pCopy;
}

// We do not directly add the vertex pointer to the graph, as each graph
// manages its set of vertices and we do not want to double-free the vertices
void Subgraph::copyVertexToSubgraph(StringGraph* pSubgraph, const Vertex* pVertex)
{
    // Make sure the vertex hasn't been added yet
    if(pSubgraph->getVertex(pVertex->getID()) == NULL)
    {
        Vertex* pCopy = new(pSubgraph->getVertexAllocator()) Vertex(pVertex->getID(), pVertex->getSeq().toString());
        pSubgraph->addVertex(pCopy);
    }
}

// Add a copy of the edge and its twin to the subgraph
void Subgraph::addEdgeToSubgraph(StringGraph* pSubGraph, const Edge* pEdge)
{
    const Vertex * pX = pEdge->getStart();
    const Vertex * pY = pEdge->getEnd();
    copyVertexToSubgraph(pSubGraph, pX);
    copyVertexToSubgraph(pSubGraph, pY);
    Overlap ovr = pEdge->getOverlap(); // contains start Id, end Id, and coordinates of the overlap
    SGAlgorithms::createEdgesFromOverlap(pSubGraph, ovr, true); // This will create two edges, one for each direction
}

// Copy edges in the eVec into the subGraph. Take care to copy an edge only once. 
void Subgraph::copyEdgesToSubgraph(StringGraph * pSubGraph, const EdgePtrVec& eVec)
{
    std::set<const Edge *> edgesAdded;
    EdgePtrVec::const_iterator i = eVec.begin();
    const EdgePtrVec::const_iterator E = eVec.end();
    for(; i != E; i++) {
        const Edge * pEdge = *i;
        const Edge * pTwinEdge = pEdge->getTwin();
        size_t pEdgeCount = edgesAdded.count(pEdge);
        size_t pTwinEdgeCount = edgesAdded.count(pTwinEdge);
        assert( pEdgeCount == pTwinEdgeCount);
        if (pEdgeCount == 0)
        {
            addEdgeToSubgraph(pSubGraph, pEdge);
            edgesAdded.insert(pEdge);
            edgesAdded.insert(pTwinEdge);
        }
    }
}
