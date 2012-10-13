#include "EdgeGenerator.h"
#include "Subgraph.h"

#include <vector>
#include <queue>
#include <set>
//#include <pair>
#include <cassert>
#include <iostream>
#include <limits>

// Perform a bounded BFS to collect edges starting from pVertex in direction dir.
EdgePtrVec boundedBFS(Vertex * pVertex, EdgeDir dir, size_t maxNodes, int maxDistance) {

    using namespace std;
    
    EdgePtrVec edges;

    typedef pair<Vertex *, bool> VBoolPair;
    typedef pair<Vertex *, int> VIntPair;
    typedef queue<VIntPair> VertexQueue;

    // Maintain a set of seen vertices and their orientation.
    set<VBoolPair> seen;

    // Maintain a queue of vertices and their starting distance.
    VertexQueue vQueue;

    // Add the first vertex to the queue
    vQueue.push(VIntPair(pVertex, -pVertex->getSeqLen()));
    bool orientation = (dir == ED_SENSE); // orientation of first vertex
    seen.insert(VBoolPair(pVertex, orientation));



    //while ( !vQueue.empty() && seen.size() < maxNodes) {
    while ( !vQueue.empty() ) {

        size_t N = vQueue.size();
        for(unsigned int i=0; i<N; i++)
        {
            // Get the next vertex
            VIntPair vi = vQueue.front();
            vQueue.pop();
            Vertex * pVertex = vi.first;
            int startPos = vi.second;
            int endPos = startPos + pVertex->getSeqLen();

            //cout << "Searching edges from vertex " << pVertex->getID() << ", startPos " << startPos << endl;

            // Get edges from the next vertex
            EdgePtrVec nextEdges = pVertex->getEdges(dir);
            EdgePtrVec::iterator iEdge = nextEdges.begin();
            const EdgePtrVec::iterator E = nextEdges.end();
            for(; iEdge != E; iEdge++) {
                Edge * pEdge = *iEdge;
                assert(pEdge->getStart() == pVertex);
                Vertex * pNextVertex = pEdge->getEnd();
                int nextStartPos = endPos - pEdge->getMatchLength();
                VBoolPair nextVertexOrientation(pNextVertex, pEdge->getTwinDir()==ED_SENSE);
                VIntPair nextVertexPosition(pNextVertex, nextStartPos);

                bool alreadySeen = (seen.count(nextVertexOrientation)>0);
                bool tooFar = (nextStartPos > maxDistance);
                if ( !alreadySeen &&
                     !tooFar ) {
                     seen.insert(nextVertexOrientation);
                     vQueue.push(nextVertexPosition);
                     edges.push_back(pEdge);
                     //cout << "Adding edge from " << pVertex->getID() << " to " << pNextVertex->getID() << ", startPos " << nextStartPos << endl;
                } else {
                    /*
                    cout << "Skipping edge from " << pVertex->getID() << " to " << pNextVertex->getID()
                         << ", startPos " << nextStartPos
                         << " AlreadySeen: " << alreadySeen << " TooFar: " << tooFar << endl;
                     */
                }
                

            }
        }
    }
    return edges;
}


// Make a subgraph of nodes that are on paths from Vertex pX to pY, distance less than maxDistance
// dX: direction of edge out of pX on walk to pY.
// dY: direction of edge out of pY on walk to pX.
// Case 1: pX Forward, pY Reverse, then dX = ED_SENSE, dY = ED_SENSE |--->.......<----|
// Case 2: pX Forward, pY Forward, then dX = ED_SENSE, dY = ED_ANTISENSE     |--->.......|---->
// Case 3: pX Reverse, pX Forward, then dX = ED_ANTISENSE, dY = ED_ANTISENSE <---|......|----->
// Case 4: pX Reverse, pY Reverse, then dX = ED_ANTISENSE, dY = ED_SENSE  <----|......<----|
StringGraph * makePathGraph(StringGraph * pGraph, Vertex * pX, EdgeDir dX, Vertex * pY, EdgeDir dY, int maxDistance)
{
    using namespace std;

    cout << "In MakePath Graph!" << endl;

    VertexID xId = pX->getID();
    VertexID yId = pY->getID();

    int maxDistance_2 = (maxDistance + 1)/2;
    //size_t maxNodes = numeric_limits<set::size_type>::max();
    size_t maxNodes = 10000;

    ///////////////////////////////////////////////////////////////
    // Create subgraph using BFS from pX and pY

    EdgePtrVec xEdges, yEdges, xyEdges;
    // Search from pX
    xyEdges = boundedBFS(pX, dX, maxNodes, maxDistance_2);
    // Search from pY
    yEdges = boundedBFS(pY, dY, maxNodes, maxDistance_2);
    xyEdges.insert(xyEdges.end(), yEdges.begin(), yEdges.end());
    // Create a subgraph with these edges
    StringGraph * pSubgraph = Subgraph::copyGraph(pGraph);
    Subgraph::copyEdgesToSubgraph(pSubgraph, pGraph, xyEdges);
    cout << "Subgraph Created:" << endl;
    pSubgraph->stats();

    ///////////////////////////////////////////////////////////////
    // Remove any nodes that are not on a path from pX to pY
    pX = pSubgraph->getVertex(xId);
    pY = pSubgraph->getVertex(yId);
    xEdges = boundedBFS(pX, dX, maxNodes, maxDistance);
    yEdges = boundedBFS(pY, dY, maxNodes, maxDistance);

    typedef set<Vertex *> VertexPtrSet;
    VertexPtrSet xVertexSet, yVertexSet, xyVertexSet;
    VertexPtrSet * pSet1, * pSet2;
    for(unsigned int i=0; i < xEdges.size(); i++)
    {
        Edge * pEdge = xEdges[i];
        xVertexSet.insert(pEdge->getStart());
        xVertexSet.insert(pEdge->getEnd());
    }
    for(unsigned int i=0; i < yEdges.size(); i++)
    {
        Edge * pEdge = yEdges[i];
        yVertexSet.insert(pEdge->getStart());
        yVertexSet.insert(pEdge->getEnd());
    }
    assert(xVertexSet.count(pX)==1);
    assert(yVertexSet.count(pY)==1);

    if (xVertexSet.size() < yVertexSet.size() ) {
        pSet1 = &xVertexSet;
        pSet2 = &yVertexSet;
    } else {
        pSet1 = &yVertexSet;
        pSet2 = &xVertexSet;
    }

    VertexPtrSet::const_iterator ib = pSet1->begin();
    VertexPtrSet::const_iterator ie = pSet1->end();
    for(; ib != ie; ib++) {
        Vertex * pVertex = *ib;
        if (pSet2->count(pVertex)==0)
        {
            pVertex->setColor(GC_BLACK);
        }
    }
    pSubgraph->sweepVertices(GC_BLACK);
    cout << "Subgraph nodes removed:" << endl;
    pSubgraph->stats();
  
    // For now, return null
    return NULL;
}


