#include "EdgeGenerator.h"
#include "Subgraph.h"

#include <vector>
#include <queue>
#include <set>
//#include <pair>
#include <cassert>
#include <iostream>
#include <limits>


#define PATHS_DEBUG 0
#define BFS_DEBUG 0

using namespace std;

class SearchEntry
{
    public:

    SearchEntry( Vertex * v, EdgeDir d, int o) :
        pVertex(v), dir(d), startPos(o) {};

    SearchEntry( const SearchEntry& se) :
        pVertex(se.pVertex), dir(se.dir), startPos(se.startPos) {};

    Vertex * pVertex;
    EdgeDir dir; // edge taken to enter the node
    int startPos; // position of the start of this node
};

ostream& operator<<(ostream& os, const SearchEntry& se)
{
    os << "Vertex: " << se.pVertex->getID()
       << " Dir: " << se.dir
       << " Start: " << se.startPos;
    return os;
}


ostream& operator<<(ostream& os, EdgePtrVec& evec)
{
    size_t N = evec.size();

    for (unsigned int i = 0; i < N; i++)
    {
        Edge * pEdge = evec[i];
        os << pEdge->getStartID() << " " << pEdge->getEndID() << " " 
           << pEdge->getDir() << " " << pEdge->getComp() << "\n";
        //os << *evec[i] << "\n";
    }
    return os;
}

// Perform a bounded BFS to collect edges starting from pVertex in direction dir.
// The end of pVertex is considered to be offset 0.
// If dir == ED_SENSE, looking for path that exits 3' end
// If dir == ED_ANTISENSE, looking ofr path that exits 5' end
EdgePtrVec boundedBFS(Vertex * pVertex, EdgeDir dir, size_t maxNodes, int maxDistance)
{
    
    EdgePtrVec edges;

    typedef pair<Vertex *, EdgeDir> VDirPair;
    typedef pair<Vertex *, int> VIntPair;
    typedef queue<SearchEntry> SearchQueue;

    // Maintain a set of seen vertices and their orientation.
    set<VDirPair> seen;

    // Maintain a queue of vertices and their starting distance and orientation.
    SearchQueue vQueue;

    // Add the first vertex to the queue
    vQueue.push(SearchEntry(pVertex, dir, -pVertex->getSeqLen()));
    seen.insert(VDirPair(pVertex, dir));

    //while ( !vQueue.empty() && seen.size() < maxNodes) {
    while ( !vQueue.empty() ) {

        size_t N = vQueue.size();
        for(unsigned int i=0; i<N; i++)
        {
            // Get the next vertex
            SearchEntry se = vQueue.front();
            vQueue.pop();
            Vertex * pVertex = se.pVertex;
            int endPos = se.startPos + pVertex->getSeqLen();

            #if BFS_DEBUG!=0
            cout << "**************************************\n"
                 << "Searching edges from: " << se << endl;
            #endif

            // Get edges from the next vertex
            EdgePtrVec nextEdges = pVertex->getEdges(se.dir);
            EdgePtrVec::iterator iEdge = nextEdges.begin();
            const EdgePtrVec::iterator E = nextEdges.end();
            for(; iEdge != E; iEdge++) {
                Edge * pEdge = *iEdge;
                assert(pEdge->getStart() == pVertex);
                Vertex * pNextVertex = pEdge->getEnd();
                int nextStartPos = endPos - pEdge->getMatchLength();

                // Example of edge direction:
                // A |------>
                //        <-----| B

                // The edge A->B has directon: ED_SENSE
                // The edge B->A has direction: ED_SENSE

                // So, if we take the edge A->B, then B is ED_ANTISENSE! We figure this out
                // by negating the direction of edge B->A
                            
                EdgeDir nextDir = !pEdge->getTwinDir();
                SearchEntry se(pNextVertex, nextDir, nextStartPos);
                VDirPair nextVertexOrientation(pNextVertex, nextDir);

                bool alreadySeen = (seen.count(nextVertexOrientation)>0);
                bool tooFar = (nextStartPos > maxDistance);
                if ( !alreadySeen && !tooFar ) {
                     seen.insert(nextVertexOrientation);
                     vQueue.push(se);
                     edges.push_back(pEdge);
                     #if BFS_DEBUG != 0
                     cout << "Adding edge from " << pVertex->getID()
                          << ": " << se  << endl;
                     #endif
                } else {
                    #if BFS_DEBUG != 0
                    cout << "Skipping edge from " << pVertex->getID()
                         << " " << se 
                         << ". AlreadySeen: " << alreadySeen << " TooFar: " << tooFar << endl;
                    #endif
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
    #if PATHS_DEBUG!=0
    cout << "X Edges: " << xyEdges.size() << endl;
    cout << xyEdges << endl;
    #endif

    // Search from pY
    yEdges = boundedBFS(pY, dY, maxNodes, maxDistance_2);
    #if PATHS_DEBUG!=0
    cout << "Y Edges: " << yEdges.size() << endl;
    cout << yEdges << endl;
    #endif

    xyEdges.insert(xyEdges.end(), yEdges.begin(), yEdges.end());
    #if PATHS_DEBUG!=0
    cout << "XY Edges: " << xyEdges.size() << endl;
    cout << xyEdges << endl;
    #endif

    // Create a subgraph with these edges
    StringGraph * pSubgraph = Subgraph::copyGraph(pGraph); // Make empty subgraph
    #if PATHS_DEBUG!=0
    cout << "Subgraph Created:" << endl;
    pSubgraph->stats();
    #endif

    // Only add unique edges
    Subgraph::copyEdgesToSubgraph(pSubgraph, pGraph, xyEdges);

    #if PATHS_DEBUG!=0
    cout << "Copied Edges:" << endl;
    pSubgraph->stats();
    #endif

    ///////////////////////////////////////////////////////////////
    // Remove any nodes that are not on a path from pX to pY
    pX = pSubgraph->getVertex(xId);
    pY = pSubgraph->getVertex(yId);
    xEdges = boundedBFS(pX, dX, maxNodes, maxDistance);
    yEdges = boundedBFS(pY, dY, maxNodes, maxDistance);
    #if PATHS_DEBUG!=0
    cout << " Searching predecessors/successors: \n"
         << " x successors: "
         << xEdges
         << "\n"
         << " y predecssors: "
         << yEdges
         << endl;
    #endif
    cout << yEdges << endl;

    typedef set<Vertex *> VertexPtrSet;
    VertexPtrSet xVertexSet, yVertexSet, xyVertexSet;
    VertexPtrSet * pSet1, * pSet2;
    for(unsigned int i=0; i < xEdges.size(); i++)
    {
        Edge * pEdge = xEdges[i];
        xVertexSet.insert(pEdge->getStart()); // Add starting node of the edge
        xVertexSet.insert(pEdge->getEnd()); // Add ending node of the edge
    }
    for(unsigned int i=0; i < yEdges.size(); i++)
    {
        Edge * pEdge = yEdges[i];
        yVertexSet.insert(pEdge->getStart()); // Add the starting node of the edge
        yVertexSet.insert(pEdge->getEnd()); // Add the ending node of the edge
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

    // Remove any nodes that are not in both the xVertexSet and yVertexSet.
    VertexPtrSet::const_iterator ib = pSet1->begin();
    VertexPtrSet::const_iterator ie = pSet1->end();
    pSubgraph->setColors(GC_BLACK);
    // Set those nodes in the intersection to WHITE
    for(; ib != ie; ib++) {
        Vertex * pVertex = *ib;
        if (pSet2->count(pVertex)==1)
        {
            pVertex->setColor(GC_WHITE);
        }
    }
    pSubgraph->sweepVertices(GC_BLACK);
    cout << "Subgraph nodes removed:" << endl;

    #if PATHS_DEBUG!=0
    pSubgraph->stats();
    #endif
  
    return pSubgraph;
}
