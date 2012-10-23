#include "EdgeGenerator.h"
#include "Subgraph.h"

#include <vector>
#include <queue>
#include <set>
//#include <pair>
#include <cassert>
#include <iostream>
#include <limits>
#include <algorithm>

#define PATHS_DEBUG 0
#define BFS_DEBUG 0

using namespace std;


// Declaration
bool pruneGraph(StringGraph * pSubgraph,
              VertexID xId, EdgeDir dirX, int maxDistanceX,
              VertexID yId, EdgeDir dirY, int maxDistanceY);

// Information about a vertex found on BFS
// Stores the orientation of the vertex, and it's starting position
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


// Perform a bounded BFS to collect edges starting from pVertex in direction dir, up to a maximum distance.
// NOTE: This search differs from the search functionality provided in SGSearch.
// Search distance is measured from the beginning of the starting vertex.
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

// If dir == ED_SENSE, looking for path that exits 3' end
// If dir == ED_ANTISENSE, looking for path that exits 5' end
// The start of pVertex is considered to be offset 0.
// If dir == ED_SENSE, then the 5' end is at offset 0.
// If dir == ED_ANTISENSE, then the 3' end is at offset 0.
EdgePtrVec boundedBFS(Vertex * pVertex, EdgeDir dir, int maxDistance)
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
    vQueue.push(SearchEntry(pVertex, dir, 0));
    seen.insert(VDirPair(pVertex, dir));

    #if BFS_DEBUG!=0
    cout << "*************************************\n"
         << "BFS: pVertex: " << pVertex->getID() << " dir: " << dir << " maxDist: " << maxDistance << endl;
    #endif
    
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
                int nextStartPos2 = se.startPos + pEdge->getTwin()->getSeqLen();
                assert(nextStartPos2 == nextStartPos);

                // Example of vertex direction inferred from edges:
                // Say A & B have this relative orientation:
                // A |------>
                //        <-----| B
                // The edge A->B has directon: ED_SENSE
                // The edge B->A has direction: ED_SENSE
                // If we take the edge A->B, then node B is ED_ANTISENSE! We figure this out
                // by negating the direction of edge B->A
                            
                EdgeDir nextDir = !pEdge->getTwinDir();
                SearchEntry se(pNextVertex, nextDir, nextStartPos);
                VDirPair nextVertexDir(pNextVertex, nextDir);

                bool alreadySeen = (seen.count(nextVertexDir)>0);
                bool tooFar = (nextStartPos > maxDistance);

                if (!tooFar)
                {
                    seen.insert(nextVertexDir);
                    edges.push_back(pEdge);
                }

                if ( !alreadySeen && !tooFar) {
                     // Explore the successors of this vertex in the next round
                     vQueue.push(se);
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


// Make a subgraph of nodes that are on gauranteed to be on paths from Vertex pX to pY
// with distance less than maxDistance.
// See description for boundedBFS for how distance is measured.

// dX: direction of edge out of pX on walk to pY.
// dY: direction of edge out of pY on walk to pX.
// Case 1: pX Forward, pY Reverse, then dX = ED_SENSE, dY = ED_SENSE |--->.......<----|
// Case 2: pX Forward, pY Forward, then dX = ED_SENSE, dY = ED_ANTISENSE     |--->.......|---->
// Case 3: pX Reverse, pX Forward, then dX = ED_ANTISENSE, dY = ED_ANTISENSE <---|......|----->
// Case 4: pX Reverse, pY Reverse, then dX = ED_ANTISENSE, dY = ED_SENSE  <----|......<----|
StringGraph * makePathGraph(StringGraph * pGraph, Vertex * pX, EdgeDir dX, Vertex * pY, EdgeDir dY, int maxDistanceX)
{
    using namespace std;

    StringGraph * pSubgraph = Subgraph::copyGraph(pGraph); // Make empty subgraph

    // |-----------> X          Y  <------------|
    // |<-------------------------------------->| startToEnd
    // |<------------------------->| maxDistanceX
    // |<------------>| halfDistanceX
    //                |<----------------------->| halfDistanceY
    //             |<---------------------------| maxDistanceY

    if (maxDistanceX <= 0)
    {
        #if PATHS_DEBUG!=0
        cout << "Warning: maxDistanceX must be >= 0. Returning empty subgraph" << endl;
        #endif
        return pSubgraph;
    }

    // Avoid a situation where we are detecting a path from X to Y which implies
    // the containment of Y by X:
    // e.g:
    // |---------------------> X
    //           <---| Y
    if (maxDistanceX + pY->getSeqLen() < pX->getSeqLen())
    {
        #if PATHS_DEBUG!=0
        cout << "Warning: PathGraph implies containment. Returning empty subgraph" << endl;
        #endif
        return pSubgraph;

    }

    int startToEnd = maxDistanceX + pY->getSeqLen();
    int maxDistanceY = startToEnd - pX->getSeqLen();
    int maxDistance_2 = (maxDistanceX+1)/2;
    int halfDistanceX = maxDistance_2; // Half the distance from start of X to start of Y
    int halfDistanceY =  startToEnd - halfDistanceX; // Half the distance from start of Y

    assert(maxDistanceY >= 0);
    assert(halfDistanceY >= 0);

    VertexID xId = pX->getID();
    VertexID yId = pY->getID();

    ///////////////////////////////////////////////////////////////
    // Create subgraph using BFS from pX and pY

    EdgePtrVec xEdges, yEdges, xyEdges;
    // Search from pX
    xEdges = boundedBFS(pX, dX, halfDistanceX);
    #if PATHS_DEBUG!=0
    cout << "X Edges: " << xEdges.size() << endl;
    cout << xEdges << endl;
    #endif

    // Search from pY
    yEdges = boundedBFS(pY, dY, halfDistanceY);
    #if PATHS_DEBUG!=0
    cout << "Y Edges: " << yEdges.size() << endl;
    cout << yEdges << endl;
    #endif

    
    if ( xEdges.size() == 0 || 
         yEdges.size() == 0 )
    {
        // This means that either pX or pY is an island.
        // Return an empty subgraph
        #if PATHS_DEBUG!=0
        cout << "Warning: X Edges or Y Edges is empty. Returning empty subgraph." << endl;
        #endif
        return pSubgraph;
    }

    // Copy the edges found into xyEdges
    xyEdges.insert(xyEdges.end(), xEdges.begin(), xEdges.end());
    xyEdges.insert(xyEdges.end(), yEdges.begin(), yEdges.end());

    #if PATHS_DEBUG!=0
    cout << "XY Edges: " << xyEdges.size() << endl;
    cout << xyEdges << endl;
    #endif

    // Create a subgraph with these edges
    Subgraph::copyEdgesToSubgraph(pSubgraph, pGraph, xyEdges);

    #if PATHS_DEBUG!=0
    cout << "Copied Edges:" << endl;
    pSubgraph->stats();
    #endif

    int numPruneRounds = 0;
    while (true)
    {
        bool graphModified = pruneGraph(pSubgraph, xId, dX, maxDistanceX, yId, dY, maxDistanceY);
        numPruneRounds++;
        if (!graphModified) break;
    }

    #if PATHS_DEBUG!=0
    cout << "Number of prune rounds: " << numPruneRounds << endl;
    cout << "Subgraph nodes removed:" << endl;
    pSubgraph->stats();
    #endif
  
    return pSubgraph;
}

// Given a StringGraph, remove any edges that are either unreachable from X within distance maxDistanceX
// or unreachable from Y within distance minDistanceY.
// Remove any vertexes that become islands as a result of the edge removal.
// Return true if the graph was modified.
bool pruneGraph(StringGraph * pSubgraph,
              VertexID xId, EdgeDir dirX, int maxDistanceX,
              VertexID yId, EdgeDir dirY, int maxDistanceY)
{

    ///////////////////////////////////////////////////////////////
    // Remove any nodes that are not on a path from pX to pY
    Vertex * pX = pSubgraph->getVertex(xId);
    Vertex * pY = pSubgraph->getVertex(yId);

    if( (pX == NULL) || (pY == NULL) )
    {
        // Remove all nodes & edges in the graph!
        int numRemoved = 0;
        pSubgraph->setColors(GC_BLACK);
        numRemoved += pSubgraph->sweepEdges(GC_BLACK);
        numRemoved += pSubgraph->sweepVertices(GC_BLACK);
        return (numRemoved > 0);
    }

    assert(pX != NULL); assert(pY != NULL);
    EdgePtrVec xEdges = boundedBFS(pX, dirX, maxDistanceX);
    EdgePtrVec yEdges = boundedBFS(pY, dirY, maxDistanceY);
    #if PATHS_DEBUG!=0
    cout << " Searching predecessors/successors: \n"
         << " x successors:\n"
         << xEdges
         << "\n"
         << " y predecssors:\n"
         << yEdges
         << endl;
    #endif

    typedef set<Edge *> EdgePtrSet;
    EdgePtrSet xEdgeSet, yEdgeSet, xyEdgeSet;
    EdgePtrSet * pSet1, * pSet2;
    for(unsigned int i=0; i < xEdges.size(); i++)
    {
        Edge * pEdge = xEdges[i];
        xEdgeSet.insert(pEdge);
        xEdgeSet.insert(pEdge->getTwin());
    }
    for(unsigned int i=0; i < yEdges.size(); i++)
    {
        Edge * pEdge = yEdges[i];
        yEdgeSet.insert(pEdge);
        yEdgeSet.insert(pEdge->getTwin());
    }

    if (xEdgeSet.size() < yEdgeSet.size() ) {
        pSet1 = &xEdgeSet;
        pSet2 = &yEdgeSet;
    } else {
        pSet1 = &yEdgeSet;
        pSet2 = &xEdgeSet;
    }

    // Remove any edges that are not in both the xEdgeSet and yEdgeSet.
    // Remove any nodes that become islands.
    EdgePtrSet::const_iterator ib = pSet1->begin();
    EdgePtrSet::const_iterator ie = pSet1->end();
    pSubgraph->setColors(GC_BLACK);
    EdgePtrVec xyEdges2;
    // Set those nodes in the intersection to WHITE
    for(; ib != ie; ib++) {
        Edge * pEdge = *ib;
        if (pSet2->count(pEdge)==1)
        {
            xyEdges2.push_back(pEdge);
            pEdge->setColor(GC_WHITE);
            pEdge->getStart()->setColor(GC_WHITE);
            pEdge->getEnd()->setColor(GC_WHITE);
        }
    }
    int numRemoved = 0;
    numRemoved  += pSubgraph->sweepEdges(GC_BLACK);
    numRemoved  += pSubgraph->sweepVertices(GC_BLACK);
    bool graphModified = (numRemoved > 0);
    #if PATHS_DEBUG!=0
    cout << " xy edges intersection: \n"
         << xyEdges2
         << endl;
    #endif
    return graphModified;
}
