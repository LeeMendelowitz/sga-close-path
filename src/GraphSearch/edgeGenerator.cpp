#include "edgeGenerator.h"
#include "subgraph.h"

#include <vector>
#include <queue>
#include <deque>
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

    SearchEntry( const Vertex * v, EdgeDir d, int o) :
        pVertex(v), dir(d), startPos(o) {};

    SearchEntry( const SearchEntry& se) :
        pVertex(se.pVertex), dir(se.dir), startPos(se.startPos) {};

    const Vertex * pVertex;
    EdgeDir dir; // edge taken to enter the node
    int startPos; // position of the start of this node
};

class SearchEntryComparison
{
    public:
    // Returns true if lhs has a smaller priority than rhs
    bool operator() (const SearchEntry& lhs, const SearchEntry& rhs)
    {
        if (rhs.startPos < lhs.startPos)
            return true;
        return false;
    }
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

// Given a list of xEdges and yEdges, take the "intersection". More specifically, return:
// xEdgeSet: The set of edges in xEdges whose twin is in yEdges.
// yEdgeSet: The set of edges in yEdges whose twin is in xEdges.
void edgeIntersection(const EdgePtrVec& xEdges, const EdgePtrVec& yEdges, EdgePtrSet& xEdgeSet, EdgePtrSet& yEdgeSet)
{
    EdgePtrSet yEdgeSetOrig = EdgePtrSet(yEdges.begin(), yEdges.end());
    yEdgeSet.clear();
    xEdgeSet.clear();

    for(EdgePtrVec::const_iterator iter = xEdges.begin(); 
        iter != xEdges.end();
        iter++)
    {
        Edge * xEdge = *iter;
        Edge * yEdge = xEdge->getTwin();
        if (yEdgeSetOrig.count(yEdge) != 0)
        {
            xEdgeSet.insert(xEdge);
            yEdgeSet.insert(yEdge);
        }
    }
}


// Perform a bounded BFS to collect edges starting from pVertex in direction dir, up to a maximum distance.
// This is done using modified Dijkstra's. We allow a node to appear once on a path in each possible orientation.
// As in Dijkstra's. we maintain a priority_queue of nodes with their distance from source (and their orientation)
//
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
EdgePtrVec boundedBFS(const Vertex * pVertex, EdgeDir dir, int maxDistance)
{
    
    EdgePtrVec edges;

    typedef pair<const Vertex *, EdgeDir> VDirPair;
    typedef priority_queue<SearchEntry, deque<SearchEntry>, SearchEntryComparison> SearchQueue;

    // Maintain a set of seen vertices and their orientation.
    set<VDirPair> seen;

    // Maintain a queue of vertices and their starting distance and orientation.
    SearchQueue vQueue;

    // Add the first vertex to the queue
    vQueue.push(SearchEntry(pVertex, dir, 0));

    #if BFS_DEBUG!=0
    cout << "*************************************\n"
         << "BFS: pVertex: " << pVertex->getID() << " dir: " << dir << " maxDist: " << maxDistance << endl;
    #endif
    
    while ( !vQueue.empty() ) {

        size_t N = vQueue.size();
        for(unsigned int i=0; i<N; i++)
        {
            // Get the next vertex. This entry will be the next closest vertex to the source.
            SearchEntry se = vQueue.top();
            vQueue.pop();

            assert(se.startPos <= maxDistance);
            const Vertex * pVertex = se.pVertex;
            VDirPair vdir(pVertex, se.dir);
            bool alreadySeen = (seen.count(vdir)>0);

            #if BFS_DEBUG!=0
            cout << "Popped " << se << ". Already Seen: " << alreadySeen << endl;
            #endif

            if (alreadySeen) continue;

            seen.insert(vdir);

            #if BFS_DEBUG!=0
            cout << "**************************************\n"
                 << "Searching edges from: " << se << endl;
            #endif

            // Get edges from the next vertex
            EdgePtrVec nextEdges = pVertex->getEdges(se.dir);
            EdgePtrVec::iterator iEdge = nextEdges.begin();
            const EdgePtrVec::iterator E = nextEdges.end();
            int endPos = se.startPos + pVertex->getSeqLen();
            for(; iEdge != E; iEdge++) {
                Edge * pEdge = *iEdge;
                assert(pEdge->getStart() == pVertex);
                const Vertex * pNextVertex = pEdge->getEnd();

                // Take the minimum value for the next starting position.
                // (There may be two possible values in the case of an inexact overlap).
                int ol = max(pEdge->getMatchLength(), pEdge->getTwin()->getMatchLength());
                int nextStartPos = endPos - ol;

                // Example of vertex direction inferred from edges:
                // Say A & B have this relative orientation:
                // A |------>
                //        <-----| B
                // The edge A->B has directon: ED_SENSE
                // The edge B->A has direction: ED_SENSE
                // If we take the edge A->B, then node B is ED_ANTISENSE! We figure this out
                // by negating the direction of edge B->A

                EdgeDir nextDir = !pEdge->getTwinDir();
                bool tooFar = (nextStartPos > maxDistance);

                if (!tooFar)
                {
                    edges.push_back(pEdge);

                    // Check if the next vertex was already seen before adding to the vQueue
                    VDirPair vdir(pVertex, nextDir);
                    bool alreadySeen = (seen.count(vdir)>0);

                    if (!alreadySeen)
                    {
                        SearchEntry se(pNextVertex, nextDir, nextStartPos);
                        #if BFS_DEBUG != 0
                        cout << "Adding edge from " << pVertex->getID()
                          << ": " << se  << endl;
                        #endif
                        vQueue.push(se);
                    }
                } else {
                    #if BFS_DEBUG != 0
                    cout << "Skipping edge from " << pVertex->getID()
                         << " " << se 
                         << " because it is too far." << endl;
                    #endif
                }
            }
        }
    }
    return edges;
}

// Perform a bounded BFS to collect edges starting from pVertex in direction dir, up to a maximum distance.
// Only use edges in the allowableEdges set.
// This is done using modified Dijkstra's. We allow a node to appear once on a path in each possible orientation.
// As in Dijkstra's. we maintain a priority_queue of nodes with their distance from source (and their orientation)
//
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
EdgePtrVec boundedBFS(const Vertex * pVertex, EdgeDir dir, int maxDistance, const EdgePtrSet& allowableEdges)
{
    
    EdgePtrVec edges;

    typedef pair<const Vertex *, EdgeDir> VDirPair;
    typedef priority_queue<SearchEntry, deque<SearchEntry>, SearchEntryComparison> SearchQueue;

    // Maintain a set of seen vertices and their orientation.
    set<VDirPair> seen;

    // Maintain a queue of vertices and their starting distance and orientation.
    SearchQueue vQueue;

    // Add the first vertex to the queue
    vQueue.push(SearchEntry(pVertex, dir, 0));

    #if BFS_DEBUG!=0
    cout << "*************************************\n"
         << "BFS: pVertex: " << pVertex->getID() << " dir: " << dir << " maxDist: " << maxDistance << endl;
    #endif
    
    while ( !vQueue.empty() ) {

        size_t N = vQueue.size();
        for(unsigned int i=0; i<N; i++)
        {
            // Get the next vertex. This entry will be the next closest vertex to the source.
            SearchEntry se = vQueue.top();
            vQueue.pop();

            assert(se.startPos <= maxDistance);
            const Vertex * pVertex = se.pVertex;
            VDirPair vdir(pVertex, se.dir);
            bool alreadySeen = (seen.count(vdir)>0);

            #if BFS_DEBUG!=0
            cout << "Popped " << se << ". Already Seen: " << alreadySeen << endl;
            #endif

            if (alreadySeen) continue;

            seen.insert(vdir);

            #if BFS_DEBUG!=0
            cout << "**************************************\n"
                 << "Searching edges from: " << se << endl;
            #endif

            // Get edges from the next vertex
            EdgePtrVec nextEdges = pVertex->getEdges(se.dir);
            EdgePtrVec::iterator iEdge = nextEdges.begin();
            const EdgePtrVec::iterator E = nextEdges.end();
            int endPos = se.startPos + pVertex->getSeqLen();
            for(; iEdge != E; iEdge++)
            {
                Edge * pEdge = *iEdge;
                // If this edge is not in the allowable edge set, do not use it.
                if (allowableEdges.count(pEdge)==0)
                {
                    continue;
                }

                assert(pEdge->getStart() == pVertex);
                const Vertex * pNextVertex = pEdge->getEnd();

                // Take the minimum value for the next starting position.
                // (There may be two possible values in the case of an inexact overlap).
                int ol = max(pEdge->getMatchLength(), pEdge->getTwin()->getMatchLength());
                int nextStartPos = endPos - ol;

                // Example of vertex direction inferred from edges:
                // Say A & B have this relative orientation:
                // A |------>
                //        <-----| B
                // The edge A->B has directon: ED_SENSE
                // The edge B->A has direction: ED_SENSE
                // If we take the edge A->B, then node B is ED_ANTISENSE! We figure this out
                // by negating the direction of edge B->A
                            
                EdgeDir nextDir = !pEdge->getTwinDir();
                bool tooFar = (nextStartPos > maxDistance);

                if (!tooFar)
                {
                    edges.push_back(pEdge);

                    // Check if the next vertex was already seen before adding to the vQueue
                    VDirPair vdir(pVertex, nextDir);
                    bool alreadySeen = (seen.count(vdir)>0);

                    if (!alreadySeen)
                    {
                        SearchEntry se(pNextVertex, nextDir, nextStartPos);
                        #if BFS_DEBUG != 0
                        cout << "Adding edge from " << pVertex->getID()
                          << ": " << se  << endl;
                        #endif
                        vQueue.push(se);
                    }

                } else {
                    #if BFS_DEBUG != 0
                    cout << "Skipping edge from " << pVertex->getID()
                         << " " << se 
                         << " because it is too far." << endl;
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
// Return a pointer to the subgraph, or NULL if the subgraph is empty.

// dX: direction of edge out of pX on walk to pY.
// dY: direction of edge out of pY on walk to pX.
// Case 1: pX Forward, pY Reverse, then dX = ED_SENSE, dY = ED_SENSE |--->.......<----|
// Case 2: pX Forward, pY Forward, then dX = ED_SENSE, dY = ED_ANTISENSE     |--->.......|---->
// Case 3: pX Reverse, pX Forward, then dX = ED_ANTISENSE, dY = ED_ANTISENSE <---|......|----->
// Case 4: pX Reverse, pY Reverse, then dX = ED_ANTISENSE, dY = ED_SENSE  <----|......<----|
StringGraph * makePathGraph(const StringGraph * pGraph, const Vertex * pX, EdgeDir dX, const Vertex * pY, EdgeDir dY, int maxDistanceX)
{
    using namespace std;


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
        return NULL;
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
        return NULL;
    }

    StringGraph * pSubgraph = Subgraph::copyGraph(pGraph); // Make empty subgraph

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

    // Copy the edges found into xyEdges
    xyEdges.insert(xyEdges.end(), xEdges.begin(), xEdges.end());
    xyEdges.insert(xyEdges.end(), yEdges.begin(), yEdges.end());

    #if PATHS_DEBUG!=0
    cout << "XY Edges: " << xyEdges.size() << endl;
    cout << xyEdges << endl;
    #endif

    // Create a subgraph with these edges
    Subgraph::copyEdgesToSubgraph(pSubgraph, xyEdges);

    // If X or Y is not in the subgraph, then return an empty graph
    // The subgraph may not have vertex X/Y because the boundedBFS may not return
    // any edges out of X/Y which satisfy the distance bound
    if (!pSubgraph->hasVertex(xId) ||
        !pSubgraph->hasVertex(yId) )
    {
        #if PATHS_DEBUG!=0
        cout << "Subgraph does not contain either X or Y, so no path satisfying path length constraints exists."
             << "Returning NULL" << endl;
        #endif
        delete pSubgraph;
        return NULL;
    }

    #if PATHS_DEBUG!=0
    cout << "Copied Edges:" << endl;
    pSubgraph->stats();
    #endif

    int numPruneRounds = 0;
    while (true)
    {
        #if PATHS_DEBUG!=0
        cout << "**********************\n"
             << "Pruning Round " << numPruneRounds << endl;
        #endif
        bool graphModified = pruneGraph(pSubgraph, xId, dX, maxDistanceX, yId, dY, maxDistanceY);
        if (!graphModified) break;
        numPruneRounds++;
    }

    #if PATHS_DEBUG!=0
    cout << "Number of prune rounds: " << numPruneRounds << endl;
    cout << "After Subgraph nodes removed:" << endl;
    pSubgraph->stats();
    #endif

    // If X or Y is not in the subgraph, then return an empty graph
    if (!pSubgraph->hasVertex(xId) ||
        !pSubgraph->hasVertex(yId) )
    {
        #if PATHS_DEBUG!=0
        cout << "Subgraph does not contain either X or Y, so no path satisfying path length constraints exists."
             << "Returning NULL" << endl;
        #endif
        delete pSubgraph;
        return NULL;
    }
  
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
    const Vertex * pX = pSubgraph->getVertex(xId);
    const Vertex * pY = pSubgraph->getVertex(yId);

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

/////////////////////////////////////////////////////////////////////////////////////////////
// Collect edges that are on gauranteed to be on paths from Vertex pX to pY
// with distance less than maxDistance.
// See description for boundedBFS for how distance is measured.

// dX: direction of edge out of pX on walk to pY.
// dY: direction of edge out of pY on walk to pX.
// Case 1: pX Forward, pY Reverse, then dX = ED_SENSE, dY = ED_SENSE |--->.......<----|
// Case 2: pX Forward, pY Forward, then dX = ED_SENSE, dY = ED_ANTISENSE     |--->.......|---->
// Case 3: pX Reverse, pX Forward, then dX = ED_ANTISENSE, dY = ED_ANTISENSE <---|......|----->
// Case 4: pX Reverse, pY Reverse, then dX = ED_ANTISENSE, dY = ED_SENSE  <----|......<----|
EdgePtrVec getPathEdges(const Vertex * pX, EdgeDir dX, const Vertex * pY, EdgeDir dY, int maxDistanceX)
{
    using namespace std;

    EdgePtrVec pathEdges;

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
        return pathEdges;
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
        return pathEdges;
    }

    int startToEnd = maxDistanceX + pY->getSeqLen();
    int maxDistanceY = startToEnd - pX->getSeqLen();
    assert(maxDistanceY >= 0);
    EdgePtrVec xEdges, yEdges;

    // Search forward from pX
    xEdges = boundedBFS(pX, dX, maxDistanceX);
    #if PATHS_DEBUG!=0
    cout << "X Edges: " << xEdges.size() << endl;
    cout << xEdges << endl;
    #endif

    // Search backwards from pY
    yEdges = boundedBFS(pY, dY, maxDistanceY);
    #if PATHS_DEBUG!=0
    cout << "Y Edges: " << yEdges.size() << endl;
    cout << yEdges << endl;
    #endif

    // Find those edges that appear in xEdges with a twin that appears in yEdges.
    EdgePtrSet xEdgeSet, yEdgeSet;
    edgeIntersection(xEdges, yEdges, xEdgeSet, yEdgeSet);



     // Iteratively refine the list of edges, until it stops changing.
     // Some of these edges are not gauranteed to be reachable from both
     // X & Y for the given distance constraints.
     // For example, an edge may be reachable from X with one path which satisfies maxDistanceX, and from Y
     // using a different path which satisfies maxDistanceY, but there is no path from X to Y which includes this edge.
     // Allowing extraneous edges to end up in the pathEdges vector can result in a repetetive graph when we search
     // for all possible paths in PCSearch.
     size_t numEdges = xEdgeSet.size();
     assert(numEdges == yEdgeSet.size());
     #if PATHS_DEBUG!=0
     cout << "XY edges in intersection: " << numEdges << std::endl;
     #endif
     if (numEdges == 0)
        return pathEdges;
     int numPruneIterations = 0;
     while (true)
     {
        EdgePtrSet xEdgeSetNew, yEdgeSetNew;

        // Search forwards from pX
        xEdges = boundedBFS(pX, dX, maxDistanceX, xEdgeSet);
        assert(xEdges.size() <= numEdges);

        // Search backwards from pY
        yEdges = boundedBFS(pY, dY, maxDistanceY, yEdgeSet);
        assert(yEdges.size() <= numEdges);

        // Take the edge intersection
        edgeIntersection(xEdges, yEdges, xEdgeSetNew, yEdgeSetNew);

        // Determine if the edge set has changed
        size_t numEdgesNew = xEdgeSetNew.size();
        assert(numEdgesNew == yEdgeSetNew.size());

        #if PATHS_DEBUG!=0
        cout << "Prune Round: " << numPruneIterations << ": Old edge set size: " << numEdges << " New size: " << numEdgesNew << std::endl;
        #endif

        if (numEdgesNew == numEdges || numEdgesNew == 0)
            break;
        xEdgeSet = xEdgeSetNew;
        yEdgeSet = yEdgeSetNew;
        numEdges = xEdgeSet.size();
        numPruneIterations++;
     }

    pathEdges = EdgePtrVec(xEdgeSet.begin(), xEdgeSet.end());
    return pathEdges;
}
