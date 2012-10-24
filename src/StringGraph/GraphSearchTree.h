//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// GraphSearchTree - Generic data structure used for implementing
// a breadth-first search of a bidirectional graph. It is designed
// to return all possible walks between the given start
// and end vertices, up to a given distance. Used to search a
// string graph or scaffold graph.
//
#ifndef GRAPHSEARCHTREE_H
#define GRAPHSEARCHTREE_H

#include "Bigraph.h"
#include "SGWalk.h"
#include <deque>
#include <queue>
#include <iostream>

#define GRAPHSEARCH_DEBUG 0
#define GRAPHDELETE_DEBUG 0

// Parameters guiding a GraphSearchTree Search
// This structure sets a default for optional search parameters
template<typename VERTEX>
class GraphSearchParams
{

    public:

    // Basic constructor allows specification of only the required parameters. All other parameters should be set by hand.
    GraphSearchParams(VERTEX * start, VERTEX * end, EdgeDir dir, int64_t maxDist) :
        pStartVertex(start),
        pEndVertex(end),
        searchDir(dir),
        goalDir(ED_SENSE),
        startDistance(0),
        maxDistance(maxDist),
        minDistance(-1),
        allowGoalRepeat(false),
        goalOriented(false),
        minDistanceEnforced(false),
        maxDistanceEnforced(false),
        nodeLimit(10000),
        selfPrune(false)
        {};

    VERTEX* pStartVertex; 
    VERTEX* pEndVertex; // The goal vertex.
    EdgeDir searchDir; // Search from the start vertex in this direction.
    EdgeDir goalDir; // If goalDir is ED_SENSE, we enter the goal on the 5' end.
    int64_t startDistance;  // The distance of the initial node.
    int64_t maxDistance; // The max distance to search. The GraphSearchTree will not expand
                         // a node with distance > maxDistance.
    int64_t minDistance; // The minimum distance to the goal vertex.
    bool allowGoalRepeat; // Allow a goal vertex to be repeated on a walk.
    bool goalOriented; // Only accept the goal if it has orientation goalDir.
    bool minDistanceEnforced; // We only accept the endVertex as a goal if it's distance is is at least minDistance.
    bool maxDistanceEnforced; // We only accept the endVertex as a goal if it's distance is is at most maxDistanceEnforced.
    size_t nodeLimit; // The maximum number of nodes allowed in the GraphSearchTree.
    bool selfPrune; // If true, prune any search nodes which do not lead to the goal.

    void print() const
    {
        std::cout << "Start: " << pStartVertex->getID() << "\n"
             << "End: "   << ((pEndVertex==NULL) ? "NULL" : pEndVertex->getID().c_str()) << "\n"
             << "SearchDir: " << searchDir << "\n"
             << "GoalDir: " << goalDir << "\n"
             << "startDistance: " << startDistance << "\n"
             << "maxDistance: " << maxDistance << "\n"
             << "minDistance: " << minDistance << "\n"
             << "allowGoalRepeat: " << allowGoalRepeat << "\n"
             << "goalOriented: " << goalOriented << "\n"
             << "minDistanceEnforced: " << minDistanceEnforced << "\n"
             << "maxDistanceEnforced: " << maxDistanceEnforced << "\n"
             << "nodeLimit: " << nodeLimit << "\n"
             << "selfPrune: " << selfPrune << std::endl;
    }
};

template<typename VERTEX, typename EDGE, typename DISTANCE>
class GraphSearchNode
{
    // Typedefs
    public:
        typedef std::deque<GraphSearchNode<VERTEX,EDGE,DISTANCE>* > GraphSearchNodePtrDeque;
        typedef std::vector<EDGE*> _EDGEPtrVector;

    public:
        GraphSearchNode(VERTEX* pVertex, EdgeDir expandDir, GraphSearchNode* pParent, EDGE* pEdgeFromParent, int distance);
        ~GraphSearchNode();

        // Reduce the child count by 1
        void decrementChildren();

        // Create the children of this node and place pointers to their nodes
        // on the queue. Returns the number of children created;
        int createChildren(GraphSearchNodePtrDeque& outQueue, const DISTANCE& distanceFunc);

        GraphSearchNode* getParent() const { return m_pParent; }
        VERTEX* getVertex() const { return m_pVertex; }
        int64_t getDistance() const { return m_distance; }
        EDGE* getEdgeFromParent() const { return m_pEdgeFromParent; }
        int getNumChildren() const { return m_numChildren; }
        bool isGoal() const { return m_isGoal; } 
        void isGoal(bool isGoal) { m_isGoal = isGoal; }


    private:

        // data
        VERTEX* m_pVertex;
        EdgeDir m_expandDir;
        GraphSearchNode* m_pParent;
        EDGE* m_pEdgeFromParent;
        int m_numChildren;
        int64_t m_distance;
        bool m_isGoal; // True if this node is a goal vertex
};

template<typename VERTEX, typename EDGE, typename DISTANCE>
class GraphSearchTree
{
    // typedefs
    typedef GraphSearchNode<VERTEX,EDGE,DISTANCE> _SearchNode;
    typedef GraphSearchParams<VERTEX> _SearchParams;
    typedef typename _SearchNode::GraphSearchNodePtrDeque _SearchNodePtrDeque;
    typedef typename std::set<_SearchNode*> _SearchNodePtrSet;
    typedef std::vector<EDGE*> WALK; // list of edges defines a walk through the graph
    typedef std::vector<WALK> WALKVector; // vector of walks
    
    typedef std::vector<VERTEX*> VertexPtrVector;
    typedef std::vector<VertexPtrVector> VertexPtrVectorVector;

    public:

        GraphSearchTree(const _SearchParams& params);

        ~GraphSearchTree();

        // Find connected components in the graph
        // Takes in a vector of all the vertices in the graph
        static void connectedComponents(VertexPtrVector allVertices, VertexPtrVectorVector& connectedComponents);

        // Returns true if the search has converged on a single vertex. In
        // other words, all walks from the start node share a common vertex,
        // which is represented by one of the nodes waiting expansion.
        // This function is the key to finding collapsed walks that represent
        // complex variation bubbles. It is a heavy operation however, as the
        // entire tree is search for each node in the expand queue. It should
        // only be used on small trees.
        // Returns true if the search converged and the pointer to the vertex
        // is return in pConvergedVertex.
        bool hasSearchConverged(VERTEX*& pConvergedVertex);

        // build walks to a given vertex
        template<typename BUILDER>
        void buildWalksContainingVertex(VERTEX* pTarget, BUILDER& walkBuilder);

        // Construct walks representing every path from the start node
        template<typename BUILDER>
        void buildWalksToAllLeaves(BUILDER& walkBuilder);

        // Construct walks representing every path from the start vertex to the goal vertex
        template<typename BUILDER>
        void buildWalksToGoal(BUILDER& walkBuilder);
                
        // Expand all nodes in the queue a single time
        // Returns false if the search has stopped
        bool stepOnce();

        // Check whether the search was aborted or not
        bool wasSearchAborted() const { return m_searchAborted; }

        // Count the number of search nodes in the goal queue. This gives the number 
        // of walks found thus far from the start node to the goal.
        size_t countWalksToGoal() const {return m_goalQueue.size();}

    private:

        // Search the branch from pNode to the root for pX.  
        bool searchBranchForVertex(_SearchNode* pNode, VERTEX* pX, _SearchNode*& pFoundNode) const;

        // Build the walks from the root to the leaves in the queue
        template<typename BUILDER>
        void _buildWalksToLeaves(const _SearchNodePtrDeque& queue, BUILDER& walkBuilder);

        //
        void addEdgesFromBranch(_SearchNode* pNode, WALK& outEdges);

        // Build a queue with all the leaves in it
        void _makeFullLeafQueue(_SearchNodePtrDeque& completeQueue) const;

        // Delete the leaf node pCurr. Follow path to root and delete any of parents that become childless.
        size_t deleteFromLeaf(_SearchNode * pCurr);

        // Delete from the leaf node pCurr. Follow path to root and delete any parents that become childless,
        // until a goal node or parent with child is reached. Set the stopNode to be the undeleted parent.
        size_t pruneFromLeaf(_SearchNode * pCurr, _SearchNode ** stopNode);

        // print the branch sequence
        void printBranch(_SearchNode* pNode) const;

        // We keep the pointers to the search nodes in queues.
        // The goal queue contains the nodes representing the vertex we are searching for.
        // The expand queue contains nodes that have not yet been explored.
        // The done queue contains nodes that will not be expanded further.
        // Together, the expand queue and done queue represent all leafs of the tree
        // NOTE: Each node in the qoal queue is duplicated in either the expand queue or
        // done queue. This depends on whether we allow a goal vertex to be repeated
        // on a walk from start to goal (controlled by m_searchParams.allowGoalRepeat).
        _SearchNodePtrDeque m_goalQueue;
        _SearchNodePtrDeque m_expandQueue;
        _SearchNodePtrDeque m_doneQueue;
        size_t m_totalNodes; // The total number of nodes in the search tree
    
        _SearchNode* m_pRootNode;
        _SearchParams m_searchParams;

        // Flag indicating the search was aborted
        bool m_searchAborted;

        // Distance functor
        DISTANCE m_distanceFunc;
};

//
template<typename VERTEX, typename EDGE, typename DISTANCE>
GraphSearchNode<VERTEX,EDGE,DISTANCE>::GraphSearchNode(VERTEX* pVertex,
                           EdgeDir expandDir,
                           GraphSearchNode<VERTEX,EDGE,DISTANCE>* pParent,
                           EDGE* pEdgeFromParent,
                           int distance) : m_pVertex(pVertex),
                                                    m_expandDir(expandDir),
                                                    m_pParent(pParent), 
                                                    m_pEdgeFromParent(pEdgeFromParent),
                                                    m_numChildren(0),
                                                    m_isGoal(false)
{
    // Set the extension distance
    if(m_pParent == NULL)
    {
        // this is the root node with distance 0
        m_distance = 0;
    }
    else
    {
        assert(m_pEdgeFromParent != NULL);
        m_distance = m_pParent->m_distance + distance;
    }
}



// Delete this node and decrement the number of children
// in the parent node. All children of a node must
// be deleted before the parent
template<typename VERTEX, typename EDGE, typename DISTANCE>
GraphSearchNode<VERTEX,EDGE,DISTANCE>::~GraphSearchNode()
{
    assert(m_numChildren == 0);
    if(m_pParent != NULL)
        m_pParent->decrementChildren();
}

//
template<typename VERTEX, typename EDGE, typename DISTANCE>
void GraphSearchNode<VERTEX,EDGE,DISTANCE>::decrementChildren()
{
    assert(m_numChildren != 0);
    m_numChildren -= 1;
}

// creates nodes for the children of this node
// and place pointers to them in the queue.
// Returns the number of nodes created
template<typename VERTEX, typename EDGE, typename DISTANCE>
int GraphSearchNode<VERTEX,EDGE,DISTANCE>::createChildren(GraphSearchNodePtrDeque& outDeque, const DISTANCE& distanceFunc)
{
    assert(m_numChildren == 0);

    _EDGEPtrVector edges = m_pVertex->getEdges(m_expandDir);

    for(size_t i = 0; i < edges.size(); ++i)
    {
        EdgeDir childExpandDir = !edges[i]->getTwin()->getDir();
        GraphSearchNode* pNode = new GraphSearchNode(edges[i]->getEnd(), childExpandDir, this, edges[i], distanceFunc(edges[i]));
        outDeque.push_back(pNode);
        m_numChildren += 1;
    }
    return edges.size();
}

//
// GraphSearchTree
//
template<typename VERTEX, typename EDGE, typename DISTANCE>
GraphSearchTree<VERTEX,EDGE,DISTANCE>::GraphSearchTree(const _SearchParams& params) :
    m_searchParams(params),
    m_searchAborted(false)
{

    // Create the root node of the search tree
    m_pRootNode = new GraphSearchNode<VERTEX, EDGE, DISTANCE>(m_searchParams.pStartVertex,
                                                              m_searchParams.searchDir, NULL, NULL,
                                                              m_searchParams.startDistance);

    // add the root to the expand queue
    m_expandQueue.push_back(m_pRootNode);

    m_totalNodes = 1;
}


template<typename VERTEX, typename EDGE, typename DISTANCE>
GraphSearchTree<VERTEX,EDGE,DISTANCE>::~GraphSearchTree()
{
    using namespace std;
    // Delete the tree
    // We delete each leaf and recurse up the tree iteratively deleting
    // parents with a single child node. This ensure that each parent is
    // deleted after all its children

    _SearchNodePtrDeque completeLeafNodes;
    _makeFullLeafQueue(completeLeafNodes);

    #if GRAPHDELETE_DEBUG > 0
    std::cout << "GraphSearchTree Destructor"
              << "\ntotal nodes: " << m_totalNodes
              << "\ndoneQueue: " << m_doneQueue.size()
              << "\nexpandQueue: " << m_expandQueue.size()
              << "\ngoalQueue: " << m_goalQueue.size()
              << "\ncompleteLeafNodes: " << completeLeafNodes.size() << endl;
    #endif 

    size_t totalDeleted = 0;
    size_t totalNodes = m_totalNodes;
    for(typename _SearchNodePtrDeque::iterator iter = completeLeafNodes.begin(); 
                                               iter != completeLeafNodes.end();
                                               ++iter)    
    {
        _SearchNode* pCurr = *iter;

        VertexID vId = pCurr->getVertex()->getID();

        #if GRAPHDELETE_DEBUG > 0
        std::cout << "Deleting from leaf: " <<  vId  << endl;
        #endif 

        size_t numDeleted = deleteFromLeaf(pCurr);

        #if GRAPHDELETE_DEBUG > 0
        std::cout << "Deleted " << numDeleted << " from leaf: " <<  vId << endl;
        #endif 

        totalDeleted += numDeleted;

    }
    assert(totalDeleted == totalNodes);
    assert(m_totalNodes == 0);
}

// Delete the leaf node pCurr. Delete any of its parents that are childless.
template<typename VERTEX, typename EDGE, typename DISTANCE>
size_t GraphSearchTree<VERTEX,EDGE,DISTANCE>::deleteFromLeaf(_SearchNode * pCurr)
{
        // loop invariant: pCurr is a deletable node
        // the loop stops when the parent is NULL or has 
        // a child other than pCurr

        #if GRAPHDELETE_DEBUG > 0
        std::cout << "Deleting from leaf: " <<  pCurr->getVertex()->getID()
                  << "\nNumChildren: " << pCurr->getNumChildren() << std::endl;
        #endif 

        size_t totalDeleted = 0;
        assert(pCurr->getNumChildren() == 0);

        while(pCurr && pCurr->getNumChildren() == 0)
        {

            assert(pCurr->getNumChildren() == 0);
            _SearchNode* pNext = pCurr->getParent();

            #if GRAPHDELETE_DEBUG > 0
            
            std::cout << "Deleting " << pCurr->getVertex()->getID()
                      << " parent: " << ( pNext ? pNext->getVertex()->getID() : "NULL" ) << std::endl;
            #endif 
            
            delete pCurr; // decrements pNext's child count
            totalDeleted += 1;

            pCurr = pNext;

            #if GRAPHDELETE_DEBUG > 0
            std::cout << " Pcurr: " << (pCurr ? pCurr->getVertex()->getID() : "NULL") << std::endl;
            #endif 
        }

        m_totalNodes -= totalDeleted;

        #if GRAPHDELETE_DEBUG > 0
        std::cout << " DONE Deleting from leaf."
                  << " totalDeleted: " << totalDeleted
                  << " m_totalNodes: " << m_totalNodes << std::endl;
        #endif 

        return totalDeleted;
}

// Delete from the leaf node pCurr. Follow path to root and delete any parents that become childless,
// until a goal node or parent with child is reached. Set the stopNode to be the undeleted parent.
template<typename VERTEX, typename EDGE, typename DISTANCE>
size_t GraphSearchTree<VERTEX,EDGE,DISTANCE>::pruneFromLeaf(_SearchNode * pCurr, _SearchNode ** stopNode)
{
        // loop invariant: pCurr is a deletable node
        // the loop stops when the parent is NULL or has 
        // a child other than pCurr

        *stopNode = NULL;

        size_t totalDeleted = 0;
        assert(pCurr->getNumChildren() == 0);

        while(pCurr && pCurr->getNumChildren() == 0)
        {
            if (pCurr->isGoal()) break;

            assert(pCurr->getNumChildren() == 0);
            _SearchNode* pNext = pCurr->getParent();
            
            delete pCurr; // decrements pNext's child count
            totalDeleted += 1;

            pCurr = pNext;
        }

        m_totalNodes -= totalDeleted;
        *stopNode = pCurr;

        return totalDeleted;
}



// Perform one step of the BFS
template<typename VERTEX, typename EDGE, typename DISTANCE>
bool GraphSearchTree<VERTEX,EDGE,DISTANCE>::stepOnce()
{
    if(m_expandQueue.empty())
        return false;

    if(m_totalNodes > m_searchParams.nodeLimit)
    {

        std::cout << "Aborting Graph Search!! "
                  << " NodeLimit: " << m_searchParams.nodeLimit
                  << " TotalNodes: " << m_totalNodes << std::endl;

        // Move all nodes in the expand queue to the done queue
        m_doneQueue.insert(m_doneQueue.end(), m_expandQueue.begin(), m_expandQueue.end());
        m_expandQueue.clear();

        // Set a flag indicating the search was aborted
        m_searchAborted = true;
        return false;
    }

    // Iterate over the expand queue. If the path to the node
    // is outside the depth limit, move that node to the done queue. It cannot
    // yield a valid path to the goal. Otherwise, add the children of the node 
    // to the incoming queue
    _SearchNodePtrDeque incomingQueue;
    while(!m_expandQueue.empty())
    {
        _SearchNode* pNode = m_expandQueue.front();
        m_expandQueue.pop_front();

        #if GRAPHSEARCH_DEBUG > 0
        std::cout << "GraphSearchTree: Touched node " << pNode->getVertex()->getID() << std::endl;
        #endif

        if(pNode->getVertex() == m_searchParams.pEndVertex)
        {
            bool minDistanceOK = true;
            bool orientationOK = true;
            bool maxDistanceOK = true;

            // Check minimum distance to goal
            if (m_searchParams.minDistanceEnforced && pNode->getDistance() < m_searchParams.minDistance) 
                minDistanceOK = false;

            // Check maximum distance to goal
            if (m_searchParams.maxDistanceEnforced && pNode->getDistance() > m_searchParams.maxDistance) 
                maxDistanceOK = false;

            // Goal Orientation check
            if (m_searchParams.goalOriented)
            {
                orientationOK = false;
                EDGE * parentEdge = pNode->getEdgeFromParent();

                // If the starting vertex and ending vertex are the same, it is possible that
                // this node is the root of the search tree and therefore there is no parent edge
                if (parentEdge != NULL)
                {
                    EdgeDir goalOrientation = !parentEdge->getTwin()->getDir();
                    orientationOK = (goalOrientation == m_searchParams.goalDir);
                }
            }

            // If we've reached the goal and pass all requirements, add it to the goalQueue
            if (minDistanceOK && orientationOK && maxDistanceOK)
            {
                #if GRAPHSEARCH_DEBUG > 0
                std::cout << "GraphSearchTree: Reached goal and satisfied requirements! node " << pNode->getVertex()->getID() << std::endl;
                #endif

                pNode->isGoal(true);
                m_goalQueue.push_back(pNode);
                if (!m_searchParams.allowGoalRepeat)
                {
                    m_doneQueue.push_back(pNode);
                    continue;
                }
            }

        }

        if(pNode->getDistance() > m_searchParams.maxDistance)
        {
            // Path to this node is too long, expand it no further
            m_doneQueue.push_back(pNode);
        }
        else
        {
            // Add the children of this node to the queue
            int numCreated = pNode->createChildren(incomingQueue, m_distanceFunc);

            #if GRAPHSEARCH_DEBUG > 0
            std::cout << "GraphSearchTree: Created " << numCreated << " children from node "
                 << pNode->getVertex()->getID() << std::endl;
            #endif

            m_totalNodes += numCreated;

            if(numCreated == 0)
            {
                // No children created, add this node to the done queue
                m_doneQueue.push_back(pNode);
            }
        }
    }

    // Check if we should prune off any leaves in the done queue
    if (m_searchParams.selfPrune)
    {
        size_t numPruned = 0;
        
        _SearchNodePtrDeque newDoneQueue; // stores the done leafs, after pruning

        // Go over through each leaf in the done queue and prune towards the root until
        // either a parent with children is found or a goal node is found.
        // Catch the node where the pruning stopped. Check if it is a leaf, and if it is,
        // add it to the done queue.
        for(typename _SearchNodePtrDeque::iterator iter = m_doneQueue.begin();
                                                  iter != m_doneQueue.end();
                                                  iter++)
        {
            _SearchNode * pNode = *iter;
            _SearchNode * pStopNode = NULL;
            numPruned +=  pruneFromLeaf(pNode, &pStopNode);

            if (pStopNode == NULL)
            {
                // This means the entire search tree was pruned!
                // Test that this is sane.
                assert(incomingQueue.size() == 0);
                assert(iter+1 == m_doneQueue.end());
            }

            if (pStopNode && pStopNode->getNumChildren() == 0)
                newDoneQueue.push_back(pStopNode);
        }

        m_doneQueue = newDoneQueue;
    }

    m_expandQueue = incomingQueue;
    return true;
}

// Return true if all the walks from the root converge
// to one vertex (ie if the search from pX converged
// to pY, then ALL paths from pX must go through pY).
template<typename VERTEX, typename EDGE, typename DISTANCE>
bool GraphSearchTree<VERTEX,EDGE,DISTANCE>::hasSearchConverged(VERTEX*& pConvergedVertex)
{
    // Construct a set of all the leaf nodes
    _SearchNodePtrDeque completeLeafNodes;
    _makeFullLeafQueue(completeLeafNodes);

    // Search all the tree for all the nodes in the expand queue
    for(typename _SearchNodePtrDeque::iterator iter = m_expandQueue.begin(); 
                                               iter != m_expandQueue.end();
                                               ++iter)
    {
        _SearchNode* pNode = *iter;
        // If this node has the same vertex as the root skip it
        // We do not want to collapse at the root
        if(pNode->getVertex() == m_pRootNode->getVertex())
            continue;

        bool isInAllBranches = true;
        for(typename _SearchNodePtrDeque::iterator leafIter = completeLeafNodes.begin();
                                                   leafIter != completeLeafNodes.end();
                                                   ++leafIter)
        {
            // Search the current branch from this leaf node to the root
            _SearchNode* pFoundNode = NULL;
            bool isInBranch = searchBranchForVertex(*leafIter, pNode->getVertex(), pFoundNode);
            if(!isInBranch)
            {
                isInAllBranches = false;
                break;
            }
        }
    
        
        if(isInAllBranches)
        {
            pConvergedVertex = pNode->getVertex();
            return true;
        }
    }

    // search has not converted
    pConvergedVertex = NULL;
    return false;
}

// Construct walks representing every path from the start node
template<typename VERTEX, typename EDGE, typename DISTANCE>
template<typename BUILDER>
void GraphSearchTree<VERTEX,EDGE,DISTANCE>::buildWalksToAllLeaves(BUILDER& walkBuilder)
{
    // Construct a queue with all leaf nodes in it
    _SearchNodePtrDeque completeLeafNodes;
    _makeFullLeafQueue(completeLeafNodes);

    _buildWalksToLeaves(completeLeafNodes, walkBuilder);
}

// Construct walks representing every path from the start vertex to the goal vertex
template<typename VERTEX, typename EDGE, typename DISTANCE>
template<typename BUILDER>
void GraphSearchTree<VERTEX,EDGE,DISTANCE>::buildWalksToGoal(BUILDER& walkBuilder)
{
    _buildWalksToLeaves(m_goalQueue, walkBuilder);
}

// Build all the walks that contain pTarget.
template<typename VERTEX, typename EDGE, typename DISTANCE>
template<typename BUILDER>
void GraphSearchTree<VERTEX,EDGE,DISTANCE>::buildWalksContainingVertex(VERTEX* pTarget, BUILDER& walkBuilder)
{
    _SearchNodePtrDeque completeLeafNodes;
    _makeFullLeafQueue(completeLeafNodes);

    // Search upwards from each leaf until pTarget is found.
    // When it is found, insert the pointer to the search node
    // in the set
    _SearchNodePtrSet leafSet;

    // Find pTarget in each branch of the graph
    for(typename _SearchNodePtrDeque::const_iterator iter = completeLeafNodes.begin();
                                                     iter != completeLeafNodes.end();
                                                     ++iter)
    {
        _SearchNode* pFoundNode = NULL;
        searchBranchForVertex(*iter, pTarget, pFoundNode);
        assert(pFoundNode != NULL);
        leafSet.insert(pFoundNode);
    }

    // Construct all the walks to the found leaves
    completeLeafNodes.clear();
    completeLeafNodes.insert(completeLeafNodes.end(), leafSet.begin(), leafSet.end());
    _buildWalksToLeaves(completeLeafNodes, walkBuilder);
}

// Main function for constructing a vector of walks from a set of leaves
template<typename VERTEX, typename EDGE, typename DISTANCE>
template<typename BUILDER>
void GraphSearchTree<VERTEX,EDGE,DISTANCE>::_buildWalksToLeaves(const _SearchNodePtrDeque& queue, BUILDER& walkBuilder)
{
    for(typename _SearchNodePtrDeque::const_iterator iter = queue.begin();
                                                     iter != queue.end();
                                                     ++iter)
    {

        // Recursively travel the tree from the leaf to the root collecting the edges in the vector
        WALK currWalk;
        addEdgesFromBranch(*iter, currWalk);

        // Reverse the walk and write it to the output structure
        walkBuilder.startNewWalk(m_pRootNode->getVertex());
        for(typename WALK::reverse_iterator iter = currWalk.rbegin(); iter != currWalk.rend(); ++iter)
            walkBuilder.addEdge(*iter);
        walkBuilder.finishCurrentWalk();
    }
}

// Return true if the vertex pX is found somewhere in the branch 
// from pNode to the root. If it is found, pFoundNode is set
// to the furtherest instance of pX from the root.
template<typename VERTEX, typename EDGE, typename DISTANCE>
bool GraphSearchTree<VERTEX,EDGE,DISTANCE>::searchBranchForVertex(_SearchNode* pNode, VERTEX* pX, _SearchNode*& pFoundNode) const
{
    if(pNode == NULL)
    {
        pFoundNode = NULL;
        return false;
    }

    if(pNode->getVertex() == pX && pNode != m_pRootNode)
    {
        pFoundNode = pNode;
        return true;
    }
    return searchBranchForVertex(pNode->getParent(), pX, pFoundNode);
}

//
template<typename VERTEX, typename EDGE, typename DISTANCE>
void GraphSearchTree<VERTEX,EDGE,DISTANCE>::addEdgesFromBranch(_SearchNode* pNode, WALK& outEdges)
{
    // Terminate the recursion at the root node and dont add an edge
    if(pNode->getParent() != NULL)
    {
        outEdges.push_back(pNode->getEdgeFromParent());
        return addEdgesFromBranch(pNode->getParent(), outEdges);
    }
}

//
template<typename VERTEX, typename EDGE, typename DISTANCE>
void GraphSearchTree<VERTEX,EDGE,DISTANCE>::_makeFullLeafQueue(_SearchNodePtrDeque& completeQueue) const
{
    completeQueue.insert(completeQueue.end(), m_expandQueue.begin(), m_expandQueue.end());
    completeQueue.insert(completeQueue.end(), m_doneQueue.begin(), m_doneQueue.end());
}

//
template<typename VERTEX, typename EDGE, typename DISTANCE>
void GraphSearchTree<VERTEX,EDGE,DISTANCE>::printBranch(_SearchNode* pNode) const
{
    if(pNode != NULL)
    {
        std::cout << pNode->getVertex()->getID() << ",";
        printBranch(pNode->getParent());
    }
}

template<typename VERTEX, typename EDGE, typename DISTANCE>
void GraphSearchTree<VERTEX,EDGE,DISTANCE>::connectedComponents(VertexPtrVector allVertices, 
                                                                VertexPtrVectorVector& connectedComponents)
{
    // Set the color of each vertex to be white signalling its not visited
    typename VertexPtrVector::iterator iter = allVertices.begin();
    for(; iter != allVertices.end(); ++iter)
    {
        assert((*iter)->getColor() == GC_WHITE);
        (*iter)->setColor(GC_WHITE);
    }

    // 
    iter = allVertices.begin();
    for(; iter != allVertices.end(); ++iter)
    {
        // Do nothing if this vertex is already part of a CC
        if((*iter)->getColor() == GC_BLACK)
            continue;

       // Start a new CC
       VertexPtrVector currComponent;

       std::queue<VERTEX*> exploreQueue;
       (*iter)->setColor(GC_GRAY); // queued color
       exploreQueue.push(*iter);

       while(!exploreQueue.empty())
       {
            VERTEX* pCurr = exploreQueue.front();
            exploreQueue.pop();

            assert(pCurr->getColor() != GC_BLACK);
            currComponent.push_back(pCurr);
            pCurr->setColor(GC_BLACK); //done with this vertex

            // Enqueue edges if they havent been visited already
            std::vector<EDGE*> edges = pCurr->getEdges();
            for(size_t i = 0; i < edges.size(); ++i)
            {
                EDGE* pEdge = edges[i];
                VERTEX* pNext = pEdge->getEnd();
                if(pNext->getColor() == GC_WHITE)
                {
                    pNext->setColor(GC_GRAY); // queued
                    exploreQueue.push(pNext);
                }
            }
       }

        connectedComponents.push_back(currComponent);
    }

    iter = allVertices.begin();
    for(; iter != allVertices.end(); ++iter)
    {
        (*iter)->setColor(GC_WHITE);
    }

    // Sanity check
    size_t totalVertices = allVertices.size();
    size_t totalInComponents = 0;

    for(size_t i = 0; i < connectedComponents.size(); ++i)
    {
        totalInComponents += connectedComponents[i].size();
    }
    assert(totalVertices == totalInComponents);
    std::cout << "[CC] total: " << totalVertices << " num components: " << connectedComponents.size() << "\n";
    std::cout << "[CC] total vertices in components: " << totalInComponents << "\n";
}

#endif
