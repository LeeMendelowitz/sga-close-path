// Depth first search
#include <stack>
#include <set>

#include "dfs.h"
#include "SGUtil.h"
#include "SGWalk.h"

using namespace DFS;

DFSearch::DFSearch(const Vertex * pStart, EdgeDir startDir, const Vertex * pEnd, EdgeDir endDir,
                   int maxDistance, int minDistance, Allocator* pAllocator, EdgePtrVec& allowableEdges) :
    pStart_(pStart),
    startDir_(startDir),
    pEnd_(pEnd),
    endDir_(endDir),
    maxDistance_(maxDistance),
    minDistance_(minDistance),
    allowableEdges_(allowableEdges),
    numSteps_(0),
    foundAll_(false),
    MAX_WALKS(20),
    pAllocator_(pAllocator)
{ 
    assert(pAllocator_);

    // Add the starting entry to the stack    
    SearchEntry * firstEntry = makeSearchEntry(pStart, startDir, 0);
    stackPush(firstEntry);
    walksFound_.reserve(MAX_WALKS);
    stack_.reserve( (size_t) maxDistance);
}

SearchEntry * DFSearch::makeSearchEntry(const Vertex * pNode, EdgeDir dir, int pos)
{
    SearchEntry * pEntry = new(pAllocator_->alloc()) SearchEntry(pNode, dir, pos);
    // Set the edges using only the allowable edges
    const EdgePtrVec::const_iterator aeB = allowableEdges_.begin();
    const EdgePtrVec::const_iterator aeE = allowableEdges_.end();
    const EdgePtrVec::const_iterator eE = pNode->getEdgesEnd();
    EdgePtrVec::const_iterator iter = pNode->getEdgesBegin();
    pEntry->edges.reserve(eE - iter);
    for( ; iter != eE; iter++)
    {
        Edge * pEdge = *iter;
        if(pEdge->getDir() != dir) continue;
        if(binary_search(aeB, aeE, pEdge))
            pEntry->edges.push_back(pEdge);
    }
    return pEntry;
}

void DFSearch::stackPush(SearchEntry * pEntry)
{
    assert(pEntry);
    stack_.push_back(pEntry);
    curNodes_.insert(VDirPair(pEntry->pNode, pEntry->dir));
    assert(stack_.size() == curNodes_.size());
}

void DFSearch::stackPop()
{
    assert(stack_.size() == curNodes_.size());
    SearchEntry * pCurEntry = stack_.back();

    // Remove from stack
    stack_.pop_back();

    // Remove node from nodes of current walk
    std::set<VDirPair>::iterator iter = curNodes_.find(VDirPair(pCurEntry->pNode, pCurEntry->dir));
    assert(iter != curNodes_.end());
    curNodes_.erase(iter);

    // Delete the current entry
    pCurEntry->~SearchEntry();
    pAllocator_->dealloc(pCurEntry);
}

// Make a walk from the stack
void DFSearch::gatherEdges(EdgePtrVec& edges)
{
    edges.reserve(edges.size() + stack_.size());
    const SearchStack::const_iterator E = stack_.end();
    for(SearchStack::const_iterator iter = stack_.begin();
        iter != E;
        iter++)
    {
        SearchEntry * pEntry = *iter;
        edges.push_back(pEntry->getCurEdge());
    }
}


// Return true if we we took a forward step;
bool DFSearch::stepOnce()
{
    assert(stack_.size() == curNodes_.size());
    bool tookStep = false;

    if (walksFound_.size() == MAX_WALKS)
    {
        foundAll_ = false;
        return false;
    }

    while (!stack_.empty())
    {
        SearchEntry * pCurEntry = stack_.back();
        Edge * pEdge = pCurEntry->getNextEdge();
        if (pEdge == NULL)
        {
            // We are done exploring here. Backtrack.
            stackPop();
            continue;
        }

        Vertex * pNext = pEdge->getEnd();
        EdgeDir nextDir = !pEdge->getTwin()->getDir();
        int curEndPos = pCurEntry->startPos + pCurEntry->pNode->getSeqLen();
        int nextStart = curEndPos - pEdge->getMatchLength();
        if (nextStart > maxDistance_) continue;

        if (  (pNext ==  pEnd_) && (nextDir == endDir_) )
        {
            // We've reached the goal. Store this walk, if it is greater than minDistance
            if (nextStart < minDistance_) continue;
            EdgePtrVec curEdges;
            gatherEdges(curEdges);
            walksFound_.push_back(SGWalk(curEdges));
            tookStep = true;
            break;
        }

        VDirPair vDir(pNext, nextDir);
        if (curNodes_.count(VDirPair(pNext, nextDir)) != 0) continue;

        // Being here means that we have not yet reached the goal and have not reached maxDistance,
        // and the next node has not been used yet on this walk.
        SearchEntry * newEntry = makeSearchEntry(pNext, nextDir, nextStart);
        stackPush(newEntry);
        tookStep = true;
        break;
    }

    if (tookStep)
        numSteps_++;

    if (!tookStep)
    {
        assert(stack_.empty());
        foundAll_ = true;
    }

    return tookStep;
}
