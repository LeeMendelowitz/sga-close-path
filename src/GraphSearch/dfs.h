// Depth first search
#ifndef DFS_H
#define DFS_H

#include <stack>
#include <set>
#include <vector>
#include <deque>

#include "SGUtil.h"
#include "SGWalk.h"

namespace DFS{

typedef std::pair<const Vertex *, EdgeDir> VDirPair;

class SearchEntry
{
    public:

    SearchEntry(const Vertex* pN, EdgeDir d, int pos) : pNode(pN), dir(d), startPos(pos), hasSolution(false), eIdx(-1) {};

    Edge * getNextEdge()
    {
        eIdx++;
        if (eIdx >= edges.size())
            return NULL;
        return edges[eIdx];
    }
    
    Edge * getCurEdge() const
    {
        if (eIdx >= edges.size())
            return NULL;
        return edges[eIdx];
    }
    const Vertex* pNode;
    const EdgeDir dir;
    EdgePtrVec edges; // possible edges that can be taken from this node
    int startPos;
    bool hasSolution;
    size_t eIdx; // index of the current edge
};

class DFSearch
{
    typedef std::deque<SearchEntry * > SearchStack;

    public:

    DFSearch(const Vertex * pStart, EdgeDir startDir, const Vertex * pEnd, EdgeDir endDir, int maxDistance,
             int minDistance, EdgePtrVec& allowableEdges);
    ~DFSearch()
    {
        for(SearchStack::iterator iter = stack_.begin();
            iter != stack_.end();
            iter++)
        {
            delete *iter;
        }
        stack_.clear();
    }

    bool stepOnce();
    bool foundAll() const { return foundAll_; }
    SGWalkVector getWalks() const { return walksFound_; }
    size_t getNumWalks() const {return walksFound_.size(); }
    size_t getNumSteps() const { return numSteps_; }
    
    private:

    SearchEntry * makeSearchEntry(const Vertex * pNode, EdgeDir dir, int pos);
    void gatherEdges(EdgePtrVec& edges);
    void stackPush(SearchEntry * pEntry);
    void stackPop();

    SearchStack stack_; 
    std::set<VDirPair> curNodes_; // nodes of the current path
    const Vertex* pStart_;
    const EdgeDir startDir_;
    const Vertex* pEnd_;
    const EdgeDir endDir_; // orientation of goal node.
    int maxDistance_;
    int minDistance_;
    EdgePtrVec allowableEdges_;
    size_t numSteps_;
    bool foundAll_;
    const size_t MAX_WALKS;
    std::vector<SGWalk> walksFound_; // number of walks found
};


}; // end DFS namespace
#endif
