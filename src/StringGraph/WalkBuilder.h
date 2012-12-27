// WalkBuilder
//
// Used with the GraphSearchTree to build walks

#ifndef WALKBUILDER_H
#define WALKBUILDER_H

#include <vector>

template <typename VERTEX, typename EDGE>
class WalkBuilder
{
    typedef std::vector<EDGE*> _EdgePtrVec;

    public:

    virtual ~WalkBuilder() {}

    // These three functions must be provided by the builder object
    // the generic graph code calls these to describe the walks through
    // the graph
    virtual void startNewWalk(VERTEX * pStartVertex) = 0;
    virtual void addEdge(EDGE* pEdge) = 0;
    virtual void finishCurrentWalk() = 0;

    // Build the walk in one shot using edges from the edgeVec
    //virtual void setEdges(_EdgePtrVec& edgeVec) = 0;
};

#endif
