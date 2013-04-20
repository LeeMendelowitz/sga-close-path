// Store walks for unique bundle closures
// Filter out closures which are contained in other closures
#ifndef CLOSURE_H
#define CLOSURE_H

#include <vector>
#include <string>
#include "SGWalk.h"
#include "SGUtil.h"

class ClosePathResult;

class Closure : public SGWalk
{
    public:
    Closure(const std::string& id, const SGWalk& w, int d1, int d2) :
        SGWalk(w), id_(id), d1max_(d1), d2max_(d2) { };

    // Return True if this walk contains the other
    bool contains(const Closure& other) const;
    
    // Define less than operator for sorting
    bool operator<(const Closure& other) const
    {
        return getNumEdges() < other.getNumEdges();
    };

    void colorInteriorEdges(GraphColor c) const;
    void getInteriorVertices(VertexPtrVec * pVec) const;

    std::string computeSeq() const;

    // Members
    std::string id_;
    int d1max_; // include the last d1max bases of the first node on walk
    int d2max_; // include the first d2max bases of the last node on walk
};


class ClosureDB
{

    public:
    typedef std::vector<Closure> ClosureVec;

    ClosureDB() {};

    // Store the closure for the ClosePathResult if it is unique
    void process(const ClosePathResult& res);

    // Filter closures which are contained
    void filterContainments();

    // Apply closures to the graph
    void addClosurePaths(StringGraph* pGraph);
    void addClosurePath(StringGraph* pGraph, const Closure& c, VertexPtrVec* pInteriorVertices);

    private:
    ClosureVec nonContained_;
    ClosureVec contained_;
};

#endif
