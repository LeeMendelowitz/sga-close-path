// Store walks for unique bundle closures
// Filter out closures which are contained in other closures
#ifndef CLOSURE_H
#define CLOSURE_H

#include <vector>
#include <string>
#include "SGWalk.h"
#include "SGUtil.h"

class ClosePathResult;
class ClosureDB;

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

    bool sfxOverlap(size_t startInd, const Closure& other, bool otherIsReverse);
    bool pfxOverlap(size_t endInd, const Closure& other, bool otherIsReverse);

    void colorInteriorEdges(GraphColor c) const;
    void getInteriorVertices(VertexPtrVec * pVec) const;

    std::string computeSeq() const;

    // Members
    std::string id_;
    int d1max_; // include the last d1max bases of the first node on walk
    int d2max_; // include the first d2max bases of the last node on walk

    friend class ClosureDB;
};


class ClosureMapEntry
{
    public:
    ClosureMapEntry(Closure * c, bool isRc) :
        c_(c), isRc_(isRc) {};
    Closure * c_;
    bool isRc_;
};

class ClosureOverlap
{
    public:
    ClosureOverlap( const Closure * c1, size_t s1, size_t e1,
                    const Closure * c2, size_t s2, size_t e2,
                    bool isRC);
    private:
    const Closure * c1_;
    const Closure * c2_;
    size_t s1_; // inclusive
    size_t e1_; // exclusive
    size_t s2_; // inclusive
    size_t e2_; // exclusive
    bool isRc_;
};            

class ClosureDB
{
    public:
    typedef std::vector<Closure *> ClosurePtrVec;
    typedef std::multimap<const Edge *, const Closure *> EdgeClosureMap;
    typedef std::vector<ClosureOverlap> ClosureOvlVec;

    ClosureDB() {};
    ~ClosureDB();

    // Store the closure for the ClosePathResult if it is unique
    void process(const ClosePathResult& res);

    // Filter closures which are contained
    void filterContainments();

    // Index nonContained closures by adding them to the firstEdgeMap_; and lastEdgeMap_;
    void indexClosures();

    void findClosureOverlaps();

    // Apply closures to the graph
    void addClosurePaths(StringGraph* pGraph);
    void addClosurePath(StringGraph* pGraph, const Closure& c, VertexPtrVec* pInteriorVertices);

    private:
    ClosurePtrVec nonContained_;
    ClosurePtrVec contained_;
    EdgeClosureMap firstEdgeMap_; // multimap from first edge in closure to the closure
    EdgeClosureMap lastEdgeMap_; // multimap from last edge in closure to the closure
    ClosureOvlVec overlaps_; // overlaps between closures
};

#endif
