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

// Reverse an EdgePtrVec (using the edge twins in reverse order)
EdgePtrVec reverse(const EdgePtrVec& vec);

class Closure : public SGWalk
{
    public:
    Closure(const std::string& id, const SGWalk& w, int d1, int d2) :
        SGWalk(w), id_(id), d1max_(d1), d2max_(d2) { };
    Closure(const std::string& id, const EdgePtrVec& edges, int d1, int d2) :
        SGWalk(edges), id_(id), d1max_(d1), d2max_(d2) { };

    // Return True if this walk contains the other
    bool contains(const Closure& other) const;
    
    // Define less than operator for sorting
    bool operator<(const Closure& other) const
    {
        return getNumEdges() < other.getNumEdges();
    };

    // Return true if the closure includes a decision node in its interior
    bool isDecisionClosure() const;

    bool sfxOverlap(size_t startInd, const Closure& other, bool otherIsReverse);
    bool pfxOverlap(size_t endInd, const Closure& other, bool otherIsReverse);

    void colorInteriorEdges(GraphColor c) const;
    void getInteriorVertices(VertexPtrVec * pVec) const;

    std::string computeSeq() const;
    std::string computeFullSeq() const;

    // Members
    std::string id_;
    int d1max_; // include the last d1max bases of the first node on walk
    int d2max_; // include the first d2max bases of the last node on walk

    friend class ClosureDB;
};

typedef std::vector<Closure *> ClosurePtrVec;
typedef std::vector<Closure> ClosureVec;

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

    void writeDecisionClosures(std::ostream& os, bool writeContained=false) const;

    private:
    ClosurePtrVec nonContained_;
    ClosurePtrVec contained_;
    EdgeClosureMap firstEdgeMap_; // multimap from first edge in closure to the closure
    EdgeClosureMap lastEdgeMap_; // multimap from last edge in closure to the closure
    ClosureOvlVec overlaps_; // overlaps between closures
};

namespace ClosureAlgorithms
{

    // Add a list of closure to the graph.
    // For each closure:
    // Replace the closure with a single node.
    // Modify edges pointing to the first/last node of the closure to point to the new node.
    // Remove the first/last node of the closure from the graph (Justification: These are presumed to be single copy,
    //   and so their sequence and overlaps are accounted for by the node created for the closure.)
    // Remove any interior edges of the closure.
    // Remove any interior vertexes of the closure if they have become islands.
    void addClosuresToGraph(StringGraph* pGraph, const ClosureVec& closures, bool removeInteriorNodes);

    // Add a single closure to the graph
    void addClosureToGraph(StringGraph* pGraph, const Closure& c);

};



#endif
