// Store walks for unique bundle closures
// Filter out closures which are contained in other closures

#include <vector>
#include "SGWalk.h"
#include "closePathProcess.h"


class Closure : public SGWalk
{
    public:
    Closure(const SGWalk& w, int d1, int d2) :
        SGWalk(w), d1max(d1), d2max(d2) { };

    // Return True if this walk contains the other
    bool contains(const Closure& other) const;
    
    // Define less than operator for sorting
    bool operator<(const Closure& other) const
    {
        return getNumEdges() < other.getNumEdges();
    };

    // Members
    int d1max; // include the last d1max bases of the first node on walk
    int d2max; // include the first d2max bases of the last node on walk
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

    private:
    ClosureVec nonContained_;
    ClosureVec contained_;

};
