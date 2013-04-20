// Class to track which edges have been covered by unique closures
#ifndef EDGETRACKER_H
#define EDGETRACKER_H

#include <map>
#include <set>
#include <ostream>

#include "SGUtil.h"

class ClosePathResult;

// Structure which maintains different measures of coverage for an edge
struct EdgeCov
{
    EdgeCov() :
        bundleCov(0),
        uniqueBundleCov(0),
        uniqueReadCov(0),
        readCovScore(0.0) {};

    int bundleCov; // Total number of bundles that have a closure with this edge
    int uniqueBundleCov; // Total number of bundles with a unique closure with this edge
    int uniqueReadCov; // Total number of read pairs with a unique closure with this edge
    float readCovScore; // Read coverage of edge. Read pairs with multiple closures can contribute fractions of a unit of coverage.
};


// Define operator which returns true if the edge should be retained, based on coverage criteria
class EdgeCovCriteria
{
    public:
    EdgeCovCriteria( int minBundleCov, int minUniqueBundleCov, int minUniqueReadCov, float minReadCovScore) :
        minBundleCov_(minBundleCov),
        minUniqueBundleCov_(minUniqueBundleCov),
        minUniqueReadCov_(minUniqueReadCov),
        minReadCovScore_(minReadCovScore) {};
  
    bool operator()(const EdgeCov& cov)
    {
        if ((cov.bundleCov >= minBundleCov_) &&
           (cov.uniqueBundleCov >= minUniqueBundleCov_) &&
           (cov.uniqueReadCov >= minUniqueReadCov_) &&
           (cov.readCovScore >= minReadCovScore_) )
           return true;
        return false;
    }

    private:
    int minBundleCov_;
    int minUniqueBundleCov_;
    int minUniqueReadCov_;
    float minReadCovScore_;
};

std::ostream& operator<<(std::ostream& os, const EdgeCov& cov);

class EdgeTracker
{
    typedef std::map<Edge *, int> EdgeIntMap;
    typedef std::map<Edge *, float> EdgeFloatMap;
    typedef std::set<Edge *> EdgeSet;
    typedef std::map<Edge *, EdgeCov> EdgeCovMap;

    public:
    void processResult(const ClosePathResult& res);
    void reset() {edgeCov_.clear();};

    // Set the color of edges for which the critera is true.
    // Criteria should have the operator()(EdgeCov) function defined
    template <class T>
    void setEdgeColors(GraphColor c, T criteria);
    void writeCoverageStats(std::ostream& os);

    private:
    void processResultReadCov(const ClosePathResult& res);
    void processResultBundleCov(const ClosePathResult& res);

    
    //EdgeIntMap edgeBundleCov_; // Coverage of edge by unique path closures

    // Store coverage of read pair paths. Note:
    // if a bundle of 10 read pairs has 5 closures, then the edges on a closure
    // get a coverage of 10.0/5.0 = 2.0 added (for each closure in which they appear)
    // If a bundle of 10 read pairs has 20 closures, then the edges each get 10/20.0 = 0.5 coverage.
    //EdgeFloatMap edgeReadCov_;

    EdgeCovMap edgeCov_;
};

template <class T>
void EdgeTracker::setEdgeColors(GraphColor c, T criteria)                                                                                                                               
{
    EdgeCovMap::const_iterator i = edgeCov_.begin();
    EdgeCovMap::const_iterator ie = edgeCov_.end();
    for(; i != ie; i++)
    {
        const EdgeCov& cov = i->second;
        if (criteria(cov))
            (i->first)->setColor(c);
    }
}



#endif
