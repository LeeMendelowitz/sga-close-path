// Class to track which edges have been covered by unique closures
#ifndef EDGETRACKER_H
#define EDGETRACKER_H

#include <map>
#include <set>
#include <ostream>

#include "SGUtil.h"

class ClosePathResult;

class EdgeTracker
{
    typedef std::map<Edge *, int> EdgeIntMap;
    typedef std::set<Edge *> EdgeSet;

    public:
    void processResult(const ClosePathResult& res);

    // Set the color of all edges to the specified color
    void setEdgeColors(GraphColor c, int minCov=0);

    void writeCoverageStats(std::ostream& os);

    private:
    EdgeIntMap edgeCov_; // Coverage of edge by unique path closures
};

#endif
