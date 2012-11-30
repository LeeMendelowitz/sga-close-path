#include "edgeTracker.h"
#include "closePathProcess.h"

using namespace std;

void EdgeTracker::processResult(const ClosePathResult& res)
{
    if (res.numClosures != 1)
         return;
    const EdgePtrVec edges = res.walks[0].getEdges();
    for(size_t i = 0; i < edges.size(); i++)
    {
        Edge * pEdge = edges[i];
        Edge * pTwin = pEdge->getTwin();
        pair<EdgeIntMap::iterator, bool> ret;

        // Add the Edge
        ret = edgeCov_.insert(EdgeIntMap::value_type(pEdge, 1));
        if (!ret.second)
        {
            ret.first->second++;
        }

        // Add the Twin
        ret = edgeCov_.insert(EdgeIntMap::value_type(pTwin, 1));
        if (!ret.second)
        {
            ret.first->second++;
        }
    }
};


void EdgeTracker::setEdgeColors(GraphColor c, int minCov)
{
    EdgeIntMap::const_iterator i = edgeCov_.begin();
    EdgeIntMap::const_iterator ie = edgeCov_.end();
    for(; i != ie; i++)
    {
        if (i->second >= minCov)
            (i->first)->setColor(c);
    }
}



void EdgeTracker::writeCoverageStats(ostream& os)
{
    EdgeIntMap::const_iterator i = edgeCov_.begin();
    EdgeIntMap::const_iterator ie = edgeCov_.end();
    for(; i != ie; i++)
    {
        os << "Edge: " << *(i->first) << " Coverage: " << i->second << endl;
    }
};
