#include "edgeTracker.h"
#include "closePathProcess.h"

using namespace std;

std::ostream& operator<<(std::ostream& os, const EdgeCov& cov)
{                                                                                                                                                                                        
    os << cov.bundleCov
    << "\t" << cov.uniqueBundleCov
    << "\t" << cov.uniqueReadCov
    << "\t" << cov.readCovScore;
    return os;
}

void EdgeTracker::processResult(const ClosePathResult& res)
{
    size_t numClosures = res.walks.size();
    if (numClosures == 0)
        return;

    const bool isUnique = (numClosures == 1);
    const Bundle* bundle = res.bundle;

    // Gather all of the edges from all of the walks, and their twins
    vector<Edge *> edgeVec;
    assert(numClosures == res.walks.size());
    for(size_t i = 0; i < numClosures; i++)
    {
        const SGWalk& walk = res.walks[i];
        const EdgePtrVec walkEdges = walk.getEdges();
        edgeVec.reserve(edgeVec.size() + 2*walkEdges.size());
        for(size_t j=0; j < walkEdges.size(); j++)
        {
            edgeVec.push_back(walkEdges[j]);
            edgeVec.push_back(walkEdges[j]->getTwin());
        }
    }
    set<Edge *> edgeSet(edgeVec.begin(), edgeVec.end());

    //  Add contributions to the read coverage score
    assert(bundle->n > 0);
    const float readCov = ((float) bundle->n)/numClosures;
    for(size_t i = 0; i < edgeVec.size(); i++)
    {
        Edge * pEdge = edgeVec[i];
        pair<EdgeCovMap::iterator, bool> ret;

        // Add the Edge
        ret = edgeCov_.insert(EdgeCovMap::value_type(pEdge, EdgeCov()));
        ret.first->second.readCovScore += readCov;
    }

    set<Edge *>::const_iterator i = edgeSet.begin();
    const set<Edge *>::const_iterator E = edgeSet.end();

    // Update attributes bundleCov, uniqueBundleCov, and uniqueReadCov
    for(; i != E; i++)
    {
        Edge * pEdge = *i;

        // Add the Edge
        EdgeCovMap::iterator iter = edgeCov_.find(pEdge);
        assert(iter != edgeCov_.end());
        iter->second.bundleCov++;
        if (isUnique)
        {
            iter->second.uniqueBundleCov++;
            iter->second.uniqueReadCov += bundle->n;
        }
    }
}

void EdgeTracker::writeCoverageStats(ostream& os)
{
    EdgeCovMap::const_iterator i = edgeCov_.begin();
    EdgeCovMap::const_iterator ie = edgeCov_.end();
    for(; i != ie; i++)
    {
        const Edge* pEdge = (i->first);
        const EdgeCov& cov = i->second;
        string edgeStr = pEdge->getStartID() + "." + pEdge->getEndID();
        os << edgeStr << "\t" << cov << "\n";
    }
};
