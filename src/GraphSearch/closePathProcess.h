// Code for sga close-path parallel processing

#ifndef CLOSEPATHPROCESS
#define CLOSEPATHPROCESS

#include <fstream>
#include "bundle.h"
#include "bundleReader.h"
#include "edgeTracker.h"
#include "overlapFinder.h"
#include "SGUtil.h"
#include "SGWalk.h"
#include "SGSearch.h"
#include "dfs.h"
#include "SimpleAllocator.h"


typedef SimpleAllocator<DFS::SearchEntry> DFSAllocator;
typedef SimpleAllocator<SGSearchNode> SGSearchNodeAllocator;

// Class for an item of work for sga close-path
// This simply wraps a pointer to a Bundle
struct ClosePathWorkItem
{
    ClosePathWorkItem() : b_(NULL) {};
    ClosePathWorkItem(Bundle * b) : b_(b) { };
    Bundle * b_;
};

// Class to generate ClosePathWorkItems
template<class INPUT>
class ClosePathWorkItemGenerator
{
    public:

    ClosePathWorkItemGenerator( BundleReader * bReader, float minStd = 0.0) : 
       bReader_(bReader),
       minStd_(minStd),
       numConsumed_(0)
    { };

    bool generate(ClosePathWorkItem& out)
    {
        Bundle * b = NULL;
        bool isValid = bReader_->getNextBundle(b);
        if (isValid)
        {
            // Adjust the bundles std estimate if it is less than minStd
            if (b->std < minStd_)
                b->std = minStd_;
            out = ClosePathWorkItem(b);
            numConsumed_++;
            return true;
        }
        else
        {
            return false;
        }
    }

    inline size_t getNumConsumed() const { return numConsumed_; }

    private:
    BundleReader * bReader_;
    float minStd_;
    size_t numConsumed_;
};


class ClosePathResult
{
    public:
    ClosePathResult(const Bundle * b) :
        bundle(b),
        tooRepetitive(false),
        overlapTooLarge(false),
        foundOverlap(false),
        shortestPath(0),
        usedFixedInterval(false)
        {};

    void setWalks(const SGWalkVector& walksIn)
    {
        walks = walksIn;
    }

    const Bundle * bundle;
    bool tooRepetitive; // Closure failed due to graph too repetitive
    bool overlapTooLarge; // Closure failed because bundle implies overlap too large
    bool foundOverlap;
    int shortestPath;
    bool usedFixedInterval;
    SGWalkVector walks;
    Overlap overlap;
};


// Worker class to find path closures for a given ClosePathWorkItem
class ClosePathProcess
{

    public:
    ClosePathProcess(StringGraph * pGraph, float numStd, int maxGap, int maxOL, int fixedIntervalWidth, bool useDFS, bool checkOverlap);
    ~ClosePathProcess();
    ClosePathResult process(const ClosePathWorkItem& item);

    private:

    // Find walks. First, try a single sided BFS. If the graph is too repetitive,
    // try a double sided BFS to collected edges which appear on a valid walk, and then
    // do a bounded DFS to find the walks.
    bool findWalks(SGSearchParams& params, ClosePathResult& res); 

    StringGraph * pGraph_;
    OverlapFinder overlapFinder_;
    const float numStd_;
    const int maxGap_;
    const int maxOL_;
    const int fixedIntervalWidth_;
    bool useDFS_;
    bool checkOverlap_;

    SGSearchNodeAllocator * pSearchNodeAllocator_;
    DFSAllocator * pDFSAllocator_;

};

// Post Processor to write results to file
class ClosePathPostProcess
{
    public:
    ClosePathPostProcess(StringGraph * pGraph, const std::string& outputPfx, float numStd, int maxGap, bool writeSubgraphs);
    ~ClosePathPostProcess();
    void process(const ClosePathWorkItem& item, const ClosePathResult& result);
    void printSummary(std::ostream& os);
    size_t addEdgesToGraph();

    // Remove edges 
    template <class T>
    void removeEdges(T criteria);
    size_t getNumProcessed() const { return numBundlesProcessed_; }

    private:
        StringGraph * pGraph_;
        std::string outputPfx_;
        float numStd_;
        int maxGap_;
        bool writeSubgraphs_; // write subgraphs for repetitive items
        std::ofstream statusFile_;
        std::ofstream statsFile_;
        std::ofstream fastaFile_;
        std::ofstream fastaFileUnique_;
        std::ofstream walksFile_;
        std::ofstream edgeCovFile_;
        EdgeTracker edgeTracker_;
        EdgePtrVec edgesToAdd_;

        // Summary of closures
        int numBundlesProcessed_;
        int numBundlesClosedUniquely_;
        int numBundlesClosed_;
        int numBundlesFailedOverlap_;
        int numBundlesFailedRepetitive_;
        int numOverlapsFound_;
        int numReadPairsProcessed_;
        int numReadPairsClosedUniquely_;
        int numReadPairsClosed_;
        int numReadPairsFailedOverlap_;
        int numReadPairsFailedRepetitive_;
        int numReadPairsOverlapFound_;

        // Write results to files
        void writeResultToStatus(const ClosePathResult & res);
        void writeResultToStats(const ClosePathResult & res);
        void writeResultToFasta(const ClosePathResult & res);
        void writeResultToWalks(const ClosePathResult & res);
        void writeStatusHeader();
        void writeStatsHeader();
        void writeSubgraphToFile(const ClosePathWorkItem& item); // write the subgraph searched to file
};

template <class T>
void ClosePathPostProcess::removeEdges(T keepCriteria)
{
   pGraph_->setColors(GC_BLACK);
   int numEdges = pGraph_->getNumEdges();
   edgeTracker_.setEdgeColors(GC_WHITE, keepCriteria);
   int numRemoved = pGraph_->sweepEdges(GC_BLACK);
   float fracRemoved = ((float) numRemoved)/numEdges;
   std::cout << "Removed " << numRemoved << " low coverage edges out of "
             << numEdges << " from the graph"
             << " (" << 100.0*fracRemoved << " %)"
             << std::endl;
}


#endif
