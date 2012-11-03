#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>

#include "bundleManager.h"
#include "PCSearch.h"

using namespace std;

#define BUNDLEMANAGER_DEBUG 0

// Class to hold the results for a single bundle closure
class CloseBundleResult
{
    public:

    CloseBundleResult(const Bundle * b) :
        bundle(b),
        tooRepetative(false),
        overlapTooLarge(false),
        numClosures(0) {};

    void setWalks(const SGWalkVector& walksIn)
    {
        walks = walksIn;
        numClosures = walks.size();
    }

    const Bundle * bundle;
    bool tooRepetative; // Closure failed due to graph too repetative
    bool overlapTooLarge; // Closure failed because bundle implies overlap too large
    size_t numClosures;
    SGWalkVector walks;
};

// Constructor
BundleManager::BundleManager(const std::string& bundleFile,
              StringGraph * pGraph,
              const std::string& outputPfx) :
              bundleFile_(bundleFile),
              pGraph_(pGraph),
              outputPfx_(outputPfx)
{
    numClosedUniquely_ = 0;
    numClosed_ = 0;
    numFailedOverlap_ = 0;
    numFailedRepetative_ = 0;

    statusFile_.open((outputPfx_ + ".status").c_str());
    statsFile_.open((outputPfx_ + ".stats").c_str());
    fastaFile_.open((outputPfx_ + ".fasta").c_str());
    fastaFileUnique_.open((outputPfx_ + ".unique.fasta").c_str());
    walksFile_.open((outputPfx_ + ".walks").c_str());
    writeStatsHeader();
    writeStatusHeader();
    readBundles();
}

// Destructor
BundleManager::~BundleManager()
{
    statusFile_.close();
    statsFile_.close();
    fastaFile_.close();
    fastaFileUnique_.close();
    walksFile_.close();

    size_t numBundles = bundles_.size();
    for( size_t i =0; i < numBundles; i++)
        delete bundles_[i];
    bundles_.clear();
}

void BundleManager::readBundles()
{
    using namespace std;

    ifstream ifs(bundleFile_.c_str());
    assert(ifs);
    string line;
    while(getline(ifs, line))
        bundles_.push_back(new Bundle(line));
    ifs.close();
}

void BundleManager::closeBundles(float maxStd, bool exhaustive)
{
    size_t numBundles = bundles_.size();

    for(size_t i = 0; i < numBundles; i++)
    {
        const Bundle * b = bundles_[i];
        CloseBundleResult res(b);
        closeBundle(b, maxStd, exhaustive, res);

        // Write to the output files
        writeResultToStatus(res);
        writeResultToStats(res);
        writeResultToFasta(res);
        writeResultToWalks(res);

        if (res.tooRepetative)
            numFailedRepetative_++;
        if (res.overlapTooLarge)
            numFailedOverlap_++;
        if (res.numClosures == 1)
            numClosedUniquely_++;
        if (res.numClosures > 0)
            numClosed_++;
    }
}


// Find the closure for a single bundle.
void BundleManager::closeBundle(const Bundle * b, float maxStd, bool exhaustive, CloseBundleResult& res)
{

    res = CloseBundleResult(b);

    Vertex * pX = pGraph_->getVertex(b->vertex1ID);
    Vertex * pY = pGraph_->getVertex(b->vertex2ID);
    assert(pX);
    assert(pY);

    int64_t maxGap = b->gap + maxStd*b->std;
    int64_t minGap = b->gap - maxStd*b->std;
    int64_t lX = pX->getSeqLen();

    #if BUNDLEMANAGER_DEBUG > 0
    int64_t lY = pY->getSeqLen();
    cout << "*****************************\n";
    cout << " V1: " << b->vertex1ID << " Length: " << lX
         << " V2: " << b->vertex2ID << " Length: " << lY
         << " Gap: " << b->gap << " std: " << b->std
         << " maxGap: " << maxGap << " minGap: " << minGap << "\n";
    #endif

    // Skip this search if the maximum allowed gap implies too large of an overlap
    if ( (maxGap < 0) && (-maxGap >= lX))
    {
        res.overlapTooLarge = true;
        return;
    }

    // Create search params for path closure search
    SGSearchParams params(pX, pY, b->dir1, 0);
    params.goalDir = !b->dir2;
    params.maxDistance = maxGap + lX;
    params.minDistance = max(minGap + lX, (int64_t) 1);
    params.allowGoalRepeat = true;
    params.goalOriented = true;
    params.minDistanceEnforced = true;
    params.maxDistanceEnforced = true;
    params.nodeLimit = 10000;
    params.selfPrune = false;
    assert(params.maxDistance > 0);
    assert(params.minDistance < params.maxDistance);
    assert(params.minDistance > 0);

    // Find paths and save results
    SGWalkVector walks;
    bool foundAll = PCSearch::findWalks(pGraph_, params, exhaustive, walks);
    res.setWalks(walks);
    res.tooRepetative = !foundAll;

    #if BUNDLEMANAGER_DEBUG > 0
    cout << "Search " << (foundAll ? "completed" : "aborted") << ". Found " << walks.size() << " walks: \n";
    #endif
}


void BundleManager::writeStatusHeader()
{
    statusFile_ << "bundle"
                << "\tNumClosures"
                << "\tTooRepetative"
                << "\tOverlapTooLarge"
                << "\n";
}

void BundleManager::writeResultToStatus(CloseBundleResult & res)
{
    statusFile_ << res.bundle->id
                << "\t" << res.numClosures
                << "\t" << res.tooRepetative
                << "\t" << res.overlapTooLarge
                << "\n";
}

void BundleManager::writeStatsHeader()
{
    statsFile_ << "bundle"
               << "\tgapEst"
               << "\tstdEst"
               << "\tnumLinks"
               << "\tnumClosures"
               << "\tclosureLengths"
               << "\n";
}

void BundleManager::writeResultToStats(CloseBundleResult & res)
{
    // Output bundle id, gap est, std est, num links, num closures, and list of gap sizes
    const Bundle * b = res.bundle;
    statsFile_ << b->id
               << "\t" << b->gap
               << "\t" << b->std
               << "\t" << b->n
               << "\t" << res.numClosures;

    // Get the gap sizes from the paths
    statsFile_ << "\t";
    for(size_t i = 0; i < res.numClosures; i++)
    {
        statsFile_ << res.walks[i].getEndToStartDistance();
        if (i != res.numClosures-1) statsFile_ << ",";
    }
    statsFile_ << "\n";
}

void BundleManager::writeResultToFasta(CloseBundleResult & res)
{
    const Bundle * b = res.bundle;

    // We must trim the entire sequence of the walk.
    // Use the read locations within the bundle to determine the trim points.
    size_t d1Max = (size_t) *max_element(b->d1List.begin(), b->d1List.end());
    size_t d2Max = (size_t) *max_element(b->d2List.begin(), b->d2List.end());
    Vertex * pX = pGraph_->getVertex(b->vertex1ID);
    Vertex * pY = pGraph_->getVertex(b->vertex2ID);
    size_t lX = pX->getSeqLen();
    size_t lY = pY->getSeqLen();
    assert(d1Max <= lX);
    assert(d2Max <= lY);
    size_t trimLeft = lX - d1Max;
    size_t trimRight = lY - d2Max;

    // Write the sequence of each closure to FASTA file.
    size_t numWalks = res.walks.size();
    for(size_t i = 0; i < numWalks; i++)
    {
        string seq = res.walks[i].getString(SGWT_START_TO_END);
        int lClosure = seq.size() - trimLeft - trimRight;
        assert(lClosure >= 0);
        string seqClosure = seq.substr(trimLeft, lClosure);
        size_t seqClosureL = seqClosure.size();
        ostringstream oss;
        oss << ">" << b->id << "-" << i << " " << seqClosureL << "\n"
           << seqClosure << "\n";

        fastaFile_ << oss.str();
        if (numWalks == 1)
        {
            fastaFileUnique_ << oss.str();
        }
    }
}

void BundleManager::writeResultToWalks(CloseBundleResult & res)
{
    // Write each walk to the walks file
    const Bundle * b = res.bundle;
    size_t numWalks = res.walks.size();
    for(size_t i = 0; i < numWalks; i++)
    {
        const SGWalk& walk = res.walks[i];
        walksFile_ << b->id << "-" << i << "\t"
                   << "NumEdges: " << walk.getNumEdges() << "\t"
                   << "Gap: " << walk.getEndToStartDistance() << "\t";
        walk.printWithOL(walksFile_);
        walksFile_ << "\n";
    }
}

void BundleManager::printSummary()
{
    int numClosedRepeat = numClosed_ - numClosedUniquely_;
    std::cout << "Num. Bundles: " <<  getNumBundles() << "\n"
              << "Num. Closed: " << numClosed_ << "\n"
              << "Num. Closed Uniquely: " << numClosedUniquely_ << "\n"
              << "Num. Closed > 1 Path: " << numClosedRepeat << "\n"
              << "Num. Failed Overlap Too Large: " << numFailedOverlap_ << "\n"
              << "Num. Failed Graph Repetative: " << numFailedRepetative_
              << std::endl;

             

}
