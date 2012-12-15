// Code for sga close-path parallel processing

#include <iostream>
#include <algorithm>
#include <sstream>
#include <string>
#include <iostream>

#include "closePathProcess.h"
#include "PCSearch.h"
#include "SGAlgorithms.h"
#include "edgeGenerator.h"
#include "subgraph.h"

using namespace std;

const int MAX_OL = 150;
const int BOUND_FUZZ = 75;

ClosePathProcess::ClosePathProcess(StringGraph * pGraph, float numStd, int maxGap, bool checkOverlap) :
    pGraph_(pGraph),
    numStd_(numStd) ,
    maxGap_(maxGap),
    checkOverlap_(checkOverlap)
    { };

ClosePathProcess::~ClosePathProcess() { };

// Given the work item, find the paths which close the bundle
ClosePathResult ClosePathProcess::process(const ClosePathWorkItem& item)
{

    Bundle * b = item.b_;

    ClosePathResult result(b);

    Vertex * pX = pGraph_->getVertex(b->vertex1ID);
    Vertex * pY = pGraph_->getVertex(b->vertex2ID);
    assert(pX);
    assert(pY);

    // Determine the upper and lower bounds for the graph search
    // Add BOUND_FUZZ to make it very likely that the interval [minGap, maxGap]
    // traps the true gap size value.
    int maxGap = b->gap + numStd_*b->std + BOUND_FUZZ;
    int minGap = b->gap - numStd_*b->std - BOUND_FUZZ;
    if (maxGap > maxGap_) maxGap = maxGap_;
    if (minGap < (-MAX_OL)) minGap = -MAX_OL;
    int lX = pX->getSeqLen();


    #if BUNDLEMANAGER_DEBUG > 0
    int lY = pY->getSeqLen();
    cout << "*****************************\n";
    cout << " V1: " << b->vertex1ID << " Length: " << lX
         << " V2: " << b->vertex2ID << " Length: " << lY
         << " Gap: " << b->gap << " std: " << b->std
         << " maxGap: " << maxGap << " minGap: " << minGap << "\n";
    #endif

    // Skip this search if the maximum allowed gap implies too large of an overlap
    if ( (maxGap < 0) && ((-maxGap >= lX) || (maxGap <= (-MAX_OL))) )
    {
        result.overlapTooLarge = true;
        return result;
    }

    assert(minGap <= maxGap);

    // Create search params for path closure search
    SGSearchParams params(pX, pY, b->dir1, 0);
    params.goalDir = !b->dir2;
    params.maxDistance = maxGap + lX;
    params.minDistance = max(minGap + lX, 0);
    params.allowGoalRepeat = true;
    params.goalOriented = true;
    params.minDistanceEnforced = true;
    params.maxDistanceEnforced = true;
    params.nodeLimit = 10000;
    //params.selfPrune = false;
    params.selfPrune = true;

    assert(params.maxDistance > 0);
    assert(params.minDistance <= params.maxDistance);
    assert(params.minDistance >= 0);

    // Find paths and save results
    SGWalkVector walks;
    bool exhaustive = true;
    //bool foundAll = PCSearch::findWalks2(pGraph_, params, exhaustive, walks);
    bool foundAll = PCSearch::findWalksDFS(pGraph_, params, exhaustive, walks, result.shortestPath);
    //bool foundAll = PCSearch::findWalks3(pGraph_, params, exhaustive, walks);
    //bool foundAll = PCSearch::findWalks(pGraph_, params, exhaustive, walks);

    // Convert the shortest distance to the gap length
    if (walks.size() > 0)
        result.shortestPath -= pX->getSeqLen();

    result.setWalks(walks);
    result.tooRepetitive = !foundAll;
    size_t numClosures = walks.size();

    // Check for overlap if there is sufficient link evidence and if no path was found.
    bool checkOverlap = checkOverlap_ && (numClosures == 0 ) && (result.bundle->n >= 2);
    result.overlap.match.numDiff = 0;
    if (checkOverlap)
    {
        Overlap overlap;
        bool foundOverlap = overlapFinder_.findOverlap(result.bundle, pGraph_, overlap);
        if (foundOverlap)
        {
            result.foundOverlap = true;
            result.overlap = overlap;
        }
    }

    #if BUNDLEMANAGER_DEBUG > 0
    cout << "Search " << (foundAll ? "completed" : "aborted") << ". Found " << walks.size() << " walks: \n";
    #endif

    return result;
};


ClosePathPostProcess::ClosePathPostProcess(StringGraph * pGraph, const std::string& outputPfx, float numStd, int maxGap, bool writeSubgraphs) :
    pGraph_(pGraph),
    outputPfx_(outputPfx),
    numStd_(numStd),
    maxGap_(maxGap),
    writeSubgraphs_(writeSubgraphs),
    numBundlesProcessed_(0),
    numBundlesClosedUniquely_(0),
    numBundlesClosed_(0),
    numBundlesFailedOverlap_(0),
    numBundlesFailedRepetitive_(0),
    numOverlapsFound_(0),
    numReadPairsProcessed_(0),
    numReadPairsClosedUniquely_(0),
    numReadPairsClosed_(0),
    numReadPairsFailedOverlap_(0),
    numReadPairsFailedRepetitive_(0),
    numReadPairsOverlapFound_(0)
{

    // Open output files
    statusFile_.open((outputPfx_ + ".status").c_str());
    statsFile_.open((outputPfx_ + ".stats").c_str());
    //fastaFile_.open((outputPfx_ + ".fasta").c_str());
    //fastaFileUnique_.open((outputPfx_ + ".unique.fasta").c_str());
    //walksFile_.open((outputPfx_ + ".walks").c_str());
    edgeCovFile_.open((outputPfx_ + ".edgeCov").c_str());

    // Write file headers
    writeStatsHeader();
    writeStatusHeader();
}

// Write the results to file
// Delete the bundle object
void ClosePathPostProcess::process(const ClosePathWorkItem& item, const ClosePathResult& result)
{
        // Write to the output files
        writeResultToStatus(result);
        writeResultToStats(result);
    //    writeResultToFasta(result);
     //   writeResultToWalks(result);
        edgeTracker_.processResult(result);

        numBundlesProcessed_++;
        numReadPairsProcessed_ += result.bundle->n;
        size_t numClosures = result.walks.size();
        if (result.tooRepetitive)
        {
            numBundlesFailedRepetitive_++;
            numReadPairsFailedRepetitive_ += result.bundle->n;

            if(writeSubgraphs_) writeSubgraphToFile(item);
        }
        else
        {
            if (numClosures == 1)
            {
                numBundlesClosedUniquely_++;
                numReadPairsClosedUniquely_ += result.bundle->n;
            }
            if (numClosures > 0)
            {
                numBundlesClosed_++;
                numReadPairsClosed_ += result.bundle->n;
            }
        }

        if (result.overlapTooLarge)
        {
            numBundlesFailedOverlap_++;
            numReadPairsFailedOverlap_ += result.bundle->n;
        }

        if (result.foundOverlap)
        {
            numOverlapsFound_++;
            numReadPairsOverlapFound_ += result.bundle->n;

            // Create the edge corresponding to this overlap and store it
            EdgePtrVec edgeVec;
            bool createdEdges = SGAlgorithms::createEdgesFromOverlap(pGraph_, result.overlap, edgeVec);
            if (createdEdges)
            {
                for (size_t i = 0; i < edgeVec.size(); i++)
                    edgesToAdd_.push_back(edgeVec[i]);
            }
        }

        // Delete the bundle object
        delete item.b_;
}

size_t ClosePathPostProcess::addEdgesToGraph()
{

    size_t N = edgesToAdd_.size();
    for (size_t i = 0; i < N; i++)
    {
        Edge * pEdge = edgesToAdd_[i];
        pGraph_->addEdge(pEdge->getStart(), pEdge);
    }
    return N;
}

void ClosePathPostProcess::writeSubgraphToFile(const ClosePathWorkItem& item)
{

    // Gather the edges for the subgraph
    Bundle * b = item.b_;
    Vertex * pX = pGraph_->getVertex(b->vertex1ID);
    Vertex * pY = pGraph_->getVertex(b->vertex2ID);
    assert(pX);
    assert(pY);


    // Determine the upper and lower bounds for the graph search
    // Add BOUND_FUZZ to make it very likely that the interval [minGap, maxGap]
    // traps the true gap size value.
    int maxGap = b->gap + numStd_*b->std + BOUND_FUZZ;
    int minGap = b->gap - numStd_*b->std - BOUND_FUZZ;
    if (maxGap > maxGap_) maxGap = maxGap_;
    if (minGap < (-MAX_OL)) minGap = -MAX_OL;
    int lX = pX->getSeqLen();
    assert(minGap <= maxGap);
    int maxDistance = maxGap + lX;
    int dummy;
    EdgePtrVec allowedEdges = getPathEdges3(pX, b->dir1, pY, b->dir2, maxDistance, dummy);
    StringGraph * pSubgraph = Subgraph::copyGraph(pGraph_);
    Subgraph::copyEdgesToSubgraph(pSubgraph, allowedEdges);


    // Write to asqg
    string subgraphFile = b->id + ".asqg.gz";
    pSubgraph->writeASQG(subgraphFile);
    delete pSubgraph;
}


// Write edge coverage statistics to file.
// Remove uncovered edges from graph
ClosePathPostProcess::~ClosePathPostProcess()
{

    // Write edge coverage statistics to file
    edgeTracker_.writeCoverageStats(edgeCovFile_);


    // Write summary to standard out
    printSummary(std::cout);

    // Close all output files
    statusFile_.close();
    statsFile_.close();
    //fastaFile_.close();
    //fastaFileUnique_.close();
    //walksFile_.close();
    edgeCovFile_.close();
}

void ClosePathPostProcess::writeStatusHeader()
{
    statusFile_ << "bundle"
                << "\tNumClosures"
                << "\tTooRepetitive"
                << "\tOverlapTooLarge"
                << "\tFoundOverlap"
                << "\n";
}

void ClosePathPostProcess::writeResultToStatus(const ClosePathResult & res)
{
    statusFile_ << res.bundle->id
                << "\t" << res.walks.size()
                << "\t" << res.tooRepetitive
                << "\t" << res.overlapTooLarge
                << "\t" << res.foundOverlap
                << "\n";
}

void ClosePathPostProcess::writeStatsHeader()
{
    statsFile_ << "bundle"
               << "\tgapEst"
               << "\tstdEst"
               << "\tnumLinks"
               << "\tnumClosures"
               << "\tfoundAll"
               << "\tshortestPath"
               << "\toverlap"
               << "\toverlapDiff"
               << "\tclosureLengths"
               << "\n";
}

void ClosePathPostProcess::writeResultToStats(const ClosePathResult & res)
{
    // Output bundle id, gap est, std est, num links, num closures, and list of gap sizes
    const Bundle * b = res.bundle;
    statsFile_ << b->id
               << "\t" << b->gap
               << "\t" << b->std
               << "\t" << b->n
               << "\t" << res.walks.size()
               << "\t" << !res.tooRepetitive
               << "\t" << res.shortestPath
               << "\t" << res.overlap.getOverlapLength(0) << "/" << res.overlap.getOverlapLength(1)
               << "\t" << res.overlap.match.numDiff;

    // Get the gap sizes from the paths
    statsFile_ << "\t";
    for(size_t i = 0; i < res.walks.size(); i++)
    {
        statsFile_ << res.walks[i].getEndToStartDistance();
        if (i != res.walks.size()-1) statsFile_ << ",";
    }
    statsFile_ << "\n";
}

void ClosePathPostProcess::writeResultToFasta(const ClosePathResult & res)
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

void ClosePathPostProcess::writeResultToWalks(const ClosePathResult & res)
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

void ClosePathPostProcess::printSummary(std::ostream& os)
{
    int numBundlesClosedRepeat = numBundlesClosed_ - numBundlesClosedUniquely_;
    int numReadPairsClosedRepeat = numReadPairsClosed_ - numReadPairsClosedUniquely_;
    int numBundlesNoPath = numBundlesProcessed_ - numBundlesClosed_ - numBundlesFailedOverlap_ - numBundlesFailedRepetitive_;
    int numReadPairsNoPath = numReadPairsProcessed_ - numReadPairsClosed_ - numReadPairsFailedOverlap_ - numReadPairsFailedRepetitive_;
    os << "\n----------------------------------------------------------\n"
       << "sga close-path summary\n\n\n"

       << "Num. Bundles Processed:  " <<  numBundlesProcessed_ << " (" << numReadPairsProcessed_ << " read pairs)\n"
       << "Num. Closed: " << numBundlesClosed_ << " (" << numReadPairsClosed_ << " read pairs)\n"
       << "Num. Closed Uniquely: " << numBundlesClosedUniquely_ <<  " (" << numReadPairsClosedUniquely_ << " read pairs)\n"
       << "Num. Closed > 1 Path: " << numBundlesClosedRepeat << " (" << numReadPairsClosedRepeat << " read pairs)\n"
       << "Num. Failed Overlap Too Large: " << numBundlesFailedOverlap_ << " (" << numReadPairsFailedOverlap_ << " read pairs)\n"
       << "Num. Failed Graph Repetitive: " << numBundlesFailedRepetitive_ << " (" << numReadPairsFailedRepetitive_ << " read pairs)\n"
       << "Num. Failed Closure (no path): " << numBundlesNoPath << " (" << numReadPairsNoPath << " read pairs)\n"
       << "Num. Overlaps Found: " << numOverlapsFound_ << " (" << numReadPairsOverlapFound_ << " read pairs)\n"
       << "\n----------------------------------------------------------" << std::endl;
}
