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

ClosePathProcess::ClosePathProcess(StringGraph * pGraph, float numStd, int maxGap, int maxOL, int fixedIntervalWidth, bool useDFS, bool checkOverlap) :
    pGraph_(pGraph),
    numStd_(numStd) ,
    maxGap_(maxGap),
    maxOL_(maxOL), 
    fixedIntervalWidth_(fixedIntervalWidth),
    useDFS_(useDFS),
    checkOverlap_(checkOverlap)
    { 
        pSearchNodeAllocator_ = new SGSearchNodeAllocator;    
        pSearchNodeAllocator_->makePool();

        pDFSAllocator_ = new DFSAllocator;
        pDFSAllocator_->makePool();
    };

ClosePathProcess::~ClosePathProcess()
{ 
    delete pSearchNodeAllocator_;    
    delete pDFSAllocator_;
};

// Given the work item, find the paths which close the bundle
void ClosePathProcess::process(const ClosePathWorkItem& item, ClosePathResult& result)
{

    Bundle * b = item.b_;
    Vertex * pX = pGraph_->getVertex(b->vertex1ID);
    Vertex * pY = pGraph_->getVertex(b->vertex2ID);
    assert(pX);
    assert(pY);
    int lX = pX->getSeqLen();

    // Set path search parameters
    SGSearchParams params(pX, pY, b->dir1, 0);
    params.goalDir = !b->dir2;
    params.allowGoalRepeat = true;
    params.goalOriented = true;
    params.minDistanceEnforced = true;
    params.maxDistanceEnforced = true;
    params.nodeLimit = 10000;
    params.selfPrune = true;
    params.pNodeAllocator = pSearchNodeAllocator_;
    params.pDFSAllocator = pDFSAllocator_;

    // Determine the upper and lower bounds for the graph search
    // First, use +/- N*b->std to set the interval
    int minStdGap, maxStdGap; // Gap bounds determined by standard deviation estimate
    int stdDelta = (int) (numStd_ * b->std + 1);
    maxStdGap = min((int) (b->gap + stdDelta), maxGap_);
    minStdGap = max((int) (b->gap - stdDelta), -maxOL_);
    params.maxDistance = lX + maxStdGap;
    params.minDistance = max(lX + minStdGap,0);
    ClosePathResult result1(b);
    result1.overlapTooLarge =  (maxStdGap < 0) && (maxStdGap <= max(-lX, -maxOL_));
    if (!result1.overlapTooLarge)
    {
        findWalks(params, result1);
    }

    // Second, use the fixed interval size: +/- fixedIntervalWidth_ (if necessary)
    int minFixedGap, maxFixedGap; // Gap bounds determined by fixed interval
    maxFixedGap = min((int) (b->gap + fixedIntervalWidth_), maxGap_);
    minFixedGap = max((int) (b->gap - fixedIntervalWidth_), -maxOL_);
    bool fixedIntervalIsLarger = (maxFixedGap > maxStdGap);
    params.maxDistance = lX + maxFixedGap;
    params.minDistance = max(lX + minFixedGap,0);
    ClosePathResult result2(b);
    result2.overlapTooLarge =  (maxFixedGap < 0) && (maxFixedGap <= max(-lX, -maxOL_));
    bool useFixedIntervalResult = false;
    if (!result2.overlapTooLarge)
    {
        if (fixedIntervalIsLarger && result1.walks.size()==0 && !result1.tooRepetitive)
        {
            // With the initial interval, no walks valid walks existed. See if any walks
            // exist in this wider interval.
            findWalks(params, result2);
            result2.usedFixedInterval = true;
            useFixedIntervalResult = true;
        }
        else if (!fixedIntervalIsLarger && result1.tooRepetitive)
        {
            // With the intitial interval, the graph was too repetitive. See if any walks
            // exist in this smaller interval.
            findWalks(params, result2);
            result2.usedFixedInterval = true;
            useFixedIntervalResult = true;
        }
    }

    ClosePathResult& myResult = (useFixedIntervalResult ? result2 : result1);

    // Check for a sequence overlap if there is sufficient link evidence and if no path was found.
    bool checkOverlap = checkOverlap_ && (myResult.walks.size() == 0 ) && (myResult.bundle->n >= 2);
    myResult.overlap.match.numDiff = 0;
    if (checkOverlap)
    {
        Overlap overlap;
        bool foundOverlap = overlapFinder_.findOverlap(myResult.bundle, pGraph_, overlap);
        if (foundOverlap)
        {
            myResult.foundOverlap = true;
            myResult.overlap = overlap;
        }
    }

    // Set all attributes of the result to be returned
    result.copyAttributesFrom(myResult);

    pSearchNodeAllocator_->reset();
    pDFSAllocator_->reset();
};


// Find walks. First, try a single sided BFS. If the graph is too repetitive,
// try a double sided BFS to collected edges which appear on a valid walk, and then
// do a bounded DFS to find the walks.
bool ClosePathProcess::findWalks(SGSearchParams& params, ClosePathResult& result)
{
    // Do one sided BFS from pX to pY
    bool foundAll = PCSearch::findWalksOneSidedBFS(pGraph_, params, true, result.walks);
    bool noPaths = result.walks.empty();
    bool tooRepetitive = !foundAll;

    // We expect regions that are tooRepetitive to have noPaths.
    if (tooRepetitive) assert(noPaths);

    // If no walks were found because graph was too repetitive, repeat the search
    // on the same interval using double sided BFS, followed by DFS.
    //
    // Note: This is more expensive than the one sided BFS, but can handle repetitive regions better,
    // which is why we try it as a second resort. The double sided search repeats some of the work done
    // in the one sided search, but it would require some code reorganization to resolve this.
    // We expect a relatively small fraction of the bundle searches to fall into this tooRepetitive category.
    if (tooRepetitive && useDFS_)
    {
        foundAll = PCSearch::findWalksDFS(pGraph_, params, true, result.walks, result.shortestPath);
        tooRepetitive = !foundAll;
    }

    result.tooRepetitive = tooRepetitive;
    
    return !tooRepetitive;
}


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

// NOTE: THIS FUNCTION MAY NEED TO BE REVISED, IF NEEDED. MAX_OL AND BOUND_FUZZ
// ARE NOW OPTIONS TO SGA CLOSE-PATH.
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
                << "\tUsedFixedInterval"
                << "\n";
}

void ClosePathPostProcess::writeResultToStatus(const ClosePathResult & res)
{
    statusFile_ << res.bundle->id
                << "\t" << res.walks.size()
                << "\t" << res.tooRepetitive
                << "\t" << res.overlapTooLarge
                << "\t" << res.foundOverlap
                << "\t" << res.usedFixedInterval
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

// Due to changes in the bundle file format, we do not read in the 
// d1List and d2List, and so we cannot compute d1Max and d2Max
void ClosePathPostProcess::writeResultToFasta(const ClosePathResult & res)
{
    std::cerr << "NOT IMPLEMENTED!" << std::endl;
    assert(false);

    /*
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
    */
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
