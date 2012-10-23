#include <iostream>
#include <string>
#include <algorithm>

//SGA Includes
#include "SGUtil.h"
#include "Bigraph.h"
#include "SGSearch.h"

//GraphSearch Includes
#include "EdgeGenerator.h"
#include "testUtils.h"
#include "Subgraph.h"
#include "PCSearch.h"
#include "bundle.h"

using namespace std;

void writeClosures(ostream& os, StringGraph * pGraph, const Bundle& b, const SGWalkVector& walks) 
{
    // We must trim the entire sequence of the walk.
    // Use the read locations within the bundle to determine the trim points.
    size_t d1Max = (size_t) *max_element(b.d1List.begin(), b.d1List.end());
    size_t d2Max = (size_t) *max_element(b.d2List.begin(), b.d2List.end());
    Vertex * pX = pGraph->getVertex(b.vertex1ID);
    Vertex * pY = pGraph->getVertex(b.vertex2ID);
    size_t lX = pX->getSeqLen();
    size_t lY = pY->getSeqLen();
    assert(d1Max <= lX);
    assert(d2Max <= lY);
    size_t trimLeft = lX - d1Max;
    size_t trimRight = lY - d2Max;

    size_t numWalks = walks.size();
    for(size_t i = 0; i < numWalks; i++)
    {
        string seq = walks[i].getString(SGWT_START_TO_END);
        int lClosure = seq.size() - trimLeft - trimRight;
        assert(lClosure >= 0);
        string seqClosure = seq.substr(trimLeft, lClosure);
        size_t seqClosureL = seqClosure.size();

        os << ">" << b.toString() << "-" << i << " " << seqClosureL << "\n"
           << seqClosure << "\n";
    }
}


void test()
{
    using namespace std;
    unsigned int minOverlap = 40;

    string asqgFile = "assemble.K27.X2.m70-graph.asqg.gz";
    //string bundleFileName = "assemble.K27.X2.m70-contigs.alignments.bundles";
    //string bundleFileName = "assemble.K27.X2.m70-contigs.alignments.bundles.same";
    //string bundleFileName = "trouble.bundles";
    string bundleFileName = "repetative.bundles";


    std::cout << "Reading Graph: " << asqgFile << std::endl;
    StringGraph * pGraph = SGUtil::loadASQG(asqgFile, minOverlap);
    pGraph->stats();

    ofstream closureFasta("closures.fasta");

    // Read the bundle file
    cout << "Reading Bundle File: " << bundleFileName << endl;
    BundleVec bundles = readBundles(bundleFileName);
    cout << "Read " << bundles.size() << " bundles!" << endl;

    const int maxStd = 5;
    bool exhaustive = true;

    // Loop over bundles and identify paths
    for (size_t i = 0; i < bundles.size(); i++)
    {

        Bundle& b = bundles[i];


        Vertex * pX = pGraph->getVertex(b.vertex1ID);
        Vertex * pY = pGraph->getVertex(b.vertex2ID);
        assert(pX);
        assert(pY);

        int64_t maxGap = b.gap + maxStd*b.std;
        int64_t minGap = b.gap - maxStd*b.std;
        int64_t lX = pX->getSeqLen();
        int64_t lY = pY->getSeqLen();

        cout << "*****************************\n";
        cout << " V1: " << b.vertex1ID << " Length: " << lX
             << " V2: " << b.vertex2ID << " Length: " << lY
             << " Gap: " << b.gap << " std: " << b.std
             << " maxGap: " << maxGap << " minGap: " << minGap << "\n";

        // Skip this search if the maximum allowed gap implies too large of an overlap
        if ( (maxGap < 0) && (-maxGap >= lX))
        {
            cout << "Skipping due to too large of an overlap!";
            continue;
        }

        // Create search params
        SGSearchParams params(pX, pY, b.dir1, 0);
        params.goalDir = !b.dir2;
        params.maxDistance = maxGap + lX;
        params.minDistance = max(minGap + lX, (int64_t) 1);
        params.allowGoalRepeat = true;
        params.goalOriented = true;
        params.minDistanceEnforced = true;
        params.maxDistanceEnforced = true;
        params.nodeLimit = 10000;
        params.print();
        assert(params.maxDistance > 0);
        assert(params.minDistance < params.maxDistance);
        assert(params.minDistance > 0);

        SGWalkVector walks;
        bool foundAll = PCSearch::findWalks(pGraph, params, exhaustive, walks);
        cout << "Search " << (foundAll ? "completed" : "aborted") << ". Found " << walks.size() << " walks: \n";
        for (size_t i=0; i < walks.size(); i++)
        {
            walks[i].print();
            cout << "\n"
                 << " End2Start Dist: " << walks[i].getEndToStartDistance() << "\n";
        }
        
        //write closures to a fasta
        writeClosures(closureFasta, pGraph, b, walks);
    }

    closureFasta.close();
}



int main()
{

    test();
    return 0;


    using namespace std;
    unsigned int minOverlap = 40;
    string asqgFile = "test.asqg";
    std::cout << "Reading Graph: " << asqgFile << std::endl;
    StringGraph * pGraph = SGUtil::loadASQG(asqgFile, minOverlap);
    pGraph->stats();


    // Do BFS for one of the vertexes
    VertexID vId = "contig-1138";
    Vertex * pX = pGraph->getVertex(vId);
    //EdgePtrVec eVec = boundedBFS(pX, ED_SENSE, 500);
    //cout << "Found " << eVec.size() << " edges." << endl;
    //printEdgePtrVec(eVec);


    // Create subgraph from edges
    VertexID yId = "contig-900";
    Vertex * pY = pGraph->getVertex(yId);

    // Find walks in the original graph:
    //cout << "*******************************************\n";
    //cout << "Finding walks in original graph: " << endl;
    SGWalkVector walks1;
    bool exhaustive = false;

    // Default settings
    SGSearchParams params1(pX, pY, ED_SENSE, 500);

    // Impose correct goal orientation
    SGSearchParams params2(pX, pY, ED_SENSE, 500);
    params2.goalDir = ED_SENSE;
    params2.maxDistance = -39;
    params2.minDistance = -40;
    params2.allowGoalRepeat = true;
    params2.goalOriented = true;
    params2.minDistanceEnforced = true;
    params2.maxDistanceEnforced = true;
    params2.startDistance = -((int64_t) pX->getSeqLen());

    // Impose incorrect goal orientation
    SGSearchParams params3(params2);
    params2.goalDir = ED_ANTISENSE;


    // Find walks in subgraph
    cout << "*******************************************\n";
    cout << "Finding walks in subgraph: " << endl;
    SGWalkVector walks2;

    cout << "\n\nparams1\n";
    bool foundAll = PCSearch::findWalks(pGraph, params1, exhaustive, walks2);
    cout << "Search " << (foundAll ? "completed" : "aborted") << ". Found " << walks2.size() << " walks: \n";
    for (size_t i=0; i < walks2.size(); i++)
    {
        walks2[i].print();
        cout << "\n"
             << " End2Start Dist: " << walks2[i].getEndToStartDistance() << "\n"
             << " Start2End Dist: " << walks2[i].getStartToEndDistance() << "\n"
             << " End2End Dist: " << walks2[i].getEndToEndDistance() << "\n";
    }

    walks2.clear();
    cout << "\n\nparams2\n";
    foundAll = PCSearch::findWalks(pGraph, params2, exhaustive, walks2);
    cout << "Search " << (foundAll ? "completed" : "aborted") << ". Found " << walks2.size() << " walks: \n";
    for (size_t i=0; i < walks2.size(); i++)
    {
        walks2[i].print();
        cout << "\n"
             << " End2Start Dist: " << walks2[i].getEndToStartDistance() << "\n"
             << " Start2End Dist: " << walks2[i].getStartToEndDistance() << "\n"
             << " End2End Dist: " << walks2[i].getEndToEndDistance() << "\n";
    }

    walks2.clear();
    cout << "\n\nparams3\n";
    foundAll = PCSearch::findWalks(pGraph, params3, exhaustive, walks2);
    cout << "Search " << (foundAll ? "completed" : "aborted") << ". Found " << walks2.size() << " walks: \n";
    for (size_t i=0; i < walks2.size(); i++)
    {
        walks2[i].print();
        cout << "\n"
             << " End2Start Dist: " << walks2[i].getEndToStartDistance() << "\n"
             << " Start2End Dist: " << walks2[i].getStartToEndDistance() << "\n"
             << " End2End Dist: " << walks2[i].getEndToEndDistance() << "\n";
    }

    //pSubgraph->stats();
    //pSubgraph->writeASQG("test.asqg.gz");

    /*
    cout << "Creating Subgraph..." << endl;
    StringGraph * pSubGraph = Subgraph::copyGraph(pGraph);
    EdgePtrVec::iterator iEdge = eVec.begin();
    EdgePtrVec::iterator eVecE = eVec.end();
    for(; iEdge != eVecE; iEdge++)
    {
        Subgraph::addEdgeToSubgraph(pSubGraph, *iEdge);
    }
    pSubGraph->stats();
    delete pSubGraph;
    */

    delete pGraph;
    //delete pSubgraph;
    return 0;
}
