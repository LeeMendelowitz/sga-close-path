#include <iostream>
#include <string>
#include <algorithm>

//SGA Includes
#include "SGUtil.h"
#include "Bigraph.h"
#include "SGSearch.h"

//GraphSearch Includes
#include "edgeGenerator.h"
#include "testUtils.h"
#include "subgraph.h"
#include "PCSearch.h"
#include "bundle.h"
#include "bundleManager.h"

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

    //string asqgFile = "assemble.K27.X2.m70-graph.asqg.gz";
    string asqgFile = "U00096.perfectCov.pp.filter.pass.fmmerged40.asqg.gz";
    //string bundleFileName = "assemble.K27.X2.m70-contigs.alignments.bundles";
    string bundleFileName = "U00096.perfectCov.pp.filter.pass.fmmerged40.alignments.bundles";
    //string bundleFileName = "assemble.K27.X2.m70-contigs.alignments.bundles.same";
    //string bundleFileName = "trouble.bundle";
    //string bundleFileName = "repetative.bundles";

    //string outputPfx = "trouble";
    string outputPfx = "bundle.out";


    std::cout << "Reading Graph: " << asqgFile << std::endl;
    StringGraph * pGraph = SGUtil::loadASQG(asqgFile, minOverlap);
    pGraph->stats();

    // Read the bundle file
    cout << "Reading Bundle File: " << bundleFileName << endl;

    BundleManager bundleManager(bundleFileName, pGraph, outputPfx);
    cout << "Read " << bundleManager.getNumBundles() << " bundles!" << endl;

    const int maxStd = 4;
    bool exhaustive = true;
    bundleManager.closeBundles(maxStd, exhaustive);
}

void test2()
{
    using namespace std;
    unsigned int minOverlap = 40;

    string asqgFile = "U00096.perfectCov.pp.filter.pass.fmmerged40.asqg.gz";
    string bundleFileName = "U00096.reads1.cov100.mean300.std30.alignments.bundles";
    string outputPfx = "bundle.out.mean300";


    std::cout << "Reading Graph: " << asqgFile << std::endl;
    StringGraph * pGraph = SGUtil::loadASQG(asqgFile, minOverlap);
    pGraph->stats();

    // Read the bundle file
    cout << "Reading Bundle File: " << bundleFileName << endl;

    BundleManager bundleManager(bundleFileName, pGraph, outputPfx);
    cout << "Read " << bundleManager.getNumBundles() << " bundles!" << endl;

    const int maxStd = 4;
    bool exhaustive = true;
    bundleManager.closeBundles(maxStd, exhaustive);
}



int main()
{
    test();
    //test2();
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
