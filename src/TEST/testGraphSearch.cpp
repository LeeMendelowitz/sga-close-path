#include <iostream>
#include <string>

//SGA Includes
#include "SGUtil.h"
#include "Bigraph.h"
#include "SGSearch.h"

//GraphSearch Includes
#include "EdgeGenerator.h"
#include "testUtils.h"
#include "Subgraph.h"
#include "PCSearch.h"

int main()
{
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

    // Search with "default" settings
    /*
    cout << "\n\nparams1\n";
    bool foundAll = SGSearch::findWalks(params1, exhaustive, walks1);
    cout << "Search " << (foundAll ? "completed" : "aborted") << ". Found " << walks1.size() << " walks: \n";
    for (size_t i=0; i < walks1.size(); i++)
        walks1[i].print();
   
    // Impose the correct orientation for goal
    walks1.clear();
    cout << "\n\nparams2\n";
    foundAll = SGSearch::findWalks(params2, exhaustive, walks1);
    cout << "Search " << (foundAll ? "completed" : "aborted") << ". Found " << walks1.size() << " walks: \n";
    cout << "Found " << walks1.size() << " walks:\n";
    for (size_t i=0; i < walks1.size(); i++)
        walks1[i].print();

    // Impose the incorrect orientation for goal
    walks1.clear();
    cout << "\n\nparams3\n";
    foundAll = SGSearch::findWalks(params3, exhaustive, walks1);
    cout << "Search " << (foundAll ? "completed" : "aborted") << ". Found " << walks1.size() << " walks: \n";
    cout << "Found " << walks1.size() << " walks:\n";
    for (size_t i=0; i < walks1.size(); i++)
        walks1[i].print();
    */

    // Create subgraph
    /*
    StringGraph * pSubgraph = makePathGraph(pGraph, pX, ED_SENSE, pY, ED_ANTISENSE, 0);
    Vertex * pX_new = pSubgraph->getVertex(vId);
    Vertex * pY_new = pSubgraph->getVertex(yId); 
    params1.pStartVertex = pX_new;
    params1.pEndVertex = pY_new;
    params2.pStartVertex = pX_new;
    params2.pEndVertex = pY_new;
    params3.pStartVertex = pX_new;
    params3.pEndVertex = pY_new;*/

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
