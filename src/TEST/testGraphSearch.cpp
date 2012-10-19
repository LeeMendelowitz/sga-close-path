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
    //EdgePtrVec eVec = boundedBFS(pX, ED_SENSE, 1000, 500);
    //cout << "Found " << eVec.size() << " edges." << endl;
    //printEdgePtrVec(eVec);


    // Create subgraph from edges
    VertexID yId = "contig-900";
    Vertex * pY = pGraph->getVertex(yId);

    // Find walks in the original graph:
    cout << "*******************************************\n";
    cout << "Finding walks in original graph: " << endl;
    SGWalkVector walks1;
    size_t maxNodes = 10000;
    bool exhaustive = true;
    // One issue:  Find walks does not impose an orientation on pY!
    SGSearch::findWalks(pX, pY, ED_SENSE, ED_SENSE, -100, 500, false, false, false, maxNodes, exhaustive, walks1);
    cout << "Found " << walks1.size() << " walks:\n";
    for (size_t i=0; i < walks1.size(); i++)
        walks1[i].print();
    SGSearch::findWalks(pX, pY, ED_SENSE, ED_ANTISENSE, 0, 500, true, true, true, maxNodes, exhaustive, walks1);
    cout << "Found " << walks1.size() << " walks:\n";
    for (size_t i=0; i < walks1.size(); i++)
        walks1[i].print();


    cout << "*******************************************"
         << "\n\n\n\n\n\n\n";
    StringGraph * pSubgraph = makePathGraph(pGraph, pX, ED_SENSE, pY, ED_ANTISENSE, 0);
    Vertex * pX_new = pSubgraph->getVertex(vId);
    Vertex * pY_new = pSubgraph->getVertex(yId); 

    // Find walks in the subgraph
    cout << "*******************************************\n";
    cout << "Finding walks in subgraph: " << endl;
    SGWalkVector walks2;
    SGSearch::findWalks(pX_new, pY_new, ED_SENSE, 500, maxNodes, exhaustive, walks2);
    cout << "Found " << walks2.size() << " walks:\n";
    for (size_t i=0; i < walks2.size(); i++)
        walks2[i].print();
    SGSearch::findWalks(pX, pY, ED_SENSE, ED_ANTISENSE, 0, 500, true, true, true, maxNodes, exhaustive, walks1);
    cout << "Found " << walks1.size() << " walks:\n";
    for (size_t i=0; i < walks1.size(); i++)
        walks1[i].print();


    pSubgraph->stats();
    pSubgraph->writeASQG("test.asqg.gz");

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
    delete pSubgraph;
    return 0;
}
