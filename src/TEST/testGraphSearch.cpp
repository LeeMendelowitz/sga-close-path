#include <iostream>
#include <string>

#include "SGUtil.h"
#include "Bigraph.h"
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
    Vertex * pVertex = pGraph->getVertex(vId);
    //EdgePtrVec eVec = boundedBFS(pVertex, ED_SENSE, 1000, 500);
    //cout << "Found " << eVec.size() << " edges." << endl;
    //printEdgePtrVec(eVec);

    // Create subgraph from edges
    VertexID yId = "contig-900";
    Vertex * pY = pGraph->getVertex(yId);

    cout << "*******************************************"
         << "\n\n\n\n\n\n\n";
    StringGraph * pSubgraph = makePathGraph(pGraph, pVertex, ED_SENSE, pY, ED_ANTISENSE, 0);
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
