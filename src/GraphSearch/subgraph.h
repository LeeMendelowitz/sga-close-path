#ifndef SUBGRAPH_H
#define SUBGRAPH_H

// Lee Mendelowitz

#include "SGUtil.h"

namespace Subgraph
{
// Code to handle creating a subgraph StringGraph from an existing StringGraph
StringGraph * copyGraph(const StringGraph * pGraph);

// We do not directly add the vertex pointer to the graph, as each graph
// manages its set of vertices and we do not want to double-free the vertices
void copyVertexToSubgraph(StringGraph* pSubgraph, const Vertex* pVertex);

void addEdgeToSubgraph(StringGraph* pSubGraph, const Edge* pEdge);

// Copy edges in the eVec into the subGraph. Take care to copy an edge only once
void copyEdgesToSubgraph(StringGraph * pSubgraph, const EdgePtrVec& eVec);

};

#endif
