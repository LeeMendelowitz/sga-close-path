#include "SGSearch.h"
#include "SGUtil.h"
#include "Bigraph.h"
#include "SimpleAllocator.h"

#include <iostream>

int main()
{
    SGSearchParams  params = SGSearchParams(0, 0, ED_SENSE, 100);
    std::cout << "size of params: " << sizeof(params) << std::endl;

    SGSearchTree myTree(params);

    std::cout << "size of tree: " << sizeof(myTree) << std::endl;

    GraphSearchNode<Vertex, Edge, SGDistanceFunction> myNode(0, ED_SENSE, 0, 0, 10);
    std::cout << "size of graph search node: " << sizeof(myNode) << std::endl;

    
    return 1;
}
