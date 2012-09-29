#include <iostream>
#include <string>


#include "SGUtil.h"
#include "Bigraph.h"

int main()
{
    using namespace std;
    std::cout << "Hello, World!" << std::endl;
    unsigned int minOverlap = 40;
    string asqgFile = "test.asqg";
    StringGraph * pGraph = SGUtil::loadASQG(asqgFile, minOverlap);
    pGraph->stats();
    delete pGraph;
    return 0;
}
