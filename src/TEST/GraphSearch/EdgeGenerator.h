// Author: Lee Mendelowitz
// Date:   9/28/2012
// Classes to generate edges for graph search


#include "SGUtil.h"


class BFSEdgeGen
{

    public:
        // Constructor
        BFSEdgeGen(const StringGraph * pGraph, EdgeDir initialDir);

        // Get the next Edge
        bool getNext(Edge * pEdge);



    private:
        const StringGraph * pGraph_;


};
