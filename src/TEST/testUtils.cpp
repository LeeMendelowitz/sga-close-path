#include "testUtils.h"

#include <iostream>

void printEdgePtrVec(const EdgePtrVec& eVec)
{
    using namespace std;
    EdgePtrVec::const_iterator i = eVec.begin();
    EdgePtrVec::const_iterator e = eVec.end();
    for(; i!= e; i++)
    {
       cout << **i << "\n";
    }
    cout << endl;
}
