#include <iostream>
#include <vector>
#include "SGUtil.h"
#include "pathCounter.h"

using namespace std;

typedef std::tr1::hash<PathMapKey> PathMapKeyHasher;
typedef DenseHashMap<PathMapKey, int, PathMapKeyHasher> PathMap;

int main()
{
    cout << "Testing HashMap!\n";

    Bigraph graph;
    SimpleAllocator<Vertex>* pAlloc = graph.getVertexAllocator();


    PathMap map_;
    map_.set_empty_key(PathMapKey());

    Vertex * v1 = new(pAlloc) Vertex("v1", "node1");
    Vertex * v2 = new(pAlloc) Vertex("v2", "node2");
    Vertex * v3 = new(pAlloc) Vertex("v3", "node3");

    PathMapKey k1(v1, ED_SENSE, v2, ED_SENSE, 0);
    PathMapKey k2(v1, ED_SENSE, v2, ED_SENSE, 1);
    PathMapKey k3(v2, ED_SENSE, v2, ED_SENSE, 0);
    PathMapKey k4(v3, ED_SENSE, v2, ED_SENSE, 0);
    PathMapKey k5(v1, ED_SENSE, v2, ED_SENSE, 0);
    PathMapKey k6(v1, ED_SENSE, v2, ED_SENSE, 0); 
    PathMapKey k7(v1, ED_SENSE, v2, ED_SENSE, 2); 

    PathMapKey keys[] = {k1, k2, k3, k4, k5, k6};
    vector<PathMapKey> keyVec(keys, keys + sizeof(keys)/sizeof(PathMapKey));

    for (size_t i = 0; i < keyVec.size(); i++)
    {
        map_[keyVec[i]] = i;
    }

    for (size_t i = 0; i < keyVec.size(); i++)
    {
        cout << "map_[" << i << "] " << map_[keyVec[i]] << "\n";
    }

    PathMap::iterator i = map_.find(k6);
    PathMap::iterator i2 = map_.find(k7);

    cout << "i is end? " << (i == map_.end()) << "\n";
    cout << "i2 is end? " << (i2 == map_.end()) << "\n";

    return 0;

}
