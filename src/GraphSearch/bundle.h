#ifndef BUNDLE_H
#define BUNDLE_H

#include <string>

#include "GraphCommon.h"

class Bundle
{
    public:
    Bundle() {};
    std::string vertex1ID;
    std::string vertex2ID;
    EdgeDir dir1; // direction of edge out of vertex1 towards vertex2
    EdgeDir dir2; // direction of edge out of vertex2 towards vertex1
    float gap; // esimate of gap between vertex1 & vertex2
    int n;
    float std;
};

typedef std::vector<Bundle> BundleVec;

BundleVec readBundles(const std::string& bundleFileName);


#endif
