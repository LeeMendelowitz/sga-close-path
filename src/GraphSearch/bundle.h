#ifndef BUNDLE_H
#define BUNDLE_H

#include <string>
#include <vector>

#include "GraphCommon.h"

class Bundle
{
    public:

    Bundle() {};
    Bundle(const std::string& line);
    std::string id;
    std::string vertex1ID;
    std::string vertex2ID;
    EdgeDir dir1; // direction of edge out of vertex1 towards vertex2
    EdgeDir dir2; // direction of edge out of vertex2 towards vertex1
    float gap; // esimate of gap between vertex1 & vertex2
    int n; // number of read pairs in this bundle
    float std;

    // The following lists measure the length of sequence in vertex
    // from a read's aligned location to the end of the vertex
    // in the direction of the pair.
    // For example:
    // |------------> X             <--------------| Y
    //        |-->    r1                 <---| r2
    //        |<--->| d1            |<------>| d2
    //std::vector<int> d1List;
    //std::vector<int> d2List;

    // Get a string representation of the bundle
    std::string toString() const;
};

typedef std::vector<Bundle> BundleVec;

BundleVec readBundles(const std::string& bundleFileName);


#endif
