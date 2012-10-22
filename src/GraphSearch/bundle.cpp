
#include "bundle.h"

#include <string>
#include <iostream>
#include <fstream>
#include <cassert>

using namespace std;

BundleVec readBundles(const std::string& bundleFileName)
{
    ifstream ifs(bundleFileName.c_str());
    assert(ifs);
    string line;
    BundleVec bundles;
    while(getline(ifs, line))
    {
        istringstream iss(line);
        Bundle b;
        char v1End, v2End;
        iss >> b.vertex1ID;
        iss >> v1End;
        iss >> b.vertex2ID;
        iss >> v2End;
        iss >> b.gap;
        iss >> b.std;
        iss >> b.n;
        assert(!iss.fail());
        b.dir1 = (v1End == 'E') ? ED_SENSE : ED_ANTISENSE;
        b.dir2 = (v2End == 'E') ? ED_SENSE : ED_ANTISENSE; 
        bundles.push_back(b);
    }

    return bundles;
}
