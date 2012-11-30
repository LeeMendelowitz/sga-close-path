// Class for reading bundles from a file
#ifndef BUNDLEREADER_H
#define BUNDLEREADER_H

#include <string>
#include <fstream>

#include "bundle.h"

class BundleReader
{
    public:

    BundleReader(const std::string& fileName);
    ~BundleReader();
    bool getNextBundle(Bundle *& b);

    private:
        std::string fileName_;
        std::ifstream ifs_;
};

#endif 
