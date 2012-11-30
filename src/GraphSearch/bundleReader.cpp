//Class for reading bundles from a file
#include "bundleReader.h"

using namespace std;

BundleReader::BundleReader(const std::string& fileName) :
    fileName_(fileName),
    ifs_(fileName_.c_str())
{
    assert(ifs_);
}

BundleReader::~BundleReader()
{
    ifs_.close();
}

bool BundleReader::getNextBundle(Bundle *& b)
{
    std::string line;
    b = NULL;
    if(getline(ifs_, line))
    {
        b = new Bundle(line);
        return true;
    }
    return false;
}
