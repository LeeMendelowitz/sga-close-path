
#include "bundle.h"

#include <string>
#include <iostream>
#include <fstream>
#include <cassert>

using namespace std;

vector<int> readDField(const string& field)
{
    size_t L = field.size();
    assert(field[0]=='(');
    assert(field[L-1]==')');
    string inner = field.substr(1,L-2);

    vector<int> d;
    const char delim = ',';
    string item;
    istringstream ss(inner);
    while(getline(ss, item, delim))
    {
        istringstream iss(item);
        int val = 0;
        iss >> val;
        assert(!iss.fail());
        d.push_back(val);
    }
    return d;
}


Bundle::Bundle(const string& line)
{
        istringstream iss(line);

        char v1End, v2End;
        iss >> vertex1ID;
        iss >> v1End;
        iss >> vertex2ID;
        iss >> v2End;
        iss >> gap;
        iss >> std;
        iss >> n;

        // Read d fields
        string d1Field, d2Field;
        iss >> d1Field >> d2Field;
        assert(!iss.fail());
        dir1 = (v1End == 'E') ? ED_SENSE : ED_ANTISENSE;
        dir2 = (v2End == 'E') ? ED_SENSE : ED_ANTISENSE; 

        // Parse d fields
        d1List = readDField(d1Field);
        d2List = readDField(d2Field);
        assert(d1List.size() == (size_t) n);
        assert(d2List.size() == (size_t) n);
}

BundleVec readBundles(const std::string& bundleFileName)
{
    ifstream ifs(bundleFileName.c_str());
    assert(ifs);
    string line;
    BundleVec bundles;
    while(getline(ifs, line))
    {
        bundles.push_back(Bundle(line));
    }

    return bundles;
}

string Bundle::toString() const
{
    ostringstream ss;
    ss << vertex1ID << "_" << dir1 << "_"
       << vertex2ID << "_" << dir2;
    return ss.str();
}
