#ifndef ASTATREADER_H
#define ASTATREADER_H
#include <string>
#include <map>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cassert>

class AstatRecord
{
    public:
    std::string contig;
    int len;
    int normalized_len;
    int numReads;
    int numAlignedBases;
    float copyNumber;
    float astat;
};

std::ostream& operator<<(std::ostream& os, AstatRecord& rec);
std::istream& operator>>(std::istream& is, AstatRecord& rec);
void operator>>(const std::string& s, AstatRecord& rec);

typedef std::map<std::string, AstatRecord> AstatMap;
class AstatReader
{
    public: 
    AstatReader(const std::string& fname, AstatMap& map) :
        fin_(fname.c_str())
    {
        map.clear();
        read(map);
        std::cout << "Read " << map.size() << " record from astat file " << fname << ".\n";
        fin_.close();
    };

    static int readAstatFile(std::string fname, AstatMap& map)
    {
        AstatReader reader(fname, map);
        return map.size();
    }

    private:

    void read(AstatMap& map)
    {
        // Skip header line
        std::string line;
        getline(fin_, line);

        AstatRecord rec; 
        while(getline(fin_, line))
        {
            line >> rec;
            assert(map.count(rec.contig) == 0);
            map.insert(AstatMap::value_type(rec.contig, rec));
        }
    };

    std::ifstream fin_;
};


#endif
