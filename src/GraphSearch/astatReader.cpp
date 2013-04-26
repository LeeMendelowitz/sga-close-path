#include "astatReader.h"

std::istream& operator>>(std::istream& is, AstatRecord& rec)
{
    is >> rec.contig >> rec.len >> rec.normalized_len
       >> rec.numReads >> rec.numAlignedBases >> rec.copyNumber
       >> rec.astat;
    //std::cout << "Read astat: " << rec << std::endl;
    assert(!is.fail());
    return is;
}

void operator>>(const std::string& s, AstatRecord& rec)
{
    std::istringstream iss(s);
    iss >> rec;
}

std::ostream& operator<<(std::ostream& os, AstatRecord& rec)
{
    os << rec.contig << ","
       << rec.len << ","
       << rec.normalized_len << ","
       << rec.numReads << ","
       << rec.numAlignedBases << ","
       << rec.copyNumber << ","
       << rec.astat;
   return os;
}
