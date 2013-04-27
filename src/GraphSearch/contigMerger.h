#ifndef CONTIGMERGER_H
#define CONTIGMERGER_H

#include <set>
#include <vector>
#include <map>
#include <string>

#include "SGUtil.h"
#include "SGWalk.h"
#include "ScaffoldRecord.h"

#include "bundle.h"
#include "closure.h"
#include "astatReader.h"


class ClosePathResult;
class ClosureData;

// For storing ClosureDatas in map
class ClosureDataKey
{
    public:
    ClosureDataKey(const ClosureData& c);
    ClosureDataKey(const std::string& root, const ScaffoldLink& link);
    bool operator<(const ClosureDataKey& other) const;
    void reorient();

    // data
    std::string v1_;
    std::string v2_;
    EdgeDir dir1_; // direction out of 1 towards 2
    EdgeDir dir2_; // direction out of 2 towards 1
};

std::ostream& operator<<(std::ostream& os, const ClosureDataKey& key);

class ClosureData
{
    public:
    ClosureData(const std::string& id, const Vertex* v1, const Vertex* v2, const SGWalkVector& walks, const Bundle& b) :
        id_(id), b_(b), v1_(v1), v2_(v2), walks_(walks){ };
    std::string id_; // bundle id
    Bundle b_;
    SGWalkVector walks_;
    const Vertex* v1_;
    const Vertex* v2_;
    const std::string& getVertex1() const {return b_.vertex1ID;}
    const std::string& getVertex2() const {return b_.vertex2ID;}
    ClosureDataKey makeKey() const { return ClosureDataKey(*this); };
};

std::ostream& operator<<(std::ostream& os, const ClosureData& data);

class ContigMerger
{
    public:
    typedef std::vector<ClosureData> ClosureVec;
    typedef std::map<ClosureDataKey, ClosureData> ClosureMap;

    ContigMerger(StringGraph * pGraph, const std::string& astatFile, const std::string& outputPfx, float astatThreshold, int minLength, int minLinksMerge,
                 bool removeInteriorNodes);

    void process(const ClosePathResult& res);
    void postProcess();
    void applyClosuresToGraph();
    Closure scaffoldToClosure(const ScaffoldRecord& rec, const std::string& scaffId, std::ostream& scaffOut, std::ostream& scaffClosureOut);

    private:
    StringGraph * pGraph_;
    const std::string outputPfx_;
    const float astatThreshold_;
    const int minLength_;
    const int minLinksMerge_;
    const bool removeInteriorNodes_;
    std::set<std::string> singleCopy_;
    ClosureMap closureMap_;
    ClosureVec closures_;
    AstatMap astatMap_;
};


#endif
