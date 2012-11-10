//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// ScaffoldRecord - A scaffold consisting of a 
// starting component and a vector of ordered
// links
//
#ifndef SCAFFOLDRECORD_H
#define SCAFFOLDRECORD_H

#include "ScaffoldLink.h"
#include "ScaffoldSequenceCollection.h"
#include "SGUtil.h"
#include "SGWalk.h"

//Placement of contig within a scaffold
struct ContigPlacement
{
    ContigPlacement (const std::string& contigId, bool isRC, int scaffStart, int scaffEnd,
                     int contigStart, int contigEnd) :
        contigId_(contigId),
        isRC_(isRC),
        scaffStart_(scaffStart),
        scaffEnd_(scaffEnd),
        contigStart_(contigStart),
        contigEnd_(contigEnd)
        {};

    std::string toString()
    {
        std::ostringstream oss;
        oss << contigId_ << ","
            << (isRC_ ? '-' : '+') << ","
            << scaffStart_ << ","
            << scaffEnd_ << ","
            << contigStart_ << ","
            << contigEnd_;
        std::string retValue = oss.str();
        return retValue;
    }

    // Reverse the start and ending scaffold coords based on the total scaffLength
    void reverseScaffCoords(int scaffLength)
    {
        int oldStart = scaffStart_;
        int oldEnd = scaffEnd_;
        scaffStart_ = scaffLength - oldEnd;
        scaffEnd_ = scaffLength - oldStart;
    }

    void reverseOrientation()
    {
        isRC_ = !isRC_;
    }

    
    std::string contigId_;
    bool isRC_;
    int scaffStart_;
    int scaffEnd_;
    int contigStart_;
    int contigEnd_;
};

// Gap resolution statistics
struct ResolveStats
{
    ResolveStats()
    {
        numGapsResolved = 0;
        numGapsAttempted = 0;
        numScaffolds = 0;
        
        graphWalkFound = 0;
        graphWalkTooMany = 0;
        graphWalkNoPath = 0;

        overlapFound = 0;
        overlapFailed = 0;
    }

    ~ResolveStats() { print(); }
    
    void print() const
    {
        printf("Num scaffolds: %d\n", numScaffolds);
        printf("Num gaps attempted: %d\n", numGapsAttempted);
        printf("Num gaps resolved: %d\n", numGapsResolved);

        printf("Num gaps resolved by graph walk: %d\n", graphWalkFound);
        printf("Num graph walks failed because of too many solutions: %d\n", graphWalkTooMany);
        printf("Num graph walks failed because of no path: %d\n", graphWalkNoPath);

        printf("Num gaps resolved by overlap: %d\n", overlapFound);
        printf("Num overlaps failed: %d\n", overlapFailed);
    }

    int numGapsResolved;
    int numGapsAttempted;
    int numScaffolds;

    // Graph resolve stats
    int graphWalkFound;
    int graphWalkTooMany;
    int graphWalkNoPath;

    // Overlap resolve stats
    int overlapFound;
    int overlapFailed;
};

// Parameter object for the scaffold record generateString function
struct ResolveParams
{
    ScaffoldSequenceCollection* pSequenceCollection;
    StringGraph* pGraph;

    int minOverlap;
    int maxOverlap;
    double maxErrorRate;
    int resolveMask;
    int minGapLength;
    double distanceFactor;
    ResolveStats* pStats;
};

// Flags indicating what level of gap resolution should be performed
const int RESOLVE_GRAPH_UNIQUE = 1;
const int RESOLVE_GRAPH_BEST = 2;
const int RESOLVE_OVERLAP = 4;

class ScaffoldRecord
{
    public:
        ScaffoldRecord();

        void setRoot(const std::string& root);
        void addLink(const ScaffoldLink& link);
        size_t getNumComponents() const;

        // Generate a sequence string representing the constructed scaffold
        std::string generateString(const ResolveParams& params, StringVector& ids) const;

        // Resolve a link by find walks through the graph
        bool graphResolve(const ResolveParams& params, const std::string& startID, 
                          const ScaffoldLink& link, std::string& extensionString,
                          SGWalk& selectedWalk) const;

        // Resolve a predicted overlap between s1 and s2 by aligning the ends of the sequences
        bool overlapResolve(const ResolveParams& params, const std::string& s1, const std::string& s2, 
                            const ScaffoldLink& link, std::string& outString) const;

        // Resolve a link by introducing a gap
        bool introduceGap(int minGapLength, const std::string& contigString, const ScaffoldLink& link, std::string& outString, int& insertedGap) const;

        void parse(const std::string& text);
        void writeScaf(std::ostream* pWriter);

    private:
        
        typedef std::vector<ScaffoldLink> LinkVector;
        
        std::string m_rootID;
        LinkVector m_links;
        ResolveStats* m_pStats;
};

#endif
