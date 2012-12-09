// Find overlaps between contigs based on read pairs.
#ifndef OVERLAPFINDER_H
#define OVERLAPFINDER_H

#include "bundle.h"

#include "SGUtil.h"

class OverlapFinder
{
    public:
    OverlapFinder() :
        maxErrorRate_(0.05),
        maxStd_(3.0)
        { };

    OverlapFinder(double maxErrorRate) : 
        maxErrorRate_(maxErrorRate),
        maxStd_(3.0)
        { }; 
    bool findOverlap(const Bundle* pBundle, const StringGraph* pGraph, Overlap& overlap);

    void setMaxErrorRate(double maxErrorRate) { maxErrorRate_ = maxErrorRate; }
    double getMaxErrorRate() {return maxErrorRate_; }

    private:
        double maxErrorRate_;
        float maxStd_;

};

#endif
