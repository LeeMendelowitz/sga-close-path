// BundleManager
//  - Reads bundles from a file and owns the Bundle objects
//  - Calls PCSearch algorithm to find path closures
//  - Write bundle status to file (number of closures, search aborted? )
//  - Write bundle stats to file (estimated gap, std, number of links, list of path lengths found)
//  - Write FASTA with path closures to file

#ifndef BUNDLE_MANAGER
#define BUNDLE_MANAGER

#include <fstream>
#include <string>
#include <iostream>
#include <vector>

//GraphSearch Includes
#include "bundle.h"
#include "SGUtil.h"
#include "SGWalk.h"

#define BUNDLE_DEBUG 0

typedef std::vector<Bundle *> BundlePtrVec;
class CloseBundleResult;

class BundleManager
{

    public:

        // Constructor
        BundleManager(const std::string& bundleFile,
                      StringGraph * pGraph,
                      const std::string& outputPfx = "bundle.out",
                      float minStd = 0.0);

        // Destructor
        ~BundleManager();

        size_t getNumBundles() const { return bundles_.size(); }

        // Close all of the bundles
        // maxStd: The maximum number of standard deviations that the path is allowed
        // to stray from the estimate of the gap size.
        void closeBundles(float maxStd, bool exhaustive);

        void printSummary();
                    
    private:
        std::string bundleFile_;
        StringGraph * pGraph_;
        std::string outputPfx_;
        std::ofstream statusFile_;
        std::ofstream statsFile_;
        std::ofstream fastaFile_;
        std::ofstream fastaFileUnique_;
        std::ofstream walksFile_;
        std::ofstream edgeCovFile_;
        BundlePtrVec bundles_;
        int numClosedUniquely_;
        int numClosed_;
        int numFailedOverlap_;
        int numFailedRepetative_;

        void readBundles();

        // Close a single bundle
        void closeBundle(const Bundle * b, float maxStd, bool exhaustive, CloseBundleResult & res);

        // Write results to files
        void writeResultToStatus(CloseBundleResult & res);
        void writeResultToStats(CloseBundleResult & res);
        void writeResultToFasta(CloseBundleResult & res);
        void writeResultToWalks(CloseBundleResult & res);
        void writeStatusHeader();
        void writeStatsHeader();
};

#endif
