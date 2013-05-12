//-----------------------------------------------
// Copyright 2011 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// closepath
//


#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>

#include "close-path.h"
#include "closePathProcess.h"
#include "Util.h"
//#include "ProcessFramework.h"
#include "ProcessFramework2.h"
#include "Timer.h"
#include "SGACommon.h"
#include "bundle.h"
#include "edgeTracker.h"
//#include "bundleManager.h"


//
// Getopt
//
#define SUBPROGRAM "close-path"
static const char *CLOSEPATH_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Lee Mendelowitz (lmendelo@umiacs.umd.edu)\n";

static const char *CLOSEPATH_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] BUNDLES.bundles GRAPH.asqg.gz\n"
" Find paths between links contigs using the overlap graph\n"
"\n"
"      --help                           display this help and exit\n"
"      -v, --verbose                    display verbose output\n"
"      -t, --threads=NUM                use NUM threads to find path closures\n"
"      -o, --output=NAME                use output prefix NAME. Defaults to bundle filename prefix.\n"
"      -n, --minLinks=N               only use bundles with at least N links.\n"
"\n\nInterval Options:\n"
"      -s, maxNumStd=FLOAT              maximum number of standard deviations allowed in path length deviation. (Default 3.0)\n"
"      --minStd=FLOAT                   minimum standard deviation to use for a bundle. All bundle standard deviations will\n"
"                                       be set to MAX(bundle.std, FLOAT)\n"
"      --intervalWidth=INT              If no paths are found in interval +/- N*std, try the fixed interval +/- F.\n"
"      --maxOL=INT                      Maximum Overlap allowed when searching for paths. (Default: 150)\n"
"      --maxGap=INT                     Upper bound on gap interval. (Default: 200)\n"
"\n\nMiscellaneous Options:\n"
"      --numRounds=NUM                  Perform NUM rounds of edge pruning.\n"
"      --writeSubgraph                  Write out the subgraph for any repetitive region.\n"
"      --minEdgeCov=FLOAT               Minimum required read coverage to keep an edge. (Default: 1.0)\n"
"      --noRemoveEdges                  Do not remove low coverage edges after path search.\n"
"      --useDFS                         Use a bounded DFS in cases where one sided BFS yields too many paths.\n"
"      -m, minOverlap                   minimum overlap used when loading the ASQG file. (Default: 0 - use all overlaps)\n"

"\n\nContig Merging Options:\n"
"      --mergeContigs                   Merge single copy contigs if there is a unique bundle between them.\n"
"      --astat=FILE                     Path to astat file.\n"
"      --minLinksMerge=INT              Minimum links to use bundles for single copy contig merging.\n"
"      --removeInteriorNodes            Removal all nodes in interior of merge contig path.\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

//static const char* PROGRAM_IDENT = PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    // Standard options
    static unsigned int verbose = 0;
    static int numThreads = 1;
    static std::string outputPfx;
    static int minLinks = -1;

    // Interval options
    static float maxNumStd = 3.0;
    static float minStd = 0.0;
    static float minEdgeCov = 1.0;
    static int intervalWidth = 50;
    static int maxOL = 150;
    static int maxGap = 200;

    // Effort
    static bool useDFS = false;

    // Miscellaneous
    static int numRounds = 1;
    static bool writeSubgraph = false;
    static int minOverlap = 0; // Minimum overlap to use when loading graph
    static bool findOverlaps = false;
    static bool removeEdges = true;

    // Merging Options
    static bool mergeContigs = false;
    static std::string astatFile;
    static float astatThreshold = 20.0;
    static int singleCopyLenThreshold = 200;
    static int minLinksMerge = 2;
    static bool removeInteriorNodes = false;
    
    // Required
    static std::string graphFile;
    static std::string bundleFile;
}

static const char* shortopts = "vm:s:t:p:o:";

enum { OPT_HELP = 1, OPT_VERSION, OPT_MINSTD, OPT_REMOVE_EDGES, OPT_NUMROUNDS, OPT_WRITESUBGRAPH, OPT_NOREMOVEEDGES, OPT_MAXOL, OPT_MAXGAP, OPT_INTERVALWIDTH,
       OPT_USEDFS, OPT_MINEDGECOV, OPT_MERGECONTIGS, OPT_ASTAT, OPT_MINLINKSMERGE, OPT_REMOVEINTERIORNODES};

static const struct option longopts[] = {
    { "verbose",       no_argument,       NULL, 'v' },
    { "threads",       required_argument, NULL, 't'},
    { "output",       required_argument, NULL, 'o' },
    { "minLinks",       required_argument, NULL, 'n' },
    { "minOverlap",       required_argument, NULL, 'm' },
    { "maxNumStd",        required_argument, NULL, 's' },
    { "minStd",        required_argument, NULL, OPT_MINSTD},
    { "numRounds",    required_argument, NULL, OPT_NUMROUNDS}, 
    { "writeSubgraph", no_argument, NULL, OPT_WRITESUBGRAPH},
    { "minEdgeCov", required_argument, NULL, OPT_MINEDGECOV},
    { "noRemoveEdges", no_argument, NULL, OPT_NOREMOVEEDGES},
    { "useDFS", no_argument, NULL, OPT_USEDFS},
    { "maxOL", required_argument, NULL, OPT_MAXOL},
    { "maxGap", required_argument, NULL, OPT_MAXGAP},
    { "intervalWidth", required_argument, NULL, OPT_INTERVALWIDTH},
    { "mergeContigs", no_argument, NULL, OPT_MERGECONTIGS},
    { "minLinksMerge", required_argument, NULL, OPT_MINLINKSMERGE},
    { "removeInteriorNodes", no_argument, NULL, OPT_REMOVEINTERIORNODES},
    { "astat", required_argument, NULL, OPT_ASTAT},
    { "help",          no_argument,       NULL, OPT_HELP },
    { "version",       no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

void printOptions();

// Main with parallel processing
int closePathMain(int argc, char** argv)
{
    using namespace std;


    /*typedef ProcessFramework<ClosePathWorkItem, 
                             ClosePathResult,
                             ClosePathWorkItemGenerator<ClosePathWorkItem>,
                             ClosePathProcess,
                             ClosePathPostProcess>  ClosePathProcessFramework;*/

    typedef ThreadScheduler<ClosePathWorkItem, 
                             ClosePathResult,
                             ClosePathWorkItemGenerator<ClosePathWorkItem>,
                             ClosePathProcess,
                             ClosePathPostProcess>  ClosePathThreadScheduler;

    typedef ClosePathWorkItemGenerator<ClosePathWorkItem> WorkGenerator;

    parseClosePathOptions(argc, argv);
    printOptions();

    // Read the graph file
    std::cout << "Reading Graph: " << opt::graphFile << std::endl;
    StringGraph * pGraph = SGUtil::loadASQG(opt::graphFile, opt::minOverlap);
    pGraph->stats();
    std::cout << "\n\n\n";

    const size_t bufferSize = 100;
    const size_t reportInterval = 2000;

    //ClosePathProcessFramework processFramework("close-path", bufferSize, reportInterval);
    ClosePathThreadScheduler  threadScheduler("close-path", bufferSize, reportInterval);

    // Edge coverage requirement
    std::vector<EdgeCovCriteria> covCriteria;
    for (int i=0; i < opt::numRounds; i++)
        covCriteria.push_back(EdgeCovCriteria(0,0,0,opt::minEdgeCov));

    //covCriteria.push_back(EdgeCovCriteria(0,1,0,0.0));  // Require 1 unique bundle closure
    //covCriteria.push_back(EdgeCovCriteria(0,1,0,0.0)); 

    //covCriteria.push_back(EdgeCovCriteria(0,0,1,0.0));  // Require the edge to be covered by 1 unique read pair closure
    //covCriteria.push_back(EdgeCovCriteria(0,0,1,0.0)); 

    size_t numRounds = covCriteria.size();

    for(size_t i = 0; i < numRounds; i ++)
    {
        size_t roundNum = i + 1;

        std::cout << "\n\n\n\n\n********************************************************\n"
                  << "sga close-path round " << roundNum << std::endl;

        std::ostringstream ssProcessName;
        ssProcessName << "close-path: round " << roundNum;
        //ClosePathProcessFramework processFramework(ssProcessName.str(), bufferSize, reportInterval);
        ClosePathThreadScheduler processFramework(ssProcessName.str(), bufferSize, reportInterval);
        BundleReader* bundleReader = new BundleReader(opt::bundleFile);
        WorkGenerator* workGenerator = new WorkGenerator(bundleReader, opt::minStd, opt::minLinks);

        std::ostringstream ssOutputPfx;
        ssOutputPfx << opt::outputPfx << ".round" << roundNum;
        std::string roundOutputPfx = ssOutputPfx.str();
        ClosePathPostProcess* postProcessor = new ClosePathPostProcess(pGraph, roundOutputPfx, opt::maxNumStd, opt::maxGap, opt::writeSubgraph,
                                                                       opt::astatFile, opt::astatThreshold, opt::singleCopyLenThreshold, opt::minLinksMerge,
                                                                       opt::removeInteriorNodes);

        // Find path closures for all bundles
        if (opt::numThreads <= 1)
        {
            // Serial Mode
            ClosePathProcess processor(pGraph, opt::maxNumStd, opt::maxGap, opt::maxOL, opt::intervalWidth, opt::useDFS, opt::findOverlaps);
            threadScheduler.processWorkSerial(*workGenerator, &processor, postProcessor);
        }
        else
        {
            // Parallel Mode
            std::vector<ClosePathProcess*> processorVector;
            for(int i = 0; i < opt::numThreads; ++i)
            {
                ClosePathProcess* pProcessor = new ClosePathProcess(pGraph, opt::maxNumStd, opt::maxGap, opt::maxOL, opt::intervalWidth, opt::useDFS, opt::findOverlaps);
                processorVector.push_back(pProcessor);
            }

            //processFramework.processWorkParallel(*workGenerator, processorVector, postProcessor);
            threadScheduler.processWorkParallel(*workGenerator, processorVector, postProcessor);

            for(int i = 0; i < opt::numThreads; ++i)
            {
                delete processorVector[i];
            }
        }

        postProcessor->writeEdgeCoverage();
        postProcessor->writeDecisionClosures();

        // Remove untrusted edges from the graph, and add missing edges
        if (opt::removeEdges)
        {
            cout << "Round " << roundNum << ": Before removing edges:\n";
            pGraph->stats();
            postProcessor->removeEdges(covCriteria[i]);
            cout << "Round " << roundNum << ": After removing edges;\n";
            pGraph->stats();

// Old code to add inexact overlaps to string graph.
// The current StringGraph implementation does not allow overlaps with indels.
// Adding edges representing overlaps of different lengths between two nodes will cause
// problems downstream.
//            size_t edgesAdded = 0;
//            if (opt::findOverlaps)
//            {
//               edgesAdded = postProcessor->addEdgesToGraph();
//            }
//            std::cout << "Added " << edgesAdded << " edges to the graph.\n";
//            std::cout << "Graph stats after round " << roundNum << " pruning:\n";

            pGraph->writeASQG(roundOutputPfx + "-pruned.asqg.gz");
        } 

        if (roundNum == numRounds)
        {
            // This is the last round. Add the unique closures to the graph as nodes.
            //std::cout << "Before adding closures:\n";
            //pGraph->stats();
            //postProcessor->overlapClosures();
            //postProcessor->addClosuresToGraph();
            //std::cout << "After adding closures:\n";
            //pGraph->stats();
            
            if (opt::mergeContigs)
            {
                cout << "Before merging contigs:\n";
                pGraph->stats();
                postProcessor->mergeContigs();
                cout << "After merging contigs:\n";
                pGraph->stats();
            }
            pGraph->writeASQG(opt::outputPfx + ".final.asqg.gz");
        }

        postProcessor->printSummary(std::cout);

        delete bundleReader;
        delete workGenerator;
        delete postProcessor;
    }

    if(opt::numThreads > 1)
        pthread_exit(NULL);

    return 0;
}

// 
// Handle command line arguments
//
void parseClosePathOptions(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'm': arg >> opt::minOverlap; break;
            case 'n': arg >> opt::minLinks; break;
            case 's': arg >> opt::maxNumStd; break;
            case 'o': arg >> opt::outputPfx; break;
            case 't': arg >> opt::numThreads; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_MINSTD: arg >> opt::minStd; break;
            case OPT_WRITESUBGRAPH: opt::writeSubgraph = true; break;
            case OPT_NOREMOVEEDGES: opt::removeEdges = false; break;
            case OPT_NUMROUNDS: arg >> opt::numRounds; break;
            case OPT_MAXOL: arg >> opt::maxOL; break;
            case OPT_MAXGAP: arg >> opt::maxGap; break;
            case OPT_INTERVALWIDTH: arg >> opt::intervalWidth; break;
            case OPT_USEDFS: opt::useDFS = true; break;
            case OPT_MINEDGECOV: arg >> opt::minEdgeCov; break;
            case OPT_MERGECONTIGS: opt::mergeContigs = true; break;
            case OPT_MINLINKSMERGE: arg >> opt::minLinksMerge; break;
            case OPT_REMOVEINTERIORNODES: opt::removeInteriorNodes = true; break;
            case OPT_ASTAT: arg >> opt::astatFile; break;
            case OPT_HELP:
                std::cout << CLOSEPATH_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << CLOSEPATH_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }

        if (arg.fail())
        {
            std::cerr << SUBPROGRAM ": malformed program option\n";
            die = true;
            break;
        }
    }

    if(opt::mergeContigs && opt::astatFile.empty())
    {
        std::cerr << SUBPROGRAM ": mergeContigs flag set true but astat file missing.\n";
        die = true;
    }
    else if (!opt::mergeContigs && !opt::astatFile.empty())
    {
        std::cerr << SUBPROGRAM " warning: Ignoring astat file since mergeContigs is false.\n";
        opt::astatFile = "";
    }

    if (!opt::removeEdges && opt::numRounds > 1)
    {
        std::cerr << SUBPROGRAM " warning: edge removal is disabled, so setting numRounds to 1.\n";
        opt::numRounds = 1;
    }


    // Validate parameters
    if (argc - optind < 2) 
    {
        std::cerr << SUBPROGRAM ": missing arguments\n";
        die = true;
    } 
    else if (argc - optind > 2) 
    {
        std::cerr << SUBPROGRAM ": too many arguments\n";
        die = true;
    }    
    
    if (opt::maxNumStd < 0)
    {
        std::cerr << SUBPROGRAM ": std must be positive\n";
        die = true;
    }

    if (opt::maxOL < 0)
    {
        std::cerr << SUBPROGRAM ": maxOL must be positive\n";
        die = true;
    }

    if (opt::minEdgeCov < 0)
    {
        std::cerr << SUBPROGRAM ": minEdgeCov must be non-negative\n";
        die = true;
    }

    if (opt::intervalWidth < 0)
    {
        std::cerr << SUBPROGRAM ": intervalWidth must be positive\n";
        die = true;
    }

    if (opt::maxGap < 0)
    {
        std::cerr << SUBPROGRAM ": maxGap must be positive\n";
        die = true;
    }

    if (opt::numRounds < 0)
    {
        std::cerr << SUBPROGRAM ": numRounds must be positive\n";
        die = true;
    }

    if (die) 
    {
        std::cout << "\n" << CLOSEPATH_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }


    opt::bundleFile = argv[optind++];
    opt::graphFile = argv[optind++];

    if (opt::outputPfx.empty())
        opt::outputPfx = stripFilename(opt::bundleFile);

}

void printOptions()
{
    if (opt::verbose > 0)
    {
        std::cerr << "Verbose: " << opt::verbose << "\n"
                  << "Threads: " << opt::numThreads << "\n"
                  << "MinLinks: " << opt::minLinks << "\n"
                  << "MinOverlap: " << opt::minOverlap << "\n"
                  << "MaxNumStd: " << opt::maxNumStd << "\n"
                  << "findOverlaps: " << opt::findOverlaps << "\n"
                  << "MinStd: " << opt::minStd << "\n"
                  << "maxOL: " << opt::maxOL << "\n"
                  << "minEdgeCov: " << opt::minEdgeCov << "\n"
                  << "maxGap: " << opt::maxGap << "\n"
                  << "intervalWidth: " << opt::intervalWidth << "\n"
                  << "useDFS: " << opt::useDFS << "\n"
                  << "mergeContigs: " << opt::mergeContigs << "\n"
                  << "astatFile: " << opt::astatFile << "\n"
                  << "minLinksMerge: " << opt::minLinksMerge << "\n"
                  << "removeInteriorNodes: " << opt::removeInteriorNodes << "\n"
                  << "graphFile: " << opt::graphFile << "\n"
                  << "bundleFile: " << opt::bundleFile << "\n"
                  << "outputPfx: " << opt::outputPfx << "\n"
                  << std::endl;
    }
}
