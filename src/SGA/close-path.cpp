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
#include "ProcessFramework.h"
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
"      --noRemoveEdges                  Do not remove low coverage edges after path search.\n"
"      -m, minOverlap                   minimum overlap used when loading the ASQG file. (Default: 0 - use all overlaps)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

//static const char* PROGRAM_IDENT = PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    // Standard options
    static unsigned int verbose = 0;
    static int numThreads = 1;
    static std::string outputPfx;

    // Interval options
    static float maxNumStd = 3.0;
    static float minStd = 0.0;
    static int intervalWidth = 50;
    static int maxOL = 150;
    static int maxGap = 200;

    // Miscellaneous
    static int numRounds = 1;
    static bool writeSubgraph = false;
    static int minOverlap = 0; // Minimum overlap to use when loading graph
    static bool findOverlaps = false;
    static bool removeEdges = true;
    
    // Required
    static std::string graphFile;
    static std::string bundleFile;
}

static const char* shortopts = "vm:s:t:p:o:";

enum { OPT_HELP = 1, OPT_VERSION, OPT_MINSTD, OPT_REMOVE_EDGES, OPT_NUMROUNDS, OPT_WRITESUBGRAPH, OPT_NOREMOVEEDGES, OPT_MAXOL, OPT_MAXGAP, OPT_INTERVALWIDTH};

static const struct option longopts[] = {
    { "verbose",       no_argument,       NULL, 'v' },
    { "threads",       required_argument, NULL, 't'},
    { "output",       required_argument, NULL, 'o' },
    { "minOverlap",       required_argument, NULL, 'm' },
    { "maxNumStd",        required_argument, NULL, 's' },
    { "minStd",        required_argument, NULL, OPT_MINSTD},
    { "numRounds",    required_argument, NULL, OPT_NUMROUNDS}, 
    { "writeSubgraph", no_argument, NULL, OPT_WRITESUBGRAPH},
    { "noRemoveEdges", no_argument, NULL, OPT_NOREMOVEEDGES},
    { "maxOL", required_argument, NULL, OPT_MAXOL},
    { "maxGap", required_argument, NULL, OPT_MAXGAP},
    { "intervalWidth", required_argument, NULL, OPT_INTERVALWIDTH},
    { "help",          no_argument,       NULL, OPT_HELP },
    { "version",       no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

void printOptions();

// Main with parallel processing
int closePathMain(int argc, char** argv)
{
    using namespace std;


    typedef ProcessFramework<ClosePathWorkItem, 
                             ClosePathResult,
                             ClosePathWorkItemGenerator<ClosePathWorkItem>,
                             ClosePathProcess,
                             ClosePathPostProcess>  ClosePathProcessFramework;

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

    const size_t bufferSize = 500;
    const size_t reportInterval = 2000;

    //ClosePathProcessFramework processFramework("close-path", bufferSize, reportInterval);
    ClosePathThreadScheduler  threadScheduler("close-path", bufferSize, reportInterval);

    std::vector<EdgeCovCriteria> covCriteria;
    for (int i=0; i < opt::numRounds; i++)
        covCriteria.push_back(EdgeCovCriteria(0,0,0,1.0)); // Require 1 read pair to cover each edge

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
        ClosePathProcessFramework processFramework(ssProcessName.str(), bufferSize, reportInterval);
        BundleReader* bundleReader = new BundleReader(opt::bundleFile);
        WorkGenerator* workGenerator = new WorkGenerator(bundleReader);

        std::ostringstream ssOutputPfx;
        ssOutputPfx << opt::outputPfx << ".round" << roundNum;
        std::string roundOutputPfx = ssOutputPfx.str();
        ClosePathPostProcess* postProcessor = new ClosePathPostProcess(pGraph, roundOutputPfx, opt::maxNumStd, opt::maxGap, opt::writeSubgraph);

        // Find path closures for all bundles
        if (opt::numThreads <= 1)
        {
            // Serial Mode
            ClosePathProcess processor(pGraph, opt::maxNumStd, opt::maxGap, opt::maxOL, opt::intervalWidth, opt::findOverlaps);
            threadScheduler.processWorkSerial(*workGenerator, &processor, postProcessor);
        }
        else
        {
            // Parallel Mode
            std::vector<ClosePathProcess*> processorVector;
            for(int i = 0; i < opt::numThreads; ++i)
            {
                ClosePathProcess* pProcessor = new ClosePathProcess(pGraph, opt::maxNumStd, opt::maxGap, opt::maxOL, opt::intervalWidth, opt::findOverlaps);
                processorVector.push_back(pProcessor);
            }

            //processFramework.processWorkParallel(*workGenerator, processorVector, postProcessor);
            threadScheduler.processWorkParallel(*workGenerator, processorVector, postProcessor);

            for(int i = 0; i < opt::numThreads; ++i)
            {
                delete processorVector[i];
            }
        }

        // Remove untrusted edges from the graph, and add missing edges
        if (opt::removeEdges)
        {
            postProcessor->removeEdges(covCriteria[i]);
            size_t edgesAdded = 0;
            if (opt::findOverlaps)
            {
               edgesAdded = postProcessor->addEdgesToGraph();
            }
            pGraph->writeASQG(roundOutputPfx + "-pruned.asqg.gz");
            std::cout << "Added " << edgesAdded << " edges to the graph.\n";
            std::cout << "Graph stats after round " << roundNum << " pruning:\n";
            pGraph->stats();
        }

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
                  << "MinOverlap: " << opt::minOverlap << "\n"
                  << "MaxNumStd: " << opt::maxNumStd << "\n"
                  << "MinStd: " << opt::minStd << "\n"
                  << "maxOL: " << opt::maxOL << "\n"
                  << "maxGap: " << opt::maxGap << "\n"
                  << "intervalWidth: " << opt::intervalWidth << "\n"
                  << "graphFile: " << opt::graphFile << "\n"
                  << "bundleFile: " << opt::bundleFile << "\n"
                  << "outputPfx: " << opt::outputPfx << "\n"
                  << std::endl;
    }
}