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
#include "Timer.h"
#include "SGACommon.h"
#include "bundle.h"
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
"      -s, maxNumStd=FLOAT              maximum number of standard deviations allowed in path length deviation. (Default 3.0)\n"
"      --minStd=FLOAT                   minimum standard deviation to use for a bundle. If a bundle has a standard deviation less than FLOAT\n"
"                                       then FLOAT will be used in its place.\n"
"      --removeEdges                    Remove zero-coverage edges from the graph.\n"
"      -m, minOverlap                   minimum overlap used when loading the ASQG file. (Default: 0 - use all overlaps)\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

//static const char* PROGRAM_IDENT = PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose = 0;
    static int minOverlap = 0;
    static int numThreads = 1;
    static float maxNumStd = 3.0;
    static float minStd = 0.0;
    static std::string graphFile;
    static std::string bundleFile;
    static std::string outputPfx;
    static bool removeEdges = false;
}

static const char* shortopts = "vm:s:t:p:o:";

enum { OPT_HELP = 1, OPT_VERSION, OPT_MINSTD, OPT_REMOVE_EDGES};

static const struct option longopts[] = {
    { "verbose",       no_argument,       NULL, 'v' },
    { "removeEdges",       no_argument,       NULL, OPT_REMOVE_EDGES},
    { "threads",       required_argument, NULL, 't'},
    { "minOverlap",       required_argument, NULL, 'm' },
    { "minStd",        required_argument, NULL, OPT_MINSTD},
    { "maxNumStd",        required_argument, NULL, 's' },
    { "output",       required_argument, NULL, 'o' },
    { "help",          no_argument,       NULL, OPT_HELP },
    { "version",       no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

void printOptions();

//
// Main
//
/*
int closePathMain(int argc, char** argv)
{
    using namespace std;

    parseClosePathOptions(argc, argv);
    printOptions();

    // Read the graph file
    std::cout << "Reading Graph: " << opt::graphFile << std::endl;
    StringGraph * pGraph = SGUtil::loadASQG(opt::graphFile, opt::minOverlap);
    pGraph->stats();

    // Read the bundle file
    cout << "Reading Bundle File: " << opt::bundleFile << endl;
    BundleManager bundleManager(opt::bundleFile, pGraph, opt::outputPfx, opt::minStd);
    cout << "Read " << bundleManager.getNumBundles() << " bundles!" << endl;

    // Close the bundles
    bool exhaustive = true;
    bundleManager.closeBundles(opt::maxNumStd, exhaustive, opt::removeEdges);
    bundleManager.printSummary();

    return 0;
}
*/

// Main with parallel processing
int closePathMain(int argc, char** argv)
{
    using namespace std;

    parseClosePathOptions(argc, argv);
    printOptions();

    // Read the graph file
    std::cout << "Reading Graph: " << opt::graphFile << std::endl;
    StringGraph * pGraph = SGUtil::loadASQG(opt::graphFile, opt::minOverlap);
    pGraph->stats();

    BundleReader bundleReader(opt::bundleFile);
    ClosePathWorkItemGenerator<ClosePathWorkItem> workGenerator(&bundleReader);
    ClosePathPostProcess postProcessor(pGraph, opt::outputPfx, opt::removeEdges);

    typedef ProcessFramework<ClosePathWorkItem, 
                             ClosePathResult,
                             ClosePathWorkItemGenerator<ClosePathWorkItem>,
                             ClosePathProcess,
                             ClosePathPostProcess>  ClosePathProcessFramework;

    const size_t bufferSize = 100;
    const size_t reportInterval = 1000;
    ClosePathProcessFramework processFramework("close-path", bufferSize, reportInterval);

    // Close bundles
    if (opt::numThreads <= 1)
    {
        // Serial Mode
        ClosePathProcess processor(pGraph, opt::maxNumStd);
        processFramework.processWorkSerial(workGenerator, &processor, &postProcessor);
    }
    else
    {
        // Parallel Mode
        std::vector<ClosePathProcess*> processorVector;
        for(int i = 0; i < opt::numThreads; ++i)
        {
            ClosePathProcess* pProcessor = new ClosePathProcess(pGraph, opt::maxNumStd);
            processorVector.push_back(pProcessor);
        }

        processFramework.processWorkParallel(workGenerator, processorVector, &postProcessor);

        for(int i = 0; i < opt::numThreads; ++i)
        {
            delete processorVector[i];
        }
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
            case OPT_REMOVE_EDGES: opt::removeEdges = true; break;
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
                  << "RemoveEdges: " << opt::removeEdges << "\n"
                  << "graphFile: " << opt::graphFile << "\n"
                  << "bundleFile: " << opt::bundleFile << "\n"
                  << "outputPfx: " << opt::outputPfx << "\n"
                  << std::endl;
    }
}
