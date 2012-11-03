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
#include "Util.h"
#include "Timer.h"
#include "SGACommon.h"
#include "bundle.h"
#include "bundleManager.h"


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
"      -o, --output=NAME                use output prefix NAME. Defaults to bundle filename prefix.\n"
"      -s, maxStd                       maximum standard deviation allowed in path length deviation.\n"
"      -m, minOverlap                   minimum overlap used when loading the ASQG file.\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

//static const char* PROGRAM_IDENT = PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static int minOverlap = 0;
    static int maxStd = 4;
    static std::string graphFile;
    static std::string bundleFile;
    static std::string outputPfx;
}

static const char* shortopts = "vm:s:p:";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "verbose",       no_argument,       NULL, 'v' },
    { "minOverlap",       required_argument, NULL, 'm' },
    { "maxStd",        required_argument, NULL, 's' },
    { "output",       required_argument, NULL, 'o' },
    { "help",          no_argument,       NULL, OPT_HELP },
    { "version",       no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int closePathMain(int argc, char** argv)
{
    using namespace std;

    parseClosePathOptions(argc, argv);

    // Read the graph file
    std::cout << "Reading Graph: " << opt::graphFile << std::endl;
    StringGraph * pGraph = SGUtil::loadASQG(opt::graphFile, opt::minOverlap);
    pGraph->stats();

    // Read the bundle file
    cout << "Reading Bundle File: " << opt::bundleFile << endl;
    BundleManager bundleManager(opt::bundleFile, pGraph, opt::outputPfx);
    cout << "Read " << bundleManager.getNumBundles() << " bundles!" << endl;

    // Close the bundles
    bool exhaustive = true;
    bundleManager.closeBundles(opt::maxStd, exhaustive);

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
            case 's': arg >> opt::maxStd; break;
            case 'o': arg >> opt::outputPfx; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
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
    
    if (opt::maxStd < 0)
    {
        std::cerr << SUBPROGRAM ": std must be positive\n";
        die = true;
    }

    if (die) 
    {
        std::cout << "\n" << CLOSEPATH_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    if (opt::outputPfx.empty())
        opt::outputPfx = stripFilename(opt::bundleFile);

    opt::bundleFile = argv[optind++];
    opt::graphFile = argv[optind++];
}
