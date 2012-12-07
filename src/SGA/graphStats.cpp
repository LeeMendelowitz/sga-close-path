//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// assemble - convert read overlaps into contigs
//
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>

#include "Util.h"
#include "graphStats.h"
#include "SGUtil.h"
#include "SGAlgorithms.h"
#include "SGPairedAlgorithms.h"
#include "SGDebugAlgorithms.h"
#include "SGVisitors.h"
#include "Timer.h"
#include "EncodedString.h"

//
// Getopt
//
#define SUBPROGRAM "graph-stats"
static const char *GRAPHSTATS_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2009 Wellcome Trust Sanger Institute\n";

static const char *GRAPHSTATS_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... ASQGFILE\n"
"Create contigs from the assembly graph ASQGFILE.\n"
"\n"
"  -v, --verbose                        display verbose output\n"
"      --help                           display this help and exit\n"
"      -o, --out-prefix=NAME            use NAME as the prefix of the output files (output files will be NAME-contigs.fa, etc)\n"
"      -m, --min-overlap=LEN            only use overlaps of at least LEN. This can be used to filter\n"
"                                       the overlap set so that the overlap step only needs to be run once.\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

namespace opt
{
    static unsigned int verbose;
    static std::string asqgFile;
    static std::string outGraphFile;
    static std::string prefix;
    static unsigned int minOverlap;
    static bool bEdgeStats = true;
    static size_t maxEdges = 128;
}

static const char* shortopts = "p:o:m:d:g:b:a:c:r:x:l:sv";

enum { OPT_HELP = 1, OPT_VERSION, OPT_VALIDATE, OPT_EDGESTATS, OPT_EXACT, OPT_MAXINDEL, OPT_TR, OPT_MAXEDGES };

static const struct option longopts[] = {
    { "verbose",               no_argument,       NULL, 'v' },
    { "out-prefix",            required_argument, NULL, 'o' },
    { "min-overlap",           required_argument, NULL, 'm' },
    { "version",               no_argument,       NULL, OPT_VERSION },
    { "validate",              no_argument,       NULL, OPT_VALIDATE},
    { NULL, 0, NULL, 0 }
};

//
// Main
//
int graphStatsMain(int argc, char** argv)
{
    Timer* pTimer = new Timer("sga assemble");
    parseGraphStatsOptions(argc, argv);
    graphStats();
    delete pTimer;

    return 0;
}

struct VertexLengthComparison
{
    bool operator()(Vertex* pX, Vertex* pY)
    {
        return (pX->getSeqLen() < pY->getSeqLen());
    }
    /*
    bool operator()(Vertex*& pX, Vertex*& pY)
    {
        return (pX->getSeqLen() < pY->getSeqLen());
    }*/
};

void graphStats()
{
    Timer t("sga graph stats");
    StringGraph* pGraph = SGUtil::loadASQG(opt::asqgFile, opt::minOverlap, true, opt::maxEdges);
    pGraph->printMemSize();

    // Visitor functors
    SGGraphStatsVisitor statsVisit;
    SGEdgeStatsVisitor edgeStatsVisit;

    SGNodeSummaryVisitor nodeSummaryVisit(opt::prefix + ".nodeSummary.txt");

    SGContainRemoveVisitor containVisit;
    SGValidateStructureVisitor validationVisit;

    // Pre-assembly graph stats
    std::cout << "[Stats] Input graph:\n";
    pGraph->visit(statsVisit);    

    // Get the nodes in the graph

    VertexPtrVec allVertices = pGraph->getAllVertices();

    // Sort the nodes in the graph in descending order of size
    std::cout << "Summarizing Vertexes...."; std::cout.flush();
    VertexLengthComparison comp;
    std::sort(allVertices.begin(), allVertices.end(), comp);
    std::reverse(allVertices.begin(), allVertices.end());
    for(size_t i = 0; i < allVertices.size(); i++)
    {
        nodeSummaryVisit.visit(pGraph, allVertices[i]);
    }
    std::cout << "done!" << std::endl;


    //std::cout << "Validating graph structure\n";
    //pGraph->visit(validationVisit);
    //pGraph->visit(edgeStatsVisit);

    delete pGraph;
}

// 
// Handle command line arguments
//
void parseGraphStatsOptions(int argc, char** argv)
{
    // Set defaults
    opt::minOverlap = 0;
    std::string prefix = "default";
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) 
        {
            case 'o': arg >> opt::prefix; break;
            case 'm': arg >> opt::minOverlap; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_EDGESTATS: opt::bEdgeStats = true; break;
            case OPT_HELP:
                std::cout << GRAPHSTATS_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << GRAPHSTATS_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
                
        }
    }



    if (argc - optind < 1) 
    {
        std::cerr << SUBPROGRAM ": missing arguments\n";
        die = true;
    } 
    else if (argc - optind > 1) 
    {
        std::cerr << SUBPROGRAM ": too many arguments\n";
        die = true;
    }

    if (die) 
    {
        std::cout << "\n" << GRAPHSTATS_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

    // Parse the input filename
    opt::asqgFile = argv[optind++];
    if (opt::prefix.empty())
    {
        opt::prefix = stripFilename(opt::asqgFile);
    }
}
