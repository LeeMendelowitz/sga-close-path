//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// assemble - Assemble reads into contigs
//
#ifndef GRAPHSTATS_H
#define GRAPHSTATS_H
#include <getopt.h>
#include "config.h"

// functions
int  graphStatsMain(int argc, char** argv);
void parseGraphStatsOptions(int argc, char** argv);
void graphStats();

#endif
