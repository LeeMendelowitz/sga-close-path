//-----------------------------------------------
// Copyright 2012
// Written by Lee Mendelowitz (lmendelo@umiacs.umd.edu)
// Released under the GPL
//-----------------------------------------------
//
// close-path: Find paths to close nodes linked by paired
// reads.
//
#ifndef CLOSEPATH_H
#define CLOSEPATH_H
#include <getopt.h>
#include "config.h"

// functions

//
int closePathMain(int argc, char** argv);

// options
void parseClosePathOptions(int argc, char** argv);

#endif
