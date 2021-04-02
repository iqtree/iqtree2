/*
 *  alisim.h
 *
 *  Created on: Mar 13, 2021
 *      Author: Nhan Ly-Trong
 */

#ifndef alisim_h
#define alisim_h

#include "utils/tools.h"
#include "alignment/alisimulatorinvar.h"
#include "alignment/alisimulatorheterogeneityinvar.h"
#include "phyloanalysis.h"
#include "utils/gzstream.h"
#include <regex>
#include <string.h>

/**
*  execute Alignment Simulator (AliSim)
*/
void runAliSim(Params &params, Checkpoint *checkpoint);

/**
*  execute AliSim without inference
*/
void runAliSimWithoutInference(Params params, IQTree *&tree);

/**
*  inferring input parameters for AliSim
*/
void inferInputParameters(Params &params, Checkpoint *checkpoint, IQTree *&tree, Alignment *&aln);

/**
*  generate a random tree
*/
void generateRandomTree(Params &params);

#endif /* alisim_h */
