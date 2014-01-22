/*
 * treeIO.h
 *
 *  Created on: Nov 22, 2012
 *      Author: tung
 */

/*
I just put some declarations of the functions that I need here.
Please extend this file. It's important to have a header file.
It make things much easier for the integration with other software.
*/

#ifndef TREEIO_H_
#define TREEIO_H_

#include "pll.h"

char *Tree2String(char *treestr, tree *tr, nodeptr p, pll_boolean printBranchLengths, pll_boolean printNames, pll_boolean printLikelihood,
		  pll_boolean rellTree, pll_boolean finalPrint, int perGene, pll_boolean branchLabelSupport, pll_boolean printSHSupport);

#endif /* TREEIO_H_ */
