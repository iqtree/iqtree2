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

char *Tree2String(char *treestr, tree *tr, nodeptr p, boolean printBranchLengths, boolean printNames, boolean printLikelihood,
		  boolean rellTree, boolean finalPrint, int perGene, boolean branchLabelSupport, boolean printSHSupport);

#endif /* TREEIO_H_ */
