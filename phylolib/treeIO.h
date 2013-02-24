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

#include "axml.h"

char *Tree2String(char *treestr, tree *tr, nodeptr p, pl_boolean printBranchLengths, pl_boolean printNames, pl_boolean printLikelihood,
		  pl_boolean rellTree, pl_boolean finalPrint, int perGene, pl_boolean branchLabelSupport, pl_boolean printSHSupport);

#endif /* TREEIO_H_ */
