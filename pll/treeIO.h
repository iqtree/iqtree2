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

char *pllTreeToNewick(char *treestr, tree *tr, nodeptr p, pllBoolean printBranchLengths, pllBoolean printNames, pllBoolean printLikelihood,
		  pllBoolean rellTree, pllBoolean finalPrint, int perGene, pllBoolean branchLabelSupport, pllBoolean printSHSupport);
double getBranchLength(pllInstance *tr, partitionList *pr, int perGene, nodeptr p);

#endif /* TREEIO_H_ */
