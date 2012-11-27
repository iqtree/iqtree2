/*
 * utils.h
 *
 *  Created on: Nov 22, 2012
 *      Author: tung
 */

/*
I just put some declarations of the functions that I need here.
Please extend this file. It's important to have a header file.
It make things much easier for the integration with other software.
*/

#ifndef UTILS_H_
#define UTILS_H_

#include "axml.h"

void read_msa(tree *tr, const char *filename);

void makeParsimonyTree(tree *tr);

#endif /* UTILS_H_ */
