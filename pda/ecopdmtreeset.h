/*
 * EcoPDmtreeset.h
 *
 *  Created on: Nov 4, 2013
 *      Author: olga
 */

#ifndef ECOPDMTREESET_H_
#define ECOPDMTREESET_H_

#include "tree/mtreeset.h"

class EcoPDmtreeset : public MTreeSet
{
public:
	EcoPDmtreeset();

	/**
		constructor, read trees from user file
		@param userTreeFile the name of the user trees
		@param is_rooted (IN/OUT) true if tree is rooted
		@param burnin the number of beginning trees to be discarded
		@param max_count max number of trees to load
	*/
	EcoPDmtreeset(const char *userTreeFile, bool &is_rooted, int burnin, int max_count,
		const char *tree_weight_file = NULL);

	void initEcoSD(const char *userTreeFile, bool &is_rooted, int burnin, int max_count,
		const char *tree_weight_file = NULL, IntVector *weights = NULL, bool compressed = false);


	virtual ~EcoPDmtreeset();
};

#endif /* ECOPDMTREESET_H_ */
