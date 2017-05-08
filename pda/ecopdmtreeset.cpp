/*
 * EcoPDmtreeset.cpp
 *
 *  Created on: Nov 4, 2013
 *      Author: olga
 */

#include "ecopdmtreeset.h"
#include "tree/mtreeset.h"

EcoPDmtreeset::EcoPDmtreeset() {
}

EcoPDmtreeset::~EcoPDmtreeset() {
}

EcoPDmtreeset::EcoPDmtreeset(const char *userTreeFile, bool &is_rooted,
	int burnin, int max_count, const char *tree_weight_file) {
	initEcoSD(userTreeFile, is_rooted, burnin, max_count, tree_weight_file);
}

void EcoPDmtreeset::initEcoSD(const char *userTreeFile, bool &is_rooted, int burnin, int max_count,
	const char *tree_weight_file, IntVector *weights, bool compressed)
{
	readTrees(userTreeFile, is_rooted, burnin, max_count, weights, compressed);
	//checkConsistency();

	if (tree_weight_file)
		readIntVector(tree_weight_file, burnin, max_count, tree_weights);
/*	else if (!weights)
		tree_weights.resize(size(), 1);*/

	if (size() != tree_weights.size())
		outError("Tree file and tree weight file have different number of entries");

}
