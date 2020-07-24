/***************************************************************************
 *   Copyright (C) 2006 by BUI Quang Minh, Steffen Klaere, Arndt von Haeseler   *
 *   minh.bui@univie.ac.at   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef SPLITGRAPH_H
#define SPLITGRAPH_H

#include <list>
#include <vector>
#include <string>
#include "split.h"
#include "ncl/ncl.h"
#include "nclextra/msplitsblock.h"
#include "nclextra/mpdablock.h"
#include "nclextra/msetsblock.h"
#include "tree/node.h"
#include "splitset.h"
#include "tree/mtree.h"
#include "utils/checkpoint.h"

class MTreeSet;

using namespace std;

/**
SplitGraph class

@author BUI Quang Minh, Steffen Klaere, Arndt von Haeseler
*/
class SplitGraph : public vector<Split*>, public CheckpointFactory
{
public:

	friend class MTree;
	friend class MTreeSet;
	friend class ECOpd;

/********************************************************
	CONSTRUCTORs, INITIALIZATION AND DESTRUCTORs
********************************************************/

	/**
		empty constructor
	*/
    SplitGraph();

	/**
		construct split graph from the parameters by calling init(params).
		@param params program parameters
	*/
    SplitGraph(Params &params);

	/**
		init split graph from the parameters
		@param params program parameters
	*/
    void init(Params &params);

    /** 
        save object into the checkpoint
    */
    virtual void saveCheckpoint();

    /** 
        restore object from the checkpoint
    */
    virtual void restoreCheckpoint();
	
	/**
		if no taxa block found, but the sets block is present, then 
		this function will be invoked. It takes the taxa names from the sets block.
	*/
	void AddTaxaFromSets();

	/**
		this function is invoked 
	*/
	void createStarTree();


	/**
		new all blocks: taxa, splits, pda
	*/
	void createBlocks();

	/** free allocated memory, called by destructor */
	void freeMem();
	
	/**
		destructor
	*/
    virtual ~SplitGraph();

	/**
		convert the collection of trees in TREES block into this split graph
		@param burnin the number of beginning trees to be discarded
		@param max_count max number of trees to consider
		@param split_threshold only keep those splits which appear more than this threshold 
		@param weight_threshold minimum weight cutoff
	*/
	void convertFromTreesBlock(int burnin, int max_count, double split_threshold, 
		int split_weight_summary, double weight_threshold, const char *tree_weight_file);

/********************************************************
	PRINT INFORMATION
********************************************************/

	/**
		print infos of split graph
		@param out the output stream
	*/	
	void report(ostream &out);

	/**
		print infos of compatibility graph of splits
		@param out the output stream
	*/	
	void reportConflict(ostream &out);

/********************************************************
	GET INFORMATION
********************************************************/

	/**
		calculate sum of weights of all splits
	*/
	double calcWeight();


	/**
		calculate sum of weights of all trivial splits
	*/
	double calcTrivialWeight();

	/**
		calculate sum of weights of preserved splits in the taxa_set
		@param taxa_set a set of taxa
	*/
	double calcWeight(Split &taxa_set);

	/**
		count how many splits are covered by the taxon set
		@param taxa_set a set of taxa
	*/
	int countSplits(Split &taxa_set);

	/**
		count how many internal splits are covered by the taxon set
		@param taxa_set a set of taxa
	*/
	int countInternalSplits(Split &taxa_set);

	/**
		generate pairs of random taxa set with overlap of taxa in common
		@param filename output file name
		@param size size of the taxa set
		@param overlap number of taxa common in both sets
		@param times number of times repeated
	*/
	void generateTaxaSet(char *filename, int size, int overlap, int times);

	/**
		scale the weight of all splits to a norm factor
		@param norm normalized factor
		@param make_int TRUE to round weights to int, FALSE otherwise
		@param precision numerical precision, default (-1) for no rounding
	*/
	void scaleWeight(double norm, bool make_int = false, int precision = -1);


	/**
		@return TRUE if split sp is contained in the split system
		@param sp target split to search for
	*/
	bool containSplit(Split &sp);

	/**
		compute the boundary length of the area set, using areas_boundary variable
		@param area a set of area ID
		@return boundary length
	*/
	double computeBoundary(Split &area);

	 /**
	  @return max split weight
	 */
	double maxWeight();
	
	/**
	 * @param name a name string
	 * @return ID of leaf corresponding to name, -1 if not found
	 */
	int findLeafName(string &name);

/********************************************************
	compatibility
********************************************************/

	/**
		find the maximum-weight set of compatible splits
		@param maxsg (OUT) set of compatible splits in a split graph class
	*/
	void findMaxCompatibleSplits(SplitGraph &maxsg);

	/**
 		check the compatibility of sp against all splits in this set
		@param sp the target split
		@return TRUE if sp is compatible with all splits here, otherwise FALSE
	*/
	bool compatible(Split *sp);

/********************************************************
	OTHER STUFFS
********************************************************/

	/**
		@return number of taxa
	*/
	int getNTaxa() {
		ASSERT(size() > 0);
		return (*begin())->ntaxa;
	}

	/**
		@return number of areas
	*/
	int getNAreas() {
		return sets->getNSets();
	}

	/**
		@return number of splits
	*/
	size_t getNSplits() {
		return size();
	}

	/**
		@return number of trivial splits
	*/
	int getNTrivialSplits();

	/**
		@return taxa block
	*/
	NxsTaxaBlock *getTaxa() {
		return taxa;
	}

	void getTaxaName(vector<string> &taxname);

	/**
		@return splits block
	*/
	MSplitsBlock *getSplitsBlock() {
		return splits;
	}

	/**
		@return PDA block
	*/
	MPdaBlock *getPdaBlock() {
		return pda;
	}

	/**
		@return SETS block
	*/
	MSetsBlock *getSetsBlock() {
		return sets;
	}

	/**
		@return TREES block
	*/
	NxsTreesBlock *getTreesBlock() {
		return trees;
	}

	MTreeSet *getMTrees() {
		return mtrees;
	}

	/**
		@return TRUE if splits graph is circular
	*/
	bool isCircular() {
		return splits->cycle.size() != 0;
	}

	/**
		@return TRUE if split system is weakly compatible
	*/
	bool isWeaklyCompatible();

	/**
		@return TRUE if it is the cost-constrained PD problem
	*/
	bool isBudgetConstraint() {
		return pda->cost_constrained;
	}

	/**
		@return TRUE if the distance matrix presents for circular splits graph
		@param mat distance matrix
	*/
	bool checkCircular(mmatrix(double) &mat);

	/**
		get the ID of the taxon around the circle in a circular splits graph
		@param i a taxon
		@return index of taxon on the circle
	*/
	int getCircleId(int i) {
		ASSERT(i >= 0 && i < getNTaxa());
		return splits->cycle[i];
	}

	/**
		generate a random circular split graph
		@param params program parameters
	*/
	void generateCircular(Params &params);

	/**
		save split systems to a file in NEXUS format
		@param out output stream
		@param omit_trivial TRUE to omit trivial splits, FALSE otherwise
	*/
	void saveFileNexus(ostream &out, bool omit_trivial = false);

	/**
		save split systems to a file in star-dot format (eg **...*)
		@param out output stream
		@param omit_trivial TRUE to omit trivial splits, FALSE otherwise
	*/
	void saveFileStarDot(ostream &out, bool omit_trivial = false);

	/**
		save split systems to a file
		@param out output file name
		@param omit_trivial TRUE to omit trivial splits, FALSE otherwise
	*/
	void saveFile(const char* out_file, InputType file_format, bool omit_trivial = false);

	/**
		calculate the distance matrix, print to file in phylip format
		@param filename output file name
	*/
	void calcDistance(char *filename);


	/**
		calculate the distance matrix
		@param dist (OUT) distance matrix
	*/
	void calcDistance(mmatrix(double) &dist);

	/**
		calculate the distance matrix, based on the taxa_order
		@param dist (OUT) distance matrix
		@param taxa_order an order of taxa
	*/
	void calcDistance(mmatrix(double) &dist, vector<int> &taxa_order);

	/**
	 * remove all trivial splits
	 * @return number of trivial splits removed
	*/
	int removeTrivialSplits();
    
protected:

	/**
		taxa block
	*/
	NxsTaxaBlock *taxa;

	/**
		splits block
	*/
	MSplitsBlock *splits;

	/**
		PDA block
	*/
	MPdaBlock *pda;

	
	/**
		SETS block
	*/
	MSetsBlock *sets;

	/**
		relationship between the sets. For example, the common boundary length between two areas.
	*/
	double *areas_boundary;

	/**
		TREES block
	*/
	NxsTreesBlock *trees;

	/**
		storing set of trees if the split graph is converted from it
	*/
	MTreeSet *mtrees;

};

#endif
