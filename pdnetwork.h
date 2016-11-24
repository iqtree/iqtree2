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
#ifndef PDNETWORK_H
#define PDNETWORK_H

#include "splitgraph.h"

/**
General Split Network for Phylogenetic Diversity Algorithm

@author BUI Quang Minh, Steffen Klaere, Arndt von Haeseler
*/
class PDNetwork : public SplitGraph
{
public:

	friend class MTree;
	friend class ECOpd;

	/**
		empty constructor
	*/
    PDNetwork();

	/**
		construct PD network from a NEXUS or NEWICK file, e.g. produced by SplitsTree
		@param params program parameters
	*/
    PDNetwork(Params &params);


	/**
		Identify the root node if specified, include it into the initial set
		@param root_name name of the root node
	*/
	void readRootNode(const char *root_name);

	/**
		read the parameter from the file and incoporate into split system
		@param params program parameters
	*/
	void readParams(Params &params);

	/**
		read the initial set of taxa to be included into PD-tree
		@param params program parameters
	*/
	void readInitialSet(Params &params);


	/**
		read the initial areas to be included into PD set
		@param params program parameters
	*/
	void readInitialAreas(Params &params);
	
	/**
		increase the weight of the split associated with initial set
	*/
	void proceedInitialSet();

	/**
		initialize when PD min specified
	*/
	void initPDMin();


	/**
		compute the minimum required costs to conserve a taxa set
		@param taxset set of taxa
		@return budget required
	*/
	int calcCost(IntVector &taxset);

	/**
		compute the minimum required costs to conserve a taxa set
		@param taxset set of taxa
		@return budget required
	*/
	int calcCost(Split &taxset);

	void printOutputSetScore(Params &params, vector<SplitSet> &pd_set);

	/**
		compute the PD score of a given taxa set in filename
		@param params program parameters
		@param taxa_set (OUT) corresponding set of taxa
		@param pd_more (OUT) more computed PD measures will be stored here
	*/
	void computePD(Params &params, SplitSet &taxa_set, PDRelatedMeasures &pd_more);

	/**
		this will be called by findPD at the beginning
		@param params program parameters
	*/
	virtual void enterFindPD(Params &params);

	/**
		main function to search for maximal phylogenetic diversity
		@param params program parameters
		@param taxa_set (OUT) the vector of set of taxa in the maximal PD set
		@param taxa_order (OUT) order of inserted taxa
	*/
	virtual void findPD(Params &params, vector<SplitSet> &taxa_set, vector<int> &taxa_order);


	/**
		this function will be called by findPD at the end
		@param taxa_set (IN/OUT) the vector of set of taxa in the maximal PD set
	*/
	virtual void leaveFindPD(vector<SplitSet> &taxa_set);

	/**
		calculate the PD gain matrix in terms of delta_k^j = pd(PD_k \/ {j}) - pd_k
		@param pd_set set of optimal PD sets
		@param delta (OUT) PD gain matrix
	*/
	void calcPDGain(vector<SplitSet> &pd_set, mmatrix(double) &delta);

	/**
		compute the PD score of a given taxa set with name in taxa_name, result is written to id_set.weight. 
		The difference from calcWeight() is that calcPD takes initialset into account
		@param id_set (IN/OUT) corresponding set of taxa
		
	*/
	void calcPD(Split &id_set);

	/**
		compute PD of a set of areas. It implicitly takes area_taxa map into account.
		@param area_id_set IDs of areas in the set
	*/
	void calcPDArea(Split &area_id_set);

	/**
		compute the EXCLUSIVE PD score of a given taxa set with name in taxa_name, result is written to id_set.weight
		@param id_set (IN/OUT) corresponding set of taxa IDs
	*/
	void calcExclusivePD(Split &id_set);

	/**
		compute the area's PD ENDEMISM of set of area
		@param area_set set of area
		@param pd_endem (OUT) corresponding PD endemism
	*/
	void calcPDEndemism(SplitSet &area_set, DoubleVector &pd_endem);

	/**
		compute the area's PD complementarity given a specific area
		@param area_set set of area
		@param area_names given area names as string separated by commas
		@param all_names all area names
		@param pd_comp (OUT) corresponding PD endemism
	*/
	void calcPDComplementarity(SplitSet &area_set, char *area_names, 
		vector<string> &all_names, DoubleVector &pd_comp);


	/**
		transform the problem into an Integer Linear Programming and write to .lp file
		@param params program parameters
		@param outfile name of output file in LP format
		@param total_size k for PD_k or total budget
		@param make_bin TRUE if creating binary programming
	*/
	void transformLP(Params &params, const char *outfile, int total_size, bool make_bin);
	void transformLP2(Params &params, const char *outfile, int total_size, bool make_bin);
	void transformEcoLP(Params &params, const char *outline, int total_size);

	/**
		transform the problem into an Integer Linear Programming and write to .lp file
		@param params program parameters
		@param outfile name of output file in LP format
		@param total_size k for PD_k or total budget
		@param make_bin TRUE if creating binary programming
	*/
	void transformLP_Area(Params &params, const char *outfile, int total_size, bool make_bin);
	void transformLP_Area2(Params &params, const char *outfile, int total_size, bool make_bin);

	/**
		transform the problem into an Integer Linear Programming and write to .lp file
		@param params program parameters
		@param outfile name of output file in LP format
		@param pd_proportion minimum PD proprotion to be conserved
		@param make_bin TRUE if creating binary programming
	*/
	void transformMinK_Area(Params &params, const char *outfile, double pd_proprotion, bool make_bin);
	void transformMinK_Area2(Params &params, const char *outfile, double pd_proportion, bool make_bin);

	/**
		transform the PD problem into linear programming and solve it
		@param params program parameters
		@param taxa_set (OUT) the vector of set of taxa in the maximal PD set
	*/
	void findPD_LP(Params &params, vector<SplitSet> &taxa_set);

	/**
		transform the PD problem into linear programming and solve it
		@param params program parameters
		@param areas_set (OUT) the vector of set of areas in the maximal PD set
	*/
	void findPDArea_LP(Params &params, vector<SplitSet> &areas_set);

	double findMinKArea_LP(Params &params, const char* filename, double pd_proportion, Split &area);

	/**
		@return TRUE if we are doing PD area optimization
	*/
	virtual bool isPDArea();

	/**
		check if all taxa are covered by the set of areas
		@return false if there exists some taxon which is not covered by any areas
	*/
	bool checkAreaCoverage();

	/**
		transform the problem into an Integer Linear Programming and write to .lp file
		@param outfile name of output file in LP format
		@param included_area (OUT) collection of areas that should always be included
	*/
	void transformLP_Area_Coverage(const char *outfile, Params &params, Split &included_area);


	/**
		@return the minimum number of areas needed to cover all taxa
		@param params program parameters
		@param area_id (OUT) minimal set of areas which cover all taxa
	*/
	int findMinAreas(Params &params, Split &area_id);



	/**
		the set of areas, each item contains the set of taxa in the area.
	*/
	SplitSet area_taxa;

	/**
	 	speciesList is used in ECOpd analysis for synchronization of species in SplitNetwork with species in FoodWeb
	 */
	void speciesList(vector<string> *speciesNames);

protected:

	/**
		extra PD when integrating initial set
	*/
	double extra_pd;

	/**
		when computing PD min (instead of PD max)
	*/
	bool min_pd;

	/**
		taxa set to be included into optimal PD set (with -i option)
	*/
	IntVector initialset;


	/**
		areas to be included into optimal PD set (with -ia option)
	*/
	IntVector initialareas;

	/**
		calculate the total maximum budget required 
		@return total maximum budget required 
	*/
	int calcMaxBudget();

/********************************************************
	hill-climbing and greedy heuristics
********************************************************/

	/**
		greedy algorithm for phylogenetic diversity of a given size 
		@param subsize the subset size
		@param taxa_set (OUT) the set of taxa in the PD-set
		@param taxa_order (OUT) order of inserted taxa
		@return the PD score of the maximal set, also returned in taxa_set.weight
	*/
	double greedyPD(int subsize, Split &taxa_set, vector<int> &taxa_order);


	/**
		local search algorithm for phylogenetic diversity of a given size 
		@param subsize the subset size
		@param taxa_set (OUT) the set of taxa in the PD-set
		@param taxa_order (IN) order of inserted taxa
		@return the PD score of the maximal set, also returned in taxa_set.weight
	*/
	double localSearchPD(int subsize, Split &taxa_set, vector<int> &taxa_order);
	
/********************************************************
	exhaustive search
********************************************************/

	/**
		exhaustive search version 2 for maximal phylogenetic diversity of a given size 
		@param subsize the subset size
		@param cur_tax current taxon
		@param curset current set
		@param find_all TRUE if wanting all max PD set
		@param best_set (OUT) the set of taxa in the maximal PD set
		@param taxa_order (OUT) order of inserted taxa
		@param rem_splits (IN) remaining splits
		@param rem_it (IN) begin of remaining iterator
		@return the PD score of the maximal set
	*/
	double exhaustPD2(int subsize, int cur_tax, Split &curset, 
		bool find_all,SplitSet &best_set, vector<int> &taxa_order, 
		IntList &rem_splits, IntList::iterator &rem_it);

	/**
		exhaustive search for maximal PD with cost-constrained
		@param cur_budget  current budget
		@param cur_tax current taxon
		@param curset current set
		@param find_all TRUE if wanting all max PD set
		@param best_set (OUT) the set of taxa in the maximal PD set
		@param taxa_order (OUT) order of inserted taxa
		@param rem_splits (IN) remaining splits
		@param rem_it (IN) begin of remaining iterator
		@return the PD score of the maximal set
	*/
	double exhaustPDBudget(int cur_budget, int cur_tax, Split &curset, 
		bool find_all,SplitSet &best_set, vector<int> &taxa_order, 
		IntList &rem_splits, IntList::iterator &rem_it);

	/**
		calculate sum of weights of preserved splits in the taxa_set
		@param taxa_set a set of taxa
		@param rem_splits remaining splits
		@param rem_it begin iterator of remaining splits
	*/
	double calcRaisedWeight(Split &taxa_set, IntList &rem_splits, IntList::iterator & rem_it);

	/**
		update the best taxa set during the search
		@param curset the current taxa set
		@param best_set the list of best taxa set so far
	*/
	void updateSplitVector(Split &curset, SplitSet &best_set);

/********************************************************
	linear programming support
********************************************************/

	/**
		y variables in the LP formulation, check if it can be dropped or equals some x variable.
		@param total_size k for PD_k or total budget
		@param y_value (OUT): vector of: -1 if cannot reduce, 1 if equals 1, or id+2 where id is the trivial split id 
	*/
	void checkYValue(int total_size, vector<int> &y_value);

	/**
		y variables in the LP formulation for PD area optimization, check if it can be dropped or equals some x variable.
		@param total_size k for PD_k or total budget
		@param y_value (OUT) vector of: -1 if cannot reduce, 1 if can be dropped, or id+2 where id is the trivial area id
		@param count1 (OUT) count of x variables in the inequality 1 
		@param count2 (OUT) count of x variables in the inequality 2 
	*/
	void checkYValue_Area(int total_size, vector<int> &y_value, vector<int> &count1, vector<int> &count2);

	/**
		check if a taxon is uniquely covered by one area
		@param taxon the taxon ID
		@param area (OUT) area the area ID that covers taxon
		@return TRUE if the 'taxon' is uniquely covered by only one area. Otherwise FALSE.
	*/
	bool isUniquelyCovered(int taxon, int &area);

	void lpObjectiveMaxSD(ostream &out, Params &params, IntVector &y_value, int total_size);

	void lpObjectiveMinK(ostream &out, Params &params);

	void lpSplitConstraint_RS(ostream &out, Params &params, IntVector &y_value, IntVector &count1, IntVector &count2, int total_size);
	void lpSplitConstraint_TS(ostream &out, Params &params, IntVector &y_value, int total_size);

	void lpK_BudgetConstraint(ostream &out, Params &params, int total_size);

	void lpMinSDConstraint(ostream &out, Params &params, IntVector &y_value, double pd_proportion);

	void lpVariableBound(ostream &out, Params &params, Split &included_vars, IntVector &y_value);

	void lpBoundaryConstraint(ostream &out, Params &params);

	void lpVariableBinary(ostream &out, Params &params, Split &included_vars);

	void lpVariableBinary(const char *outfile, Params &params, Split &included_vars);
	void lpVariableBinary(const char *outfile, Params &params, IntVector &initialset);
	void lpInitialArea(ostream &out, Params &params);

	void computeFeasibleBudget(Params &params, IntVector &list_k);

};

#endif
