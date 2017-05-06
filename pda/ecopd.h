/*
 * ecopd.h
 *
 *  Created on: Oct 30, 2013
 *      Author: Olga
 */

#ifndef ECOPD_H
#define ECOPD_H

#include "tree/mtree.h"
#include "tree/mtreeset.h"
#include "ecopdmtreeset.h"
#include "pdnetwork.h"

/* ===============================================================================
 *	Class for processing IP problem - PD/SD with ecological constraints
 * ===============================================================================*/

class ECOpd : public MTree
{
public:

	/**
		CONSTRUCTORs, INITIALIZATION and DESTRUCTOR
	*/
    ECOpd(const char *userTreeFile, bool &is_rooted);
    ECOpd();
   ~ECOpd();

   void initializeEcoPD();
   void initializeEcoPD(Params &params);

   /*
    * Checks whether taxon with id i ("real" id, that is when calling call with i+1 ) is present on the tree/split network
    */
   bool OUT_tree(int i);

   /*
    * Reading and processing Diet Composition matrix
    */
   void readDAG(const char *infile);
   void readDAG(istream &in);

   /*
    * Transform problem into IP problem and print it to .lp file, for rooted trees
    */
   void printECOlpRooted(const char* fileOUT,ECOpd &tree);

   /*
    * Transform problem into IP problem and print it to .lp file, for UNrooted trees
    */
   void printECOlpUnrooted(const char* fileOUT,ECOpd &tree);

   /*
    * Transform problem into IP problem and print it to .lp file, for split system
    */
	void printInfDAG (const char* fileOUT,PDNetwork &splitsys,Params &params);

	/*
	 * Synchronization of species in the food web with species on the tree
	 */
	void synchTreeDAG(ECOpd &tree);

	/*
	 * Synchronization of species in the food web with species in the split network
	 */
	void synchSplitDAG(PDNetwork &system);

	/*
	 * some left_overs from mtree class, function which is not there anymore..
	 */
	void getBranchOrdered(NodeVector &nodes, NodeVector &nodes2,Node *node = NULL, Node *dad = NULL);

	/*
	 * Find the id of the species on tree by name
	 */
	int findPhyloID(string name);

	/*
	 * Find the id of the species in the food web by their phylo id (id of this species in the tree)
	 */
	int findFoodWebID(int id);

	/*
	 * checks whether there are some species present on the tree/splitSys, but not in the foodWeb, and wise versa.
	 */
	void detectMissingSpecies();

	/*
	 * List of species missing either on the tree/splitSys (missInPhylo), or in the food web (missInDAG)
	 */
	vector<string> missInPhylo,missInDAG;

	/*
	 * Finding Taxon from tree/splitSys among DAG Species
	 */
	bool findTaxaDAG(int i);

	/*
	 * Finding Species from DAG among Taxa on tree/splitSys
	 */
	bool findSpeciesPhylo(int i);

	/*
	 * the number of links in the food web
	 */
	int linksNUM;

	/*
	 * synchronization of species on Tree/SplitSys and in Food Web
	 */
	void synchronizeSpecies();

	/*
	 * Reading taxa to be included in the final optimal subset
	 */
	void readInitialTaxa(const char *infile);
	void readInitialTaxa(istream &in);

	/*
	 * list of taxa (names) to be included in the final optimal set
	 */
	vector<string> initialTaxa;

	/*
	 * Check if the species in InitialTaxa are actually present either on tree/network or in the food web
	 */
	void checkInitialTaxa();

	/*
	 * find an id (among nvar) of a given species by name
	 */
	int findSpeciesIDname(string *name);

	/*
	 * Define the subset size and check if it's >1 and <nvar (#of all species in the analysis = (TaxaNUM > SpeciesNUM) ? TaxaNUM : SpeciesNUM)
	 */
	void defineK(Params &params);

	/*
	 * Check whether the food web is acyclic or not
	 */
	void checkGraph();

	/*
	 * Diet Composition Matrix (entries either 0/1 or [0,100] )
	 */
	vector<double*> DAG;

	/*
	 * Prints the sub food web corresponding to the optimal subset
	 */
	void printSubFoodWeb(char* fileOUT, double* variables);

	/*
	 * t for tree or n for networks
	 */
	string phyloType;

	/*
	 * Structure of the DAG: taxa with neighbors being their preys
	 */
	NodeVector taxaDAG;

	/*
	 * two vectors of nodes, corresponding to ends of branches
	 */
	NodeVector nodes1,nodes2;

	/*
	 * Ids of species not present on tree/split network
	 */
	vector<int> OUTtreeTaxa;

	/*
	 * Ids of species not present in the food web
	 */
	vector<int> OUTdagTaxa;

	/*
	 * contains the ids of species based on tree/split network (phylo_oder[i] = j species with id=j in the food web has id=i on tree/split network)
	 */
	vector<int> phylo_order;

	/*
	 * for each species contains information about its longest food chain (excluding species itself)
	 */
	vector<int> levelDAG;

	/*
	 * the names of species present in the food web and on the tree/split network respectively
	 */
	vector<string> dagNames,phyloNames;

	/*
	 * names of all species: union of species on the tree/network and in the food web
	 */
	vector<string*> names;

	/*
	 * flag for whether to treat the food web as weighted or not weighted
	 */
	bool weighted;

	/*
	 * the size of an optimal subset to be chosen
	 */
	int k;

	/*
	 * the diet portion to be conserved for each predator, when equals to 0 corresponds to a naive viability
	 */
	double T;

	/*
	 * the number of species in the food web (if rooted tree counts also the root, technical)
	 */
	int SpeciesNUM;

	/*
	 * the number of species on the tree/split network (if rooted tree counts also the root)
	 */
	int TaxaNUM;

	/*
	 * the number of all species: union of species on the tree/network and in the food web
	 */
	int nvar;

	/*
	 * calculates for each predator the diet proportional conserved
	 */
	void dietConserved(double *variables);

	/*
	 * for each predator the diet proportional conserved
	 */
	vector<double> dietVAL;

	/*
	 * print the results
	 */

	void printResults(char* fileOUT,double* variables, double score, Params &params);

	/*
	 * Splits number and total SD
	 */
	int splitsNUM;
	double totalSD;

	/**************************
	 * Miscellaneous
	 **************************/

	/*
	 * These function were used when we analyzed results from LP problems, to set the fractional values to be integers
	 * now, it is not used
	 */
	void readREC(const char *infile);
	void readREC(istream &in);
	void generateFirstMultinorm(vector<int> &x, int n, int k);
	bool generateNextMultinorm(vector<int> &x);

	vector<string> fractVAR;
	vector<int> dvec,hvec;

	/*
	 * Assigns to a given tree topology random branch lengths
	 */
	void randomBranLenTrees(Params &params);

};

#endif
