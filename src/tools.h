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

#ifndef TOOLS_H
#define TOOLS_H

#include <vector>
#include <string>
#include <set>
#include <map>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

//#include <sys/time.h>
//#include <time.h>
#include <sys/stat.h>
#include <math.h>
#include "ncl/ncl.h"
#include "msetsblock.h"


#define USE_HASH_MAP

#ifdef USE_HASH_MAP
	#if defined(WIN32)
	#	include <hash_map>
	using namespace stdext;
	#else
	#	include <ext/hash_map>
	using namespace __gnu_cxx;
	#endif
#else
	#include <map>
#endif


using namespace std;


/**
	vector of double number
*/
typedef vector<double> DoubleVector;

/**
	vector of int
*/
typedef vector<int> IntList;


/**
	vector of int
*/
typedef vector<int> IntVector;

/**
	vector of bool
*/
typedef vector<bool> BoolVector;


/**
	vector of char
*/
typedef vector<char> CharVector;

/**
	vector of string
*/
typedef vector<string> StrVector;


/**
	matrix of double number
*/
#define matrix(T) vector<vector<T> >

/**
	matrix of double
*/
/*
class DoubleMatrix {
public:
	double *value;
	int rows, cols, size;
	DoubleMatrix(int arows, int acols);
	//inline double operator() (int i, int j);
	inline double &operator() (int i, int j) { return value[i * cols + j]; };
	inline double *operator[] (int i) {	return value + (i*cols); };
	virtual ~DoubleMatrix();
	void setZero();
};
*/
typedef matrix(double) DoubleMatrix;

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/**
	run mode of program
*/
enum RunMode {DETECTED, GREEDY, PRUNING, BOTH_ALG, EXHAUSTIVE, DYNAMIC_PROGRAMMING, CALC_DIST, PD_USER_SET, PRINT_TAXA, PRINT_AREA, SCALE_BRANCH_LEN, SCALE_NODE_NAME, PD_DISTRIBUTION, LINEAR_PROGRAMMING};

/**
	type of generating trees or splits graphs
*/
enum TreeGenType {NONE, YULE_HARDING, UNIFORM, CATERPILLAR, BALANCED,
		CIRCULAR_SPLIT_GRAPH, TAXA_SET};

/**
	when writing tree:
		BR_LEN - output branch length
		BR_NONE - don't output branch length
		BR_CLADE - put branch length into internal node name
*/
const int WT_BR_LEN    = 1;
const int WT_BR_CLADE  = 2;
const int WT_TAXON_ID  = 4;

/**
	search mode
*/
//enum SearchMode {EXHAUSTIVE, EXHAUSTIVE_CIRCULAR};

/**
	input type, tree or splits graph
*/
enum InputType {IN_NEWICK, IN_NEXUS, IN_OTHER};

/**
	verbose mode, determine how verbose should the screen be printed.
*/
enum VerboseMode {VB_MIN, VB_MED, VB_MAX, VB_DEBUG};

/**
	verbose level on the screen
*/
extern VerboseMode verbose_mode;

/**
	consensus reconstruction type
*/
enum ConsensusType {CONSENSUS_TREE, CONSENSUS_NETWORK, ASSIGN_BOOTSTRAP};

enum TestType {TEST_NONE, TEST_COMPATIBLE, TEST_CIRCULAR, TEST_WEAKLY_COMPATIBLE, TEST_K_COMPATIBLE};

/**
	State frequency type
*/
enum StateFreqType {FREQ_UNKNOWN, FREQ_USER_DEFINED, FREQ_EQUAL, FREQ_EMPIRICAL, FREQ_ESTIMATE};


const double MAX_GENETIC_DIST = 100.0;

extern bool simple_nni;

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/**
	program parameters, everything is specified here
*/
struct Params {

	/**
		 input file name
	*/
	char *user_file;

	/**
		alignment file name
	*/
	char *aln_file;

	/**
		compute parsimony score on trees
	*/
	bool parsimony;

	/**
		 output file name
	*/
	char *out_file;

	/**
		 size of the maximal PD-tree
	*/
	int sub_size;

	/**
		 min size of the maximal PD-tree
		 used to calculate all PD-k trees from min_size to sub_size
	*/
	int min_size;

	/**
		step_size when running from min_size to sub_size
	*/
	int step_size;

	/**
		sample size for computing PD distribution
	*/
	int sample_size;


	/**
		TRUE if want to find all optimal PD-k set
		with the same maximal PD score
	*/
	bool find_all;

	/**
		 type of random tree to be generated
	*/
	TreeGenType tree_gen;

	/**
		when generating random split graph, specify the number of
		splits here!
	*/
	int num_splits;

	/**
		 running mode: which algorithms to be applied
	*/
	RunMode run_mode;

	/**
		 real running mode if run_mode == DETECTED
	*/
	RunMode detected_mode;

	/**
		 parameter file
	*/
	char *param_file;

	/**
		file containing taxa names to be included into the PD-tree
	*/
	char *initial_file;

	/**
		file containing area names to be included into the PD set
	*/
	char *initial_area_file;

	/**
		file containing a list of specific taxa sets which user wants
		to compute PD score on these sets only
	*/
	char *pdtaxa_file;

	/**
		output file to store the distance matrix
	*/
	char *dist_file;

	/**
		file containing budget information
	*/
	char *budget_file;

	/**
		used when generating pair of taxa set with overlapping
	*/
	int overlap;

	// private use
	/**
		 number of times to repeat the algorithms
	*/
	int repeated_time;

	/**
		 print no tree to output
	*/
	int nr_output;

	/**
		input type, tree or splits graph
	*/
	InputType intype;

	/**
		total budget, for cost constrained PD problem
	*/
	int budget;

	/**
		minimum budget, for cost constrained PD problem
	*/
	int min_budget;

	/**
		step_budget when running from min_budget to budget
	*/
	int step_budget;

	/**
		name of the root taxon
	*/
	char *root;

	/**
		true if tree is forced to be rooted
	*/
	bool is_rooted;


	/**
		min branch length, used to create random tree/network
	*/
	double min_len;

	/**
		mean branch length, used to create random tree/network
	*/
	double mean_len;

	/**
		max branch length, used to create random tree/network
	*/
	double max_len;

	/**
		random number seed
	*/
	unsigned int ran_seed;

	/**
		run time of the algorithm
	*/
	long run_time;

	/**
		limit on the number of optimal PD sets
	*/
	int pd_limit;

	/**
		TRUE if one wants to calculate the PD gain matrix in terms of delta_k^j = pd(PD_k \/ {j}) - pd_k
	*/
	bool calc_pdgain;

	/**
		TRUE if tree file contains more than 1 tree
	*/
	bool multi_tree;

	/**
		file name containing all trees from bootstrap analysis
	*/
	char *boot_trees;

	/**
		type of consensus building
	*/
	ConsensusType calc_consensus;

	/**
		set the TRUE if want to find the minimal PD set, instead of the default maximal PD set
	*/
	bool find_pd_min;

	/**
		set TRUE to find area's endemic PD instead of regular PD
	*/
	bool endemic_pd;

	/**
		set TRUE to find exclusive PD instead of regular PD
	*/
	bool exclusive_pd;

	/**
		to find PD complementarity given this area
	*/
	char *complement_area;

	/**
		used for likelihood mapping: for each branch, print the four cluster
	*/
	int branch_cluster;

	/**
		file containing taxa order
	*/
	char *taxa_order_file;

	/**
		to scale branch length or clade support with a factor
	*/
	double scaling_factor;

	/**
		TRUE if always use binary linear programming
	*/
	bool binary_programming;

	/**
		test the input split system in one of the TestType
	*/
	TestType test_input;

	/**
		burnin value: number of beginning trees to be discarded
	*/
	int tree_burnin;

	/**
		threshold of split frequency, splits appear less than threshold will be discarded
	*/
	double split_threshold;

	/**
		true if one wants to optimize tree by subtree pruning and regrafting
	*/
	bool tree_spr;

	/**
		true if printing out of optimal sets in NEXUS format
	*/
	bool nexus_output;

	/**
		k-representative parameter, used for IQP algorithm
	*/
	int k_representative;

	/**
		probability of deleting a leaf, used for IQP algorithm
	*/
	double p_delete;

	/**
		number of iqpnni iterations
	*/
	int iqpnni_iterations;

	/**
		name of the substitution model (e.g., HKY, GTR, TN+I+G, JC+G, etc.)
	*/
	string model_name;

	/**
		state frequency type
	*/
	StateFreqType freq_type;


	/**
		the number of rate categories
	*/
	int num_rate_cats;

	/**
		TRUE if you want to optimize branch lengths by Newton-Raphson method
	*/
	bool optimize_by_newton;

};

/**
	related measures for PD
*/
struct PDRelatedMeasures {
	/**
		names of areas
	*/
	vector<NxsString> setName;

	/**
		pd scores of areas
	*/
	DoubleVector PDScore;

	/**
		exclusive PD scores
	*/
	DoubleVector exclusivePD;

	/**
		endemic pd scores of an area given all other areas
	*/
	DoubleVector PDEndemism;

	/**
		pd-complementarity scores of an area given some provided area
	*/
	DoubleVector PDComplementarity;
};



/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/


/**
	@return TRUE of ch is a control character (ascii <= 32)
*/
inline bool controlchar(char ch) {
	return ch <= 32;
}

inline bool is_newick_token(char ch) {
	return ch == ':' || ch == ';' || ch == ',' || ch == ')' || ch == '(';
}

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/**
	print error message then exit program
*/
//void outError(char *error);

/**
	print error message then exit program
*/
void outError(const char *error);

/**
	print error message then exit program
*/
void outError(string error);


/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/**
	print double error messages then exit program
*/
void outError(const char *error, const char *msg);

/**
	print double error messages then exit program
*/
void outError(const char *error, string msg);

/**
	Output a warning message to screen
	@param error warning message
*/
void outWarning(const char *warn);
void outWarning(string warn);


/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/


/**
	generate a random branch length under an exponential distribution
	with mean params.mean_len. Also make sure that the resulting
	length is in the range (params.min_len, params.max_len)
	@return the random branch length
*/
double randomLen(Params &params);

/**
	convert string to int, with error checking
	@param str original string
	@return the integer
*/



/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*
	Error messages
*/
const char ERR_NO_TAXON[]          = "Find no taxon with name ";
const char ERR_NO_AREA[]          = "Find no area with name ";
const char ERR_NO_ROOT[]           = "Root taxon not found: ";
const char ERR_ROOT_NET[]          = "-root option is not available for network";
const char ERR_CONFLICT_ROOT[]     = "Tree is already rooted, -o <taxon> is not allowed.";
const char ERR_DUPLICATED_TAXA[]   = "Duplicated taxa name in the tree.";
const char ERR_FEW_TAXA[]          = "Number of taxa must be greater than 2.";
const char ERR_NO_SPLITS[]         = "No splits found!";
const char ERR_FEW_SPLITS[]        = "Number of splits must be at least equal to the number of taxa";
const char ERR_NEG_BRANCH[]        = "Negative branch length not allowed.";
const char ERR_NO_MEMORY[]         = "Not enough memory!";

const char ERR_READ_INPUT[]        = "File not found or incorrect input, pls check it again.";
const char ERR_UNEXPECTED_EOF[]    = "Unexpected end of file.";
const char ERR_READ_ANY[]          = "Unidentified error while reading file, pls check it carefully again.";
const char ERR_WRITE_OUTPUT[]      = "Cannot write to file ";

const char ERR_NO_K[]              = "You must specify the number of taxa in the PD set.";
const char ERR_TOO_SMALL_K[]       = "Size of PD-set must be at least the size of initial set.";
const char ERR_NO_BUDGET[]         = "Total budget is not specified or less than zero.";
const char ERR_TOO_SMALL_BUDGET[]  = "Not enough budget to conserve the inital set of taxa.";

const char ERR_INTERNAL[]          = "Internal error, pls contact authors!";

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/**
 * convert int to string
 * @param int
 * @return string
 */
string convertIntToString(int number);

/**
 * Check if the file exists
 * @param strFilename
 * @return
 */
bool fileExists(string strFilename);

int convert_int(const char *str) throw (string);

/**
	convert string to double, with error checking
	@param str original string
	@return the double
*/
double convert_double(const char *str) throw (string);


/**
	convert a string to to range lower:upper:step_size with error checking
	@param str original string
	@param lower (OUT) lower bound of the range
	@param upper (OUT) upper bound of the range
	@param step_size (OUT) step size of the range
*/
void convert_range(const char *str, int &lower, int &upper, int &step_size) throw (string);

/**
	convert a string to to range lower:upper:step_size with error checking
	@param str original string
	@param lower (OUT) lower bound of the range
	@param upper (OUT) upper bound of the range
	@param step_size (OUT) step size of the range
*/
void convert_range(const char *str, double &lower, double &upper, double &step_size) throw (string);


/**
	read the file containing branch/split scaling factor and taxa weights
	@param params program parameters
	@param ntaxa total number of taxa
	@param scale (OUT) scaling factor
	@param tax_name (OUT) vector of taxa names
	@param tax_weight (OUT) vector of corresponding taxa weights
*/
void readWeightFile(Params &params, int ntaxa, double &scale, StrVector &tax_name, DoubleVector &tax_weight);

/**
	read the initial taxa set from the file
	@param params program parameters
	@param ntaxa number of taxa
	@param tax_name (OUT) vector of taxa names
*/
void readInitTaxaFile(Params &params, int ntaxa, StrVector &tax_name);

/**
	read the initial area set from the file
	@param params program parameters
	@param nareas number of areas
	@param area_name (OUT) vector of area names
*/
void readInitAreaFile(Params &params, int nareas, StrVector &area_name);


/**
	read a list of taxa set from a file, not in nexus format but as follows:
	n1
	tax-name-1
	...
	tax-name-n1

	n2
	tax-name-1
	...
	tax-name-n2
	....

	@param filename file name
	@param sets (OUT) the returned sets of taxa
*/
void readTaxaSets(char *filename, MSetsBlock *sets);


/**
	parse program argument into params
	@param argc number of arguments
	@param argv list of arguments
	@param params (OUT) program parameters
*/
void parseArg(int argc, char *argv[], Params &params);

/**
	detect the format of input file
	@param input_file file name
	@return
		IN_NEWICK if file in newick format,
		IN_NEXUS if in nexus format,
		IN_OTHER if file format unknown.
*/
InputType detectInputFile(char *input_file);

/**
	if file exists, ask user to overwrite it or not
	@param filename file name
	@return TRUE if agree to overwrite an existing file, or simply file does not exist
*/
bool overwriteFile(char *filename);

/**
	print usage information
	@param argv program arguments list
	@param full_command TRUE to print all available commands, FALSE to print normal usage dialog
*/
void usage(char* argv[], bool full_command);

/**
	parse area name string, where names are separated by commas
	@param area_names a string of name
	@param areas (OUT) a set of name string
*/
void parseAreaName(char *area_names, set<string> &areas);

#endif
