/*
 * modelcodon.h
 *
 *  Created on: May 24, 2013
 *      Author: minh
 */

#ifndef MODELCODON_H_
#define MODELCODON_H_

#include "gtrmodel.h"

/**
 * parameter constraint
 */
struct ParamConstraint {
	bool fixed; // TRUE if this parameter is fixed
	// minimum, initial, and maximum value
	double min_value, init_value, max_value;
	char opr; // operator: '*', '/', or 0, to force this parameter
	int param1, param2; // index of 2 parameters for operator
	double opr_value; // instead of multiplying 2 parameters, one can multiply a parameter with this constant
};

/**
 * Codon substitution models
 */
class ModelCodon: public GTRModel {
public:
	/**
		constructor
		@param model_name model name, e.g., GY,YN
		@param freq state frequency type
		@param tree associated phylogenetic tree
	*/
	ModelCodon(const char *model_name, string model_params, StateFreqType freq, string freq_params,
    		PhyloTree *tree, bool count_rates = true);

	/**
	 * destructor
	 */
	virtual ~ModelCodon();

	StateFreqType initCodon(const char *model_name);

	/**
		initialization, called automatically by the constructor, no need to call it
		@param model_name model name, e.g., JC, HKY.
		@param freq state frequency type
	*/
	virtual void init(const char *model_name, string model_params, StateFreqType freq, string freq_params);

	/**
	 * @return model name with parameters in form of e.g. GTR{a,b,c,d,e,f}
	 */
	virtual string getNameParams() { return name; }

	/**
	 * set rates into groups, rates within a group are equal
	 * @param group assignment of each rate into group
	 */
	void setRateGroup(IntVector group);

	/**
	 * set rates into groups, rates within a group are equal
	 * @param group assignment of each rate into group
	 */
	void setRateGroup(const char *group);

	/**
	 * Set constraints for rate-groups, a comma-separated string of constraints.
	 * Each constraint has the following format:
	 *   xi=?       : rate of group x_i will be estimated from data
	 *   xi=?value  : rate of group x_i will be initialized at value and then estimated from data
	 *   xi=value   : rate of group x_i is fixed at a specific floating-point value
	 *   xi>value   : rate of group x_i must be > value
	 *   xi<value   : rate of group x_i must be < value
	 *   xi=xj*xk   : rate of group x_i is constrained to equal to group x_j * group x_k
	 *   xi=xj/xk   : rate of group x_i is constrained to equal to group x_j / group x_k
	 *   @param constraint comma-separated string of constraints
	 */
	void setRateGroupConstraint(string constraint);

	/**
		Read the rate parameters from a comma-separated string
		It will throw error messages if failed
		@param in input stream
	*/
	virtual void readRates(string str) throw(const char*);

	/**
	 * @return true if codon1<->codon2 involves more than 1 nucleotide
	 */
	bool isMultipleSubst(int state1, int state2);

	/**
	 * @return if single nucleotide substitution i (0<=i<=3) is involved from state1->state2, return
	 * j*4+i, where j is the codon position of substitution. Otherwise return -1.
	 */
	int targetNucleotide(int state1, int state2);

	/**
	 * @return true if codon1<->codon2 is a synonymous substitution
	 */
	bool isSynonymous(int state1, int state2);

	/**
	 * @return true if codon1<->codon2 involves exactly one nucleotide transversion
	 */
	bool isTransversion(int state1, int state2);


	/** 3x4 matrix of nucleotide frequencies at 1st,2nd,3rd codon position */
	double *ntfreq;


	/**
	 * read codon model from a stream, modying rates and state_freq accordingly
	 * @param in input stream containing lower triangular matrix of rates, frequencies and list of codons
	 */
	void readCodonModel(istream &in);

	/**
	 * read codon model from a string, modying rates and state_freq accordingly
	 * @param str input string containing lower triangular matrix of rates, frequencies and list of codons
	 */
	void readCodonModel(string &str);

	/**
		write information
		@param out output stream
	*/
	virtual void writeInfo(ostream &out);

protected:

	/** initialize Muse-Gaut 1994 model */
	void initMG94();

	/** initialize Goldman-Yang 1994 model (simplified version with 2 parameters omega and kappa */
	void initGY94();

	/**
		this function is served for the multi-dimension optimization. It should pack the model parameters
		into a vector that is index from 1 (NOTE: not from 0)
		@param variables (OUT) vector of variables, indexed from 1
	*/
	virtual void setVariables(double *variables);

	/**
		this function is served for the multi-dimension optimization. It should assign the model parameters
		from a vector of variables that is index from 1 (NOTE: not from 0)
		@param variables vector of variables, indexed from 1
	*/
	virtual void getVariables(double *variables);

	/** assignment of each rate into group */
	IntVector rate_group;

	/** constraint for each rate group */
	vector<ParamConstraint> rate_constraints;

	/** empirical rates for empirical codon model or parametric+empirical codon model */
	double *empirical_rates;

	/** extra rate multiplier (e.g., frequency of target nucleotide for MG model) */
	double *extra_rates;
};

#endif /* MODELCODON_H_ */
