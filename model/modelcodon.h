/*
 * modelcodon.h
 *
 *  Created on: May 24, 2013
 *      Author: minh
 */

#ifndef MODELCODON_H_
#define MODELCODON_H_

#include "modelgtr.h"

/**
 * parameter constraint
 */
struct ParamConstraint {
	bool fixed; // TRUE if this parameter is fixed
	// minimum, initial, and maximum value
	double min_value, init_value, max_value;
	char opr, opr2; // operator: '*', '/', or 0, to force this parameter
	int param1, param2, param3; // index of 2 parameters for operator
	double opr_value, opr_value2; // instead of multiplying 2 parameters, one can multiply a parameter with this constant
};

/** CF_TARGET_NT: frequency of target nucleotide is multiplied with the rate entry (Muse and Gaut 1994)
    CF_TARGET_CODON: frequency of target codon is multiplied with the rate entry (Goldman Yang 1994)
    */
enum CodonFreqStyle {CF_TARGET_NT, CF_TARGET_CODON};

enum CodonKappaStyle {CK_ONE_KAPPA, CK_ONE_KAPPA_TS, CK_ONE_KAPPA_TV, CK_TWO_KAPPA};

const int CA_STOP_CODON   = 1; // stop codon substitution
const int CA_MULTI_NT     = 2; // codon substitution involves > 1 NT
const int CA_SYNONYMOUS   = 4; // synonymous codon substitution
const int CA_NONSYNONYMOUS= 8; // synonymous codon substitution
const int CA_TRANSVERSION = 16; // codon substitution involves 1 NT transversion
const int CA_TRANSITION   = 32; // codon substitution involves 1 NT transition
const int CA_TRANSVERSION_1NT = 64; // codon substitution involve the 1st NT which is also a transversion
const int CA_TRANSVERSION_2NT = 128; // codon substitution involve the 2nd NT which is also a transversion
const int CA_TRANSVERSION_3NT = 256; // codon substitution involve the 3rd NT which is also a transversion
const int CA_TRANSITION_1NT   = 512; // codon substitution involve the 1st NT which is also a transversion
const int CA_TRANSITION_2NT   = 1024; // codon substitution involve the 2nd NT which is also a transversion
const int CA_TRANSITION_3NT   = 2048; // codon substitution involve the 3rd NT which is also a transversion

/**
 * Codon substitution models
 */
class ModelCodon: public ModelGTR {
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

	/**
		@return the number of rate entries, equal to the number of non-diagonal elements of the rate matrix
        since we store full matrix here
	*/
	virtual int getNumRateEntries() { return num_states*(num_states); }

	StateFreqType initCodon(const char *model_name, StateFreqType freq);

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

    /** main function to compute rate matrix */
    void computeCodonRateMatrix();

	/**
		decompose the rate matrix into eigenvalues and eigenvectors
	*/
	virtual void decomposeRateMatrix();

	/**
	 * set rates into groups, rates within a group are equal
	 * @param group assignment of each rate into group
	 */
//	void setRateGroup(IntVector &group);

	/**
	 * set rates into groups, rates within a group are equal
	 * @param group assignment of each rate into group
	 */
//	void setRateGroup(const char *group);

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
//	void setRateGroupConstraint(string constraint);

	/**
		Read the rate parameters from a comma-separated string
		It will throw error messages if failed
		@param in input stream
	*/
//	virtual void readRates(string str) throw(const char*);

	/**
	 * @return true if codon1<->codon2 involves more than 1 nucleotide
	 */
//	bool isMultipleSubst(int state1, int state2);

	/**
	 * @return if single nucleotide substitution i (0<=i<=3) is involved from state1->state2, return
	 * j*4+i, where j is the codon position of substitution. Otherwise return -1.
	 */
//	int targetNucleotide(int state1, int state2);

	/**
	 * @return true if codon1<->codon2 is a synonymous substitution
	 */
//	bool isSynonymous(int state1, int state2);

	/**
	 * @return true if codon1<->codon2 involves exactly one nucleotide transversion
	 */
//	bool isTransversion(int state1, int state2);

    /**
        count number of transition and transversion between 2 codons
        @param state1 codon 1
        @param state2 codon 2
        @param ts (OUT) number of transitions (between 0 and 3)
        @param tv (OUT) number of transversions (between 0 and 3)
    */
//    void countTsTv(int state1, int state2, int &ts, int &tv);

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

    /** dn/ds rate ratio */
    double omega;
    
    /** TRUE to fix omega, default: FALSE */
    bool fix_omega; 

    /** style for kappa */
    CodonKappaStyle codon_kappa_style;

    /** ts/tv rate ratio */
    double kappa;
    
    /** TRUE to fix kappa, default: FALSE */
    bool fix_kappa;

    /** ts/tv rate ratio for 2-kappa model (Kosiol et al 2007) */
    double kappa2;
    
    /** TRUE to fix kappa2, default: FALSE */
    bool fix_kappa2;
    
    /** GY- or MG-style codon frequencies */
    CodonFreqStyle codon_freq_style;
    
    /** rate atrributes */
    int *rate_attr;
    
    /** compute rate_attr for all codoni->codoni substitution */
    void computeRateAttributes();
    
//    bool combined_rate_ntfreq;
    
    /** combine rates with target nucleotide frequency (ntfreq) for MG-style model */
    void combineRateNTFreq();

    /** separate rates from target nucleotide frequency (ntfreq) for MG-style model */
//    void separateRateNTFreq();
    
protected:

    void computeCodonRateMatrix_1KAPPA();
    void computeCodonRateMatrix_1KAPPATS();
    void computeCodonRateMatrix_1KAPPATV();
    void computeCodonRateMatrix_2KAPPA();

	/** initialize Muse-Gaut 1994 model 
        @param fix_kappa whether or not to fix kappa
        @param freq input frequency
        @return default frequency type
    */
	StateFreqType initMG94(bool fix_kappa, StateFreqType freq, CodonKappaStyle kappa_style);

	/** initialize Goldman-Yang 1994 model (simplified version with 2 parameters omega and kappa 
        @param fix_kappa whether or not to fix kappa
        @param kappa_style: CK_ONE_KAPPA for traditional GY model, others follow Kosiol et al 2007
        @return default frequency type
    */
	StateFreqType initGY94(bool fix_kappa, CodonKappaStyle kappa_style);
//
//    /**
//        @param kappa_ts true for Kappa(ts) model, false for Kappa(tv) (Kosiol et al 2007)
//    */
//	void initGY94_1K(bool kappa_ts);
//
//    /**
//        Kappa(ts)Kappa(tv) model of Kosiol et al 2007
//    */
//	void initGY94_2K();

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
//	IntVector rate_group;

	/** constraint for each rate group */
//	vector<ParamConstraint> rate_constraints;

	/** empirical rates for empirical codon model or parametric+empirical codon model */
	double *empirical_rates;

	/** extra rate multiplier (e.g., frequency of target nucleotide for MG model) */
//	double *extra_rates;
};

#endif /* MODELCODON_H_ */
