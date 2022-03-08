/*
 * modelcodon.h
 *
 *  Created on: May 24, 2013
 *      Author: minh
 */

#ifndef MODELCODON_H_
#define MODELCODON_H_

#include "modelmarkov.h"

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
class ModelCodon: public ModelMarkov {
public:
	/**
		constructor
		@param model_name model name, e.g., GY,YN
		@param freq state frequency type
		@param tree associated phylogenetic tree
	*/
	ModelCodon(const char *model_name, string model_params, StateFreqType freq, string freq_params,
    		PhyloTree *tree);

	/**
	 * destructor
	 */
	virtual ~ModelCodon();

    /**
        start structure for checkpointing
    */
    virtual void startCheckpoint();

    /**
        save object into the checkpoint
    */
    virtual void saveCheckpoint();

    /**
        restore object from the checkpoint
    */
    virtual void restoreCheckpoint();

	/**
		@return the number of rate entries, equal to the number of non-diagonal elements of the rate matrix
        since we store full matrix here
	*/
	virtual int getNumRateEntries() { return num_states*(num_states); }

	/**
		initialization, called automatically by the constructor, no need to call it
		@param model_name model name, e.g., JC, HKY.
		@param freq state frequency type
	*/
	virtual void init(const char *model_name, string model_params, StateFreqType freq, string freq_params);

	StateFreqType initCodon(const char *model_name, StateFreqType freq, bool reset_params, string freq_params="");


	/**
	 * @return model name with parameters in form of e.g. GTR{a,b,c,d,e,f}
	 */
	virtual string getNameParams(bool show_fixed_params = false) { return name; }

    /** main function to compute rate matrix */
    void computeCodonRateMatrix();

	/**
		decompose the rate matrix into eigenvalues and eigenvectors
	*/
	virtual void decomposeRateMatrix();

	/**
	 * read codon model from a stream, modying rates and state_freq accordingly
	 * @param in input stream containing lower triangular matrix of rates, frequencies and list of codons
     * @reset_params true to reset parameters, false otherwise
	 */
	void readCodonModel(istream &in, bool reset_params);

	/**
	 * read codon model from a string, modying rates and state_freq accordingly
	 * @param str input string containing lower triangular matrix of rates, frequencies and list of codons
     * @reset_params true to reset parameters, false otherwise
	 */
	void readCodonModel(string &str, bool reset_params);

	/**
	 * read codon model from a file in PAML format
	 * @param filename input file containing lower triangular matrix of rates, frequencies and list of codons
     * @reset_params true to reset parameters, false otherwise
	 */
    void readCodonModelFile(const char *filename, bool reset_params);

	/**
		write information
		@param out output stream
	*/
	virtual void writeInfo(ostream &out);

    /** compute rate_attr for all codoni->codoni substitution */
    void computeRateAttributes();
    
    /** combine rates with target nucleotide frequency (ntfreq) for MG-style model */
    void combineRateNTFreq();

    /** compute the corrected empirical omega (Kosiol et al 2007) */
    double computeEmpiricalOmega();
    
	/**
	 * setup the bounds for joint optimization with BFGS
	 */
	virtual void setBounds(double *lower_bound, double *upper_bound, bool *bound_check);

	/**
		optimize model parameters
		@return the best likelihood
	*/
	virtual double optimizeParameters(double gradient_epsilon);

	/** 3x4 matrix of nucleotide frequencies at 1st,2nd,3rd codon position */
	double *ntfreq;

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
    
	/** empirical rates for empirical codon model or parametric+empirical codon model */
	double *empirical_rates;
    
protected:

    void computeCodonRateMatrix_1KAPPA();
    void computeCodonRateMatrix_1KAPPATS();
    void computeCodonRateMatrix_1KAPPATV();
    void computeCodonRateMatrix_2KAPPA();

	/** initialize Muse-Gaut 1994 model 
        @param fix_kappa whether or not to fix kappa
        @param freq input frequency
        @param freq_params is user-specified frequency params (assuming 4 params for F1x4 and 12 params for F3x4).
        @return default frequency type
    */
	StateFreqType initMG94(bool fix_kappa, StateFreqType freq, CodonKappaStyle kappa_style, string freq_params="");

	/** initialize Goldman-Yang 1994 model (simplified version with 2 parameters omega and kappa 
        @param fix_kappa whether or not to fix kappa
        @param kappa_style: CK_ONE_KAPPA for traditional GY model, others follow Kosiol et al 2007
        @return default frequency type
    */
	StateFreqType initGY94(bool fix_kappa, CodonKappaStyle kappa_style);

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
		@return TRUE if parameters are changed, FALSE otherwise (2015-10-20)
	*/
	virtual bool getVariables(double *variables);

};

#endif /* MODELCODON_H_ */
