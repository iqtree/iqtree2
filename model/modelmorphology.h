/*
 * modelmorphology.h
 *
 *  Created on: Apr 15, 2014
 *      Author: minh
 */

#ifndef MODELMORPHOLOGY_H_
#define MODELMORPHOLOGY_H_

#include "modelmarkov.h"

/**
 * This class implement ML model for morphological data. Such models are:
 * - Mk (Lewis 2001) a JC-type model
 * - ORDERED: allowing only transition from state i to i-1 and i+1
 * TODO: Mkv to account for absence of constant sites
 */
class ModelMorphology: public ModelMarkov {
public:
	/**
		constructor
		@param model_name model name, e.g., JC, HKY.
		@param freq state frequency type
		@param tree associated phylogenetic tree
	*/
	ModelMorphology(const char *model_name, string model_params, StateFreqType freq, string freq_params, PhyloTree *tree);


	/**
		initialization, called automatically by the constructor, no need to call it
		@param model_name model name, e.g., JC, HKY.
		@param freq state frequency type
	*/
	virtual void init(const char *model_name, string model_params, StateFreqType freq, string freq_params);

    /**
     return the number of dimensions
     */
    virtual int getNDim();

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
     * @return model name with parameters in form of e.g. GTR{a,b,c,d,e,f}
     */
    virtual string getNameParams(bool show_fixed_params = false);

    /**
     write information
     @param out output stream
     */
    virtual void writeInfo(ostream &out);

    /**
     write parameters, used with modeltest
     @param out output stream
     */
    virtual void writeParameters(ostream &out);

	/**
		read the rates from an input stream. it will throw error messages if failed
		@param in input stream
	*/
	virtual void readRates(istream &in) throw(const char*, string);

    virtual ~ModelMorphology();
};

#endif /* MODELMORPHOLOGY_H_ */
