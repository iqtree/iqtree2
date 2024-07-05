/***************************************************************************
 *   Copyright (C) 2009 by BUI Quang Minh   *
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
#ifndef MODELPROTEIN_H
#define MODELPROTEIN_H

#include "modelmarkov.h"

extern const char* builtin_prot_models;

/**
Substitution models for protein sequences

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class ModelProtein : public ModelMarkov
{
public:
	/**
		constructor
		@param model_name model name, e.g., JTT, WAG.
		@param freq state frequency type
		@param tree associated phylogenetic tree
	*/
    ModelProtein(const char *model_name, string model_params, StateFreqType freq, string freq_params, PhyloTree *tree, ModelsBlock *models_block);

	/**
		initialization, called automatically by the constructor, no need to call it
		@param model_name model name, e.g., JTT, WAG.
		@param freq state frequency type
	*/
	virtual void init(const char *model_name, string model_params, StateFreqType freq, string freq_params);

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
		read the rates from an input stream. it will throw error messages if failed
		@param in input stream
	*/
    virtual void readRates(istream &in) noexcept(false);

	/**
	 * @return model name with parameters in form of e.g. GTR{a,b,c,d,e,f}
	 */
    virtual string getNameParams(bool show_fixed_params = false);

    /** compute the tip likelihood vector of a state for Felsenstein's pruning algorithm
     @param state character state
     @param[out] state_lk state likehood vector of size num_states
     */
    virtual void computeTipLikelihood(PML::StateType state, double *state_lk);

    /**
     * Print the model information in a format that can be accepted by MrBayes, using lset and prset.<br>
     * By default, it simply prints a warning to the log and to the stream, stating that this model is not supported by MrBayes.
     * @param rate the rate information
     * @param out the ofstream to print the result to
     * @param partition the partition to apply lset and prset to
     * @param charset the current partition's charset. Useful for getting information from the checkpoint file
     * @param isSuperTree whether the tree is a super tree. Useful for retrieving information from the checkpoint file, which has different locations for PhyloTree and PhyloSuperTree
     * @param inclParams whether to include IQTree optimized parameters for the model
     */
    virtual void printMrBayesModelText(RateHeterogeneity* rate, ofstream& out, string partition, string charset, bool isSuperTree, bool inclParams);

private:
    ModelsBlock *models_block;

};

#endif
