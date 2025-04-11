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
#ifndef MODELBIN_H
#define MODELBIN_H

#include "modelmarkov.h"

/**
Model for Binary data

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class ModelBIN : public ModelMarkov
{
public:
	/**
		constructor
		@param model_name model name, e.g., JC, HKY.
		@param freq state frequency type
		@param tree associated phylogenetic tree
	*/
    ModelBIN(const char *model_name, string model_params, StateFreqType freq, string freq_params, PhyloTree *tree);

	/**
		initialization, called automatically by the constructor, no need to call it
		@param model_name model name, e.g., JC, HKY.
		@param freq state frequency type
	*/
	virtual void init(const char *model_name, string model_params, StateFreqType freq, string freq_params);

	/**
	 * @return model name with parameters in form of e.g. GTR{a,b,c,d,e,f}
	 */
    virtual string getNameParams(bool show_fixed_params = false);

    /**
        start structure for checkpointing
    */
    virtual void startCheckpoint();

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
    virtual void printMrBayesModelText(ofstream& out, string partition, string charset, bool isSuperTree, bool inclParams);

};

#endif
