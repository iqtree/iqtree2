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
#ifndef PARTITIONMODEL_H
#define PARTITIONMODEL_H

#include "tree/phylosupertree.h"
#include "modelfactory.h"

/**
Partition model

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class PartitionModel : public ModelFactory
{
public:
    PartitionModel();
	/**
		constructor
		create partition model with possible rate heterogeneity. Create proper class objects
		for two variables: model and site_rate. It takes the following field of params into account:
			model_name, num_rate_cats, freq_type, store_trans_matrix
		@param params program parameters
		@param tree associated phylogenetic super-tree
	*/
	PartitionModel(Params &params, PhyloSuperTree *tree, ModelsBlock *models_block);

    ~PartitionModel();

    /**
        set checkpoint object
        @param checkpoint
    */
    virtual void setCheckpoint(Checkpoint *checkpoint);

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
     * @param brlen_type either BRLEN_OPTIMIZE, BRLEN_FIX or BRLEN_SCALE
     * @return #parameters of the model + # branches
     */
    virtual int getNParameters(int brlen_type);

	/**
		optimize model parameters and tree branch lengths
        NOTE 2016-08-20: refactor the semantic of fixed_len
		@param fixed_len 0: optimize branch lengths, 1: fix branch lengths, 2: scale branch lengths
        @param write_info TRUE to write model parameters every optimization step, FALSE to only print at the end
        @param logl_epsilon log-likelihood epsilon to stop
        @param gradient_epsilon gradient (derivative) epsilon to stop
		@return the best likelihood 
	*/
	virtual double optimizeParameters(int fixed_len = BRLEN_OPTIMIZE, bool write_info = true,
                                      double logl_epsilon = 0.1, double gradient_epsilon = 0.0001);

	/**
	 *  optimize model parameters and tree branch lengths for the +I+G model
	 *  using restart strategy.
	 * 	@param fixed_len TRUE to fix branch lengths, default is false
	 *	@return the best likelihood
	 */
	virtual double optimizeParametersGammaInvar(int fixed_len = BRLEN_OPTIMIZE, bool write_info = true, double logl_epsilon = 0.1, double gradient_epsilon = 0.0001);

	/**
	 * @return TRUE if parameters are at the boundary that may cause numerical unstability
	 */
	virtual bool isUnstableParameters();

	/** optimize linked alpha parameter of over all partitions with Gamma rate */
	double optimizeLinkedAlpha(bool write_info, double gradient_epsilon);

	/**
		override function from Optimization class, used by the minimizeOneDimen() to optimize
		gamma shape parameter
	*/
	virtual double computeFunction(double shape);


protected:

	/** linked Gamma shape alpha between partitions */
	double linked_alpha;

};

#endif
