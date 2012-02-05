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
#include "partitionmodel.h"
#include "superalignment.h"

PartitionModel::PartitionModel()
 : ModelFactory()
{
}

PartitionModel::PartitionModel(Params &params, PhyloSuperTree *tree) 
 : ModelFactory(params, tree) 
{
	string model_name = params.model_name; 	
	PhyloSuperTree::iterator it;
	int part;
	for (it = tree->begin(), part = 0; it != tree->end(); it++, part++) {
		assert(!((*it)->model_factory));
		params.model_name = tree->part_info[part].model_name;
		(*it)->setModelFactory(new ModelFactory(params, (*it)));
		(*it)->setModel((*it)->getModelFactory()->model);
		(*it)->setRate((*it)->getModelFactory()->site_rate);
		params.model_name = model_name;
		string taxa_set = ((SuperAlignment*)tree->aln)->getPattern(part);
		(*it)->copyTree(tree, taxa_set);
		//(*it)->drawTree(cout);
	}
}

double PartitionModel::optimizeParameters(bool fixed_len, bool write_info) {
	PhyloSuperTree *tree = (PhyloSuperTree*)site_rate->getTree();
	double tree_lh = 0.0;
	int part = 0;
	for (PhyloSuperTree::iterator it = tree->begin(); it != tree->end(); it++, part++) {
		cout << "Optimizing " << (*it)->getModelName() << " parameters for partition " << tree->part_info[part].name << endl;
		tree_lh += (*it)->getModelFactory()->optimizeParameters(fixed_len, write_info);
	}
	//return ModelFactory::optimizeParameters(fixed_len, write_info);
	return tree_lh;
}

PartitionModel::~PartitionModel()
{
}


