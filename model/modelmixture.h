/*
 * modelmixture.h
 *
 *  Created on: Nov 29, 2014
 *      Author: minh
 */

#ifndef MODELMIXTURE_H_
#define MODELMIXTURE_H_

#include "phylotree.h"
#include "modelsubst.h"


/**
 * mixture model
 */
class ModelMixture: public ModelSubst, vector<ModelSubst*> {
public:
	ModelMixture(PhyloTree *tree);
	virtual ~ModelMixture();
};

#endif /* MODELMIXTURE_H_ */
