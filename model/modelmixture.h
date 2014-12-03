/*
 * modelmixture.h
 *
 *  Created on: Nov 29, 2014
 *      Author: minh
 */

#ifndef MODELMIXTURE_H_
#define MODELMIXTURE_H_

#include "phylotree.h"
#include "modelgtr.h"


/**
 * mixture model
 */
class ModelMixture: public ModelGTR, vector<ModelGTR*> {
public:
	ModelMixture(PhyloTree *tree);
	virtual ~ModelMixture();

	/**
	 * @return the number of mixture model components
	 */
	virtual int getMixtureComponents() {return size(); }
};

#endif /* MODELMIXTURE_H_ */
