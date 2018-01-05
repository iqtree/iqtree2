/*
 * modelsblock.h
 *
 *  Created on: Jan 9, 2015
 *      Author: minh
 */

#ifndef MODELSBLOCK_H_
#define MODELSBLOCK_H_

#include "ncl/ncl.h"
#include "utils/tools.h"

const int NM_ATOMIC  = 1; // NxsModel is not mixture or +G etc. model
const int NM_FREQ    = 2;   // NxsModel contains state frequency
const int NM_PROTEIN = 4;   // NxsModel contains emprical protein model

class NxsModel {
public:
	/* model name */
	string name;

	/* model description */
	string description;

	/* flag as NM_ATOMIC or NM_FREQ or both */
	int flag;

    NxsModel() {
        flag = 0;
    }

    NxsModel(string name) {
        this->name = name;
        flag = 0;
    }

	virtual ~NxsModel() {}
};

/**
 * Class to parse MODELS block in NEXUS file
 */
class ModelsBlock: public NxsBlock, public unordered_map<string, NxsModel> {
public:
	/** constructor */
	ModelsBlock();
	/** destructor */
	virtual ~ModelsBlock();

    /**
        @param name model name
        @return pointer to model with the name or NULL if not found
    */
	NxsModel *findModel(string name);

    /**
        @param name model name
        @return pointer to a mixed model with the name or NULL if not found
    */
	NxsModel *findMixModel(string name);


protected:

	/**
		main method to read block from file
		@param token a token reader
	*/
	virtual void Read(NxsToken &token);

};

#endif /* MODELSBLOCK_H_ */
