/*
 * searchinfo.cpp
 *
 *  Created on: Nov 18, 2014
 *      Author: tung
 */

#include "searchinfo.h"

void SearchInfo::init(bool nni5, double initPS, int maxNNISteps, bool reduction, bool speedNNI) {
	this->nni5 = nni5;
	this->initPS = initPS;
	this->maxNNISteps = maxNNISteps;
	this->reduction = reduction;
	this->speedNNI = speedNNI;
	this->curPS = initPS;
}

SearchInfo::SearchInfo() {
	this->nni5 = true;
	this->initPS = 0.5;
	this->maxNNISteps = 10000;
	this->reduction = false;
	this->speedNNI = false;
	this->curPS = this->initPS;
	this->nniOptimal = false;
	this->numDup = 0;
}

SearchInfo::~SearchInfo() {
}

double SearchInfo::getCurPs() const {
	return curPS;
}

void SearchInfo::setCurPs(double curPs) {
	curPS = curPs;
}

double SearchInfo::getInitPs() const {
	return initPS;
}

void SearchInfo::setInitPs(double initPs) {
	initPS = initPs;
}

int SearchInfo::getMaxNniSteps() const {
	return maxNNISteps;
}

void SearchInfo::setMaxNniSteps(int maxNniSteps) {
	maxNNISteps = maxNniSteps;
}

bool SearchInfo::isNni5() const {
	return nni5;
}

void SearchInfo::setNni5(bool nni5) {
	this->nni5 = nni5;
}

bool SearchInfo::isReduction() const {
	return reduction;
}

void SearchInfo::setReduction(bool reduction) {
	this->reduction = reduction;
}

bool SearchInfo::isSpeedNni() const {
	return speedNNI;
}

void SearchInfo::setSpeedNni(bool speedNni) {
	speedNNI = speedNni;
}
