/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2012  BUI Quang Minh <email>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include "modelset.h"

ModelSet::ModelSet(PhyloTree *tree) : GTRModel(tree)
{

}

void ModelSet::computeTransMatrix(double time, double* trans_matrix)
{
	for (iterator it = begin(); it != end(); it++) {
		(*it)->computeTransMatrix(time, trans_matrix);
		trans_matrix += (num_states * num_states);
	}
}

void ModelSet::computeTransMatrixFreq(double time, double* trans_matrix)
{
	for (iterator it = begin(); it != end(); it++) {
		(*it)->computeTransMatrixFreq(time, trans_matrix);
		trans_matrix += (num_states * num_states);
	}
}

void ModelSet::computeTransDerv(double time, double* trans_matrix, double* trans_derv1, double* trans_derv2)
{
	for (iterator it = begin(); it != end(); it++) {
		(*it)->computeTransDerv(time, trans_matrix, trans_derv1, trans_derv2);
		trans_matrix += (num_states * num_states);
		trans_derv1 += (num_states * num_states);
		trans_derv2 += (num_states * num_states);
	}
}

void ModelSet::computeTransDervFreq(double time, double rate_val, double* trans_matrix, double* trans_derv1, double* trans_derv2)
{
	for (iterator it = begin(); it != end(); it++) {
		(*it)->computeTransDervFreq(time, rate_val, trans_matrix, trans_derv1, trans_derv2);
		trans_matrix += (num_states * num_states);
		trans_derv1 += (num_states * num_states);
		trans_derv2 += (num_states * num_states);
	}
}

int ModelSet::getNDim()
{
	assert(size());
    return front()->getNDim();
}

void ModelSet::writeInfo(ostream& out)
{
	assert(size());
	if (verbose_mode >= VB_MAX) {
		for (iterator it = begin(); it != end(); it++)
			(*it)->writeInfo(out);
	} else {
		front()->writeInfo(out);
	}
}


void ModelSet::getVariables(double* variables)
{
	assert(size());
	front()->getVariables(variables);
}

void ModelSet::setVariables(double* variables)
{
	assert(size());
	for (iterator it = begin(); it != end(); it++) {
		(*it)->setVariables(variables);
	}
}


ModelSet::~ModelSet()
{

}

