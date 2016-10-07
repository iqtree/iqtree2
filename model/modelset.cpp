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

ModelSet::ModelSet(const char *model_name, PhyloTree *tree) : ModelGTR(tree)
{
	name = full_name = model_name;
	name += "+SSF";
	full_name += "+site-specific state-frequency model (unpublished)";
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

int ModelSet::getPtnModelID(int ptn)
{
	assert(ptn >= 0 && ptn < pattern_model_map.size());
	assert(pattern_model_map[ptn] >= 0 && pattern_model_map[ptn] < size());
    return pattern_model_map[ptn];
}


double ModelSet::computeTrans(double time, int model_id, int state1, int state2) {
	return at(model_id)->computeTrans(time, state1, state2);
}

double ModelSet::computeTrans(double time, int model_id, int state1, int state2, double &derv1, double &derv2) {
	return at(model_id)->computeTrans(time, state1, state2, derv1, derv2);
	
}

int ModelSet::getNDim()
{
	assert(size());
    return front()->getNDim();
}

void ModelSet::writeInfo(ostream& out)
{
    if (empty())
        return;
	if (verbose_mode >= VB_DEBUG) {
		int i = 1;
		for (iterator it = begin(); it != end(); it++, i++) {
			out << "Partition " << i << ":" << endl;
			(*it)->writeInfo(out);
		}
	} else {
		front()->writeInfo(out);
	}
}

void ModelSet::decomposeRateMatrix()
{
    if (empty())
        return;
	for (iterator it = begin(); it != end(); it++)
		(*it)->decomposeRateMatrix();
	if (phylo_tree->vector_size == 1)
		return;
	// rearrange eigen to obey vector_size
	size_t vsize = phylo_tree->vector_size;
	size_t states2 = num_states*num_states;
	size_t ptn, i, x;
    double new_eval[num_states*vsize];
    double new_evec[states2*vsize];
    double new_inv_evec[states2*vsize];

	for (ptn = 0; ptn < size(); ptn += vsize) {
		double *eval_ptr = &eigenvalues[ptn*num_states];
		double *evec_ptr = &eigenvectors[ptn*states2];
		double *inv_evec_ptr = &inv_eigenvectors[ptn*states2];
		for (i = 0; i < vsize; i++) {
			for (x = 0; x < num_states; x++)
				new_eval[x*vsize+i] = eval_ptr[x];
			for (x = 0; x < states2; x++) {
				new_evec[x*vsize+i] = evec_ptr[x];
				new_inv_evec[x*vsize+i] = inv_evec_ptr[x];
			}
			eval_ptr += num_states;
			evec_ptr += states2;
			inv_evec_ptr += states2;
		}
		// copy new values
        memcpy(&eigenvalues[ptn*num_states], new_eval, sizeof(double)*num_states*vsize);
        memcpy(&eigenvectors[ptn*states2], new_evec, sizeof(double)*states2*vsize);
        memcpy(&inv_eigenvectors[ptn*states2], new_inv_evec, sizeof(double)*states2*vsize);
	}
}


bool ModelSet::getVariables(double* variables)
{
	assert(size());
    bool changed = false;
	for (iterator it = begin(); it != end(); it++)
		changed |= (*it)->getVariables(variables);
    return changed;
}

void ModelSet::setVariables(double* variables)
{
	assert(size());
	front()->setVariables(variables);
}


ModelSet::~ModelSet()
{
	for (reverse_iterator rit = rbegin(); rit != rend(); rit++) {
		(*rit)->eigenvalues = NULL;
		(*rit)->eigenvectors = NULL;
		(*rit)->inv_eigenvectors = NULL;
		delete (*rit);
	}
}

void ModelSet::joinEigenMemory() {
    size_t nmixtures = get_safe_upper_limit(size());
	if (eigenvalues) aligned_free(eigenvalues);
	if (eigenvectors) aligned_free(eigenvectors);
	if (inv_eigenvectors) aligned_free(inv_eigenvectors);

    size_t states2 = num_states*num_states;

	eigenvalues = aligned_alloc<double>(num_states*nmixtures);
	eigenvectors = aligned_alloc<double>(states2*nmixtures);
	inv_eigenvectors = aligned_alloc<double>(states2*nmixtures);

	// assigning memory for individual models
	size_t m = 0;
	for (iterator it = begin(); it != end(); it++, m++) {
        // first copy memory for eigen stuffs
        memcpy(&eigenvalues[m*num_states], (*it)->eigenvalues, num_states*sizeof(double));
        memcpy(&eigenvectors[m*states2], (*it)->eigenvectors, states2*sizeof(double));
        memcpy(&inv_eigenvectors[m*states2], (*it)->inv_eigenvectors, states2*sizeof(double));
        // then delete
		if ((*it)->eigenvalues) aligned_free((*it)->eigenvalues);
		if ((*it)->eigenvectors) aligned_free((*it)->eigenvectors);
		if ((*it)->inv_eigenvectors) aligned_free((*it)->inv_eigenvectors);
//		if ((*it)->eigen_coeff) aligned_free((*it)->eigen_coeff);

        // and assign new memory
		(*it)->eigenvalues = &eigenvalues[m*num_states];
		(*it)->eigenvectors = &eigenvectors[m*states2];
		(*it)->inv_eigenvectors = &inv_eigenvectors[m*states2];
	}

    // copy dummy values
    for (m = size(); m < nmixtures; m++) {
        memcpy(&eigenvalues[m*num_states], &eigenvalues[(m-1)*num_states], sizeof(double)*num_states);
        memcpy(&eigenvectors[m*states2], &eigenvectors[(m-1)*states2], sizeof(double)*states2);
        memcpy(&inv_eigenvectors[m*states2], &inv_eigenvectors[(m-1)*states2], sizeof(double)*states2);
    }
}
