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

ModelSet::ModelSet(const char *model_name, PhyloTree *tree) : ModelMarkov(tree)
{
	name = full_name = model_name;
	name += "+SSF";
	full_name += "+site-specific state-frequency model (unpublished)";
}

void ModelSet::computeTransMatrix(double time, double* trans_matrix, int mixture, int selected_row)
{
    // TODO not working with vectorization
    ASSERT(0);
	for (iterator it = begin(); it != end(); it++) {
		(*it)->computeTransMatrix(time, trans_matrix, mixture, selected_row);
		trans_matrix += (num_states * num_states);
	}
}

void ModelSet::computeTransDerv(double time, double* trans_matrix, double* trans_derv1, double* trans_derv2, int mixture)
{
    // TODO not working with vectorization
    ASSERT(0);
	for (iterator it = begin(); it != end(); it++) {
		(*it)->computeTransDerv(time, trans_matrix, trans_derv1, trans_derv2, mixture);
		trans_matrix += (num_states * num_states);
		trans_derv1 += (num_states * num_states);
		trans_derv2 += (num_states * num_states);
	}
}

int ModelSet::getPtnModelID(int ptn)
{
	ASSERT(ptn >= 0 && ptn < pattern_model_map.size());
	ASSERT(pattern_model_map[ptn] >= 0 && pattern_model_map[ptn] < size());
    return pattern_model_map[ptn];
}


double ModelSet::computeTrans(double time, int model_id, int state1, int state2) {
    if (phylo_tree->vector_size == 1) {
        return at(model_id)->computeTrans(time, state1, state2);
    }
	// temporary fix problem with vectorized eigenvectors
	int i;
    int vsize = phylo_tree->vector_size;
    int states_vsize = num_states*vsize;
    int model_vec_id = model_id % vsize;
    int start_ptn = model_id - model_vec_id;
    double *evec = &eigenvectors[start_ptn*num_states*num_states + model_vec_id + state1*num_states*vsize];
    double *inv_evec = &inv_eigenvectors[start_ptn*num_states*num_states + model_vec_id + state2*vsize];
    double *eval = &eigenvalues[start_ptn*num_states + model_vec_id];
	double trans_prob = 0.0;
	for (i = 0; i < states_vsize; i+=vsize) {
        double val = eval[i];
		double trans = evec[i] * inv_evec[i*num_states] * exp(time * val);
		trans_prob += trans;
	}
	return trans_prob;
}

double ModelSet::computeTrans(double time, int model_id, int state1, int state2, double &derv1, double &derv2) {
    if (phylo_tree->vector_size == 1) {
        return at(model_id)->computeTrans(time, state1, state2, derv1, derv2);
    }
	// temporary fix problem with vectorized eigenvectors
	int i;
    int vsize = phylo_tree->vector_size;
    int states_vsize = num_states*vsize;
    int model_vec_id = model_id % vsize;
    int start_ptn = model_id - model_vec_id;
    double *evec = &eigenvectors[start_ptn*num_states*num_states + model_vec_id + state1*num_states*vsize];
    double *inv_evec = &inv_eigenvectors[start_ptn*num_states*num_states + model_vec_id + state2*vsize];
    double *eval = &eigenvalues[start_ptn*num_states + model_vec_id];
	double trans_prob = 0.0;
	derv1 = derv2 = 0.0;
	for (i = 0; i < states_vsize; i+=vsize) {
        double val = eval[i];
		double trans = evec[i] * inv_evec[i*num_states] * exp(time * val);
		double trans2 = trans * val;
		trans_prob += trans;
		derv1 += trans2;
		derv2 += trans2 * val;
	}
	return trans_prob;
}

int ModelSet::getNDim()
{
	ASSERT(size());
    return front()->getNDim();
}

void ModelSet::writeInfo(ostream& out)
{
    if (empty()) {
        return;
    }
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
    if (empty()) {
        return;
    }
    for (iterator it = begin(); it != end(); it++) {
        (*it)->decomposeRateMatrix();
    }
	if (phylo_tree->vector_size == 1)
		return;
	// rearrange eigen to obey vector_size
	size_t vsize = phylo_tree->vector_size;
	size_t states2 = num_states*num_states;

    size_t max_size = get_safe_upper_limit(size());

    // copy dummy values
    for (size_t m = size(); m < max_size; m++) {
        memcpy(&eigenvalues[m*num_states], &eigenvalues[(m-1)*num_states], sizeof(double)*num_states);
        memcpy(&eigenvectors[m*states2], &eigenvectors[(m-1)*states2], sizeof(double)*states2);
        memcpy(&inv_eigenvectors[m*states2], &inv_eigenvectors[(m-1)*states2], sizeof(double)*states2);
        memcpy(&inv_eigenvectors_transposed[m*states2], &inv_eigenvectors_transposed[(m-1)*states2], sizeof(double)*states2);
    }

    double new_eval[num_states*vsize];
    double new_evec[states2*vsize];
    double new_inv_evec[states2*vsize];

    for (size_t ptn = 0; ptn < size(); ptn += vsize) {
        double *eval_ptr = &eigenvalues[ptn*num_states];
        double *evec_ptr = &eigenvectors[ptn*states2];
        double *inv_evec_ptr = &inv_eigenvectors[ptn*states2];
        for (size_t i = 0; i < vsize; i++) {
            for (size_t x = 0; x < num_states; x++)
                new_eval[x*vsize+i] = eval_ptr[x];
            for (size_t x = 0; x < states2; x++) {
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
        calculateSquareMatrixTranspose(new_inv_evec, num_states
                                       , &inv_eigenvectors_transposed[ptn*states2]);
    }
}

bool ModelSet::getVariables(double* variables)
{
	ASSERT(size());
    bool changed = false;
	for (iterator it = begin(); it != end(); it++)
		changed |= (*it)->getVariables(variables);
    return changed;
}

void ModelSet::setVariables(double* variables)
{
	ASSERT(size());
	front()->setVariables(variables);
}

ModelSet::~ModelSet()
{
    for (reverse_iterator rit = rbegin(); rit != rend(); rit++) {
        (*rit)->eigenvalues = nullptr;
        (*rit)->eigenvectors = nullptr;
        (*rit)->inv_eigenvectors = nullptr;
        (*rit)->inv_eigenvectors_transposed = nullptr;
        delete (*rit);
    }
}

void ModelSet::joinEigenMemory() {
    size_t nmixtures = get_safe_upper_limit(size());
    aligned_free(eigenvalues);
    aligned_free(eigenvectors);
    aligned_free(inv_eigenvectors);
    aligned_free(inv_eigenvectors_transposed);
    
    size_t states2 = num_states*num_states;
    
    eigenvalues = aligned_alloc<double>(num_states*nmixtures);
    eigenvectors = aligned_alloc<double>(states2*nmixtures);
    inv_eigenvectors = aligned_alloc<double>(states2*nmixtures);
    inv_eigenvectors_transposed = aligned_alloc<double>(states2*nmixtures);
    
    // assigning memory for individual models
    size_t m = 0;
    for (iterator it = begin(); it != end(); it++, m++) {
        // first copy memory for eigen stuffs
        memcpy(&eigenvalues[m*num_states], (*it)->eigenvalues, num_states*sizeof(double));
        memcpy(&eigenvectors[m*states2], (*it)->eigenvectors, states2*sizeof(double));
        memcpy(&inv_eigenvectors[m*states2], (*it)->inv_eigenvectors, states2*sizeof(double));
        memcpy(&inv_eigenvectors_transposed[m*states2], (*it)->inv_eigenvectors_transposed, states2*sizeof(double));
        // then delete
        aligned_free((*it)->eigenvalues);
        aligned_free((*it)->eigenvectors);
        aligned_free((*it)->inv_eigenvectors);
        aligned_free((*it)->inv_eigenvectors_transposed);
        
        // and assign new memory
        (*it)->eigenvalues = &eigenvalues[m*num_states];
        (*it)->eigenvectors = &eigenvectors[m*states2];
        (*it)->inv_eigenvectors = &inv_eigenvectors[m*states2];
        (*it)->inv_eigenvectors_transposed = &inv_eigenvectors_transposed[m*states2];
    }
    
    // copy dummy values
    for (m = size(); m < nmixtures; m++) {
        memcpy(&eigenvalues[m*num_states], &eigenvalues[(m-1)*num_states], sizeof(double)*num_states);
        memcpy(&eigenvectors[m*states2], &eigenvectors[(m-1)*states2], sizeof(double)*states2);
        memcpy(&inv_eigenvectors[m*states2], &inv_eigenvectors[(m-1)*states2], sizeof(double)*states2);
        memcpy(&inv_eigenvectors_transposed[m*states2], &inv_eigenvectors_transposed[(m-1)*states2], sizeof(double)*states2);
    }
}
