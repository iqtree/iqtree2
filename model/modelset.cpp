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


ModelSet::ModelSet(const char* model_name, PhyloTree* tree) : ModelMarkov(tree)
{
	name = full_name = model_name;
	name += "+SSF";
	full_name += "+site-specific state-frequency model (unpublished)";
}

void ModelSet::computeTransMatrix(double time, double* trans_matrix, 
                                  int mixture) const
{
    // TODO not working with vectorization
    ASSERT(0);
	for (auto it = models.begin(); it != models.end(); it++) {
		(*it)->computeTransMatrix(time, trans_matrix, mixture);
		trans_matrix += (num_states * num_states);
	}
}

void ModelSet::computeTransDerv(double time, double* trans_matrix,
                                double* trans_derv1, double* trans_derv2,
                                int mixture) const
{
    // TODO not working with vectorization
    ASSERT(0);
	for (auto it = models.begin(); it != models.end(); it++) {
		(*it)->computeTransDerv(time, trans_matrix, trans_derv1,
                                trans_derv2, mixture);
		trans_matrix += (num_states * num_states);
		trans_derv1 += (num_states * num_states);
		trans_derv2 += (num_states * num_states);
	}
}

int ModelSet::getPtnModelID(int ptn) const
{
	ASSERT(ptn >= 0 && ptn < pattern_model_map.size());
	ASSERT(pattern_model_map[ptn] >= 0 && pattern_model_map[ptn] < models.size());
    return pattern_model_map[ptn];
}


double ModelSet::computeTrans(double time, int model_id,
                              int state1, int state2) const {
    if (phylo_tree->getVectorSize() == 1) {
        return models.at(model_id)->computeTrans(time, state1, state2);
    }
	// temporary fix problem with vectorized eigenvectors
	int i;
    int vsize        = phylo_tree->getVectorSizeAsInt();
    int states_vsize = num_states*vsize;
    int model_vec_id = model_id % vsize;
    int start_ptn    = model_id - model_vec_id;
    int square       = num_states * num_states;
    double *evec     = &eigenvectors[start_ptn*square + model_vec_id + state1*states_vsize];
    double *inv_evec = &inv_eigenvectors[start_ptn*square + model_vec_id + state2*vsize];
    double *eval     = &eigenvalues[start_ptn*num_states + model_vec_id];
	double trans_prob = 0.0;
	for (i = 0; i < states_vsize; i+=vsize) {
        double val   = eval[i];
		double trans = evec[i] * inv_evec[i*num_states] * exp(time * val);
		trans_prob  += trans;
	}
	return trans_prob;
}

double ModelSet::computeTrans(double time, int model_id, 
                              int state1, int state2,
                              double &derv1, double &derv2) const {
    if (phylo_tree->getVectorSize() == 1) {
        return models.at(model_id)->computeTrans(time, state1, state2, derv1, derv2);
    }
	// temporary fix problem with vectorized eigenvectors
    int     vsize        = phylo_tree->getVectorSizeAsInt();
    int     states_vsize = num_states*vsize;
    int     model_vec_id = model_id % vsize;
    int     start_ptn    = model_id - model_vec_id;
    int     square       = num_states * num_states;
    double* evec         = &eigenvectors[start_ptn*square + model_vec_id + state1*states_vsize];
    double* inv_evec     = &inv_eigenvectors[start_ptn*square + model_vec_id + state2*vsize];
    double* eval         = &eigenvalues[start_ptn*num_states + model_vec_id];
	double  trans_prob   = 0.0;
	derv1 = derv2 = 0.0;
	for (int i = 0; i < states_vsize; i+=vsize) {
        double val    = eval[i];
		double trans  = evec[i] * inv_evec[i*num_states] * exp(time * val);
		double trans2 = trans * val;
		trans_prob   += trans;
		derv1        += trans2;
		derv2        += trans2 * val;
	}
	return trans_prob;
}

int ModelSet::getNDim() const
{
	ASSERT(0<models.size());
    return models.front()->getNDim();
}

void ModelSet::writeInfo(ostream& out)
{
    if (models.empty()) {
        return;
    }
	if (verbose_mode >= VerboseMode::VB_DEBUG) {
		int i = 1;
		for (auto it = models.begin(); it != models.end(); it++, i++) {
			out << "Partition " << i << ":" << endl;
			(*it)->writeInfo(out);
		}
	} else {
		models.front()->writeInfo(out);
	}
}

void ModelSet::decomposeRateMatrix()
{
    if (models.empty()) {
        return;
    }
    for (auto it = models.begin(); it != models.end(); it++) {
        (*it)->decomposeRateMatrix();
    }
    if (phylo_tree->getVectorSize() == 1) {
        return;
    }
	// rearrange eigen to obey vector_size
	int  vsize    = phylo_tree->getVectorSizeAsInt();
	int  square   = num_states*num_states;
    int  max_size = static_cast<int>(get_safe_upper_limit(models.size()));

    // copy dummy values
    for (size_t m = models.size(); m < max_size; m++) {
        memcpy(&eigenvalues[m*num_states], &eigenvalues[(m-1)*num_states], sizeof(double)*num_states);
        memcpy(&eigenvectors[m*square], &eigenvectors[(m-1)*square], sizeof(double)*square);
        memcpy(&inv_eigenvectors[m*square], &inv_eigenvectors[(m-1)*square], sizeof(double)*square);
        memcpy(&inv_eigenvectors_transposed[m*square], &inv_eigenvectors_transposed[(m-1)*square], sizeof(double)*square);
    }

    std::vector<double> eval_vector(num_states*vsize);
    std::vector<double> evec_vector(square*vsize);
    std::vector<double> inv_vector(square*vsize);

    double* new_eval     = eval_vector.data();
    double* new_evec     = evec_vector.data();
    double* new_inv_evec = inv_vector.data();

    intptr_t num_models = static_cast<intptr_t>(models.size());
    for (intptr_t ptn = 0; ptn < num_models; ptn += vsize) {
        double* eval_ptr     = &eigenvalues[ptn*num_states];
        double* evec_ptr     = &eigenvectors[ptn*square];
        double* inv_evec_ptr = &inv_eigenvectors[ptn*square];
        for (int i = 0; i < vsize; i++) {
            for (int x = 0; x < num_states; x++) {
                new_eval[x * vsize + i] = eval_ptr[x];
            }
            for (int x = 0; x < square; x++) {
                new_evec[x*vsize+i] = evec_ptr[x];
                new_inv_evec[x*vsize+i] = inv_evec_ptr[x];
            }
            eval_ptr     += num_states;
            evec_ptr     += square;
            inv_evec_ptr += square;
        }
        // copy new values
        memcpy(&eigenvalues[ptn*num_states],  &new_eval[0],     sizeof(double)*num_states*vsize);
        memcpy(&eigenvectors[ptn*square],     &new_evec[0],     sizeof(double)*square*vsize);
        memcpy(&inv_eigenvectors[ptn*square], &new_inv_evec[0], sizeof(double)*square*vsize);
        calculateSquareMatrixTranspose(&new_inv_evec[0], num_states
                                       , &inv_eigenvectors_transposed[ptn*square]);
    }
}

bool ModelSet::getVariables(const double* variables)
{
	ASSERT(models.size());
    bool changed = false;
	for (auto it = models.begin(); it != models.end(); ++it) {
		changed |= (*it)->getVariables(variables);
    }
    return changed;
}

void ModelSet::logVariablesTo(std::stringstream& var_list) const {
    var_list << " set of " << models.size() << "models";
    for (auto model : models) {
        var_list << "\n  ";
        model->logVariablesTo(var_list);
    }
}


void ModelSet::setVariables(double* variables)
{
	ASSERT(models.size());
	models.front()->setVariables(variables);
}

ModelSet::~ModelSet()
{
    for (auto rit = models.rbegin(); rit != models.rend(); ++rit) {
        (*rit)->eigenvalues = nullptr;
        (*rit)->eigenvectors = nullptr;
        (*rit)->inv_eigenvectors = nullptr;
        (*rit)->inv_eigenvectors_transposed = nullptr;
        delete (*rit);
        (*rit) = nullptr;
    }
    models.clear();
}

void ModelSet::joinEigenMemory() {
    size_t nmixtures = get_safe_upper_limit(models.size());
    aligned_free(eigenvalues);
    aligned_free(eigenvectors);
    aligned_free(inv_eigenvectors);
    aligned_free(inv_eigenvectors_transposed);
    
    size_t square = num_states*num_states;
    
    eigenvalues      = aligned_alloc<double>(num_states*nmixtures);
    eigenvectors     = aligned_alloc<double>(square*nmixtures);
    inv_eigenvectors = aligned_alloc<double>(square*nmixtures);
    inv_eigenvectors_transposed = aligned_alloc<double>(square*nmixtures);
    
    // assigning memory for individual models
    size_t m = 0;
    for (auto it = models.begin(); it != models.end(); it++, m++) {
        // first copy memory for eigen stuffs
        memcpy(&eigenvalues[m*num_states],  (*it)->eigenvalues,      num_states*sizeof(double));
        memcpy(&eigenvectors[m*square],     (*it)->eigenvectors,     square*sizeof(double));
        memcpy(&inv_eigenvectors[m*square], (*it)->inv_eigenvectors, square*sizeof(double));
        memcpy(&inv_eigenvectors_transposed[m*square],
               (*it)->inv_eigenvectors_transposed, square*sizeof(double));
        // then delete
        aligned_free((*it)->eigenvalues);
        aligned_free((*it)->eigenvectors);
        aligned_free((*it)->inv_eigenvectors);
        aligned_free((*it)->inv_eigenvectors_transposed);
        
        // and assign new memory
        (*it)->eigenvalues      = &eigenvalues[m*num_states];
        (*it)->eigenvectors     = &eigenvectors[m*square];
        (*it)->inv_eigenvectors = &inv_eigenvectors[m*square];
        (*it)->inv_eigenvectors_transposed = &inv_eigenvectors_transposed[m*square];
    }
    
    // copy dummy values
    for (m = models.size(); m < nmixtures; m++) {
        memcpy(&eigenvalues[m*num_states], &eigenvalues[(m-1)*num_states], sizeof(double)*num_states);
        memcpy(&eigenvectors[m*square], &eigenvectors[(m-1)*square], sizeof(double)*square);
        memcpy(&inv_eigenvectors[m*square], &inv_eigenvectors[(m-1)*square], sizeof(double)*square);
        memcpy(&inv_eigenvectors_transposed[m*square], &inv_eigenvectors_transposed[(m-1)*square], sizeof(double)*square);
    }
}

ModelMarkov* ModelSet::firstModel() const {
    return models.front();
}

void ModelSet::addModel(ModelMarkov* model_to_add) {

}

void ModelSet::getDivergentModels
        (DivergentModels& div_models) {
    for (auto model: models) {
        model->getDivergentModels(div_models);
    }
}

void ModelSet::setPatternInvar
        (double* ptn_invar, bool take_ownership) {
    super::setPatternInvar(ptn_invar, take_ownership);
    for (auto model: models) {
        model->setPatternInvar(ptn_invar, false);
    }
}

