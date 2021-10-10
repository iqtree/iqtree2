#include "modeldivergent.h"
#include "variablebounds.h" //for VariableBounds class
#include "modelmarkov.h"    //for RATE_TOL
#include <tree/phylotree.h> //for PhyloTree::vector_size

ModelDivergent::ModelDivergent(): super(), optimize_steps(0) {
    //If it's still zero when optimizeParameters() is called
}

ModelDivergent::ModelDivergent(PhyloTree* tree, 
                               PhyloTree* report_to_tree): 
    super(tree, report_to_tree),
    optimize_steps(0) {
    //If it's still zero when optimizeParameters() is called
}


ModelDivergent::~ModelDivergent() {
    for (ModelMarkov* subtree_model : subtree_models) {
        delete subtree_model;
    }
    subtree_models.clear();
}

bool ModelDivergent::isDivergentModel() const {
    return true;
}

void ModelDivergent::setTree(PhyloTree *tree) {
    for (ModelMarkov* subtree_model : subtree_models) {
        subtree_model->setTree(tree);
    }
    phylo_tree = tree;
}

void ModelDivergent::checkModelReversibility() {
    // check reversibility
    int subtree_count = static_cast<int>(subtree_models.size());
    ASSERT(0<subtree_count);
    bool rev      = subtree_models.front()->isReversible();
    bool err      = false;
    for (int i = 1; i < subtree_count; ++i) {
        if (subtree_models[i]->isReversible() != rev) {
            cerr << "ERROR: Model " << subtree_models[i]->name 
                 << " has different reversible property" << endl;
            err = true;
        }
    }
    if (err) {
        outError("Model reversibility is not consistent");
    }
    if (rev != isReversible()) {
        setReversible(rev);
    }
}

void ModelDivergent::setNumberOfVariableRates(int param_count) {
    num_params = param_count;
}

void ModelDivergent::setNumberOfStates(int states) {
    for (ModelMarkov* subtree_model : subtree_models) {
        subtree_model->setNumberOfStates(states);
    }
}

bool ModelDivergent::getSpecifiedAscertainmentBiasCorrection
        (ASCType& asc_type) {
    ASSERT(!subtree_models.empty());
    return subtree_models[0]->getSpecifiedAscertainmentBiasCorrection(asc_type);
}

int  ModelDivergent::getNDim()     const {
    int parameters = static_cast<int>(own_parameters.size());
    for (ModelMarkov* subtree_model : subtree_models) {
        parameters += subtree_model->getNDim();
    }
    return parameters;
}

int  ModelDivergent::getNDimFreq() const {
    int frequency_parameters = 0; 
    //Note: Can't have its own frequency parameters... yet.
    for (ModelMarkov* subtree_model : subtree_models) {
        frequency_parameters += subtree_model->getNDimFreq();
    }
    return frequency_parameters;
}

bool ModelDivergent::isReversible() const {
    ASSERT(!subtree_models.empty());
    return subtree_models[0]->isReversible();
}

int  ModelDivergent::getNumRateEntries() const {
    ASSERT(!subtree_models.empty());
    return subtree_models[0]->getNumRateEntries();
}

int  ModelDivergent::getNumberOfRates() const {
    int num_rates = 0;
    for (ModelMarkov* subtree_model : subtree_models) {
        num_rates += subtree_model->getNumberOfRates();
    }
}

int  ModelDivergent::getTransMatrixSize() const  { 
    ASSERT(!subtree_models.empty());
    return num_states * num_states 
            * static_cast<int>(subtree_models.size()); 
}

double ModelDivergent::computeTrans
        (double time, int model_id, 
         int state1, int state2) const {
    ASSERT(0<=model_id && model_id<subtree_models.size());
    return subtree_models[model_id]->computeTrans
            (time, state1, state2);
}

double ModelDivergent::computeTrans
        (double time, int model_id, 
         int state1, int state2, 
         double &derv1, double &derv2) const {
    ASSERT(0<=model_id && model_id<subtree_models.size());
    return subtree_models[model_id]->computeTrans
            (time, state1, state2, derv1, derv2);
}

void ModelDivergent::getStateFrequency(double *state_freq, 
                                       int model_id) const {
    ASSERT(0<=model_id && model_id<subtree_models.size());
    subtree_models[model_id]->getStateFrequency(state_freq, 0);
}

void ModelDivergent::setStateFrequency(double *state_freq) {
    //Todo: But... wouldn't it be better for this to ASSERT(0)
    for (ModelMarkov* subtree_model : subtree_models) {
        subtree_model->setStateFrequency(state_freq);
    }
}

bool ModelDivergent::scaleStateFreq() {
    return true;
}

void ModelDivergent::setRateMatrix(double* rate_mat) {
    ASSERT(0 && "Cannot set RateMatrix of ModelDivergent");
}

void ModelDivergent::setRateMatrixFromModel() {
     for (ModelMarkov* subtree_model : subtree_models) {
        subtree_model->setRateMatrixFromModel();
    }   
}

void ModelDivergent::decomposeRateMatrix() {
    if (subtree_models.empty()) {
        return;
    }
    for (ModelMarkov* subtree_model : subtree_models) {
        subtree_model->decomposeRateMatrix();
    }
    ASSERT(phylo_tree != nullptr);
    if (phylo_tree->vector_size == 1) {
        return;
    }
    ASSERT(0 && "ModelDivergent::decomposeRateMatrix"
           "doesn't work with vectorized kernels, as yet");
    //Todo: Can probably adapt a lot of what is done in
    //ModelSet's version of decompareRateMatrix
}

void ModelDivergent::setOptimizeSteps(int steps) {
    optimize_steps = steps;
}

double ModelDivergent::optimizeParameters
        (double gradient_epsilon, PhyloTree* report_to_tree) {
    if (phylo_tree!=nullptr) {
        if (!phylo_tree->getModelFactory()->unobserved_ptns.empty()) {
            outError("Divergent model +ASC is not supported yet." 
                    " Contact author if needed.");
        }
    }
	int ndim = getNDim();	
    if (ndim == 0) {
        //nothing to be optimized!  Out we go!
        return 0.0;
    }    
    TREE_LOG_LINE(*report_to_tree, VerboseMode::VB_MAX, 
                  "Optimizing " << name << " model parameters...");
    VariableBounds vb(ndim+1);
	double score;

	setVariables(vb.variables);
	setBounds(vb.lower_bound, vb.upper_bound, vb.bound_check);
    //todo: optimization algorithm as a member variable instad.
    if (phylo_tree->params->optimize_alg_freerate.find("BFGS-B") == string::npos) {
        score = -minimizeMultiDimen(vb.variables, ndim, 
                                    vb.lower_bound, vb.upper_bound, 
                                    vb.bound_check, max(gradient_epsilon, 
                                    TOL_RATE));
    }
    else {
        score = -L_BFGS_B(ndim, vb.variables + 1, 
                          vb.lower_bound + 1, vb.upper_bound + 1, 
                          max(gradient_epsilon, TOL_RATE));
    }
	bool changed = getVariables(vb.variables);
    if (changed) {
        decomposeRateMatrix();
        phylo_tree->clearAllPartialLH();
        score = phylo_tree->computeLikelihood();
    }
	return score;
}

void ModelDivergent::setBounds(double* lower_bound, double* upper_bound, 
                               bool*   bound_check) {
    int i = 1; //Numbering starts at 1
    for (ModelVariable& p : own_parameters) {
        lower_bound[i] = p.getRange().first;
        upper_bound[i] = p.getRange().second;
        bound_check[i] = p.isFixed();
        ++i;
    }
    int offset = static_cast<int>(own_parameters.size());
    for (ModelMarkov* subtree_model : subtree_models) {
        subtree_model->setBounds(lower_bound+offset, 
                                 upper_bound+offset,
                                 bound_check+offset);
        offset += subtree_model->getNDim();
    }
}

void ModelDivergent::setVariables(double *variables) {
    int offset = 0; 
    for (ModelVariable& p : own_parameters) {
        variables[offset+1] = p.getValue();
        ++offset;
    }
    for (ModelMarkov* subtree_model : subtree_models) {
        subtree_model->setVariables(variables+offset);
        offset+= subtree_model->getNDim();
    }
}

bool ModelDivergent::getVariables(const double *variables) {
    bool changed = false;
    int  offset  = 0;

    for (ModelVariable& p : own_parameters) {
        double newValue = variables[offset+1];
        double delta    = fabs(p.getValue() - newValue);
        if (delta >= TOL_RATE) {
            changed = true;
        }
        p.setValue(newValue);
        ++offset;
    }
    for (ModelMarkov* subtree_model : subtree_models) {
        changed |= subtree_model->getVariables(variables+offset);
        offset  += subtree_model->getNDim();
    }
    return changed;
}


void ModelDivergent::afterVariablesChanged() {
    for (ModelMarkov* subtree_model : subtree_models) {
        subtree_model->afterVariablesChanged();
    }
}

bool ModelDivergent::isUnstableParameters() const {
    bool unstable = false;
    for (ModelMarkov* subtree_model : subtree_models) {
        if (subtree_model->isUnstableParameters()) {
            unstable = true;
        }
    }
    return unstable;
}

void ModelDivergent::writeInfo(ostream &out) {
    ASSERT(0 && "ModelDivergent:writeInfo not implemented");
}

void ModelDivergent::report(ostream &out) {
    ASSERT(0 && "ModelDivergent::report not implemented");
}

uint64_t ModelDivergent::getMemoryRequired() const {
    uint64_t mem_needed = own_parameters.size() * sizeof(double);
    for (ModelMarkov* subtree_model : subtree_models) {
        mem_needed += subtree_model->getMemoryRequired();
    }
    return mem_needed;
}

void ModelDivergent::startCheckpoint() {
    auto model_count = own_parameters.size();
    auto nmod_string = convertIntToString(model_count);
    checkpoint->startStruct("ModelDivergent" + nmod_string);
}

void ModelDivergent::saveCheckpoint() {
    startCheckpoint();
    int model_num = 1;
    for (ModelMarkov* subtree_model : subtree_models) {
        auto model_num_string = convertIntToString(model_num);
        checkpoint->startStruct("Model" + model_num_string);
        subtree_model->saveCheckpoint();
        checkpoint->endStruct();
        ++model_num;
    }
    endCheckpoint();
}

void ModelDivergent::restoreCheckpoint() {
    startCheckpoint();
    int model_num = 1;
    for (ModelMarkov* subtree_model : subtree_models) {
        auto model_num_string = convertIntToString(model_num);
        checkpoint->startStruct("Model" + model_num_string);
        subtree_model->restoreCheckpoint();
        checkpoint->endStruct();
        ++model_num;
    }
    endCheckpoint();
    decomposeRateMatrix();
    if (phylo_tree!=nullptr) {
        phylo_tree->clearAllPartialLH();
    }
}

void ModelDivergent::setMixtureWeight(int cat, double weight) {
    ASSERT(0 && "Cannot set mixture weight on ModelDivergent");
}

void ModelDivergent::setFixMixtureWeight(bool fix_prop) {
    ASSERT(0 && "Cannot fix mixture weight on ModelDivergent");
}

void ModelDivergent::setMixtureClass(int cat, ModelSubst* m) {
    ASSERT(0 && "Cannot set mixture class on ModelDivergent");
}

void ModelDivergent::computeTransMatrix
        (double time, double *trans_matrix, 
         int subtree_number /* = 0 */) const {
    ASSERT(0 && "Cannot compute transition matrix of ModelDivergent");
    //Because the transition matrices will be different, 
    //in different regions of the phylogenetic tree.
}

double ModelDivergent::computeTrans
        (double time, int state1, int state2) const {
    ASSERT(0 && "Cannot compute transition probabilities for ModelDivergent");
    return 0;
}

double ModelDivergent::computeTrans
        (double time, int state1, int state2, 
         double &derv1, double &derv2) const {
    ASSERT(0 && "Cannot compute transition probabilities for ModelDivergent");                            
    return 0;
}

void ModelDivergent::getRateMatrix(double *rate_mat) const {
    ASSERT(0 && "Cannot compute transition probabilities for ModelDivergent");
}

void ModelDivergent::getQMatrix(double *q_mat) const {
    ASSERT(0 && "Cannot compute transition probabilities for ModelDivergent");
}

void ModelDivergent::computeTransDerv
        (double time, double *trans_matrix, 
         double *trans_derv1, double *trans_derv2, 
         int mixture) const {
    ASSERT(0 && "Cannot compute derivative" 
           " of transition matrix for ModelDivergent");
}

double* ModelDivergent::getEigenvalues()                   const {
    ASSERT(0 && "Cannot compute eigen values" 
           " for ModelDivergent");
    return nullptr;
}

double* ModelDivergent::getEigenvectors()                  const {
    ASSERT(0 && "Cannot compute eigen vectors" 
           " for ModelDivergent");
    return nullptr;
}

double* ModelDivergent::getInverseEigenvectors()           const {
    ASSERT(0 && "Cannot compute inverse" 
           " eigen vectors for ModelDivergent");
    return nullptr;
}

double* ModelDivergent::getInverseEigenvectorsTransposed() const {
    ASSERT(0 && "Cannot compute transpose of"
           " inverse eigen vectors for ModelDivergent");
    return nullptr;
}

