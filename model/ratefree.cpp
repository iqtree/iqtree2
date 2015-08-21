/*
 * ratefree.cpp
 *
 *  Created on: Nov 3, 2014
 *      Author: minh
 */

#include "phylotree.h"
#include "ratefree.h"

#include "model/modelfactory.h"
#include "model/modelmixture.h"


const double MIN_FREE_RATE = 0.001;
const double MAX_FREE_RATE = 1000.0;
const double TOL_FREE_RATE = 0.0001;

// Modified by Thomas on 13 May 2015
//const double MIN_FREE_RATE_PROP = 0.0001;
//const double MAX_FREE_RATE_PROP = 0.9999;
const double MIN_FREE_RATE_PROP = 0.001;
const double MAX_FREE_RATE_PROP = 1000;

RateFree::RateFree(int ncat, double start_alpha, string params, bool sorted_rates, string opt_alg, PhyloTree *tree) : RateGamma(ncat, start_alpha, false, tree) {
	fix_params = false;
	prop = NULL;
    this->sorted_rates = sorted_rates;
    optimizing_params = 0;
    this->optimize_alg = opt_alg;
	setNCategory(ncat);

	if (params.empty()) return;
	DoubleVector params_vec;
	try {
		convert_double_vec(params.c_str(), params_vec);
		if (params_vec.size() != ncategory*2)
			outError("Number of parameters for FreeRate model must be twice number of categories");
		int i;
		double sum, sum_prop;
		for (i = 0, sum = 0.0, sum_prop = 0.0; i < ncategory; i++) {
			prop[i] = params_vec[i*2];
			rates[i] = params_vec[i*2+1];
			sum += prop[i]*rates[i];
			sum_prop += prop[i];
		}
		if (fabs(sum_prop-1.0) > 1e-5)
			outError("Sum of category proportions not equal to 1");
		for (i = 0; i < ncategory; i++)
			rates[i] /= sum;
		fix_params = true;
	} catch (string &str) {
		outError(str);
	}
}

void RateFree::setNCategory(int ncat) {

    // initialize with gamma rates
    RateGamma::setNCategory(ncat);
	if (prop) delete [] prop;
	prop  = new double[ncategory];

    int i;
	for (i = 0; i < ncategory; i++)
        prop[i] = 1.0/ncategory;
    
//	double sum_prop = (ncategory)*(ncategory+1)/2.0;
//	double sum = 0.0;
//	int i;
	// initialize rates as increasing
//	for (i = 0; i < ncategory; i++) {
//		prop[i] = (double)(ncategory-i) / sum_prop;
//        prop[i] = 1.0 / ncategory;
//		rates[i] = (double)(i+1);
//		sum += prop[i]*rates[i];
//	}
//	for (i = 0; i < ncategory; i++)
//		rates[i] /= sum;

	name = "+R";
	name += convertIntToString(ncategory);
	full_name = "FreeRate";
	full_name += " with " + convertIntToString(ncategory) + " categories";
}

void RateFree::setRateAndProp(RateFree *input) {
    assert(input->ncategory == ncategory-1);
    int k = 0, i;
    // get the category k with largest proportion
    for (i = 1; i < ncategory-1; i++)
        if (input->prop[i] > input->prop[k]) k = i;

    memcpy(rates, input->rates, (k+1)*sizeof(double));
    memcpy(prop, input->prop, (k+1)*sizeof(double));
    rates[k+1] = 1.414*rates[k]; // sqrt(2)
    prop[k+1] = prop[k]/2;
    rates[k] = 0.586*rates[k];
    prop[k] = prop[k]/2;
    if (k < ncategory-2) {
        memcpy(&rates[k+2], &input->rates[k+1], (ncategory-2-k)*sizeof(double));
        memcpy(&prop[k+2], &input->prop[k+1], (ncategory-2-k)*sizeof(double));
    }
    // copy half of k to the last category


//    rates[ncategory-1] = rates[k];
//    prop[ncategory-1] = prop[k] / 2;
//    prop[k] = prop[k] / 2;
    // sort the rates in increasing order
    if (sorted_rates)
        quicksort(rates, 0, ncategory-1, prop);
    phylo_tree->clearAllPartialLH();
}


RateFree::~RateFree() {
	if (prop) delete [] prop;
	prop = NULL;
}

string RateFree::getNameParams() {
	stringstream str;
	str << "+R" << ncategory << "{";
	for (int i = 0; i < ncategory; i++) {
		if (i > 0) str << ",";
		str << prop[i]<< "," << rates[i];
	}
	str << "}";
	return str.str();
}

int RateFree::getNDim() { 
    if (fix_params) return 0;
    if (optimizing_params == 0) return (2*ncategory-2); 
    if (optimizing_params == 1) // rates
        return ncategory;
    if (optimizing_params == 2) // proportions
        return ncategory-1;
    return 0;
}

double RateFree::targetFunk(double x[]) {
	getVariables(x);
    if (optimizing_params == 1)
        // only clear partial_lh if optimizing rates
        phylo_tree->clearAllPartialLH();
	return -phylo_tree->computeLikelihood();
}



/**
	optimize parameters. Default is to optimize gamma shape
	@return the best likelihood
*/
double RateFree::optimizeParameters(double gradient_epsilon) {

	int ndim = getNDim();

	// return if nothing to be optimized
	if (ndim == 0)
		return phylo_tree->computeLikelihood();

    if (optimize_alg == "EM")
        return optimizeWithEM();

	if (verbose_mode >= VB_MAX)
		cout << "Optimizing " << name << " model parameters by BFGS..." << endl;

	//if (freq_type == FREQ_ESTIMATE) scaleStateFreq(false);

	double *variables = new double[ndim+1];
	double *upper_bound = new double[ndim+1];
	double *lower_bound = new double[ndim+1];
	bool *bound_check = new bool[ndim+1];
	double score;

//    score = optimizeWeights();

    for (optimizing_params = 1; optimizing_params <= 2; optimizing_params++) {
    ndim = getNDim();
	// by BFGS algorithm
	setVariables(variables);
	setBounds(lower_bound, upper_bound, bound_check);

//	score = -minimizeMultiDimen(variables, ndim, lower_bound, upper_bound, bound_check, max(gradient_epsilon, TOL_FREE_RATE));
    score = -L_BFGS_B(ndim, variables+1, lower_bound+1, upper_bound+1, max(gradient_epsilon, TOL_FREE_RATE));

	getVariables(variables);
    // sort the rates in increasing order
    if (sorted_rates)
        quicksort(rates, 0, ncategory-1, prop);
	phylo_tree->clearAllPartialLH();
    }
    optimizing_params = 0;

	delete [] bound_check;
	delete [] lower_bound;
	delete [] upper_bound;
	delete [] variables;

	return score;
}

void RateFree::setBounds(double *lower_bound, double *upper_bound, bool *bound_check) {
	if (getNDim() == 0) return;
	int i;
    if (optimizing_params == 2) {
        // proportions
        for (i = 1; i < ncategory; i++) {
            lower_bound[i] = MIN_FREE_RATE_PROP;
            upper_bound[i] = MAX_FREE_RATE_PROP;
            bound_check[i] = false;
        }
    } else {
        // rates
        for (i = 1; i <= ncategory; i++) {
            lower_bound[i] = MIN_FREE_RATE;
            upper_bound[i] = MAX_FREE_RATE;
            bound_check[i] = false;
        }
    }
//	for (i = ncategory; i <= 2*ncategory-2; i++) {
//		lower_bound[i] = MIN_FREE_RATE;
//		upper_bound[i] = MAX_FREE_RATE;
//		bound_check[i] = false;
//	}
}


void RateFree::setVariables(double *variables) {
	if (getNDim() == 0) return;
	int i;

	// Modified by Thomas on 13 May 2015
	// --start--
	/*
	variables[1] = prop[0];
	for (i = 2; i < ncategory; i++)
		variables[i] = variables[i-1] + prop[i-1];
	*/
    
    if (optimizing_params == 2) {    
        // proportions
        for (i = 0; i < ncategory-1; i++)
            variables[i+1] = prop[i] / prop[ncategory-1];
    } else {
        // rates
        for (i = 0; i < ncategory; i++)
            variables[i+1] = rates[i];
    }

//	for (i = 0; i < ncategory-1; i++)
//		variables[i+ncategory] = rates[i] / rates[ncategory-1];
}

void RateFree::getVariables(double *variables) {
	if (getNDim() == 0) return;
	int i;

	// Modified by Thomas on 13 May 2015
	// --start--
	/*
	double *y = new double[2*ncategory+1];
	double *z = y+ncategory+1;
	//  site proportions: y[0..c] <-> (0.0, variables[1..c-1], 1.0)
	y[0] = 0; y[ncategory] = 1.0;
	memcpy(y+1, variables+1, (ncategory-1) * sizeof(double));
	std::sort(y+1, y+ncategory);

	// category rates: z[0..c-1] <-> (variables[c..2*c-2], 1.0)
	memcpy(z, variables+ncategory, (ncategory-1) * sizeof(double));
	z[ncategory-1] = 1.0;
	//std::sort(z, z+ncategory-1);

	double sum = 0.0;
	for (i = 0; i < ncategory; i++) {
		prop[i] = (y[i+1]-y[i]);
		sum += prop[i] * z[i];
	}
	for (i = 0; i < ncategory; i++) {
		rates[i] = z[i] / sum;
	}

	delete [] y;
	*/

	double sum = 1.0;
    if (optimizing_params == 2) {
        // proportions
        for (i = 0; i < ncategory-1; i++) {
            sum += variables[i+1];
        }
        for (i = 0; i < ncategory-1; i++) {
            prop[i] = variables[i+1] / sum;
        }
        prop[ncategory-1] = 1.0 / sum;
    } else {
        // rates
        for (i = 0; i < ncategory; i++)
            rates[i] = variables[i+1];
    }
//	sum = prop[ncategory-1];
//	for (i = 0; i < ncategory-1; i++) {
//		sum += prop[i] * variables[i+ncategory];
//	}
//	for (i = 0; i < ncategory-1; i++) {
//		rates[i] = variables[i+ncategory] / sum;
//	}
//	rates[ncategory-1] = 1.0 / sum;
	// --end--

}

/**
	write information
	@param out output stream
*/
void RateFree::writeInfo(ostream &out) {
	out << "Site proportion and rates: ";
	for (int i = 0; i < ncategory; i++)
		out << " (" << prop[i] << "," << rates[i] << ")";
	out << endl;
}

/**
	write parameters, used with modeltest
	@param out output stream
*/
void RateFree::writeParameters(ostream &out) {
	for (int i = 0; i < ncategory; i++)
		out << "\t" << prop[i] << "\t" << rates[i];

}

double RateFree::optimizeWithEM() {
    size_t ptn, c;
    size_t nptn = phylo_tree->aln->getNPattern();
    size_t nmix = ncategory;
    
//    double *lk_ptn = aligned_alloc<double>(nptn);
    double *new_prop = aligned_alloc<double>(nmix);
    PhyloTree *tree = new PhyloTree;
    tree->copyPhyloTree(phylo_tree);
    tree->optimize_by_newton = phylo_tree->optimize_by_newton;
    tree->setLikelihoodKernel(phylo_tree->sse);
    // initialize model
    ModelFactory *model_fac = new ModelFactory();
    model_fac->joint_optimize = phylo_tree->params->optimize_model_rate_joint;

    RateHeterogeneity *site_rate = new RateHeterogeneity; 
    tree->setRate(site_rate);
    site_rate->setTree(tree);
            
    model_fac->site_rate = site_rate;
    tree->model_factory = model_fac;
    tree->setParams(phylo_tree->params);
        
    // EM algorithm loop described in Wang, Li, Susko, and Roger (2008)
    for (int step = 0; step < 10; step++) {
        // first compute _pattern_lh_cat
        double score;
        if (!phylo_tree->getModel()->isMixture())
            score = phylo_tree->computeLikelihoodBranchEigen((PhyloNeighbor*)phylo_tree->root->neighbors[0], (PhyloNode*)phylo_tree->root); 
        else if (phylo_tree->getModelFactory()->fused_mix_rate) {
            score = phylo_tree->computeMixrateLikelihoodBranchEigen((PhyloNeighbor*)phylo_tree->root->neighbors[0], (PhyloNode*)phylo_tree->root); 
        } else {
            outError("Mixture model does not work with FreeRate model!");
            score = phylo_tree->computeMixtureLikelihoodBranchEigen((PhyloNeighbor*)phylo_tree->root->neighbors[0], (PhyloNode*)phylo_tree->root); 
        }
        memset(new_prop, 0, nmix*sizeof(double));
                
        // E-step
        // decoupled weights (prop) from _pattern_lh_cat to obtain L_ci and compute pattern likelihood L_i
        for (ptn = 0; ptn < nptn; ptn++) {
            double *this_lk_cat = phylo_tree->_pattern_lh_cat + ptn*nmix;
            double lk_ptn = 0.0;
            for (c = 0; c < nmix; c++) {
                lk_ptn += this_lk_cat[c];
            }
            lk_ptn = phylo_tree->ptn_freq[ptn] / lk_ptn;
            
            // transform _pattern_lh_cat into posterior probabilities of each category
            for (c = 0; c < nmix; c++) {
                this_lk_cat[c] *= lk_ptn;
                new_prop[c] += this_lk_cat[c];
            }
            
        } 
        
        // M-step, update weights according to (*)        
        
        bool converged = true;
        for (c = 0; c < nmix; c++) {
            new_prop[c] = new_prop[c] / phylo_tree->getAlnNSite();
            // check for convergence
            converged = converged && (fabs(prop[c]-new_prop[c]) < 1e-4);
            prop[c] = new_prop[c];
        }
        
        // now optimize rates one by one
        double sum = 0.0;
        for (c = 0; c < nmix; c++) {
            tree->copyPhyloTree(phylo_tree);
            ModelGTR *subst_model;
            if (phylo_tree->getModel()->isMixture())
                subst_model = ((ModelMixture*)phylo_tree->getModel())->at(c);
            else
                subst_model = (ModelGTR*)phylo_tree->getModel();
            tree->setModel(subst_model);
            subst_model->setTree(tree);
            model_fac->model = subst_model;
                        
            // initialize likelihood
            tree->initializeAllPartialLh();
            // copy posterior probability into ptn_freq
            tree->computePtnFreq();
            double *this_lk_cat = phylo_tree->_pattern_lh_cat+c;
            for (ptn = 0; ptn < nptn; ptn++)
                tree->ptn_freq[ptn] = this_lk_cat[ptn*nmix];
            double scaling = rates[c];
            tree->scaleLength(scaling);
            tree->optimizeTreeLengthScaling(scaling, 0.001);
            converged = converged && (fabs(rates[c] - scaling) < 1e-4);
            rates[c] = scaling;
            sum += prop[c] * rates[c];
            // reset subst model
            tree->setModel(NULL);
            subst_model->setTree(phylo_tree);
            
        }
        
        phylo_tree->clearAllPartialLH();
        if (converged) break;
    }
    
    delete tree;
    aligned_free(new_prop);
    return phylo_tree->computeLikelihood();
}
