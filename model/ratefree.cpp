/*
 * ratefree.cpp
 *
 *  Created on: Nov 3, 2014
 *      Author: minh
 */

#include "tree/phylotree.h"
#include "ratefree.h"
#include "rateinvar.h"

#include "model/modelfactory.h"
#include "model/modelmixture.h"
#include "utils/timeutil.h" //temporary : for time log-lining

const double MIN_FREE_RATE = 0.001;
const double MAX_FREE_RATE = 1000.0;
const double TOL_FREE_RATE = 0.0001;

// Modified by Thomas on 13 May 2015
//const double MIN_FREE_RATE_PROP = 0.0001;
//const double MAX_FREE_RATE_PROP = 0.9999;
const double MIN_FREE_RATE_PROP = 0.001;
const double MAX_FREE_RATE_PROP = 1000;

RateFree::RateFree(int ncat, double start_alpha, string params, bool sorted_rates, string opt_alg, PhyloTree *tree) : RateGamma(ncat, start_alpha, false, tree) {
	fix_params = 0;
	prop = NULL;
    this->sorted_rates = sorted_rates;
    optimizing_params = 0;
    this->optimize_alg = opt_alg;
	setNCategory(ncat);

	if (params.empty()) return;
	DoubleVector params_vec;
	try {
		convert_double_vec(params.c_str(), params_vec);
		int i;
		double sum, sum_prop;
        if (params_vec.size() == ncategory) {
            // only inputing prop
            for (i = 0, sum_prop = 0.0; i < ncategory; i++) {
                prop[i] = params_vec[i];
                rates[i] = 1.0;
                sum_prop += prop[i];
            }
            fix_params = (Params::getInstance().optimize_from_given_params) ? 0 : 1;
        } else {
            if (params_vec.size() != ncategory*2)
                outError("Number of parameters for FreeRate model must be twice number of categories");
            for (i = 0, sum = 0.0, sum_prop = 0.0; i < ncategory; i++) {
                prop[i] = params_vec[i*2];
                rates[i] = params_vec[i*2+1];
                sum += prop[i]*rates[i];
                sum_prop += prop[i];
            }
            for (i = 0; i < ncategory; i++)
                rates[i] /= sum;
            fix_params = (Params::getInstance().optimize_from_given_params) ? 0 : 2;
        }
		if (fabs(sum_prop-1.0) > 1e-5)
			outError("Sum of category proportions not equal to 1");
	} catch (string &str) {
		outError(str);
	}
}

void RateFree::startCheckpoint() {
    checkpoint->startStruct("RateFree" + convertIntToString(ncategory));
}

void RateFree::saveCheckpoint() {
    startCheckpoint();
//    CKP_SAVE(fix_params);
//    CKP_SAVE(sorted_rates);
//    CKP_SAVE(optimize_alg);
    CKP_ARRAY_SAVE(ncategory, prop);
    CKP_ARRAY_SAVE(ncategory, rates);
    endCheckpoint();
//    RateGamma::saveCheckpoint();
}

void RateFree::restoreCheckpoint() {
//    RateGamma::restoreCheckpoint();
    startCheckpoint();
//    CKP_RESTORE(fix_params);
//    CKP_RESTORE(sorted_rates);
//    CKP_RESTORE(optimize_alg);
    CKP_ARRAY_RESTORE(ncategory, prop);
    CKP_ARRAY_RESTORE(ncategory, rates);
    endCheckpoint();

//	setNCategory(ncategory);
}

void RateFree::setNCategory(int ncat) {

    // initialize with gamma rates
    RateGamma::setNCategory(ncat);
	if (prop) delete [] prop;
	prop  = new double[ncategory];

    for (int i = 0; i < ncategory; i++) {
        prop[i] = (1.0-getPInvar())/ncategory;
    }
    
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

void RateFree::initFromCatMinusOne() {
    ncategory--;
    restoreCheckpoint();
    ncategory++;

    int first = 0;
    // get the category k with largest proportion
    for (int i = 1; i < ncategory-1; i++) {
        if (prop[i] > prop[first]) {
            first = i;
        }
    }
    int second = (first == 0) ? 1 : 0;
    for (int i = 0; i < ncategory-1; i++)
        if (prop[i] > prop[second] && i != first)
            second = i;

//    memmove(rates, input->rates, (k+1)*sizeof(double));
//    memmove(prop, input->prop, (k+1)*sizeof(double));

    // divide highest category into 2 of the same prop
    // 2018-06-12: fix bug negative rates
    if (-rates[second] + 3*rates[first] > 0.0) {
        rates[ncategory-1] = (-rates[second] + 3*rates[first])/2.0;
        rates[first] = (rates[second]+rates[first])/2.0;
    } else {
        rates[ncategory-1] = (3*rates[first])/2.0;
        rates[first] = (rates[first])/2.0;
    }
    prop[ncategory-1] = prop[first]/2;
    prop[first] = prop[first]/2;
//    if (k < ncategory-2) {
//        memcpy(&rates[k+2], &input->rates[k+1], (ncategory-2-k)*sizeof(double));
//        memcpy(&prop[k+2], &input->prop[k+1], (ncategory-2-k)*sizeof(double));
//    }
    // copy half of k to the last category


//    rates[ncategory-1] = rates[k];
//    prop[ncategory-1] = prop[k] / 2;
//    prop[k] = prop[k] / 2;
    // sort the rates in increasing order
    if (sorted_rates) {
        quicksort(rates, 0, ncategory-1, prop);
    }
    phylo_tree->clearAllPartialLH();
}


RateFree::~RateFree() {
    delete [] prop;
    prop = nullptr;
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

double RateFree::meanRates() {
	double ret = 0.0;
	for (int i = 0; i < ncategory; i++)
		ret += prop[i] * rates[i];
	return ret;
}

/**
 * rescale rates s.t. mean rate is equal to 1, useful for FreeRate model
 * @return rescaling factor
 */
double RateFree::rescaleRates() {
	double norm = meanRates();
	for (int i = 0; i < ncategory; i++)
		rates[i] /= norm;
	return norm;
}

int RateFree::getNDim() { 
    if (fix_params == 2) return 0;
    if (fix_params == 1) // only fix prop
        return (ncategory-1);
    if (optimizing_params == 0) return (2*ncategory-2); 
    if (optimizing_params == 1) // rates
        return ncategory-1;
    if (optimizing_params == 2) // proportions
        return ncategory-1;
    return 0;
}

double RateFree::targetFunk(double x[]) {
	getVariables(x);
    if (optimizing_params != 2) {
        // only clear partial_lh if optimizing rates
        phylo_tree->clearAllPartialLH();
    }
	return -phylo_tree->computeLikelihood();
}

/**
	optimize parameters. Default is to optimize gamma shape
	@return the best likelihood
*/
double RateFree::optimizeParameters(double gradient_epsilon) {

	int ndim = getNDim();

	// return if nothing to be optimized
    if (ndim == 0) {
        return phylo_tree->computeLikelihood();
    }
    TREE_LOG_LINE(*phylo_tree, VB_MED, "Optimizing " << name << " model parameters by " << optimize_alg << " algorithm...");
    // TODO: turn off EM algorithm for +ASC model
    if ((optimize_alg.find("EM") != string::npos && phylo_tree->getModelFactory()->unobserved_ptns.empty())) {
        if (fix_params == 0) {
            return optimizeWithEM();
        }
    }
	//if (freq_type == FREQ_ESTIMATE) scaleStateFreq(false);

	double *variables = new double[ndim+1];
	double *upper_bound = new double[ndim+1];
	double *lower_bound = new double[ndim+1];
	bool *bound_check = new bool[ndim+1];
	double score;

//    score = optimizeWeights();

    int left = 1, right = 2;
    if (fix_params == 1) // fix proportions
        right = 1;
    if (optimize_alg.find("1-BFGS") != string::npos) {
        left = 0; 
        right = 0;
    }

    // changed to Wi -> Ri by Thomas on Sept 11, 15
    for (optimizing_params = right; optimizing_params >= left; optimizing_params--) {
    
        ndim = getNDim();
        // by BFGS algorithm
        setVariables(variables);
        setBounds(lower_bound, upper_bound, bound_check);

//        if (optimizing_params == 2 && optimize_alg.find("-EM") != string::npos)
//            score = optimizeWeights();
//        else 
        if (optimize_alg.find("BFGS-B") != string::npos)
            score = -L_BFGS_B(ndim, variables+1, lower_bound+1, upper_bound+1, max(gradient_epsilon, TOL_FREE_RATE));
        else
            score = -minimizeMultiDimen(variables, ndim, lower_bound, upper_bound, bound_check, max(gradient_epsilon, TOL_FREE_RATE));

        getVariables(variables);
        // sort the rates in increasing order
        if (sorted_rates)
            quicksort(rates, 0, ncategory-1, prop);
        phylo_tree->clearAllPartialLH();
        score = phylo_tree->computeLikelihood();
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
    } else if (optimizing_params == 1){
        // rates
        for (i = 1; i < ncategory; i++) {
            lower_bound[i] = MIN_FREE_RATE;
            upper_bound[i] = MAX_FREE_RATE;
            bound_check[i] = false;
        }
    } else {
        // both weights and rates
        for (i = 1; i < ncategory; i++) {
            lower_bound[i] = MIN_FREE_RATE_PROP;
            upper_bound[i] = MAX_FREE_RATE_PROP;
            bound_check[i] = false;
        }
        for (i = 1; i < ncategory; i++) {
            lower_bound[i+ncategory-1] = MIN_FREE_RATE;
            upper_bound[i+ncategory-1] = MAX_FREE_RATE;
            bound_check[i+ncategory-1] = false;
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
    } else if (optimizing_params == 1) {
        // rates
        for (i = 0; i < ncategory-1; i++)
            variables[i+1] = rates[i];
    } else {
        // both rates and weights
        for (i = 0; i < ncategory-1; i++) {
            variables[i+1] = prop[i] / prop[ncategory-1];
        }
        for (i = 0; i < ncategory-1; i++) {
            variables[i+ncategory] = rates[i] / rates[ncategory-1];
        }
    }
}

bool RateFree::getVariables(double *variables) {
	if (getNDim() == 0) return false;
	int i;
    bool changed = false;
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
            changed |= (prop[i] != variables[i+1] / sum);
            prop[i] = variables[i+1] / sum;
        }
        changed |= (prop[ncategory-1] != 1.0 / sum);
        prop[ncategory-1] = 1.0 / sum;
        // added by Thomas on Sept 10, 15
        // update the values of rates, in order to
        // maintain the sum of prop[i]*rates[i] = 1
//        sum = 0;
//        for (i = 0; i < ncategory; i++) {
//            sum += prop[i] * rates[i];
//        }
//        for (i = 0; i < ncategory; i++) {
//            rates[i] = rates[i] / sum;
//        }
    } else if (optimizing_params == 1) {
        // rates
        for (i = 0; i < ncategory-1; i++) {
            changed |= (rates[i] != variables[i+1]);
            rates[i] = variables[i+1];
        }
        // added by Thomas on Sept 10, 15
        // need to normalize the values of rates, in order to
        // maintain the sum of prop[i]*rates[i] = 1
//        sum = 0;
//        for (i = 0; i < ncategory; i++) {
//            sum += prop[i] * rates[i];
//        }
//        for (i = 0; i < ncategory; i++) {
//            rates[i] = rates[i] / sum;
//        }
    } else {
        // both weights and rates
        for (i = 0; i < ncategory-1; i++) {
            sum += variables[i+1];
        }
        for (i = 0; i < ncategory-1; i++) {
            changed |= (prop[i] != variables[i+1] / sum);
            prop[i] = variables[i+1] / sum;
        }
        changed |= (prop[ncategory-1] != 1.0 / sum);
        prop[ncategory-1] = 1.0 / sum;
        
        // then rates
    	sum = prop[ncategory-1];
    	for (i = 0; i < ncategory-1; i++) {
    		sum += prop[i] * variables[i+ncategory];
    	}
    	for (i = 0; i < ncategory-1; i++) {
            changed |= (rates[i] != variables[i+ncategory] / sum);
    		rates[i] = variables[i+ncategory] / sum;
    	}
        changed |= (rates[ncategory-1] != 1.0 / sum);
    	rates[ncategory-1] = 1.0 / sum;
    }
	// --end--
    return changed;
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
    const double MIN_PROP = 1e-4;
    
//    double *lk_ptn = aligned_alloc<double>(nptn);
    double *new_prop = aligned_alloc<double>(nmix);
    PhyloTree *tree = new PhyloTree;

    // attach memory to save space
//    tree->central_partial_lh = phylo_tree->central_partial_lh;
//    tree->central_scale_num = phylo_tree->central_scale_num;
//    tree->central_partial_pars = phylo_tree->central_partial_pars;

    tree->copyPhyloTree(phylo_tree, true);
    tree->optimize_by_newton = phylo_tree->optimize_by_newton;
    tree->setParams(phylo_tree->params);
    tree->setLikelihoodKernel(phylo_tree->sse);
    tree->setNumThreads(phylo_tree->num_threads);

    // initialize model
    ModelFactory *model_fac = new ModelFactory();
    model_fac->joint_optimize = phylo_tree->params->optimize_model_rate_joint;
//    model_fac->unobserved_ptns = phylo_tree->getModelFactory()->unobserved_ptns;

    RateHeterogeneity *site_rate = new RateHeterogeneity; 
    tree->setRate(site_rate);
    site_rate->setTree(tree);
            
    model_fac->site_rate = site_rate;
    tree->model_factory = model_fac;
    tree->setParams(phylo_tree->params);
    double old_score = 0.0;
    // EM algorithm loop described in Wang, Li, Susko, and Roger (2008)
    for (int step = 0; step < ncategory; step++) {
        // first compute _pattern_lh_cat
        double score;
        score = phylo_tree->computePatternLhCat(WSL_RATECAT);
        if (score > 0.0) {
            phylo_tree->printTree(cout, WT_BR_LEN+WT_NEWLINE);
            writeInfo(cout);
        }
        ASSERT(score < 0);
        
        if (step > 0) {
            if (score <= old_score-0.1) {
                phylo_tree->printTree(cout, WT_BR_LEN+WT_NEWLINE);
                writeInfo(cout);
                cout << "Partition " << phylo_tree->aln->name << endl;
                cout << "score: " << score << "  old_score: " << old_score << endl;
            }
            if (!Params::getInstance().ignore_any_errors) {
                ASSERT(score > old_score-0.1);
            }
        }
        old_score = score;
        
                
        // E-step
        // decoupled weights (prop) from _pattern_lh_cat to obtain L_ci and compute pattern likelihood L_i
        memset(new_prop, 0, nmix*sizeof(double));
        for (ptn = 0; ptn < nptn; ptn++) {
            double *this_lk_cat = phylo_tree->tree_buffers._pattern_lh_cat + ptn*nmix;
            double lk_ptn = phylo_tree->ptn_invar[ptn];
            for (c = 0; c < nmix; c++) {
                lk_ptn += this_lk_cat[c];
            }
            ASSERT(lk_ptn != 0.0);
            lk_ptn = phylo_tree->ptn_freq[ptn] / lk_ptn;
            
            // transform _pattern_lh_cat into posterior probabilities of each category
            for (c = 0; c < nmix; c++) {
                this_lk_cat[c] *= lk_ptn;
                new_prop[c] += this_lk_cat[c];
            }
        }
        
        // M-step, update weights according to (*)
        int maxpropid = 0;
        double new_pinvar = 0.0;
        for (c = 0; c < nmix; c++) {
            new_prop[c] = new_prop[c] / phylo_tree->getAlnNSite();
            if (new_prop[c] > new_prop[maxpropid])
                maxpropid = c;
        }
        // regularize prop
        bool zero_prop = false;
        for (c = 0; c < nmix; c++) {
            if (new_prop[c] < MIN_PROP) {
                new_prop[maxpropid] -= (MIN_PROP - new_prop[c]);
                new_prop[c] = MIN_PROP;
                zero_prop = true;
            }
        }
        // break if some probabilities too small
        if (zero_prop) break;

        bool converged = true;
        double sum_prop = 0.0;
        for (c = 0; c < nmix; c++) {
//            new_prop[c] = new_prop[c] / phylo_tree->getAlnNSite();
            // check for convergence
            sum_prop += new_prop[c];
            converged = converged && (fabs(prop[c]-new_prop[c]) < 1e-4);
            prop[c] = new_prop[c];
            new_pinvar += new_prop[c];
        }

        new_pinvar = 1.0 - new_pinvar;

        if (new_pinvar > 1e-4 && getPInvar() != 0.0) {
            converged = converged && (fabs(getPInvar()-new_pinvar) < 1e-4);
            if (isFixPInvar())
                outError("Fixed given p-invar is not supported");
            setPInvar(new_pinvar);
//            setOptimizePInvar(false);
            phylo_tree->computePtnInvar();
        }
        
        ASSERT(fabs(sum_prop+new_pinvar-1.0) < MIN_PROP);
        
        // now optimize rates one by one
        double sum = 0.0;
        for (c = 0; c < nmix; c++) {
            tree->copyPhyloTree(phylo_tree, true);
            ModelMarkov *subst_model;
            if (phylo_tree->getModel()->isMixture() && phylo_tree->getModelFactory()->fused_mix_rate)
                subst_model = (ModelMarkov*)phylo_tree->getModel()->getMixtureClass(c);
            else
                subst_model = (ModelMarkov*)phylo_tree->getModel();
            tree->setModel(subst_model);
            subst_model->setTree(tree);
            model_fac->model = subst_model;
            if (subst_model->isMixture() || subst_model->isSiteSpecificModel() || !subst_model->isReversible())
                tree->setLikelihoodKernel(phylo_tree->sse);

            // initialize likelihood
            tree->initializeAllPartialLh();
            // copy posterior probability into ptn_freq
            tree->computePtnFreq();
            double *this_lk_cat = phylo_tree->tree_buffers._pattern_lh_cat+c;
            for (ptn = 0; ptn < nptn; ptn++) {
                tree->ptn_freq[ptn] = this_lk_cat[ptn*nmix];
            }
            double scaling = rates[c];
            tree->scaleLength(scaling);
            tree->optimizeTreeLengthScaling(MIN_PROP, scaling, 1.0/prop[c], 0.001);
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
    
    // sort the rates in increasing order
    if (sorted_rates) {
        quicksort(rates, 0, ncategory-1, prop);
    }
    
    // deattach memory
//    tree->central_partial_lh = NULL;
//    tree->central_scale_num = NULL;
//    tree->central_partial_pars = NULL;

    delete tree;
    aligned_free(new_prop);
    return phylo_tree->computeLikelihood();
}
