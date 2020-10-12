//
//  rateheterotachy.cpp
//  iqtree
//
//  Created by Minh Bui on 11/8/16.
//
//

#include "tree/phylotree.h"
#include "rateheterotachy.h"

RateHeterotachy::RateHeterotachy(int ncat, string params, PhyloTree *tree) : RateHeterogeneity() {
    phylo_tree = tree;
    prop = NULL;
    fix_params = 0;
    optimize_steps = 0;
    setNCategory(ncat);

	if (params.empty()) return;
	DoubleVector params_vec;
	try {
		convert_double_vec(params.c_str(), params_vec);
		if (params_vec.size() != ncategory)
			outError("Number of parameters for rate heterotachy model must equal number of categories");
		int i;
		double sum_prop;
		for (i = 0, sum_prop = 0.0; i < ncategory; i++) {
			prop[i] = params_vec[i];
			sum_prop += prop[i];
		}
		if (fabs(sum_prop-1.0) > 1e-5)
			outError("Sum of category proportions not equal to 1");
		// Minh: Please double check this one. It isn't quite so
		// clear what fix_params is doing, as it seems to take values
		// 0, 1 or 2.  -- MDW
        //BQM: that OK

		if (!(tree->params->optimize_from_given_params)) {
		        fix_params = 1;
		} // else fix_params == 0 still. 
	} catch (string &str) {
		outError(str);
	}
}

/**
    destructor
*/
RateHeterotachy::~RateHeterotachy() {
    delete [] prop;
    prop = nullptr;
}

void RateHeterotachy::setNCategory(int ncat) {
    ncategory = ncat;
    if (optimize_steps == 0) {
        optimize_steps = ncat * 100;
    }
    // initialize with gamma rates
	delete [] prop;
	prop  = new double[ncategory];

    int i;
	for (i = 0; i < ncategory; i++)
        prop[i] = (1.0-getPInvar())/ncategory;
    
	name = "+H";
	name += convertIntToString(ncategory);
	full_name = "Rate heterotachy";
	full_name += " with " + convertIntToString(ncategory) + " categories";
}

void RateHeterotachy::startCheckpoint() {
    checkpoint->startStruct("RateHeterotachy" + convertIntToString(ncategory));
}

/**
    save object into the checkpoint
*/
void RateHeterotachy::saveCheckpoint() {
    startCheckpoint();
    CKP_ARRAY_SAVE(ncategory, prop);
    endCheckpoint();
    RateHeterogeneity::saveCheckpoint();
}

/**
    restore object from the checkpoint
*/
void RateHeterotachy::restoreCheckpoint() {
    RateHeterogeneity::restoreCheckpoint();
    startCheckpoint();
    CKP_ARRAY_RESTORE(ncategory, prop);
    endCheckpoint();
}

int RateHeterotachy::getNDim() {
    if (fix_params) return 0;
    return ncategory-1;
}

/**
 * @return model name with parameters in form of e.g. GTR{a,b,c,d,e,f}
 */
string RateHeterotachy::getNameParams() {
	stringstream str;
	str << "+H" << ncategory << "{";
	for (int i = 0; i < ncategory; i++) {
		if (i > 0) str << ",";
		str << prop[i];
	}
	str << "}";
	return str.str();
}

void RateHeterotachy::writeInfo(ostream &out) {
    if (fix_params != 2) {
        out << "Heterotachy weights:     ";
        for (int i = 0; i < ncategory; i++)
            out << " " << prop[i];
        out << endl;
    }
    DoubleVector lenvec;
    phylo_tree->treeLengths(lenvec);
    out << "Heterotachy tree lengths:";
    for (int j = 0; j < lenvec.size(); j++)
        out << " " << lenvec[j];
    out << endl;
}

void RateHeterotachy::writeParameters(ostream &out) {
	for (int i = 0; i < ncategory; i++)
		out << "\t" << prop[i];
}

/**
    optimize parameters. Default is to optimize gamma shape
    @return the best likelihood
*/
double RateHeterotachy::optimizeParameters(double gradient_epsilon) {
    if (fix_params) {
        return phylo_tree->computeLikelihood();
    }
    if (verbose_mode >= VB_MED) {
        cout << "Optimizing " << name << " model parameters by EM algorithm..." << endl;
    }
    return optimizeWithEM();
}

double RateHeterotachy::optimizeWithEM() {

    // first compute _pattern_lh_cat
    phylo_tree->computePatternLhCat(WSL_RATECAT);
    size_t nptn = phylo_tree->aln->getNPattern();
    size_t nmix = ncategory;
    
    double *new_prop = aligned_alloc<double>(nmix);
    double *ratio_prop = aligned_alloc<double>(nmix);

    // EM algorithm loop described in Wang, Li, Susko, and Roger (2008)

    for (int step = 0; step < optimize_steps; step++) {
        // E-step

        if (step > 0) {
            // convert _pattern_lh_cat taking into account new weights
            for (size_t ptn = 0; ptn < nptn; ptn++) {
                double *this_lk_cat = phylo_tree->tree_buffers._pattern_lh_cat + ptn*nmix;
                for (size_t c = 0; c < nmix; c++) {
                    this_lk_cat[c] *= ratio_prop[c];
                }
            } 
        }
        memset(new_prop, 0, nmix*sizeof(double));
        for (size_t ptn = 0; ptn < nptn; ptn++) {
            double *this_lk_cat = phylo_tree->tree_buffers._pattern_lh_cat + ptn*nmix;
            double lk_ptn = phylo_tree->ptn_invar[ptn];
            for (size_t c = 0; c < nmix; c++) {
                lk_ptn += this_lk_cat[c];
            }
            ASSERT(lk_ptn != 0.0);
            lk_ptn = phylo_tree->ptn_freq[ptn] / lk_ptn;
            for (size_t c = 0; c < nmix; c++) {
                new_prop[c] += this_lk_cat[c] * lk_ptn;
            }
        } 
        bool converged = true;
        double new_pinvar = 0.0;    
        for (size_t c = 0; c < nmix; c++) {
            new_prop[c] /= phylo_tree->getAlnNSite();
            // Make sure that probabilities do not get zero
            if (new_prop[c] < 1e-10) new_prop[c] = 1e-10;
            // check for convergence
            converged = converged && (fabs(prop[c]-new_prop[c]) < 1e-4);
            ratio_prop[c] = new_prop[c] / prop[c];
            if (std::isnan(ratio_prop[c])) {
                cerr << "BUG: " << new_prop[c] << " " << prop[c] << " " << ratio_prop[c] << endl;
            }
            prop[c] = new_prop[c];
            new_pinvar += prop[c];
        }
        new_pinvar = fabs(1.0 - new_pinvar);
        if (new_pinvar > 1e-6) {
            converged = converged && (fabs(getPInvar()-new_pinvar) < 1e-4);
            // TODO fix p_pinvar
            setPInvar(new_pinvar);
//            phylo_tree->getRate()->setOptimizePInvar(false);
            phylo_tree->computePtnInvar();
            phylo_tree->clearAllPartialLH();
            
        }
        if (converged) break;

    }
    
    aligned_free(ratio_prop);
    aligned_free(new_prop);
//    aligned_free(lk_ptn);
    return phylo_tree->computeLikelihood();

/*
    size_t ptn, c;
    size_t nptn = phylo_tree->aln->getNPattern();
    size_t nmix = ncategory;
    const double MIN_PROP = 1e-4;
    
    double new_prop[nmix];

    // EM algorithm loop described in Wang, Li, Susko, and Roger (2008)

    // first compute _pattern_lh_cat
    double score;
    score = phylo_tree->computePatternLhCat(WSL_RATECAT);

    memset(new_prop, 0, nmix*sizeof(double));
            
    // E-step
    // decoupled weights (prop) from _pattern_lh_cat to obtain L_ci and compute pattern likelihood L_i
    for (ptn = 0; ptn < nptn; ptn++) {
        double *this_lk_cat = phylo_tree->_pattern_lh_cat + ptn*nmix;
        double lk_ptn = phylo_tree->ptn_invar[ptn];
        for (c = 0; c < nmix; c++) {
            lk_ptn += this_lk_cat[c];
        }
        assert(lk_ptn != 0.0);
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

    double sum_prop = 0.0;
    for (c = 0; c < nmix; c++) {
        // check for convergence
        sum_prop += new_prop[c];
        prop[c] = new_prop[c];
        new_pinvar += new_prop[c];
    }

    new_pinvar = 1.0 - new_pinvar;

    if (new_pinvar > 1e-4 && getPInvar() != 0.0) {
        setPInvar(new_pinvar);
//            setOptimizePInvar(false);
        phylo_tree->computePtnInvar();
    }
    
    assert(fabs(sum_prop+new_pinvar-1.0) < MIN_PROP);
    
    // sort the rates in increasing order
    if (sorted_rates) {
        // TODO sort tree lengths per category
    }

    return phylo_tree->computeLikelihood();
*/
}
