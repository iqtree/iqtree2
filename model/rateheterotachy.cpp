//
//  rateheterotachy.cpp
//  iqtree
//
//  Created by Minh Bui on 11/8/16.
//
//

#include "phylotree.h"
#include "rateheterotachy.h"

RateHeterotachy::RateHeterotachy(int ncat, string params, bool sorted_rates, PhyloTree *tree) : RateHeterogeneity() {
    phylo_tree = tree;
    prop = NULL;
    fix_params = false;
    this->sorted_rates = sorted_rates;
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
		fix_params = true;
	} catch (string &str) {
		outError(str);
	}
}

/**
    destructor
*/
RateHeterotachy::~RateHeterotachy() {
    if (prop)
        delete [] prop;
    prop = NULL;
}

void RateHeterotachy::setNCategory(int ncat) {
    ncategory = ncat;
    // initialize with gamma rates
	if (prop) delete [] prop;
	prop  = new double[ncategory];

    int i;
	for (i = 0; i < ncategory; i++)
        prop[i] = (1.0-getPInvar())/ncategory;
    
	name = "+H";
	name += convertIntToString(ncategory);
	full_name = "Rate heterotachy";
	full_name += " with " + convertIntToString(ncategory) + " categories";
}


/**
    save object into the checkpoint
*/
void RateHeterotachy::saveCheckpoint() {
    checkpoint->startStruct("RateHeterotachy");
    CKP_ARRAY_SAVE(ncategory, prop);
    checkpoint->endStruct();
    RateHeterogeneity::saveCheckpoint();
}

/**
    restore object from the checkpoint
*/
void RateHeterotachy::restoreCheckpoint() {
    RateHeterogeneity::restoreCheckpoint();
    checkpoint->startStruct("RateHeterotachy");
    CKP_ARRAY_RESTORE(ncategory, prop);
    checkpoint->endStruct();
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
	out << "Heterotachy weights: ";
	for (int i = 0; i < ncategory; i++)
		out << " " << prop[i];
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
    if (fix_params)
        return phylo_tree->computeLikelihood();
	if (verbose_mode >= VB_MED)
		cout << "Optimizing " << name << " model parameters by EM algorithm..." << endl;

    return optimizeWithEM();
}

double RateHeterotachy::optimizeWithEM() {

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

}
