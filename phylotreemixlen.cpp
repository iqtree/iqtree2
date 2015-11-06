//
//  phylotreemixlen.cpp
//  iqtree
//
//  Created by Minh Bui on 24/08/15.
//
//

#include "phylotreemixlen.h"
#include "phylonodemixlen.h"
#include "model/modelmixture.h"
#include "model/ratefree.h"

PhyloTreeMixlen::PhyloTreeMixlen() : IQTree() {
	mixlen = 1;
    cur_mixture = -1;
    relative_rate = NULL;
}

PhyloTreeMixlen::PhyloTreeMixlen(Alignment *aln, int mixlen) : IQTree(aln) {
	cout << "Initializing heterotachy model with " << mixlen << " mixture branch lengths" << endl;
    cur_mixture = -1;
    relative_rate = NULL;
    setMixlen(mixlen);
}

PhyloTreeMixlen::~PhyloTreeMixlen() {
    if (relative_rate)
        delete relative_rate;
}

Node* PhyloTreeMixlen::newNode(int node_id, const char* node_name) {
    return (Node*) (new PhyloNodeMixlen(node_id, node_name));
}

Node* PhyloTreeMixlen::newNode(int node_id, int node_name) {
    return (Node*) (new PhyloNodeMixlen(node_id, node_name));
}

void PhyloTreeMixlen::setMixlen(int mixlen) {
	this->mixlen = mixlen;
}


void PhyloTreeMixlen::initializeMixBranches(PhyloNode *node, PhyloNode *dad) {
    if (!node) {
        node = (PhyloNode*)root;
        // exit if already initialized
        if (!((PhyloNeighborMixlen*)root->neighbors[0])->lengths.empty())
            return;
    }
    int i;
    FOR_NEIGHBOR_IT(node, dad, it) {
        // assign length of left branch
        PhyloNeighborMixlen *nei = (PhyloNeighborMixlen*)(*it);
        nei->lengths.resize(mixlen, nei->length);
        assert(nei->length >= 0);
        for (i = 0; i < mixlen; i++)
            nei->lengths[i] = max(MIN_BRANCH_LEN, nei->length * relative_rate->getRate(i));

        // assign length of right branch
        nei = (PhyloNeighborMixlen*)((*it)->node->findNeighbor(node));
        nei->lengths.resize(mixlen, nei->length);
        for (i = 0; i < mixlen; i++)
            nei->lengths[i] = max(MIN_BRANCH_LEN, nei->length * relative_rate->getRate(i));
            
        // recursive call
        initializeMixBranches((PhyloNode*)(*it)->node, node);
    }
}

void PhyloTreeMixlen::assignMeanMixBranches(Node *node, Node *dad) {
    if (!node) node = root;
    FOR_NEIGHBOR_IT(node, dad, it) {
        PhyloNeighborMixlen *nei = (PhyloNeighborMixlen*)(*it);
        double mean_len = 0.0;
        for (int i = 0; i < nei->lengths.size(); i++)
            mean_len += nei->lengths[i];
        mean_len /= nei->lengths.size();
        nei->length = mean_len;
        
        nei = (PhyloNeighborMixlen*)(*it)->node->findNeighbor(node);
        nei->length = mean_len;
        assignMeanMixBranches((*it)->node, node);
    }
}

void PhyloTreeMixlen::initializeMixlen(double tolerance) {
    if (((PhyloNeighborMixlen*)root->neighbors[0])->lengths.empty()) {
        // initialize mixture branch lengths if empty

        if (!relative_rate) {
            RateHeterogeneity *saved_rate = getRate();
            bool saved_fused_mix_rate = model_factory->fused_mix_rate;

            // create new rate model
            // random alpha
//            relative_rate = new RateGamma(mixlen, 0.0, params->gamma_median, this);
            relative_rate = new RateFree(mixlen, params->gamma_shape, "", true, params->optimize_alg, this);
            relative_rate->setTree(this);
            
            // setup new rate model
            setRate(relative_rate);
            model_factory->site_rate = relative_rate;
            if (getModel()->isMixture()) {
                model_factory->fused_mix_rate = true;
                setLikelihoodKernel(sse);
            }

            // optimize rate model
            double tree_lh = relative_rate->optimizeParameters(tolerance);
            cout << "tree_lh = " << tree_lh << endl;
            // restore rate model
            setRate(saved_rate);
            model_factory->site_rate = saved_rate;
            model_factory->fused_mix_rate = saved_fused_mix_rate;
            if (getModel()->isMixture()) {
                setLikelihoodKernel(sse);
                ModelMixture *mm = (ModelMixture*)getModel();
                double pinvar = site_rate->getPInvar();
                if (!mm->fix_prop)
                    for (int i = 0; i < mm->getNMixtures(); i++)
                        mm->prop[i] = relative_rate->getProp(i)*(1.0-pinvar);
            }
            
        }
        
        // assign branch length from rate model
        initializeMixBranches();
        clearAllPartialLH();
    }

}

void PhyloTreeMixlen::optimizeOneBranch(PhyloNode *node1, PhyloNode *node2, bool clearLH, int maxNRStep) {
    current_it = (PhyloNeighbor*) node1->findNeighbor(node2);
    assert(current_it);
    current_it_back = (PhyloNeighbor*) node2->findNeighbor(node1);
    assert(current_it_back);

    size_t ptn, c;
    size_t nptn = aln->getNPattern();
    size_t nmix = model->getNMixtures();
    assert(nmix == mixlen);

    // first compute _pattern_lh_cat
    double tree_lh;
    if (!getModel()->isMixture())
        tree_lh = computeLikelihoodBranchEigen(current_it, node1); 
    else if (getModelFactory()->fused_mix_rate) {
        outError("Heterotachy with fused mixture not supported");
        tree_lh = computeMixrateLikelihoodBranchEigen(current_it, node1); 
    } else {
        tree_lh = computeMixtureLikelihoodBranchEigen(current_it, node1); 
    }
//    cout << "Init LnL = " << tree_lh << endl;

    // E-step
    // decoupled weights (prop) from _pattern_lh_cat to obtain L_ci and compute pattern likelihood L_i
    for (ptn = 0; ptn < nptn; ptn++) {
        double *this_lk_cat = _pattern_lh_cat + ptn*nmix;
        double lk_ptn = ptn_invar[ptn];
        for (c = 0; c < nmix; c++) {
            lk_ptn += this_lk_cat[c];
        }
        assert(lk_ptn != 0.0);
        lk_ptn = ptn_freq[ptn] / lk_ptn;
        // transform _pattern_lh_cat into posterior probabilities of each category
        for (c = 0; c < nmix; c++) {
            this_lk_cat[c] *= lk_ptn;
        }
        
    } 
 
    double negative_lh;
    double optx;
    theta_computed = false;
    computePtnFreq();
    
    for (cur_mixture = 0; cur_mixture < mixlen; cur_mixture++) {

        double *this_lk_cat = _pattern_lh_cat+cur_mixture;
        for (ptn = 0; ptn < nptn; ptn++)
            ptn_freq[ptn] = this_lk_cat[ptn*nmix];
        
        double current_len = current_it->getLength(cur_mixture);
        assert(current_len >= 0.0);
        // Newton-Raphson method
        optx = minimizeNewton(MIN_BRANCH_LEN, current_len, MAX_BRANCH_LEN, TOL_BRANCH_LEN, negative_lh, maxNRStep);

        current_it->setLength(cur_mixture, optx);
        current_it_back->setLength(cur_mixture, optx);
    
    }
    
    cur_mixture = -1;
    // reset ptn_freq
    ptn_freq_computed = false;
    computePtnFreq();

    if (clearLH) {
        node1->clearReversePartialLh(node2);
        node2->clearReversePartialLh(node1);
    }

}


double PhyloTreeMixlen::optimizeAllBranches(int my_iterations, double tolerance, int maxNRStep) {

	initializeMixlen(tolerance);    
    clearAllPartialLH();
    
    double tree_lh = PhyloTree::optimizeAllBranches(my_iterations, tolerance, maxNRStep);    
    assignMeanMixBranches();
    
    return tree_lh;
}

void PhyloTreeMixlen::printBranchLength(ostream &out, int brtype, bool print_slash, Neighbor *length_nei) {
    if (((PhyloNeighborMixlen*)length_nei)->lengths.empty())
        return PhyloTree::printBranchLength(out, brtype, print_slash, length_nei);

    PhyloNeighborMixlen *nei = (PhyloNeighborMixlen*) length_nei;
    if (brtype & WT_BR_LEN) 
        out << ":";
    else if ((brtype & WT_BR_CLADE) && print_slash)
        out << "/";
        
    for (int i = 0; i < mixlen; i++) {
        if (i > 0) out << BRANCH_LENGTH_SEPARATOR;
        double length = nei->lengths[i];
        if (brtype & WT_BR_SCALE) length *= len_scale;
        if (brtype & WT_BR_LEN_ROUNDING) length = round(length);
        if (brtype & WT_BR_LEN) {
            if (brtype & WT_BR_LEN_FIXED_WIDTH)
                out << fixed << length;
            else
                out << length;
        } else if (brtype & WT_BR_CLADE) {
            out << length;
        }
    }
}

void PhyloTreeMixlen::computeFuncDerv(double value, double &df, double &ddf) {
    current_it->setLength(cur_mixture, value);
    current_it_back->setLength(cur_mixture, value);

    PhyloNeighbor* dad_branch = current_it;
    PhyloNode *dad = (PhyloNode*) current_it_back->node;

    PhyloNode *node = (PhyloNode*) dad_branch->node;
    PhyloNeighbor *node_branch = (PhyloNeighbor*) node->findNeighbor(dad);
    if (!central_partial_lh)
        initializeAllPartialLh();
    if (node->isLeaf()) {
    	PhyloNode *tmp_node = dad;
    	dad = node;
    	node = tmp_node;
    	PhyloNeighbor *tmp_nei = dad_branch;
    	dad_branch = node_branch;
    	node_branch = tmp_nei;
    }
    
    assert(dad_branch->partial_lh_computed & 1);
    assert(node_branch->partial_lh_computed & 1);
//    if ((dad_branch->partial_lh_computed & 1) == 0)
//        computePartialLikelihood(dad_branch, dad);
//    if ((node_branch->partial_lh_computed & 1) == 0)
//        computePartialLikelihood(node_branch, node);
    size_t nstates = aln->num_states;
    size_t ncat = site_rate->getNRate();
    size_t nmixture = model->getNMixtures();

    size_t block = ncat * nstates * nmixture;
    size_t statemix = nstates * nmixture;
    size_t statecat = nstates * ncat;
    size_t ptn; // for big data size > 4GB memory required
    size_t c, i, m = cur_mixture;
    size_t orig_nptn = aln->size();
    size_t nptn = aln->size()+model_factory->unobserved_ptns.size();
    double *eval = model->getEigenvalues();
    assert(eval);

	assert(theta_all);
	if (!theta_computed) {
		// precompute theta for fast branch length optimization

	    if (dad->isLeaf()) {
	    	// special treatment for TIP-INTERNAL NODE case
#ifdef _OPENMP
#pragma omp parallel for private(ptn, i, m)
#endif
	    	for (ptn = 0; ptn < nptn; ptn++) {
				double *partial_lh_dad = dad_branch->partial_lh + ptn*block;
				double *theta = theta_all + ptn*block;
				double *lh_tip = tip_partial_lh +
						((int)((ptn < orig_nptn) ? (aln->at(ptn))[dad->id] :  model_factory->unobserved_ptns[ptn-orig_nptn]))*statemix;
				for (m = 0; m < nmixture; m++) {
					for (i = 0; i < statecat; i++) {
						theta[m*statecat+i] = lh_tip[m*nstates + i%nstates] * partial_lh_dad[m*statecat+i];
					}
				}

			}
			// ascertainment bias correction
	    } else {
	    	// both dad and node are internal nodes
		    double *partial_lh_node = node_branch->partial_lh;
		    double *partial_lh_dad = dad_branch->partial_lh;

	    	size_t all_entries = nptn*block;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	    	for (i = 0; i < all_entries; i++) {
				theta_all[i] = partial_lh_node[i] * partial_lh_dad[i];
			}
	    }
		theta_computed = true;
	}

    double *val0 = new double[statecat];
    double *val1 = new double[statecat];
    double *val2 = new double[statecat];
	for (c = 0; c < ncat; c++) {
		double prop = site_rate->getProp(c);
        for (i = 0; i < nstates; i++) {
            double cof = eval[cur_mixture*nstates+i]*site_rate->getRate(c);
            // length for heterotachy model
            double val = exp(cof*dad_branch->getLength(cur_mixture)) * prop * ((ModelMixture*)model)->prop[cur_mixture];
            double val1_ = cof*val;
            val0[(c)*nstates+i] = val;
            val1[(c)*nstates+i] = val1_;
            val2[(c)*nstates+i] = cof*val1_;
		}
	}


    double my_df = 0.0, my_ddf = 0.0, prob_const = 0.0, df_const = 0.0, ddf_const = 0.0;

#ifdef _OPENMP
#pragma omp parallel for reduction(+: my_df, my_ddf, prob_const, df_const, ddf_const) private(ptn, i)
#endif
    for (ptn = 0; ptn < nptn; ptn++) {
		double lh_ptn = ptn_invar[ptn], df_ptn = 0.0, ddf_ptn = 0.0;
		double *theta = theta_all + ptn*block + cur_mixture*statecat;
		for (i = 0; i < statecat; i++) {
			lh_ptn += val0[i] * theta[i];
			df_ptn += val1[i] * theta[i];
			ddf_ptn += val2[i] * theta[i];
		}

//        assert(lh_ptn > 0.0);
        lh_ptn = fabs(lh_ptn);

        if (ptn < orig_nptn) {
			double df_frac = df_ptn / lh_ptn;
			double ddf_frac = ddf_ptn / lh_ptn;
			double freq = ptn_freq[ptn];
			double tmp1 = df_frac * freq;
			double tmp2 = ddf_frac * freq;
			my_df += tmp1;
			my_ddf += tmp2 - tmp1 * df_frac;
		} else {
			// ascertainment bias correction
			prob_const += lh_ptn;
			df_const += df_ptn;
			ddf_const += ddf_ptn;
		}
    }
	df = my_df;
	ddf = my_ddf;
    if (isnan(df) || isinf(df)) {
        df = 0.0;
        ddf = 0.0;
//        outWarning("Numerical instability (some site-likelihood = 0)");
    }

	if (orig_nptn < nptn) {
    	// ascertainment bias correction
    	prob_const = 1.0 - prob_const;
    	double df_frac = df_const / prob_const;
    	double ddf_frac = ddf_const / prob_const;
    	int nsites = aln->getNSite();
    	df += nsites * df_frac;
    	ddf += nsites *(ddf_frac + df_frac*df_frac);
    }


    delete [] val2;
    delete [] val1;
    delete [] val0;


	df = -df;
    ddf = -ddf;

//    return lh;
}
