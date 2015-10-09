/*
 * modelfactorymixlen.cpp
 *
 *  Created on: Sep 2, 2015
 *      Author: minh
 */

#include "phylotreemixlen.h"
#include "timeutil.h"
#include "model/modelfactorymixlen.h"
#include "model/modelgtr.h"
#include "model/modelmixture.h"

ModelFactoryMixlen::ModelFactoryMixlen(Params &params, PhyloTree *tree, ModelsBlock *models_block) :   
    ModelFactory(params, tree, models_block) {
    if (!model->isMixture())
        outError("Model is not mixture");
    if (((PhyloTreeMixlen*)tree)->mixlen != model->getNMixtures())
        outError("#mixture categories and #mixture branch lengths do not match");
}

double ModelFactoryMixlen::optimizeParameters(bool fixed_len, bool write_info, double logl_epsilon, double gradient_epsilon) {

	PhyloTreeMixlen *tree = (PhyloTreeMixlen*)site_rate->getTree();
	assert(tree);
    
    tree->initializeMixlen(logl_epsilon);

    return ModelFactory::optimizeParameters(fixed_len, write_info, logl_epsilon, gradient_epsilon);

	assert(model);
	assert(site_rate);

    double defaultEpsilon = logl_epsilon;

	double begin_time = getRealTime();
	double cur_lh;

	stopStoringTransMatrix();

    // EM algorithm
    size_t ptn, c;
    size_t nptn = tree->aln->getNPattern();
    size_t nmix = model->getNMixtures();
    assert(nmix == tree->mixlen);
    tree->print_mix_brlen = false;
    double *new_prop = aligned_alloc<double>(nmix);

	int i;

	for (i = 1; i < tree->params->num_param_iterations; i++) {

        if (!tree->getModel()->isMixture()) {
            outError("Heterotachy must be used with mixture model");
            cur_lh = tree->computeLikelihoodBranchEigen((PhyloNeighbor*)tree->root->neighbors[0], (PhyloNode*)tree->root); 
        } else if (tree->getModelFactory()->fused_mix_rate) {
            outError("Heterotachy with fused mixture not supported");
            cur_lh = tree->computeMixrateLikelihoodBranchEigen((PhyloNeighbor*)tree->root->neighbors[0], (PhyloNode*)tree->root); 
        } else {
            cur_lh = tree->computeMixtureLikelihoodBranchEigen((PhyloNeighbor*)tree->root->neighbors[0], (PhyloNode*)tree->root); 
        }

        tree->setCurScore(cur_lh);
        if (verbose_mode >= VB_MED || write_info) 
            cout << i << ". Log-likelihood: " << cur_lh << endl;

        memset(new_prop, 0, nmix*sizeof(double));

        // E-step
        // decoupled weights (prop) from _pattern_lh_cat to obtain L_ci and compute pattern likelihood L_i
        for (ptn = 0; ptn < nptn; ptn++) {
            double *this_lk_cat = tree->_pattern_lh_cat + ptn*nmix;
            double lk_ptn = 0.0;
            for (c = 0; c < nmix; c++) {
                lk_ptn += this_lk_cat[c];
            }
            assert(lk_ptn != 0.0);
            lk_ptn = tree->ptn_freq[ptn] / lk_ptn;
            // transform _pattern_lh_cat into posterior probabilities of each category
            for (c = 0; c < nmix; c++) {
                this_lk_cat[c] *= lk_ptn;
                new_prop[c] += this_lk_cat[c];
            }
            
        } 

        for (c = 0; c < nmix; c++) {
            new_prop[c] = new_prop[c] / tree->getAlnNSite();
            ((ModelMixture*)tree->getModel())->prop[c] = new_prop[c];
        }

        double new_lh;
        // EM algorithm
        // now optimize categories one by one
        for (c = 0; c < nmix; c++) {
            tree->assignMixBranches(c);
            tree->cat_tree->copyPhyloTree(tree);
            
            ModelGTR *subst_model;
            if (tree->getModel()->isMixture())
                subst_model = ((ModelMixture*)tree->getModel())->at(c);
            else
                subst_model = (ModelGTR*)tree->getModel();
            tree->cat_tree->setModel(subst_model);
            subst_model->setTree(tree->cat_tree);
            tree->cat_tree->getModelFactory()->model = subst_model;
                        
            // initialize likelihood
            tree->cat_tree->initializeAllPartialLh();
            // copy posterior probability into ptn_freq
            tree->cat_tree->computePtnFreq();
            double *this_lk_cat = tree->_pattern_lh_cat+c;
            for (ptn = 0; ptn < nptn; ptn++)
                tree->cat_tree->ptn_freq[ptn] = this_lk_cat[ptn*nmix];
            
            // optimize branch lengths of mixture category
            if (!fixed_len)
                tree->cat_tree->optimizeAllBranches(min(i,3), logl_epsilon/nmix);  // loop only 3 times in total (previously in v0.9.6 5 times)
            if (tree->getModel()->isMixture())
                tree->cat_tree->getModelFactory()->optimizeParametersOnly(gradient_epsilon/nmix);
            
            // copy optimized branch lengths
            tree->copyMixBranches(tree->cat_tree, c);
            
            // reset subst model
            tree->cat_tree->setModel(NULL);
            subst_model->setTree(tree);
        }
        
        tree->assignMeanMixBranches();
        
        tree->clearAllPartialLH();
    
        new_lh = tree->computeLikelihood();
    
		if (verbose_mode >= VB_MED) {
			model->writeInfo(cout);
			site_rate->writeInfo(cout);
		}
		if (new_lh > cur_lh + logl_epsilon) {
            if (Params::getInstance().testAlpha && i == 3) {
                double newEpsilon = (new_lh - cur_lh) * 0.01;
                if (newEpsilon > defaultEpsilon) {
                    logl_epsilon = newEpsilon;
                    cout << "Estimate model parameters with new epsilon = " << logl_epsilon << endl;
                }
            }
			cur_lh = new_lh;
//			if (verbose_mode >= VB_MED || write_info)
//				cout << i << ". Current log-likelihood: " << cur_lh << endl;
		} else {
			site_rate->classifyRates(new_lh);
			if (!fixed_len) cur_lh = tree->optimizeAllBranches(100, logl_epsilon);
				break;
		}
	}

    tree->setCurScore(cur_lh);

	// normalize rates s.t. branch lengths are #subst per site
//    double mean_rate = site_rate->rescaleRates();
//    if (mean_rate != 1.0) {
//		tree->scaleLength(mean_rate);
//		tree->clearAllPartialLH();
//    }
    
	if (verbose_mode >= VB_MED || write_info)
		cout << "Optimal log-likelihood: " << cur_lh << endl;

	if (verbose_mode <= VB_MIN && write_info) {
		model->writeInfo(cout);
		site_rate->writeInfo(cout);
	}
	double elapsed_secs = getRealTime() - begin_time;
	if (write_info)
		cout << "Parameters optimization took " << i-1 << " rounds (" << elapsed_secs << " sec)" << endl << endl;
	startStoringTransMatrix();

    tree->print_mix_brlen = true;
    aligned_free(new_prop);

	return cur_lh;
}
