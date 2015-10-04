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
}

PhyloTreeMixlen::PhyloTreeMixlen(Alignment *aln, int mixlen) : IQTree(aln) {
	cout << "Initializing heterotachy model with " << mixlen << " mixture branch lengths" << endl;
    setMixlen(mixlen);
    cur_mixture = -1;
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
        for (i = 0; i < mixlen; i++)
            nei->lengths[i] = nei->length * getRate()->getRate(i);

        // assign length of right branch
        nei = (PhyloNeighborMixlen*)((*it)->node->findNeighbor(node));
        nei->lengths.resize(mixlen, nei->length);
        for (i = 0; i < mixlen; i++)
            nei->lengths[i] = nei->length * getRate()->getRate(i);
            
        // recursive call
        initializeMixBranches((PhyloNode*)(*it)->node, node);
    }
}

void PhyloTreeMixlen::assignMixBranches(int category, Node *node, Node *dad) {
    if (!node) node = root;
    FOR_NEIGHBOR_IT(node, dad, it) {
        PhyloNeighborMixlen *nei = (PhyloNeighborMixlen*)(*it);
        assert(category < nei->lengths.size());
        nei->length = nei->lengths[category];
        nei = (PhyloNeighborMixlen*)(*it)->node->findNeighbor(node);
        assert(category < nei->lengths.size());
        nei->length = nei->lengths[category];
        assignMixBranches(category, (*it)->node, node);
    }
}

void PhyloTreeMixlen::copyMixBranches(PhyloTree *tree, int category) {
    NodeVector this_nodes1, this_nodes2, tree_nodes1, tree_nodes2;
    tree->getBranches(tree_nodes1, tree_nodes2, tree->root->neighbors[0]->node);
    getBranches(this_nodes1, this_nodes2, root->neighbors[0]->node);
    for (int i = 0; i < this_nodes1.size(); i++) {
        PhyloNeighborMixlen *nei = (PhyloNeighborMixlen*)this_nodes1[i]->findNeighbor(this_nodes2[i]);
        assert(category < nei->lengths.size());
        nei->lengths[category] = tree_nodes1[i]->findNeighbor(tree_nodes2[i])->length;
        
        nei = (PhyloNeighborMixlen*)this_nodes2[i]->findNeighbor(this_nodes1[i]);
        assert(category < nei->lengths.size());
        nei->lengths[category] = tree_nodes2[i]->findNeighbor(tree_nodes1[i])->length;
    }
}


double PhyloTreeMixlen::optimizeAllBranches(int my_iterations, double tolerance, int maxNRStep) {
    if (((PhyloNeighborMixlen*)root->neighbors[0])->lengths.empty()) {
        // initialize mixture branch lengths if empty

        RateHeterogeneity *saved_rate = getRate();
        bool saved_fused_mix_rate = model_factory->fused_mix_rate;

        // create new rate model
        RateGamma *tmp_rate = new RateGamma(mixlen, -1.0, params->gamma_median, this);
        tmp_rate->setTree(this);
        
        // setup new rate model
        setRate(tmp_rate);
        model_factory->site_rate = tmp_rate;
        if (getModel()->isMixture()) {
            model_factory->fused_mix_rate = true;
            setLikelihoodKernel(sse);
        }

        // optimize rate model
        tmp_rate->optimizeParameters(tolerance);
        
        // assign branch length from rate model
        initializeMixBranches();
        
        // restore rate model
        setRate(saved_rate);
        model_factory->site_rate = saved_rate;
        model_factory->fused_mix_rate = saved_fused_mix_rate;
        if (getModel()->isMixture()) {
            setLikelihoodKernel(sse);
        }
        
        clearAllPartialLH();
    }
    // EM algorithm
    size_t ptn, c;
    size_t nptn = aln->getNPattern();
    size_t nmix = model->getNMixtures();
    assert(nmix == mixlen);

    PhyloTree *tree = new PhyloTree;
    tree->copyPhyloTree(this);
    tree->optimize_by_newton = optimize_by_newton;
    tree->setLikelihoodKernel(sse);
    // initialize model
    ModelFactory *model_fac = new ModelFactory();
    model_fac->joint_optimize = params->optimize_model_rate_joint;

    RateHeterogeneity *site_rate = new RateHeterogeneity; 
    tree->setRate(site_rate);
    site_rate->setTree(tree);
            
    model_fac->site_rate = site_rate;
    tree->model_factory = model_fac;
    tree->setParams(params);
    // first compute _pattern_lh_cat
    double tree_lh;

    if (!getModel()->isMixture())
        tree_lh = computeLikelihoodBranchEigen((PhyloNeighbor*)root->neighbors[0], (PhyloNode*)root); 
    else if (getModelFactory()->fused_mix_rate) {
        outError("Heterotachy with fused mixture not supported");
        tree_lh = computeMixrateLikelihoodBranchEigen((PhyloNeighbor*)root->neighbors[0], (PhyloNode*)root); 
    } else {
        tree_lh = computeMixtureLikelihoodBranchEigen((PhyloNeighbor*)root->neighbors[0], (PhyloNode*)root); 
    }
    cout << "Init LnL = " << tree_lh << endl;

    // E-step
    // decoupled weights (prop) from _pattern_lh_cat to obtain L_ci and compute pattern likelihood L_i
    for (ptn = 0; ptn < nptn; ptn++) {
        double *this_lk_cat = _pattern_lh_cat + ptn*nmix;
        double lk_ptn = 0.0;
        for (c = 0; c < nmix; c++) {
            lk_ptn += this_lk_cat[c];
        }
        lk_ptn = ptn_freq[ptn] / lk_ptn;
        // transform _pattern_lh_cat into posterior probabilities of each category
        for (c = 0; c < nmix; c++) {
            this_lk_cat[c] *= lk_ptn;
        }
        
    } 
    
    double new_tree_lh = 0.0;
    
    // now optimize categories one by one
    for (c = 0; c < nmix; c++) {
        assignMixBranches(c);
        tree->copyPhyloTree(this);
        
        ModelGTR *subst_model;
        if (getModel()->isMixture())
            subst_model = ((ModelMixture*)getModel())->at(c);
        else
            subst_model = (ModelGTR*)getModel();
        tree->setModel(subst_model);
        subst_model->setTree(tree);
        model_fac->model = subst_model;
                    
        // initialize likelihood
        tree->initializeAllPartialLh();
        // copy posterior probability into ptn_freq
        tree->computePtnFreq();
        double *this_lk_cat = _pattern_lh_cat+c;
        for (ptn = 0; ptn < nptn; ptn++)
            tree->ptn_freq[ptn] = this_lk_cat[ptn*nmix];
        
        // optimize branch lengths of mixture category
        tree->optimizeAllBranches(my_iterations, tolerance, maxNRStep);
        
        // copy optimized branch lengths
        copyMixBranches(tree, c);
        
        // reset subst model
        tree->setModel(NULL);
        subst_model->setTree(this);
    }
    
    clearAllPartialLH();
    new_tree_lh = computeLikelihood();
    cout << "Optimized LnL = " << new_tree_lh << endl;
    assert(new_tree_lh >= tree_lh - 0.1);
    
    delete tree;
    return new_tree_lh;
}
