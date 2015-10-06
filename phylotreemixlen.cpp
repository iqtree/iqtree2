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
    print_mix_brlen = false;
    relative_rate = NULL;
    cat_tree = NULL;
}

PhyloTreeMixlen::PhyloTreeMixlen(Alignment *aln, int mixlen) : IQTree(aln) {
	cout << "Initializing heterotachy model with " << mixlen << " mixture branch lengths" << endl;
    cur_mixture = -1;
    print_mix_brlen = false;
    relative_rate = NULL;
    cat_tree = NULL;
    setMixlen(mixlen);
}

PhyloTreeMixlen::~PhyloTreeMixlen() {
    if (relative_rate)
        delete relative_rate;
    if (cat_tree)
        delete cat_tree;
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

void PhyloTreeMixlen::optimizeOneBranch(PhyloNode *node1, PhyloNode *node2, bool clearLH, int maxNRStep) {
    return PhyloTree::optimizeOneBranch(node1, node2, clearLH, maxNRStep);
    size_t ptn, c;
    size_t nptn = aln->getNPattern();
    size_t nmix = model->getNMixtures();
    assert(nmix == mixlen);

    assert(cat_tree);
    // first compute _pattern_lh_cat
    double tree_lh;

    if (!getModel()->isMixture())
        tree_lh = computeLikelihoodBranchEigen((PhyloNeighbor*)node1->findNeighbor(node2), node1); 
    else if (getModelFactory()->fused_mix_rate) {
        outError("Heterotachy with fused mixture not supported");
        tree_lh = computeMixrateLikelihoodBranchEigen((PhyloNeighbor*)node1->findNeighbor(node2), node1); 
    } else {
        tree_lh = computeMixtureLikelihoodBranchEigen((PhyloNeighbor*)node1->findNeighbor(node2), node1); 
    }
   
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
        cat_tree->copyPhyloTree(this);
        
        ModelGTR *subst_model;
        if (getModel()->isMixture())
            subst_model = ((ModelMixture*)getModel())->at(c);
        else
            subst_model = (ModelGTR*)getModel();
        cat_tree->setModel(subst_model);
        subst_model->setTree(cat_tree);
        cat_tree->getModelFactory()->model = subst_model;
                    
        // initialize likelihood
        cat_tree->initializeAllPartialLh();
        // copy posterior probability into ptn_freq
        cat_tree->computePtnFreq();
        double *this_lk_cat = _pattern_lh_cat+c;
        for (ptn = 0; ptn < nptn; ptn++)
            cat_tree->ptn_freq[ptn] = this_lk_cat[ptn*nmix];
        
        // TODO optimize branch lengths of mixture category
        cat_tree->optimizeOneBranch(NULL, NULL, clearLH, maxNRStep);
        
        // copy optimized branch lengths
        copyMixBranches(cat_tree, c);
        
        // reset subst model
        cat_tree->setModel(NULL);
        subst_model->setTree(this);
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
            relative_rate = new RateFree(mixlen, 0.0, "", true, params->optimize_alg, this);
            relative_rate->setTree(this);
            
            // setup new rate model
            setRate(relative_rate);
            model_factory->site_rate = relative_rate;
            if (getModel()->isMixture()) {
                model_factory->fused_mix_rate = true;
                setLikelihoodKernel(sse);
            }

            // optimize rate model
            relative_rate->optimizeParameters(tolerance);
            // restore rate model
            setRate(saved_rate);
            model_factory->site_rate = saved_rate;
            model_factory->fused_mix_rate = saved_fused_mix_rate;
            if (getModel()->isMixture()) {
                setLikelihoodKernel(sse);
                ModelMixture *mm = (ModelMixture*)getModel();
                for (int i = 0; i < mm->getNMixtures(); i++)
                    mm->prop[i] = relative_rate->getProp(i);
            }
            
        }
        
        // assign branch length from rate model
        initializeMixBranches();
        clearAllPartialLH();
    }

    if (!cat_tree) {
        // set up category tree for EM algorithm
        cat_tree = new PhyloTree;
        cat_tree->copyPhyloTree(this);
        cat_tree->optimize_by_newton = optimize_by_newton;
        cat_tree->setLikelihoodKernel(sse);
        // initialize model
        ModelFactory *model_fac = new ModelFactory();
        model_fac->joint_optimize = params->optimize_model_rate_joint;
        RateHeterogeneity *site_rate = new RateHeterogeneity; 
        cat_tree->setRate(site_rate);
        site_rate->setTree(cat_tree);
        model_fac->site_rate = site_rate;
        cat_tree->model_factory = model_fac;
        cat_tree->setParams(params);
    }

}

double PhyloTreeMixlen::optimizeAllBranches(int my_iterations, double tolerance, int maxNRStep) {

	initializeMixlen(tolerance);

    // EM algorithm
    size_t ptn, c;
    size_t nptn = aln->getNPattern();
    size_t nmix = model->getNMixtures();
    assert(nmix == mixlen);

    print_mix_brlen = false;

            
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
//    cout << "Init LnL = " << tree_lh << endl;

    // E-step
    // decoupled weights (prop) from _pattern_lh_cat to obtain L_ci and compute pattern likelihood L_i
    for (ptn = 0; ptn < nptn; ptn++) {
        double *this_lk_cat = _pattern_lh_cat + ptn*nmix;
        double lk_ptn = 0.0;
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
    
    double new_tree_lh = 0.0;
    
    // now optimize categories one by one
    for (c = 0; c < nmix; c++) {
        assignMixBranches(c);
        cat_tree->copyPhyloTree(this);
        
        ModelGTR *subst_model;
        if (getModel()->isMixture())
            subst_model = ((ModelMixture*)getModel())->at(c);
        else
            subst_model = (ModelGTR*)getModel();
        cat_tree->setModel(subst_model);
        subst_model->setTree(cat_tree);
        cat_tree->getModelFactory()->model = subst_model;
                    
        // initialize likelihood
        cat_tree->initializeAllPartialLh();
        // copy posterior probability into ptn_freq
        cat_tree->computePtnFreq();
        double *this_lk_cat = _pattern_lh_cat+c;
        for (ptn = 0; ptn < nptn; ptn++)
            cat_tree->ptn_freq[ptn] = this_lk_cat[ptn*nmix];
        
        // optimize branch lengths of mixture category
        cat_tree->optimizeAllBranches(my_iterations, tolerance, maxNRStep);
        
        // copy optimized branch lengths
        copyMixBranches(cat_tree, c);
        
        // reset subst model
        cat_tree->setModel(NULL);
        subst_model->setTree(this);
    }
    
    assignMeanMixBranches();
    
    clearAllPartialLH();
    new_tree_lh = computeLikelihood();
//    cout << "Optimized LnL = " << new_tree_lh << endl;
    assert(new_tree_lh >= tree_lh - 0.1);
    
//    delete tree;
    
    print_mix_brlen = true;
    
    return new_tree_lh;
}

void PhyloTreeMixlen::printBranchLength(ostream &out, int brtype, bool print_slash, Neighbor *length_nei) {
    if (!print_mix_brlen || ((PhyloNeighborMixlen*)length_nei)->lengths.empty())
        return PhyloTree::printBranchLength(out, brtype, print_slash, length_nei);

    PhyloNeighborMixlen *nei = (PhyloNeighborMixlen*) length_nei;
    if (brtype & WT_BR_LEN) 
        out << ":";
    else if ((brtype & WT_BR_CLADE) && print_slash)
        out << "/";
        
    for (int i = 0; i < mixlen; i++) {
        if (i > 0) out << "_";
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
