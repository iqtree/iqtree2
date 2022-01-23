//
//  phylotreemixlen.cpp
//  iqtree
//
//  Created by Minh Bui on 24/08/15.
//
//

#include "phylotreemixlen.h"
#include "phylonodemixlen.h"

#include <model/modelfactorymixlen.h>
#include <model/modeldivergent.h>
#include <model/modelmixture.h>
#include <model/ratefree.h>
#include "utils/MPIHelper.h"

#ifdef _MSC_VER
#include <boost/scoped_array.hpp>
#endif

#ifdef USE_CPPOPTLIB
#include "cppoptlib/solver/newtondescentsolver.h"
#include "cppoptlib/solver/lbfgsbsolver.h"
#endif

PhyloTreeMixlen::PhyloTreeMixlen() : IQTree()
#ifdef USE_CPPOPTLIB
, cppoptlib::BoundedProblem<double>()
#endif
    ,mixlen(1), cur_mixture(-1), initializing_mixlen(false)
{
}

PhyloTreeMixlen::PhyloTreeMixlen(Alignment *aln, int mixlen) 
    : IQTree(aln)
#ifdef USE_CPPOPTLIB
    , cppoptlib::BoundedProblem<double>(mixlen)
#endif
    , cur_mixture(-1), initializing_mixlen(false)
{
    setMixlen(mixlen);
}

PhyloTreeMixlen::~PhyloTreeMixlen() {
}

void PhyloTreeMixlen::startCheckpoint() {
    if (mixlen > 0)
        checkpoint->startStruct("PhyloTreeMixlen" + convertIntToString(getMixlen()));
    else
        PhyloTree::startCheckpoint();
}

void PhyloTreeMixlen::saveCheckpoint() {
    if (mixlen > 0) {
        startCheckpoint();
        if (this->relative_treelen.size() > 0) {
            ASSERT(mixlen == this->relative_treelen.size());
#ifndef _MSC_VER
            double relative_treelen[mixlen];
#else
            boost::scoped_array<double> relative_treelen(new double[mixlen]);
#endif
            for (int i = 0; i < mixlen; i++) {
                relative_treelen[i] = this->relative_treelen[i];
            }
            CKP_ARRAY_SAVE(mixlen, &relative_treelen[0]);
        }
        endCheckpoint();
    }
    IQTree::saveCheckpoint();
}

/** 
    restore object from the checkpoint
*/
void PhyloTreeMixlen::restoreCheckpoint() {
    if (mixlen > 0) {
        startCheckpoint();
#ifndef _MSC_VER
        double relative_treelen[mixlen];
#else
        boost::scoped_array<double> relative_treelen(new double[mixlen]);
#endif
        if (CKP_ARRAY_RESTORE(mixlen, &relative_treelen[0])) {
            this->relative_treelen.resize(mixlen);
            for (int i = 0; i < mixlen; i++)
                this->relative_treelen[i] = relative_treelen[i];
        }
        endCheckpoint();
    }
    IQTree::restoreCheckpoint();
    if (!root) {
        // if not success, try to restore from PhyloTree
        int orig_mixlen = mixlen;
        mixlen = 0;
        PhyloTree::restoreCheckpoint();
        mixlen = orig_mixlen;
    }
}

PhyloNodeMixlen* PhyloTreeMixlen::getRoot() const {
    return dynamic_cast<PhyloNodeMixlen*>(root);
}

PhyloNodeMixlen* PhyloTreeMixlen::newNode(int node_id, const char* node_name) {
    return new PhyloNodeMixlen(node_id, node_name);
}

PhyloNodeMixlen* PhyloTreeMixlen::newNode(int node_id, int node_name) {
    return new PhyloNodeMixlen(node_id, node_name);
}

void PhyloTreeMixlen::setMixlen(int mixlen) {
	this->mixlen = mixlen;
}

void PhyloTreeMixlen::readTreeString(const string &tree_string, 
                                     bool nodes_have_name) {
    IQTree::readTreeString(tree_string, nodes_have_name);
    treeLengths(relative_treelen);
    if (mixlen > 0 && relative_treelen[0] == 0.0) {
        relative_treelen.clear();
    }
}

void PhyloTreeMixlen::initializeModel(Params &params, string model_name,
                                      ModelsBlock *models_block,
                                      PhyloTree* report_to_tree) {
    try {
        if (!getModelFactory()) {
            setModelFactory(new ModelFactoryMixlen(params, model_name, this,
                                                   models_block, report_to_tree));
        }
    } catch (string & str) {
        outError(str);
    }
    IQTree::initializeModel(params, model_name,
                            models_block, report_to_tree);
}

void PhyloTreeMixlen::treeLengths(DoubleVector &lenvec, 
                                  Node *node, Node *dad) {
    if (lenvec.empty()) {
        lenvec.resize(mixlen, 0.0);
    }
    if (!node) {
        node = root;
    }
    FOR_NEIGHBOR_IT(node, dad, it) {
        treeLengths(lenvec, (*it)->node, node);
        for (int i = 0; i < mixlen; i++) {
            lenvec[i] += (*it)->getLength(i);
        }
    }
}

#define FOR_EACH_PHYLO_NEIGHBOR_MIXLEN(mynode, mydad, it, nei) \
    for (PhyloNeighborMixlen* nei=nullptr, *nei2x=dynamic_cast<PhyloNodeMixlen*>(mynode)->firstNeighbor(); nei2x!=nullptr ; nei2x=nullptr) \
        for (NeighborVec::iterator it = (mynode)->neighbors.begin(); it != (mynode)->neighbors.end(); ++it) \
                if ((nei = dynamic_cast<PhyloNeighborMixlen*>(*it)) && nei->getNode() != (mydad) )




void PhyloTreeMixlen::initializeMixBranches(PhyloNodeMixlen* node, 
                                            PhyloNodeMixlen* dad) {
    if (!node) {
        node = getRoot();
        // exit if already initialized
        // if (! getRoot()->firstNeighbor()->lengths.empty() )
        //     return;
    }
    FOR_EACH_PHYLO_NEIGHBOR_MIXLEN(node, dad, it, nei) {
        // assign length of left branch
        PhyloNeighborMixlen *back_nei = nei->getNode()->findNeighbor(node);
        if (nei->lengths.empty()) {
            // no branch lengths, initialize with relative rates
            ASSERT(nei->length >= 0);
            nei->lengths.resize(mixlen, nei->length);
            back_nei->lengths.resize(mixlen, back_nei->length);
            for (int i = 0; i < mixlen; i++) {
                auto max_length = max(params->min_branch_length, 
                                      nei->length * relative_treelen[i]);
                nei->lengths[i] = back_nei->lengths[i] = max_length;
            }
        } else if (nei->lengths.size() > mixlen) {
            // too many lengths, cut down
            nei->lengths.resize(mixlen);
            back_nei->lengths.resize(mixlen);
        } else {
            // too few lengths, add more
            int cur = static_cast<int>(nei->lengths.size());
            nei->lengths.resize(mixlen, nei->length);
            back_nei->lengths.resize(mixlen, back_nei->length);
            double avglen = 0.0;
            for (int i = 0; i < cur; i++) {
                avglen += nei->lengths[i];
            }
            avglen /= cur;
            for (int i = cur; i < mixlen; i++) {
                auto max_length = max(params->min_branch_length, 
                                      avglen * relative_treelen[i]);
                nei->lengths[i] = back_nei->lengths[i] = max_length;
            }
        }

        double mean_len = 0.0;
        for (int i = 0; i < mixlen; i++) {
            mean_len += nei->lengths[i] 
                      * site_rate->getProp(i);
        }
        // mean_len /= mixlen;
        nei->length = back_nei->length = mean_len;

        // recursive call
        initializeMixBranches(nei->getNode(), node);
    }
}

void PhyloTreeMixlen::assignMeanMixBranches(PhyloNodeMixlen* node, 
                                            PhyloNodeMixlen* dad) {
    if (!node) {
        node = getRoot();
    }
    RateHeterogeneity* other_rate = nullptr;
    auto rate_model = getRateModelForBranch(node, dad, other_rate);
    FOR_EACH_PHYLO_NEIGHBOR_MIXLEN(node, dad, it, nei) {
        double mean_len = 0.0;
        for (int i = 0; i < nei->lengths.size(); ++i) {
            mean_len += nei->lengths[i] * rate_model->getProp(i);
        }
        //mean_len /= nei->lengths.size();
        nei->length      = mean_len;
        auto back_nei    = nei->getNode()->findNeighbor(node);
        back_nei->length = mean_len;
        assignMeanMixBranches(nei->getNode(), node);
    }
}

void PhyloTreeMixlen::initializeMixlen(double tolerance, bool write_info,
                                       PhyloTree* report_to_tree) {
    // initialize mixture branch lengths if empty
    if (initializing_mixlen) {
        return;
    }
    initializing_mixlen = true;

    if (relative_treelen.empty()) {
        RateHeterogeneity* saved_rate = getRate();
        bool saved_fused_mix_rate = model_factory->fused_mix_rate;

        string param;
        if (getRate()->getFixParams()) {
            stringstream ss;
            const char* separator = "";
            for (int i = 0; i < mixlen; i++) {
                ss << separator << getRate()->getProp(i);
                separator = ",";
            }
            param = ss.str();
        }
        RateFree relative_rate(mixlen, params->gamma_shape, param, 
                               false, params->optimize_alg_freerate, 
                               this);
        relative_rate.setTree(this);
        
        // setup new rate model
        setRate(&relative_rate);
        model_factory->site_rate = &relative_rate;
        if (getModel()->isMixture()) {
            setLikelihoodKernel(sse);
        }

        // optimize rate model
        /*double tree_lh =*/ (void) relative_rate.optimizeParameters
                                    (tolerance, report_to_tree);

        //model_factory->optimizeParameters(params->fixed_branch_length,
        //                                  false, tolerance, PhyloTree* report_to_tree);

        //optimizeModelParameters();

        // 2016-07-22: BUGFIX should rescale rates
        double mean_rate = relative_rate.rescaleRates();
        if (fabs(mean_rate-1.0) > 1e-6 && 
            params->fixed_branch_length != BRLEN_FIX) {
            // scaleLength(mean_rate);
        }
        if (write_info) {
            cout << "Initial LogL: " << curScore << ", ";
            relative_rate.writeInfo(cout);
        }
        // make the rates more distinct
        auto first_rate = relative_rate.getRate(0);
        auto last_rate  = relative_rate.getRate(mixlen-1);
        if (mixlen > 1 && first_rate / last_rate > 0.9) {
            cout << "Making the rates more distinct..." << endl;
            relative_rate.setRate(0, first_rate*0.95);
            relative_rate.setRate(mixlen-1, last_rate*1.05);
        }
        double treelen = treeLength();
        relative_treelen.resize(mixlen);
        // 2017-12-21: BUG: moved this out of write_info if
        // otherwise, relative_treelen is not initialized
        for (int i = 0; i < mixlen; i++) {
            relative_treelen[i] = treelen * relative_rate.getRate(i);
        }
        if (write_info) {
            cout << "relative_treelen:";
            for (int i = 0; i < mixlen; i++) {
                cout << " " << relative_treelen[i];
            }
            cout << endl;
        }
        // restore rate model
        setRate(saved_rate);
        model_factory->site_rate      = saved_rate;
        model_factory->fused_mix_rate = saved_fused_mix_rate;
        setLikelihoodKernel(sse);

        // set the weights of heterotachy model
        double pinvar = site_rate->getPInvar();
        if (!site_rate->getFixParams()) {
            for (int i = 0; i < mixlen; ++i) {
                site_rate->setProp(i, relative_rate.getProp(i)*(1.0-pinvar));
            }
        }       
        clearAllPartialLH();
        clearAllPartialParsimony(false);
    }

    if (getRoot()->firstNeighbor()->lengths.size() != mixlen) {
        // assign branch length from rate model
        DoubleVector saved_treelen = relative_treelen;
        DoubleVector lenvec;
        treeLengths(lenvec);
        for (int i = 0; i < mixlen; i++) {
            relative_treelen[i] = relative_treelen[i] / lenvec[i];
        }
        if (verbose_mode >= VerboseMode::VB_MED) {
            cout << "relative_ratio:";
            for (int i = 0; i < mixlen; i++) {
                cout << " " << relative_treelen[i];
            }
            cout << endl;
        }
        initializeMixBranches();
        clearAllPartialLH();
        clearAllPartialParsimony(false);
        relative_treelen = saved_treelen;
    }
    initializing_mixlen = false;
}

void PhyloTreeMixlen::fixOneNegativeBranch
        ( double               branch_length, 
          Neighbor* dad_branch, 
          Node*     dad) {
    PhyloTree::fixOneNegativeBranch(branch_length, dad_branch, dad);
    /*
    if (dad_branch->lengths.empty())
        return;
    int i;
    for (i = 0; i < dad_branch->lengths.size(); i++) {
        dad_branch->lengths[i] = branch_length * relative_rate.getRate(i);
    }
    PhyloNeighborMixlen *br = dad_branch->getNode()->findNeighbor(dad);
    for (i = 0; i < br->lengths.size(); i++)
        br->lengths[i] = branch_length * relative_rate.getRate(i);
    */
}

void PhyloTreeMixlen::optimizeOneBranch
        ( PhyloNode* node1, PhyloNode* node2, 
          bool clearLH, int maxNRStep) {
    if (initializing_mixlen) {
        return PhyloTree::optimizeOneBranch(node1, node2, clearLH, maxNRStep);
    }
    current_it =  node1->findNeighbor(node2);
    ASSERT(current_it);
    current_it_back = node2->findNeighbor(node1);
    ASSERT(current_it_back);
    tree_buffers.theta_computed = false;

#ifdef USE_CPPOPTLIB
    auto alg = params->optimize_alg_mixlen;
    if (contains(alg,"cppopt")) {
        //----- using cppoptlib ------//

        TVector lower_bound(mixlen), upper_bound(mixlen), variables(mixlen);

    //    variables.resize(mixlen);
        for (i = 0; i < mixlen; i++) {
            lower_bound[i] = params->min_branch_length;
            variables[i] = current_it->getLength(i);
            upper_bound[i] = params->max_branch_length;
        }

        setBoxConstraint(lower_bound, upper_bound);

        cppoptlib::NewtonDescentSolver<PhyloTreeMixlen> solver;
    //    cppoptlib::LbfgsbSolver<PhyloTreeMixlen> solver;
        solver.minimize(*this, variables);
        for (i = 0; i < mixlen; i++) {
            current_it->setLength(i, variables[i]);
            current_it_back->setLength(i, variables[i]);
        }
    } else
#endif

    auto alg = params->optimize_alg_mixlen;
    if (contains(alg,"newton")) {

        //----- Newton-Raphson -----//
#ifndef _MSC_VER
        double lower_bound[mixlen];
        double upper_bound[mixlen];
        double variables[mixlen];
#else
        boost::scoped_array<double> lower_bound(new double[mixlen]);
        boost::scoped_array<double> upper_bound(new double[mixlen]);
        boost::scoped_array<double> variables(new double[mixlen]);
#endif
        for (int i = 0; i < mixlen; i++) {
            lower_bound[i] = params->min_branch_length;
            variables[i]   = current_it->getLength(i);
            upper_bound[i] = params->max_branch_length;
        }
        minimizeNewtonMulti(&lower_bound[0], &variables[0], 
                            &upper_bound[0], params->min_branch_length, 
                            mixlen);
        for (int i = 0; i < mixlen; i++) {
            current_it->setLength(i, variables[i]);
            current_it_back->setLength(i, variables[i]);
        }
    } else if (params->optimize_alg_mixlen.find("BFGS") != string::npos) {

        // BFGS method to simultaneously 
        // optimize all lengths per branch
        // It is often better than the true Newton method 
        // (Numerical Recipes in C++, chap. 10.7)
#ifndef _MSC_VER
        double variables  [mixlen+1];
        double upper_bound[mixlen+1];
        double lower_bound[mixlen+1];
        bool   bound_check[mixlen+1];
#else 
        boost::scoped_array<double> variables ( new double [mixlen + 1]);
        boost::scoped_array<double> upper_bound ( new double [mixlen + 1]);
        boost::scoped_array<double> lower_bound ( new double [mixlen + 1]);
        boost::scoped_array<bool> bound_check( new bool [mixlen + 1]);
#endif
        for (int i = 0; i < mixlen; i++) {
            lower_bound[i+1] = params->min_branch_length;
            variables[i+1] = current_it->getLength(i);
            upper_bound[i+1] = params->max_branch_length;
            bound_check[i+1] = false;
        }

#ifndef _MSC_VER
        double grad[mixlen + 1];
        double hessian[mixlen * mixlen];
#else
        boost::scoped_array<double> grad( new double [mixlen + 1]);
        boost::scoped_array<double> hessian( new double [mixlen * mixlen]);
#endif 
        computeFuncDervMulti(&variables[1], &grad[0], &hessian[0]);
        double score;
        if (params->optimize_alg_mixlen.find("BFGS-B") != string::npos)
            score = -L_BFGS_B(mixlen, &variables[1], &lower_bound[1], 
                              &upper_bound[1], params->min_branch_length);
        else
            score = -minimizeMultiDimen(&variables[0], mixlen, 
                                        &lower_bound[0], &upper_bound[0], 
                                        &bound_check[0], params->min_branch_length, 
                                        &hessian[0]);

        for (int i = 0; i < mixlen; i++) {
            current_it->setLength(i, variables[i+1]);
            current_it_back->setLength(i, variables[i+1]);
        }
        if (verbose_mode >= VerboseMode::VB_DEBUG) {
            cout << "Mixlen-LnL: " << score << endl;
        }
    } else {
        if (!model_factory->fused_mix_rate && getModel()->isMixture()) {
            outError("Please use option -optlen BFGS to disable EM algorithm");
        }
        // EM algorithm
        intptr_t nptn = aln->getNPattern();
        size_t   nmix = site_rate->getNRate();
        ASSERT(nmix == mixlen);

        // first compute _pattern_lh_cat
        double tree_lh = -DBL_MAX;

        // 2 steps are empirically determined to be best!
        for (int EM_step = 0; EM_step < 2; EM_step++) {
            double new_tree_lh = computePatternLhCat(WSL_RATECAT);
            if (new_tree_lh+1.0 < tree_lh) {
                cout << "WARN: at EM step " << EM_step 
                     << " new_tree_lh " << new_tree_lh 
                     << " worse than tree_lh " << tree_lh << endl;
            }
            if (new_tree_lh-params->min_branch_length < tree_lh) {
                break;
            }
            tree_lh = new_tree_lh;
            // E-step
            // decoupled weights (prop) from _pattern_lh_cat 
            // to obtain L_ci and compute pattern likelihood L_i
            for (intptr_t ptn = 0; ptn < nptn; ptn++) {
                double* this_lk_cat = tree_buffers._pattern_lh_cat 
                                    + ptn*nmix;
                double  lk_ptn = ptn_invar[ptn];
                for (size_t c = 0; c < nmix; c++) {
                    lk_ptn += this_lk_cat[c];
                }
                ASSERT(lk_ptn != 0.0);
                lk_ptn = ptn_freq[ptn] / lk_ptn;
                // transform _pattern_lh_cat into 
                // posterior probabilities of each category
                for (size_t c = 0; c < nmix; c++) {
                    this_lk_cat[c] *= lk_ptn;
                }
            } 
         
            double negative_lh;
            double optx;
            tree_buffers.theta_computed = false;
            computePtnFreq();
            
            for (cur_mixture = 0; cur_mixture < mixlen; cur_mixture++) {
                double *this_lk_cat = tree_buffers._pattern_lh_cat+cur_mixture;
                for (intptr_t ptn = 0; ptn < nptn; ptn++) {
                    ptn_freq[ptn] = this_lk_cat[ptn*nmix];
                }                
                double current_len = current_it->getLength(cur_mixture);
                ASSERT(current_len >= 0.0);
                // Newton-Raphson method
                optx = minimizeNewton(params->min_branch_length, 
                                      current_len, params->max_branch_length, 
                                      params->min_branch_length, negative_lh, 
                                      maxNRStep);
                current_it->setLength(cur_mixture, optx);
                current_it_back->setLength(cur_mixture, optx);
            }
            cur_mixture = -1;
            // reset ptn_freq
            ptn_freq_computed = false;
            computePtnFreq();
        } // for EM_step
    }
    if (clearLH) {
        node1->clearReversePartialLh(node2);
        node2->clearReversePartialLh(node1);
    }
}

/**
    return the number of dimensions
*/
int PhyloTreeMixlen::getNDim() const {
    return mixlen;
}

/**
    the target function which needs to be optimized
    @param x the input vector x
    @return the function value at x
*/
double PhyloTreeMixlen::targetFunk(double x[]) {
    for (int i = 0; i < mixlen; i++) {
        current_it->setLength(i, x[i+1]);
        current_it_back->setLength(i, x[i+1]);
    }
    if (tree_buffers.theta_computed) {
        return -computeLikelihoodFromBuffer();
    }
    else {
        return -computeLikelihoodBranch(current_it, current_it_back->getNode(),
                                        tree_buffers);
    }
}

double PhyloTreeMixlen::derivativeFunk(double x[], double dfx[]) {
    for (int i = 0; i < mixlen; ++i) {
        ASSERT(!std::isnan(x[i+1]));
        current_it->setLength(i, x[i+1]);
        current_it_back->setLength(i, x[i+1]);
    }
#ifndef _MSC_VER
    double df[mixlen + 1];
    double ddf[mixlen*mixlen];
#else
    boost::scoped_array<double> df  ( new double [mixlen + 1] );
    boost::scoped_array<double> ddf ( new double [static_cast<size_t>(mixlen) * static_cast<size_t>(mixlen)]);
#endif
    computeLikelihoodDerv(current_it, current_it_back->getNode(), 
                          &df[0], &ddf[0], tree_buffers);
    for (int i = 0; i < mixlen; ++i) {
        df[i] = -df[i];
    }
    memcpy(dfx+1, &df[0], sizeof(double)*mixlen);
    return -df[mixlen];
}

void PhyloTreeMixlen::computeFuncDervMulti
        (double *value, double *df, double *ddf) {
    for (int i = 0; i < mixlen; i++) {
        current_it->setLength(i, value[i]);
        current_it_back->setLength(i, value[i]);
    }
    computeLikelihoodDerv(current_it, current_it_back->getNode(), 
                          df, ddf, tree_buffers);

    // last element of df is the tree log-ikelihood
    for (int i = 0; i <= mixlen; i++) {
        df[i] = -df[i];
    }
    int mixlen2 = mixlen * mixlen;
    for (int i = 0; i < mixlen2; i++)
        ddf[i] = -ddf[i];
}

double PhyloTreeMixlen::optimizeAllBranches(int my_iterations, double tolerance,
                                            int maxNRStep, bool were_lengths_consistent,
                                            PhyloTree* report_to_tree) {
    if (report_to_tree==nullptr) {
        report_to_tree = this;
    }
	initializeMixlen(tolerance, false, report_to_tree);
    clearAllPartialLH();
    clearAllPartialParsimony(false);
    
    double tree_lh = PhyloTree::optimizeAllBranches(my_iterations, tolerance,
                                                    maxNRStep, were_lengths_consistent,
                                                    report_to_tree);
    if (!initializing_mixlen) {
        assignMeanMixBranches();
    }
    return tree_lh;
}

pair<int, int> PhyloTreeMixlen::optimizeNNI(bool speedNNI, const char* context) {
    DoubleVector meanlenvec;
    treeLengths(meanlenvec);
    // compute mean branch length
    for (int j = 0; j < mixlen; j++) {
        meanlenvec[j] /= (branchNum);
    }
    // scan over all branches and fix short/long branches
    PhyloNodeMixlenVector nodes1, nodes2;
    getBranches(nodes1, nodes2);
    int num_fixed = 0;
    for (int i = 0; i < nodes1.size(); i++) {
        PhyloNeighborMixlen* nei      = nodes1[i]->findNeighbor(nodes2[i]);
        PhyloNeighborMixlen* nei_back = nodes2[i]->findNeighbor(nodes1[i]);
        for (int j = 0; j < mixlen; j++) {
            if (nei->lengths[j] > params->max_branch_length*0.9) {
                // if too long or too short branch, assign with mean branch length
                nei->lengths[j] = nei_back->lengths[j] = meanlenvec[j];
                num_fixed++;
            }
        }
    }
    if (num_fixed > 0) {
        optimizeBranches(num_fixed);
    }
    return IQTree::optimizeNNI(speedNNI, context);
}

void PhyloTreeMixlen::printBranchLength(ostream &out, int brtype, 
                                        bool print_slash, 
                                        Neighbor *length_nei) {
    PhyloNeighborMixlen* nei = dynamic_cast<PhyloNeighborMixlen*>(length_nei);
    if (nei->lengths.empty()) {
        return PhyloTree::printBranchLength(out, brtype, print_slash, length_nei);
    }
    if ((brtype & (WT_BR_LEN+WT_BR_SCALE)) == 0) {
        return;
    }
    if (cur_mixture == -1) {
        // print mixture branch lengths
        out << "[";
        for (int i = 0; i < mixlen; i++) {
            if (i > 0) {
                out << BRANCH_LENGTH_SEPARATOR;
            }
            double length = nei->lengths[i];
            if (brtype & WT_BR_SCALE) length *= len_scale;
            if (brtype & WT_BR_LEN_ROUNDING) length = round(length);
            if (brtype & WT_BR_LEN) {
                if (brtype & WT_BR_LEN_FIXED_WIDTH) {
                    out << fixed << length;
                }
                else {
                    out << length;
                }
            } else if (brtype & WT_BR_CLADE && 
                       length_nei->node->name != ROOT_NAME) {
                out << length;
            }
        }
        out << "]";
    }

    if (brtype & WT_BR_LEN) {
        out << ":";
    }
    else if ((brtype & WT_BR_CLADE) && print_slash && 
             length_nei->node->name != ROOT_NAME) {
        out << "/";
    }    
    double length = length_nei->length;

    if (cur_mixture >= 0) {
        // print branch length of a mixture component only!
        length = nei->lengths[cur_mixture];
    }
    if (brtype & WT_BR_SCALE) {
        length *= len_scale;
    }
    if (brtype & WT_BR_LEN_ROUNDING) {
        length = round(length);
    }
    if (brtype & WT_BR_LEN) {
        if (brtype & WT_BR_LEN_FIXED_WIDTH) {
            out << fixed << length;
        }
        else {
            out << length;
        }
    } else if (brtype & WT_BR_CLADE && 
               length_nei->node->name != ROOT_NAME) {
        out << length;
    }
}

void PhyloTreeMixlen::printResultTree(string suffix) {
    if (MPIHelper::getInstance().isWorker()) {
        return;
    }
    if (params->suppress_output_flags & OUT_TREEFILE) {
        return;
    }
    setRootNode(params->root);
    string tree_file_name = params->out_prefix;
    tree_file_name += ".treefile";
    if (suffix.compare("") != 0) {
        tree_file_name += "." + suffix;
    }
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(tree_file_name.c_str());
        cur_mixture = -1;
        printTree(out, WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | 
                       WT_SORT_TAXA | WT_NEWLINE);
        for (cur_mixture = 0; cur_mixture < mixlen; cur_mixture++) {
            //out << "[Heterotachy class " << cur_mixture+1 << "]" << endl;
            printTree(out, WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | 
                           WT_SORT_TAXA | WT_NEWLINE);
        }
        cur_mixture = -1;
        out.close();
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, tree_file_name);
    }
    if (verbose_mode >= VerboseMode::VB_MED) {
        cout << "Best tree printed to " << tree_file_name << endl;
    }
}


/*************** Using cppoptlib for branch length optimization ***********/

#ifdef USE_CPPOPTLIB
double PhyloTreeMixlen::value(const TVector &x) {
    double xx[mixlen+1];
    for (int i = 0; i < mixlen; i++)
        xx[i+1] = x(i);
    return targetFunk(xx);
}

void PhyloTreeMixlen::gradient(const TVector &x, TVector &grad) {
    int i;
    double xx[mixlen];
    for (i = 0; i < mixlen; i++)
        xx[i] = x(i);
    double df[mixlen+1], ddf[mixlen*mixlen];
    computeFuncDervMulti(xx, df, ddf);
    for (i = 0; i < mixlen; i++)
        grad(i) = df[i];
}

void PhyloTreeMixlen::hessian(const TVector &x, THessian &hessian) {
    int i, j;
    double xx[mixlen];
    for (i = 0; i < mixlen; i++)
        xx[i] = x(i);
    int mixlen2 = mixlen*mixlen;
    double df[mixlen+1], ddf[mixlen2];
    computeFuncDervMulti(xx, df, ddf);

    for (i = 0; i < mixlen; i++)
        for (j = 0; j < mixlen; j++)
            hessian(i, j) = ddf[i*mixlen+j];
}
#endif

// defining log-likelihood derivative function for EM algorithm
void PhyloTreeMixlen::computeFuncDerv(double value, double &df, double &ddf) {

    if (cur_mixture < 0) {
        return PhyloTree::computeFuncDerv(value, df, ddf);
    }
    current_it->setLength(cur_mixture, value);
    current_it_back->setLength(cur_mixture, value);
    (this->*computeLikelihoodDervMixlenPointer)
        (current_it,  current_it_back->getNode(), 
         df, ddf, tree_buffers);
	df  = -df;
    ddf = -ddf;
    return;

    PhyloNeighbor* dad_branch  = current_it;
    PhyloNode*     dad         = current_it_back->getNode();
    PhyloNode*     node        = dad_branch->getNode();
    PhyloNeighbor* node_branch = node->findNeighbor(dad);
    if (!central_partial_lh) {
        initializeAllPartialLh();
    }
    if (node->isLeaf()) {
        std::swap(dad, node);
        std::swap(dad_branch, node_branch);
    }
    ASSERT(dad_branch->isLikelihoodComputed()  || node->isLeaf());
    ASSERT(node_branch->isLikelihoodComputed() || dad->isLeaf());

    size_t nstates   = aln->num_states;
    size_t ncat      = site_rate->getNRate();
    size_t nmixture  = model->getNMixtures();
    size_t block     = ncat * nstates * nmixture;
    size_t statemix  = nstates * nmixture;
    size_t statecat  = nstates * ncat;
    intptr_t orig_nptn = aln->size();
    intptr_t nptn    = aln->size()+model_factory->unobserved_ptns.size();
    intptr_t maxptn  = get_safe_upper_limit(nptn);

    ModelSubst*        model_to_use = nullptr;
    RateHeterogeneity* rate_model   = nullptr;
    ModelSubst*        other_model  = nullptr;
    RateHeterogeneity* other_rate   = nullptr;
    double*            tip_lh       = tip_partial_lh;
    getModelAndTipLikelihood(dad, node, model_to_use, other_model,
                             rate_model, other_rate, tip_lh);

    double*            eval         = model_to_use->getEigenvalues();
    ASSERT(eval);

    ASSERT(tree_buffers.theta_all);
    if (!tree_buffers.theta_computed) {
		// precompute theta for fast branch length optimization

	    if (dad->isLeaf() && 0 <= dad->id) {
	    	// special treatment for TIP-INTERNAL NODE case
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (intptr_t ptn = 0; ptn < nptn; ptn++) {
                double* partial_lh_dad = dad_branch->partial_lh + ptn*block;
                double* theta          = tree_buffers.theta_all + ptn*block;
                
                // TODO: check with vectorclass!
                double *lh_tip         = tip_lh +
                ((int)((ptn < orig_nptn) 
                    ? (aln->at(ptn))[dad->id] 
                    :  model_factory->unobserved_ptns[ptn-orig_nptn][dad->id]))*statemix;
                for (size_t m = 0; m < nmixture; m++) {
                    for (size_t i = 0; i < statecat; i++) {
                        theta[m*statecat+i] = lh_tip[m*nstates 
                                            + i%nstates] * partial_lh_dad[m*statecat+i];
                    }
                }
            }
			// ascertainment bias correction
	    } else {
	    	// both dad and node are internal nodes
		    double* partial_lh_node = node_branch->partial_lh;
		    double* partial_lh_dad  = dad_branch->partial_lh;

            intptr_t all_entries = nptn*block;
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (intptr_t i = 0; i < all_entries; i++) {
                tree_buffers.theta_all[i] = partial_lh_node[i] * partial_lh_dad[i];
            }
        }
        if (nptn < maxptn) {
            // copy dummy values
            for (intptr_t ptn = nptn; ptn < maxptn; ptn++)
                memcpy(&tree_buffers.theta_all[ptn*block],
                       &tree_buffers.theta_all[(ptn-1)*block], block*sizeof(double));
        }
        tree_buffers.theta_computed = true;
    }
    std::vector<double> val_vector(statecat*3);
    double* val0 = val_vector.data();
    double* val1 = val0 + statecat;
    double* val2 = val1 + statecat;
    for (int c = 0; c < ncat; c++) {
        double prop = rate_model->getProp(c);
        for (size_t i = 0; i < nstates; i++) {
            double cof   = eval[cur_mixture*nstates+i]*rate_model->getRate(c);
                           // length for heterotachy model
            double val   = exp(cof*dad_branch->getLength(cur_mixture))
                               * prop * model_to_use->getMixtureWeight(cur_mixture);
            double val1_ = cof * val;
            val0[(c)*nstates+i] = val;
            val1[(c)*nstates+i] = val1_;
            val2[(c)*nstates+i] = cof*val1_;
        }
    }
    double my_df = 0.0, my_ddf = 0.0, prob_const = 0.0;
    double df_const = 0.0, ddf_const = 0.0;

#ifdef _OPENMP
#pragma omp parallel for reduction(+:my_df,my_ddf,prob_const,df_const,ddf_const)
#endif
    for (intptr_t ptn = 0; ptn < nptn; ptn++) {
        double  lh_ptn  = ptn_invar[ptn];
        double  df_ptn  = 0.0;
        double  ddf_ptn = 0.0;
        double* theta   = tree_buffers.theta_all + ptn*block 
                        + cur_mixture*statecat;
        for (size_t i = 0; i < statecat; i++) {
            lh_ptn  += val0[i] * theta[i];
            df_ptn  += val1[i] * theta[i];
            ddf_ptn += val2[i] * theta[i];
        }
        lh_ptn = fabs(lh_ptn);
        if (ptn < orig_nptn) {
            double df_frac = df_ptn / lh_ptn;
            double ddf_frac = ddf_ptn / lh_ptn;
            double freq = ptn_freq[ptn];
            double tmp1 = df_frac * freq;
            double tmp2 = ddf_frac * freq;
            my_df  += tmp1;
            my_ddf += tmp2 - tmp1 * df_frac;
        } else {
            // ascertainment bias correction
            prob_const += lh_ptn;
            df_const   += df_ptn;
            ddf_const  += ddf_ptn;
        }
    }
    df  = my_df;
    ddf = my_ddf;
    if (std::isnan(df) || std::isinf(df)) {
        df = 0.0;
        ddf = 0.0;
    }
    if (orig_nptn < nptn) {
        // ascertainment bias correction
        prob_const = 1.0 - prob_const;
        double df_frac = df_const / prob_const;
        double ddf_frac = ddf_const / prob_const;
        size_t nsites = aln->getNSite();
        df += nsites * df_frac;
        ddf += nsites *(ddf_frac + df_frac*df_frac);
    }

    df  = -df;
    ddf = -ddf;
}
