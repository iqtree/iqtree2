/*
 * phylotesting.cpp
 *
 *  Created on: Sep 21, 2019
 *      Author: minh
 */



#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iqtree_config.h>

#include "treetesting.h"
#include "tree/phylotree.h"
#include "tree/phylosupertree.h"
#include "tree/iqtreemix.h"
#include "gsl/mygsl.h"
#include "utils/timeutil.h"


void printSiteLh(const char*filename, PhyloTree *tree, double *ptn_lh,
                 bool append, const char *linename) {
    double *pattern_lh;

    if (!tree->isTreeMix()) {
        if (!ptn_lh) {
            pattern_lh = new double[tree->getAlnNPattern()];
            tree->computePatternLikelihood(pattern_lh);
        } else
            pattern_lh = ptn_lh;
    }

    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        if (append) {
            out.open(filename, ios::out | ios::app);
        } else {
            out.open(filename);
        }
        
        if (tree->isTreeMix()) {
            // tree mixture model
            IQTreeMix* treeMix = (IQTreeMix*) tree;
            treeMix->showLhProb(out);
            /*
            // also print out the pattern's likelihoods
            string pfilename = filename;
            pfilename += ".ptn";
            ofstream pout;
            pout.exceptions(ios::failbit | ios::badbit);
            if (append) {
                pout.open(pfilename.c_str(), ios::out | ios::app);
            } else {
                pout.open(pfilename.c_str());
            }
            treeMix->showPatternLhProb(pout);
            pout.close();
            if (!append)
                cout << "pattern log-likelihoods printed to " << pfilename << endl;
            */
        } else {
            if (!append)
                out << 1 << " " << tree->getAlnNSite() << endl;
            if (!linename)
                out << "Site_Lh   ";
            else {
                out.width(10);
                out << left << linename;
            }
            IntVector pattern_index;
            tree->aln->getSitePatternIndex(pattern_index);
            for (size_t i = 0; i < tree->getAlnNSite(); i++)
                out << " " << pattern_lh[pattern_index[i]];
            out << endl;
        }
        out.close();
        if (!append)
            cout << "Site log-likelihoods printed to " << filename << endl;
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, filename);
    }
    
    if (!tree->isTreeMix() && !ptn_lh)
        delete[] pattern_lh;
}

void printPartitionLh(const char*filename, PhyloTree *tree, double *ptn_lh,
                      bool append, const char *linename) {
    
    ASSERT(tree->isSuperTree());
    PhyloSuperTree *stree = (PhyloSuperTree*)tree;
    double *pattern_lh;
    if (!ptn_lh) {
        pattern_lh = new double[tree->getAlnNPattern()];
        tree->computePatternLikelihood(pattern_lh);
    } else
        pattern_lh = ptn_lh;
    
    double partition_lh[stree->size()];
    double *pattern_lh_ptr = pattern_lh;
    for (int part = 0; part < stree->size(); part++) {
        size_t nptn = stree->at(part)->getAlnNPattern();
        partition_lh[part] = 0.0;
        for (int i = 0; i < nptn; i++)
            partition_lh[part] += pattern_lh_ptr[i] * stree->at(part)->ptn_freq[i];
        pattern_lh_ptr += nptn;
    }
    
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        if (append) {
            out.open(filename, ios::out | ios::app);
        } else {
            out.open(filename);
            out << 1 << " " << stree->size() << endl;
        }
        if (!linename)
            out << "Part_Lh   ";
        else {
            out.width(10);
            out << left << linename;
        }
        for (int i = 0; i < stree->size(); ++i) {
            out << " " << partition_lh[i];
        }
        out << endl;
        out.close();
        if (!append)
            cout << "Partition log-likelihoods printed to " << filename << endl;
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, filename);
    }    
    delete[] pattern_lh;
}

void printSiteLhCategory(const char*filename, PhyloTree *tree, SiteLoglType wsl) {
    
    if (wsl == WSL_NONE || wsl == WSL_SITE)
        return;
    int ncat = tree->getNumLhCat(wsl);
    if (tree->isSuperTree()) {
        PhyloSuperTree *stree = (PhyloSuperTree*)tree;
        for (auto it = stree->begin(); it != stree->end(); it++) {
            int part_ncat = (*it)->getNumLhCat(wsl);
            if (part_ncat > ncat)
                ncat = part_ncat;
        }
    }
    int i;
    
    
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(filename);
        out << "# Site likelihood per rate/mixture category" << endl
        << "# This file can be read in MS Excel or in R with command:" << endl
        << "#   tab=read.table('" <<  filename << "',header=TRUE,fill=TRUE)" << endl
        << "# Columns are tab-separated with following meaning:" << endl;
        if (tree->isSuperTree()) {
            out << "#   Part:   Partition ID (1=" << ((PhyloSuperTree*)tree)->at(0)->aln->name << ", etc)" << endl
            << "#   Site:   Site ID within partition (starting from 1 for each partition)" << endl;
        } else
            out << "#   Site:   Alignment site ID" << endl;
        
        out << "#   LnL:    Logarithm of site likelihood" << endl
        << "#           Thus, sum of LnL is equal to tree log-likelihood" << endl
        << "#   LnLW_k: Logarithm of (category-k site likelihood times category-k weight)" << endl
        << "#           Thus, sum of exp(LnLW_k) is equal to exp(LnL)" << endl;
        
        if (tree->isSuperTree()) {
            out << "Part\tSite\tLnL";
        } else
            out << "Site\tLnL";
        for (i = 0; i < ncat; i++)
            out << "\tLnLW_" << i+1;
        out << endl;
        out.precision(4);
        out.setf(ios::fixed);
        
        tree->writeSiteLh(out, wsl);
        
        out.close();
        cout << "Site log-likelihoods per category printed to " << filename << endl;
        /*
         if (!tree->isSuperTree()) {
         cout << "Log-likelihood of constant sites: " << endl;
         double const_prob = 0.0;
         for (i = 0; i < tree->aln->getNPattern(); i++)
         if (tree->aln->at(i).isConst()) {
         Pattern pat = tree->aln->at(i);
         for (Pattern::iterator it = pat.begin(); it != pat.end(); it++)
         cout << tree->aln->convertStateBackStr(*it);
         cout << ": " << pattern_lh[i] << endl;
         const_prob += exp(pattern_lh[i]);
         }
         cout << "Probability of const sites: " << const_prob << endl;
         }
         */
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, filename);
    }
    
}

void printAncestralSequences(const char *out_prefix, PhyloTree *tree, AncestralSeqType ast) {
    
    //    int *joint_ancestral = NULL;
    //
    //    if (tree->params->print_ancestral_sequence == AST_JOINT) {
    //        joint_ancestral = new int[nptn*tree->leafNum];
    //        tree->computeJointAncestralSequences(joint_ancestral);
    //    }
    
    string filename = (string)out_prefix + ".state";
    //    string filenameseq = (string)out_prefix + ".stateseq";
    
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(filename.c_str());
        out.setf(ios::fixed, ios::floatfield);
        out.precision(5);
        
        //        ofstream outseq;
        //        outseq.exceptions(ios::failbit | ios::badbit);
        //        outseq.open(filenameseq.c_str());
        
        NodeVector nodes;
        tree->getInternalNodes(nodes);
        
        double *marginal_ancestral_prob;
        int *marginal_ancestral_seq;
        
        //        if (tree->params->print_ancestral_sequence == AST_JOINT)
        //            outseq << 2*(tree->nodeNum-tree->leafNum) << " " << nsites << endl;
        //        else
        //            outseq << (tree->nodeNum-tree->leafNum) << " " << nsites << endl;
        //
        //        int name_width = max(tree->aln->getMaxSeqNameLength(),6)+10;
        
        out << "# Ancestral state reconstruction for all nodes in " << tree->params->out_prefix << ".treefile" << endl
        << "# This file can be read in MS Excel or in R with command:" << endl
        << "#   tab=read.table('" <<  tree->params->out_prefix << ".state',header=TRUE)" << endl
        << "# Columns are tab-separated with following meaning:" << endl
        << "#   Node:  Node name in the tree" << endl;
        if (tree->isSuperTree()) {
            PhyloSuperTree *stree = (PhyloSuperTree*)tree;
            out << "#   Part:  Partition ID (1=" << stree->at(0)->aln->name << ", etc)" << endl
            << "#   Site:  Site ID within partition (starting from 1 for each partition)" << endl;
        } else
            out << "#   Site:  Alignment site ID" << endl;
        
        out << "#   State: Most likely state assignment" << endl
        << "#   p_X:   Posterior probability for state X (empirical Bayesian method)" << endl;
        
        if (tree->isSuperTree()) {
            PhyloSuperTree *stree = (PhyloSuperTree*)tree;
            out << "Node\tPart\tSite\tState";
            for (size_t i = 0; i < stree->front()->aln->num_states; i++)
                out << "\tp_" << stree->front()->aln->convertStateBackStr(i);
        } else {
            out << "Node\tSite\tState";
            for (size_t i = 0; i < tree->aln->num_states; i++)
                out << "\tp_" << tree->aln->convertStateBackStr(i);
        }
        out << endl;
        
        
        bool orig_kernel_nonrev;
        tree->initMarginalAncestralState(out, orig_kernel_nonrev, marginal_ancestral_prob, marginal_ancestral_seq);
        
        for (NodeVector::iterator it = nodes.begin(); it != nodes.end(); it++) {
            PhyloNode *node = (PhyloNode*)(*it);
            PhyloNode *dad = (PhyloNode*)node->neighbors[0]->node;
            
            tree->computeMarginalAncestralState((PhyloNeighbor*)dad->findNeighbor(node), dad,
                                                marginal_ancestral_prob, marginal_ancestral_seq);
            
            //            int *joint_ancestral_node = joint_ancestral + (node->id - tree->leafNum)*nptn;
            
            // set node name if neccessary
            if (node->name.empty() || !isalpha(node->name[0])) {
                node->name = "Node" + convertIntToString(node->id-tree->leafNum+1);
            }
            
            // print ancestral state probabilities
            tree->writeMarginalAncestralState(out, node, marginal_ancestral_prob, marginal_ancestral_seq);
            
            // print ancestral sequences
            //            outseq.width(name_width);
            //            outseq << left << node->name << " ";
            //            for (i = 0; i < nsites; i++)
            //                outseq << tree->aln->convertStateBackStr(marginal_ancestral_seq[pattern_index[i]]);
            //            outseq << endl;
            //
            //            if (tree->params->print_ancestral_sequence == AST_JOINT) {
            //                outseq.width(name_width);
            //                outseq << left << (node->name+"_joint") << " ";
            //                for (i = 0; i < nsites; i++)
            //                    outseq << tree->aln->convertStateBackStr(joint_ancestral_node[pattern_index[i]]);
            //                outseq << endl;
            //            }
        }
        
        tree->endMarginalAncestralState(orig_kernel_nonrev, marginal_ancestral_prob, marginal_ancestral_seq);
        
        out.close();
        //        outseq.close();
        cout << "Ancestral state probabilities printed to " << filename << endl;
        //        cout << "Ancestral sequences printed to " << filenameseq << endl;
        
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, filename);
    }
    
    //    if (joint_ancestral)
    //        delete[] joint_ancestral;
    
}

void printSiteProbCategory(const char*filename, PhyloTree *tree, SiteLoglType wsl) {
    
    if (wsl == WSL_NONE || wsl == WSL_SITE)
        return;
    // error checking
    if (!tree->getModel()->isMixture()) {
        if (wsl != WSL_RATECAT) {
            outWarning("Switch now to '-wspr' as it is the only option for non-mixture model");
            wsl = WSL_RATECAT;
        }
    } else {
        // mixture model
        if (wsl == WSL_MIXTURE_RATECAT && tree->getModelFactory()->fused_mix_rate) {
            outWarning("-wspmr is not suitable for fused mixture model, switch now to -wspm");
            wsl = WSL_MIXTURE;
        }
    }
    size_t cat, ncat = tree->getNumLhCat(wsl);
    double *ptn_prob_cat = new double[((size_t)tree->getAlnNPattern())*ncat];
    tree->computePatternProbabilityCategory(ptn_prob_cat, wsl);
    
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(filename);
        if (tree->isSuperTree())
            out << "Set\t";
        out << "Site";
        for (cat = 0; cat < ncat; cat++)
            out << "\tp" << cat+1;
        out << endl;
        IntVector pattern_index;
        if (tree->isSuperTree()) {
            PhyloSuperTree *super_tree = (PhyloSuperTree*)tree;
            size_t offset = 0;
            for (PhyloSuperTree::iterator it = super_tree->begin(); it != super_tree->end(); it++) {
                size_t part_ncat = (*it)->getNumLhCat(wsl);
                (*it)->aln->getSitePatternIndex(pattern_index);
                size_t nsite = (*it)->aln->getNSite();
                for (size_t site = 0; site < nsite; ++site) {
                    out << (it-super_tree->begin())+1 << "\t" << site+1;
                    double *prob_cat = ptn_prob_cat + (offset+pattern_index[site]*part_ncat);
                    for (cat = 0; cat < part_ncat; cat++)
                        out << "\t" << prob_cat[cat];
                    out << endl;
                }
                offset += (*it)->aln->getNPattern()*(*it)->getNumLhCat(wsl);
            }
        } else {
            tree->aln->getSitePatternIndex(pattern_index);
            size_t nsite = tree->getAlnNSite();
            for (size_t site = 0; site < nsite; ++site) {
                out << site+1;
                double *prob_cat = ptn_prob_cat + pattern_index[site]*ncat;
                for (cat = 0; cat < ncat; cat++) {
                    out << "\t" << prob_cat[cat];
                }
                out << endl;
            }
        }
        out.close();
        cout << "Site probabilities per category printed to " << filename << endl;
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, filename);
    }
    
}


void printSiteStateFreq(const char*filename, PhyloTree *tree, double *state_freqs) {
    size_t nsites = tree->getAlnNSite();
    size_t nstates = tree->aln->num_states;
    double *ptn_state_freq;
    if (state_freqs) {
        ptn_state_freq = state_freqs;
    } else {
        ptn_state_freq = new double[((size_t)tree->getAlnNPattern()) * nstates];
        tree->computePatternStateFreq(ptn_state_freq);
    }
    
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(filename);
        IntVector pattern_index;
        tree->aln->getSitePatternIndex(pattern_index);
        for (size_t i = 0; i < nsites; ++i) {
            out.width(6);
            out << left << i+1 << " ";
            double *state_freq = &ptn_state_freq[pattern_index[i]*nstates];
            for (size_t j = 0; j < nstates; ++j) {
                out.width(15);
                out << state_freq[j] << " ";
            }
            out << endl;
        }
        out.close();
        cout << "Site state frequency vectors printed to " << filename << endl;
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, filename);
    }
    if (!state_freqs)
        delete [] ptn_state_freq;
}

void printSiteStateFreq(const char* filename, Alignment *aln) {
    if (aln->site_state_freq.empty())
        return;
    size_t nsites  = aln->getNSite();
    int    nstates = aln->num_states;
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(filename);
        IntVector pattern_index;
        aln->getSitePatternIndex(pattern_index);
        for (size_t i = 0; i < nsites; ++i) {
            out.width(6);
            out << left << i+1 << " ";
            double *state_freq = aln->site_state_freq[pattern_index[i]];
            for (size_t j = 0; j < nstates; ++j) {
                out.width(15);
                out << state_freq[j] << " ";
            }
            out << endl;
        }
        out.close();
        cout << "Site state frequency vectors printed to " << filename << endl;
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, filename);
    }
}

int countDistinctTrees(istream &in, bool rooted, IQTree *tree, IntVector &distinct_ids, bool exclude_duplicate) {
    StringIntMap treels;
    int tree_id;
    for (tree_id = 0; !in.eof(); tree_id++) {
        if (exclude_duplicate) {
            tree->freeNode();
            tree->readTree(in, rooted);
            tree->setAlignment(tree->aln);
            tree->setRootNode(tree->params->root);
            StringIntMap::iterator it = treels.end();
            ostringstream ostr;
            tree->printTree(ostr, WT_TAXON_ID | WT_SORT_TAXA);
            it = treels.find(ostr.str());
            if (it != treels.end()) { // already in treels
                distinct_ids.push_back(it->second);
            } else {
                distinct_ids.push_back(-1);
                treels[ostr.str()] = tree_id;
            }
        } else {
            // ignore tree
            char ch;
            do {
                in >> ch;
            } while (!in.eof() && ch != ';');
            distinct_ids.push_back(-1);
        }
        char ch;
        in.exceptions(ios::goodbit);
        (in) >> ch;
        if (in.eof()) break;
        in.unget();
        in.exceptions(ios::failbit | ios::badbit);
    }
    in.clear();
    if (exclude_duplicate)
        return treels.size();
    else
        return distinct_ids.size();
}

int countDistinctTrees(const char *filename, bool rooted, IQTree *tree, IntVector &distinct_ids, bool exclude_duplicate)
{
    try {
        ifstream in;
        in.exceptions(ios::failbit | ios::badbit);
        in.open(filename);
        // remove the failbit
        in.exceptions(ios::badbit);
        
        int res = countDistinctTrees(in, rooted, tree, distinct_ids, exclude_duplicate);
        in.close();
        return res;
    } catch (ios::failure) {
        outError("Cannot read file ", filename);
        return 0;
    }
}
//const double TOL_RELL_SCORE = 0.01;

/*
 Problem: solve the following linear system equation:
 a_1*x + b_1*y = c_1
 a_2*x + b_2*y = c_2
 ....
 a_n*x + b_n*y = c_n
 
 becomes minimizing weighted least square:
 
 sum_k { w_k*[ c_k - (a_k*x + b_k*y) ]^2 }
 
 
 the solution is:
 
 x = [(sum_k w_k*b_k*c_k)*(sum_k w_k*a_k*b_k) - (sum_k w_k*a_k*c_k)(sum_k w_k*b_k^2)] /
 [ (sum_k w_k*a_k*b_k)^2 - (sum_k w_k*a_k^2)*(sum_k w_k*b_k^2) ]
 
 y = [(sum_k w_k*a_k*c_k)*(sum_k w_k*a_k*b_k) - (sum_k w_k*b_k*c_k)(sum_k w_k*a_k^2)] /
 [ (sum_k w_k*a_k*b_k)^2 - (sum_k w_k*a_k^2)*(sum_k w*k*b_k^2) ]
 
 @param n number of data points
 @param w weight vector of length n
 @param a a value vector of length n
 @param b b value vector of length n
 @param c c value vector of length n
 @param[out] x x-value
 @param[out] y y-value
 @return least square value
 */
void doWeightedLeastSquare(int n, double *w, double *a, double *b, double *c, double &x, double &y, double &se) {
    int k;
    double BC = 0.0, AB = 0.0, AC = 0.0, A2 = 0.0, B2 = 0.0;
    double denom;
    for (k = 0; k < n; k++) {
        double wa = w[k]*a[k];
        double wb = w[k]*b[k];
        AB += wa*b[k];
        BC += wb*c[k];
        AC += wa*c[k];
        A2 += wa*a[k];
        B2 += wb*b[k];
    }
    denom = 1.0/(AB*AB - A2*B2);
    x = (BC*AB - AC*B2) * denom;
    y = (AC*AB - BC*A2) * denom;
    
    se = -denom*(B2+A2+2*AB);
    ASSERT(se >= 0.0);
}

/**
 MLE estimates for AU test
 */
class OptimizationAUTest : public Optimization {
    
public:
    
    OptimizationAUTest(double d, double c, int nscales, double *bp, double *rr, double *rr_inv) {
        this->d = d;
        this->c = c;
        this->bp = bp;
        this->rr = rr;
        this->rr_inv = rr_inv;
        this->nscales = nscales;
        
    }
    
    /**
     return the number of dimensions
     */
    virtual int getNDim() { return 2; }
    
    
    /**
     the target function which needs to be optimized
     @param x the input vector x
     @return the function value at x
     */
    virtual double targetFunk(double x[]) {
        d = x[1];
        c = x[2];
        double res = 0.0;
        for (int k = 0; k < nscales; k++) {
            double cdf = gsl_cdf_ugaussian_P(d*rr[k] + c*rr_inv[k]);
            res += bp[k] * log(1.0 - cdf) + (1.0-bp[k])*log(cdf);
        }
        return res;
    }
    
    void optimizeDC() {
        double x[3], lower[3], upper[3];
        bool bound_check[3];
        x[1] = d;
        x[2] = c;
        lower[1] = lower[2] = 1e-4;
        upper[1] = upper[2] = 100.0;
        bound_check[1] = bound_check[2] = false;
        minimizeMultiDimen(x, 2, lower, upper, bound_check, 1e-4);
        d = x[1];
        c = x[2];
    }
    
    double d, c;
    int nscales;
    double *bp;
    double *rr;
    double *rr_inv;
};


/* BEGIN CODE WAS TAKEN FROM CONSEL PROGRAM */

/* binary search for a sorted vector
 find k s.t. vec[k-1] <= t < vec[k]
 */
int cntdist2(double *vec, int bb, double t)
{
    int i,i0,i1;
    
    i0=0; i1=bb-1;
    if(t < vec[0]) return 0;
    else if(vec[bb-1] <= t) return bb;
    
    while(i1-i0>1) {
        i=(i0+i1)/2;
        if(vec[i] <= t) i0=i;
        else i1=i;
    }
    
    return i1;
}

/*
 smoothing the counting for a sorted vector
 the piecewise linear function connecting
 F(v[i]) =  1/(2n) + i/n, for i=0,...,n-1
 F(1.5v[0]-0.5v[1]) = 0
 F(1.5v[n-1]-0.5v[n-2]) = 1.
 
 1. F(x)=0 for x<=1.5v[0]-0.5v[1]
 
 2. F(x)=1/(2n) + (1/n)*(x-v[0])/(v[1]-v[0])
 for 1.5v[0]-0.5v[1] < x <= v[0]
 
 3. F(x)=1/(2n) + i/n + (1/n)*(x-v[i])/(v[i]-v[i+1])
 for v[i] < x <= v[i+1], i=0,...,
 
 4. F(x)=1-(1/2n) + (1/n)*(x-v[n-1])/(v[n-1]-v[n-2])
 for v[n-1] < x <= 1.5v[n-1]-0.5v[n-2]
 
 5. F(x)=1 for x > 1.5v[n-1]-0.5v[n-2]
 */
double cntdist3(double *vec, int bb, double t)
{
    double p,n;
    int i;
    i=cntdist2(vec,bb,t)-1; /* to find vec[i] <= t < vec[i+1] */
    n=(double)bb;
    if(i<0) {
        if(vec[1]>vec[0]) p=0.5+(t-vec[0])/(vec[1]-vec[0]);
        else p=0.0;
    } else if(i<bb-1) {
        if(vec[i+1]>vec[i]) p=0.5+(double)i+(t-vec[i])/(vec[i+1]-vec[i]);
        else p=0.5+(double)i; /* <- should never happen */
    } else {
        if(vec[bb-1]-vec[bb-2]>0) p=n-0.5+(t-vec[bb-1])/(vec[bb-1]-vec[bb-2]);
        else p=n;
    }
    if(p>n) p=n; else if(p<0.0) p=0.0;
    return p;
}

double log3(double x)
{
    double y,z1,z2,z3,z4,z5;
    if(fabs(x)>1.0e-3) {
        y=-log(1.0-x);
    } else {
        z1=x; z2=z1*x; z3=z2*x; z4=z3*x; z5=z4*x;
        y=((((z5/5.0)+z4/4.0)+z3/3.0)+z2/2.0)+z1;
    }
    return y;
}

int mleloopmax=30;
double mleeps=1e-10;
int mlecoef(double *cnts, double *rr, double bb, int kk,
            double *coef0, /* set initinal value (size=2) */
            double *lrt, int *df, /* LRT statistic */
            double *se
            )
{
    int i,m,loop;
    double coef[2], update[2];
    double d1f, d2f, d11f, d12f, d22f; /* derivatives */
    double v11, v12, v22; /* inverse of -d??f */
    double a,e;
    double s[kk], r[kk],c[kk], b[kk],z[kk],p[kk],d[kk],g[kk],h[kk];
    
    m=0;
    for(i=0;i<kk;i++)
        {
        r[m]=rr[i]; s[m]=sqrt(rr[i]); c[m]=cnts[i]*bb; b[m]=bb;
        m++;
        }
    if(m<2) return 1;
    
    coef[0]=coef0[0]; /* signed distance */
    coef[1]=coef0[1]; /* curvature */
    
    for(loop=0;loop<mleloopmax;loop++) {
        d1f=d2f=d11f=d12f=d22f=0.0;
        for(i=0;i<m;i++) {
            z[i]=coef[0]*s[i]+coef[1]/s[i];
            p[i]=gsl_cdf_ugaussian_P(-z[i]);
            d[i]=gsl_ran_ugaussian_pdf(z[i]);
            if(p[i]>0.0 && p[i]<1.0) {
                g[i]=d[i]*( d[i]*(-c[i]+2.0*c[i]*p[i]-b[i]*p[i]*p[i])/
                           (p[i]*p[i]*(1.0-p[i])*(1.0-p[i]))
                           + z[i]*(c[i]-b[i]*p[i])/(p[i]*(1.0-p[i])) );
                h[i]=d[i]*(c[i]-b[i]*p[i])/(p[i]*(1.0-p[i]));
            } else { g[i]=h[i]=0.0; }
            d1f+= -h[i]*s[i]; d2f+= -h[i]/s[i];
            d11f+= g[i]*r[i]; d12f+= g[i]; d22f+= g[i]/r[i];
        }
        
        a=d11f*d22f-d12f*d12f;
        if(a==0.0) {
            return 2;
        }
        v11=-d22f/a; v12=d12f/a; v22=-d11f/a;
        
        /* Newton-Raphson update */
        update[0]=v11*d1f+v12*d2f; update[1]=v12*d1f+v22*d2f;
        coef[0]+=update[0]; coef[1]+=update[1];
        
        /* check convergence */
        e=-d11f*update[0]*update[0]-2.0*d12f*update[0]*update[1]
        -d22f*update[1]*update[1];
        
        if(e<mleeps) break;
    }
    
    /* calc log-likelihood */
    *lrt=0.0; *df=0;
    for(i=0;i<m;i++) {
        if(p[i]>0.0 && p[i]<1.0) {
            *df+=1;
            if(c[i]>0.0) a=c[i]*log(c[i]/b[i]/p[i]); else a=0.0;
            if(c[i]<b[i]) a+=(b[i]-c[i])*(log3(p[i])-log3(c[i]/b[i]));
            *lrt += a;
        }
    }
    *lrt *= 2.0; *df -= 2;
    
    /* write back the results */
    coef0[0]=coef[0]; coef0[1]=coef[1];
    *se = v11 + v22 - 2*v12;
    //  vmat[0][0]=v11;vmat[0][1]=vmat[1][0]=v12;vmat[1][1]=v22;
    if(loop==mleloopmax || *df< 0) i=1; else i=0;
    return i;
}

/* END CODE WAS TAKEN FROM CONSEL PROGRAM */

/**
 @param tree_lhs RELL score matrix of size #trees x #replicates
 */
void performAUTest(Params &params, PhyloTree *tree, double *pattern_lhs, vector<TreeInfo> &info) {
    
    if (params.topotest_replicates < 10000)
        outWarning("Too few replicates for AU test. At least -zb 10000 for reliable results!");
    
    /* STEP 1: specify scale factors */
    size_t nscales = 10;
    double r[] = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4};
    double rr[] = {sqrt(0.5), sqrt(0.6), sqrt(0.7), sqrt(0.8), sqrt(0.9), 1.0,
        sqrt(1.1), sqrt(1.2), sqrt(1.3), sqrt(1.4)};
    double rr_inv[] = {sqrt(1/0.5), sqrt(1/0.6), sqrt(1/0.7), sqrt(1/0.8), sqrt(1/0.9), 1.0,
        sqrt(1/1.1), sqrt(1/1.2), sqrt(1/1.3), sqrt(1/1.4)};
    
    /* STEP 2: compute bootstrap proportion */
    size_t ntrees = info.size();
    size_t nboot = params.topotest_replicates;
    //    double nboot_inv = 1.0 / nboot;
    
    size_t nptn = tree->getAlnNPattern();
    size_t maxnptn = get_safe_upper_limit(nptn);
    
    //    double *bp = new double[ntrees*nscales];
    //    memset(bp, 0, sizeof(double)*ntrees*nscales);
    
    double *treelhs;
    cout << (ntrees*nscales*nboot*sizeof(double) >> 20) << " MB required for AU test" << endl;
    treelhs = new double[ntrees*nscales*nboot];
    if (!treelhs)
        outError("Not enough memory to perform AU test!");
    
    size_t k, tid, ptn;
    
    double start_time = getRealTime();
    
    cout << "Generating " << nscales << " x " << nboot << " multiscale bootstrap replicates... ";
    
#ifdef _OPENMP
#pragma omp parallel private(k, tid, ptn)
    {
    int *rstream;
    init_random(params.ran_seed + omp_get_thread_num(), false, &rstream);
#else
    int *rstream = randstream;
#endif
    size_t boot;
    int *boot_sample = aligned_alloc<int>(maxnptn);
    memset(boot_sample, 0, maxnptn*sizeof(int));
    
    double *boot_sample_dbl = aligned_alloc<double>(maxnptn);
    
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (int k = 0; k < nscales; ++k) {
        string str = "SCALE=" + convertDoubleToString(r[k]);
        for (boot = 0; boot < nboot; boot++) {
            if (r[k] == 1.0 && boot == 0)
                // 2018-10-23: get one of the bootstrap sample as the original alignment
                tree->aln->getPatternFreq(boot_sample);
            else
                tree->aln->createBootstrapAlignment(boot_sample, str.c_str(), rstream);
            for (ptn = 0; ptn < maxnptn; ptn++)
                boot_sample_dbl[ptn] = boot_sample[ptn];
            double max_lh = -DBL_MAX, second_max_lh = -DBL_MAX;
            int max_tid = -1;
            for (tid = 0; tid < ntrees; tid++) {
                double *pattern_lh = pattern_lhs + (tid*maxnptn);
                double tree_lh;
                if (params.SSE == LK_386) {
                    tree_lh = 0.0;
                    for (ptn = 0; ptn < nptn; ptn++)
                        tree_lh += pattern_lh[ptn] * boot_sample_dbl[ptn];
                } else {
                    tree_lh = tree->dotProductDoubleCall(pattern_lh, boot_sample_dbl, nptn);
                }
                // rescale lh
                tree_lh /= r[k];
                
                // find the max and second max
                if (tree_lh > max_lh) {
                    second_max_lh = max_lh;
                    max_lh = tree_lh;
                    max_tid = tid;
                } else if (tree_lh > second_max_lh)
                    second_max_lh = tree_lh;
                
                treelhs[(tid*nscales+k)*nboot + boot] = tree_lh;
            }
            
            // compute difference from max_lh
            for (tid = 0; tid < ntrees; tid++)
                if (tid != max_tid)
                    treelhs[(tid*nscales+k)*nboot + boot] = max_lh - treelhs[(tid*nscales+k)*nboot + boot];
                else
                    treelhs[(tid*nscales+k)*nboot + boot] = second_max_lh - max_lh;
            //            bp[k*ntrees+max_tid] += nboot_inv;
        } // for boot
        
        // sort the replicates
        for (tid = 0; tid < ntrees; tid++) {
            quicksort<double,int>(treelhs + (tid*nscales+k)*nboot, 0, nboot-1);
        }
        
    } // for scale
    
    aligned_free(boot_sample_dbl);
    aligned_free(boot_sample);
    
#ifdef _OPENMP
    finish_random(rstream);
    }
#endif
    
    //    if (verbose_mode >= VB_MED) {
    //        cout << "scale";
    //        for (k = 0; k < nscales; k++)
    //            cout << "\t" << r[k];
    //        cout << endl;
    //        for (tid = 0; tid < ntrees; tid++) {
    //            cout << tid;
    //            for (k = 0; k < nscales; k++) {
    //                cout << "\t" << bp[tid+k*ntrees];
    //            }
    //            cout << endl;
    //        }
    //    }
    
    cout << getRealTime() - start_time << " seconds" << endl;
    
    /* STEP 3: weighted least square fit */
    
    double *cc = new double[nscales];
    double *w = new double[nscales];
    double *this_bp = new double[nscales];
    cout << "TreeID\tAU\tRSS\td\tc" << endl;
    for (tid = 0; tid < ntrees; tid++) {
        double *this_stat = treelhs + tid*nscales*nboot;
        double xn = this_stat[(nscales/2)*nboot + nboot/2], x;
        double c, d; // c, d in original paper
        int idf0 = -2;
        double z = 0.0, z0 = 0.0, thp = 0.0, th = 0.0, ze = 0.0, ze0 = 0.0;
        double pval, se;
        int df;
        double rss = 0.0;
        int step;
        const int max_step = 30;
        bool failed = false;
        for (step = 0; step < max_step; step++) {
            x = xn;
            int num_k = 0;
            for (k = 0; k < nscales; k++) {
                this_bp[k] = cntdist3(this_stat + k*nboot, nboot, x) / nboot;
                if (this_bp[k] <= 0 || this_bp[k] >= 1) {
                    cc[k] = w[k] = 0.0;
                } else {
                    double bp_val = this_bp[k];
                    cc[k] = -gsl_cdf_ugaussian_Pinv(bp_val);
                    double bp_pdf = gsl_ran_ugaussian_pdf(cc[k]);
                    w[k] = bp_pdf*bp_pdf*nboot / (bp_val*(1.0-bp_val));
                    num_k++;
                }
            }
            df = num_k-2;
            if (num_k >= 2) {
                // first obtain d and c by weighted least square
                doWeightedLeastSquare(nscales, w, rr, rr_inv, cc, d, c, se);
                
                // maximum likelhood fit
                double coef0[2] = {d, c};
                int mlefail = mlecoef(this_bp, r, nboot, nscales, coef0, &rss, &df, &se);
                
                if (!mlefail) {
                    d = coef0[0];
                    c = coef0[1];
                }
                
                se = gsl_ran_ugaussian_pdf(d-c)*sqrt(se);
                
                // second, perform MLE estimate of d and c
                //            OptimizationAUTest mle(d, c, nscales, this_bp, rr, rr_inv);
                //            mle.optimizeDC();
                //            d = mle.d;
                //            c = mle.c;
                
                /* STEP 4: compute p-value according to Eq. 11 */
                pval = gsl_cdf_ugaussian_Q(d-c);
                z = -pval;
                ze = se;
                // compute sum of squared difference
                rss = 0.0;
                for (k = 0; k < nscales; k++) {
                    double diff = cc[k] - (rr[k]*d + rr_inv[k]*c);
                    rss += w[k] * diff * diff;
                }
                
            } else {
                // not enough data for WLS
                int num0 = 0;
                for (k = 0; k < nscales; k++)
                    if (this_bp[k] <= 0.0) num0++;
                if (num0 > nscales/2)
                    pval = 0.0;
                else
                    pval = 1.0;
                se = 0.0;
                d = c = 0.0;
                rss = 0.0;
                if (verbose_mode >= VB_MED)
                    cout << "   error in wls" << endl;
                //info[tid].au_pvalue = pval;
                //break;
            }
            
            
            if (verbose_mode >= VB_MED) {
                cout.unsetf(ios::fixed);
                cout << "\t" << step << "\t" << th << "\t" << x << "\t" << pval << "\t" << se << "\t" << nscales-2 << "\t" << d << "\t" << c << "\t" << z << "\t" << ze << "\t" << rss << endl;
            }
            
            if(df < 0 && idf0 < 0) { failed = true; break;} /* degenerated */
            
            if ((df < 0) || (idf0 >= 0 && (z-z0)*(x-thp) > 0.0 && fabs(z-z0)>0.1*ze0)) {
                if (verbose_mode >= VB_MED)
                    cout << "   non-monotone" << endl;
                th=x;
                xn=0.5*x+0.5*thp;
                continue;
            }
            if(idf0 >= 0 && (fabs(z-z0)<0.01*ze0)) {
                if(fabs(th)<1e-10)
                    xn=th;
                else th=x;
            } else
                xn=0.5*th+0.5*x;
            info[tid].au_pvalue = pval;
            thp=x;
            z0=z;
            ze0=ze;
            idf0 = df;
            if(fabs(x-th)<1e-10) break;
        } // for step
        
        if (failed && verbose_mode >= VB_MED)
            cout << "   degenerated" << endl;
        
        if (step == max_step) {
            if (verbose_mode >= VB_MED)
                cout << "   non-convergence" << endl;
            failed = true;
        }
        
        double pchi2 = (failed) ? 0.0 : computePValueChiSquare(rss, df);
        cout << tid+1 << "\t" << info[tid].au_pvalue << "\t" << rss << "\t" << d << "\t" << c;
        
        // warning if p-value of chi-square < 0.01 (rss too high)
        if (pchi2 < 0.01)
            cout << " !!!";
        cout << endl;
    }
    
    delete [] this_bp;
    delete [] w;
    delete [] cc;
    
    cout << "Time for AU test: " << getRealTime() - start_time << " seconds" << endl;
    //    delete [] bp;
}


void evaluateTrees(istream &in, Params &params, IQTree *tree, vector<TreeInfo> &info, IntVector &distinct_ids)
{
    cout << endl;
    //MTreeSet trees(treeset_file, params.is_rooted, params.tree_burnin, params.tree_max_count);
    size_t ntrees = countDistinctTrees(in, params.is_rooted, tree, distinct_ids, params.distinct_trees);
    in.seekg(0);
    if (ntrees < distinct_ids.size()) {
        cout << "WARNING: " << distinct_ids.size() << " trees detected but only " << ntrees << " distinct trees will be evaluated" << endl;
    } else {
        cout << ntrees << (params.distinct_trees ? " distinct" : "") << " trees detected" << endl;
    }
    if (ntrees == 0) return;
    
    //if (trees.size() == 1) return;
    //string tree_file = treeset_file;
    string tree_file = params.out_prefix;
    tree_file += ".trees";
    ofstream treeout;
    //if (!params.fixed_branch_length) {
    treeout.open(tree_file.c_str());
    //}
    string score_file = params.out_prefix;
    score_file += ".treelh";
    ofstream scoreout;
    if (params.print_tree_lh)
        scoreout.open(score_file.c_str());
    string site_lh_file = params.out_prefix;
    site_lh_file += ".sitelh";
    if (params.print_site_lh) {
        ofstream site_lh_out(site_lh_file.c_str());
        site_lh_out << ntrees << " " << tree->getAlnNSite() << endl;
        site_lh_out.close();
    }
    
    if (params.print_partition_lh && !tree->isSuperTree()) {
        outWarning("-wpl does not work with non-partition model");
        params.print_partition_lh = false;
    }
    string part_lh_file = params.out_prefix;
    part_lh_file += ".partlh";
    if (params.print_partition_lh) {
        ofstream part_lh_out(part_lh_file.c_str());
        part_lh_out << ntrees << " " << ((PhyloSuperTree*)tree)->size() << endl;
        part_lh_out.close();
    }
    
    double time_start = getRealTime();
    
    int *boot_samples = NULL;
    //double *saved_tree_lhs = NULL;
    double *tree_lhs = NULL; // RELL score matrix of size #trees x #replicates
    double *pattern_lh = NULL;
    double *pattern_lhs = NULL;
    double *orig_tree_lh = NULL; // Original tree log-likelihoods
    double *max_lh = NULL;
    double *lhdiff_weights = NULL;
    size_t nptn = tree->getAlnNPattern();
    size_t maxnptn = get_safe_upper_limit(nptn);
    
    if (params.topotest_replicates && ntrees > 1) {
        size_t mem_size = (size_t)params.topotest_replicates*nptn*sizeof(int) +
        ntrees*params.topotest_replicates*sizeof(double) +
        (nptn + ntrees*3 + params.topotest_replicates*2)*sizeof(double) +
        ntrees*sizeof(TreeInfo) +
        params.do_weighted_test*(ntrees * nptn * sizeof(double) + ntrees*ntrees*sizeof(double));
        cout << "Note: " << ((double)mem_size/1024)/1024 << " MB of RAM required!" << endl;
        if (mem_size > getMemorySize()-100000)
            outWarning("The required memory does not fit in RAM!");
        cout << "Creating " << params.topotest_replicates << " bootstrap replicates..." << endl;
        if (!(boot_samples = new int [params.topotest_replicates*nptn]))
            outError(ERR_NO_MEMORY);
#ifdef _OPENMP
#pragma omp parallel if(nptn > 10000)
        {
        int *rstream;
        init_random(params.ran_seed + omp_get_thread_num(), false, &rstream);
#pragma omp for schedule(static)
#else
        int *rstream = randstream;
#endif
        for (size_t boot = 0; boot < params.topotest_replicates; boot++)
            if (boot == 0)
                tree->aln->getPatternFreq(boot_samples + (boot*nptn));
            else
                tree->aln->createBootstrapAlignment(boot_samples + (boot*nptn), params.bootstrap_spec, rstream);
#ifdef _OPENMP
        finish_random(rstream);
        }
#endif
        cout << "done" << endl;
        //if (!(saved_tree_lhs = new double [ntrees * params.topotest_replicates]))
        //    outError(ERR_NO_MEMORY);
        if (!(tree_lhs = new double [ntrees * params.topotest_replicates]))
            outError(ERR_NO_MEMORY);
        if (params.do_weighted_test || params.do_au_test) {
            if (!(lhdiff_weights = new double [ntrees * ntrees]))
                outError(ERR_NO_MEMORY);
            pattern_lhs = aligned_alloc<double>(ntrees*maxnptn);
            //            if (!(pattern_lhs = new double[ntrees* nptn]))
            //                outError(ERR_NO_MEMORY);
        }
        pattern_lh = aligned_alloc<double>(maxnptn);
        //        if (!(pattern_lh = new double[nptn]))
        //            outError(ERR_NO_MEMORY);
        if (!(orig_tree_lh = new double[ntrees]))
            outError(ERR_NO_MEMORY);
        if (!(max_lh = new double[params.topotest_replicates]))
            outError(ERR_NO_MEMORY);
    }
    int tree_index, tid, tid2;
    info.resize(ntrees);
    string saved_tree;
    saved_tree = tree->getTreeString();
    //for (MTreeSet::iterator it = trees.begin(); it != trees.end(); it++, tree_index++) {
    for (tree_index = 0, tid = 0; tree_index < distinct_ids.size(); tree_index++) {
        
        cout << "Tree " << tree_index + 1;
        if (distinct_ids[tree_index] >= 0) {
            cout << " / identical to tree " << distinct_ids[tree_index]+1 << endl;
            // ignore tree
            char ch;
            do {
                in >> ch;
            } while (!in.eof() && ch != ';');
            continue;
        }
        tree->freeNode();
        tree->readTree(in, tree->rooted);
        if (!tree->findNodeName(tree->aln->getSeqName(0))) {
            outError("Taxon " + tree->aln->getSeqName(0) + " not found in tree");
        }
        
        if (tree->rooted && tree->getModelFactory()->isReversible()) {
            if (tree->leafNum != tree->aln->getNSeq()+1)
                outError("Tree does not have same number of taxa as alignment");
            tree->convertToUnrooted();
//            cout << "convertToUnrooted" << endl;
        } else if (!tree->rooted && !tree->getModelFactory()->isReversible()) {
            if (tree->leafNum != tree->aln->getNSeq())
                outError("Tree does not have same number of taxa as alignment");
            tree->convertToRooted();
//            cout << "convertToRooted" << endl;
        }
        tree->setAlignment(tree->aln);
        tree->setRootNode(params.root);
        if (tree->isSuperTree())
            ((PhyloSuperTree*) tree)->mapTrees();
        
        tree->initializeAllPartialLh();
        tree->fixNegativeBranch(false);
        if (params.fixed_branch_length) {
            tree->setCurScore(tree->computeLikelihood());
        } else if (params.topotest_optimize_model) {
            tree->getModelFactory()->optimizeParameters(BRLEN_OPTIMIZE, false, params.modelEps);
            tree->setCurScore(tree->computeLikelihood());
        } else {
            tree->setCurScore(tree->optimizeAllBranches(100, 0.001));
        }
        treeout << "[ tree " << tree_index+1 << " lh=" << tree->getCurScore() << " ]";
        tree->printTree(treeout);
        treeout << endl;
        if (params.print_tree_lh)
            scoreout << tree->getCurScore() << endl;
        
        cout << " / LogL: " << tree->getCurScore() << endl;
        
        if (pattern_lh) {
            double curScore = tree->getCurScore();
            memset(pattern_lh, 0, maxnptn*sizeof(double));
            tree->computePatternLikelihood(pattern_lh, &curScore);
            if (params.do_weighted_test || params.do_au_test)
                memcpy(pattern_lhs + tid*maxnptn, pattern_lh, maxnptn*sizeof(double));
        }
        if (params.print_site_lh) {
            string tree_name = "Tree" + convertIntToString(tree_index+1);
            printSiteLh(site_lh_file.c_str(), tree, pattern_lh, true, tree_name.c_str());
        }
        if (params.print_partition_lh) {
            string tree_name = "Tree" + convertIntToString(tree_index+1);
            printPartitionLh(part_lh_file.c_str(), tree, pattern_lh, true, tree_name.c_str());
        }
        info[tid].logl = tree->getCurScore();
        
        if (!params.topotest_replicates || ntrees <= 1) {
            tid++;
            continue;
        }
        // now compute RELL scores
        orig_tree_lh[tid] = tree->getCurScore();
        double *tree_lhs_offset = tree_lhs + (tid*params.topotest_replicates);
        for (size_t boot = 0; boot < params.topotest_replicates; boot++) {
            double lh = 0.0;
            int *this_boot_sample = boot_samples + (boot*nptn);
            for (size_t ptn = 0; ptn < nptn; ptn++)
                lh += pattern_lh[ptn] * this_boot_sample[ptn];
            tree_lhs_offset[boot] = lh;
        }
        tid++;
    }
    
    ASSERT(tid == ntrees);
    
    if (params.topotest_replicates && ntrees > 1) {
        double *tree_probs = new double[ntrees];
        memset(tree_probs, 0, ntrees*sizeof(double));
        int *tree_ranks = new int[ntrees];
        
        /* perform RELL BP method */
        cout << "Performing RELL-BP test..." << endl;
        int *maxtid = new int[params.topotest_replicates];
        double *maxL = new double[params.topotest_replicates];
        int *maxcount = new int[params.topotest_replicates];
        memset(maxtid, 0, params.topotest_replicates*sizeof(int));
        memcpy(maxL, tree_lhs, params.topotest_replicates*sizeof(double));
        for (size_t boot = 0; boot < params.topotest_replicates; ++boot)
            maxcount[boot] = 1;
        for (tid = 1; tid < ntrees; tid++) {
            double *tree_lhs_offset = tree_lhs + (tid * params.topotest_replicates);
            for (size_t boot = 0; boot < params.topotest_replicates; ++boot)
                if (tree_lhs_offset[boot] > maxL[boot] + params.ufboot_epsilon) {
                    maxL[boot] = tree_lhs_offset[boot];
                    maxtid[boot] = tid;
                    maxcount[boot] = 1;
                } else if (tree_lhs_offset[boot] > maxL[boot] - params.ufboot_epsilon &&
                           random_double() <= 1.0/(maxcount[boot]+1)) {
                    maxL[boot] = max(maxL[boot],tree_lhs_offset[boot]);
                    maxtid[boot] = tid;
                    maxcount[boot]++;
                }
        }
        for ( size_t boot = 0; boot < params.topotest_replicates; ++boot)
            tree_probs[maxtid[boot]] += 1.0;
        for (tid = 0; tid < ntrees; tid++) {
            tree_probs[tid] /= params.topotest_replicates;
            info[tid].rell_confident = false;
            info[tid].rell_bp = tree_probs[tid];
        }
        sort_index(tree_probs, tree_probs + ntrees, tree_ranks);
        double prob_sum = 0.0;
        // obtain the confidence set
        for (tid = ntrees-1; tid >= 0; tid--) {
            info[tree_ranks[tid]].rell_confident = true;
            prob_sum += tree_probs[tree_ranks[tid]];
            if (prob_sum > 0.95) break;
        }
        
        // sanity check
        for (tid = 0, prob_sum = 0.0; tid < ntrees; tid++)
            prob_sum += tree_probs[tid];
        if (fabs(prob_sum-1.0) > 0.01)
            outError("Internal error: Wrong ", __func__);
        
        delete [] maxcount;
        delete [] maxL;
        delete [] maxtid;
        
        /* now do the SH test */
        cout << "Performing KH and SH test..." << endl;
        // SH centering step
        for (size_t boot = 0; boot < params.topotest_replicates; ++boot)
            max_lh[boot] = -DBL_MAX;
        double *avg_lh = new double[ntrees];
        for (tid = 0; tid < ntrees; tid++) {
            avg_lh[tid] = 0.0;
            double *tree_lhs_offset = tree_lhs + (tid * params.topotest_replicates);
            for (size_t boot = 0; boot < params.topotest_replicates; ++boot)
                avg_lh[tid] += tree_lhs_offset[boot];
            avg_lh[tid] /= params.topotest_replicates;
            for (size_t boot = 0; boot < params.topotest_replicates; ++boot) {
                max_lh[boot] = max(max_lh[boot], tree_lhs_offset[boot] - avg_lh[tid]);
            }
        }
        
        double orig_max_lh = orig_tree_lh[0];
        size_t orig_max_id = 0;
        double orig_2ndmax_lh = -DBL_MAX;
        size_t orig_2ndmax_id = -1;
        // find the max tree ID
        for (tid = 1; tid < ntrees; tid++)
            if (orig_max_lh < orig_tree_lh[tid]) {
                orig_max_lh = orig_tree_lh[tid];
                orig_max_id = tid;
            }
        // find the 2nd max tree ID
        for (tid = 0; tid < ntrees; tid++)
            if (tid != orig_max_id && orig_2ndmax_lh < orig_tree_lh[tid]) {
                orig_2ndmax_lh = orig_tree_lh[tid];
                orig_2ndmax_id = tid;
            }
        
        
        // SH compute p-value
        for (tid = 0; tid < ntrees; tid++) {
            double *tree_lhs_offset = tree_lhs + (tid * params.topotest_replicates);
            // SH compute original deviation from max_lh
            info[tid].kh_pvalue = 0.0;
            info[tid].sh_pvalue = 0.0;
            size_t max_id = (tid != orig_max_id) ? orig_max_id : orig_2ndmax_id;
            double orig_diff = orig_tree_lh[max_id] - orig_tree_lh[tid] - avg_lh[tid];
            double *max_kh = tree_lhs + (max_id * params.topotest_replicates);
            for (size_t boot = 0; boot < params.topotest_replicates; ++boot) {
                if (max_lh[boot] - tree_lhs_offset[boot] > orig_diff)
                    info[tid].sh_pvalue += 1.0;
                //double max_kh_here = max(max_kh[boot]-avg_lh[max_id], tree_lhs_offset[boot]-avg_lh[tid]);
                double max_kh_here = (max_kh[boot]-avg_lh[max_id]);
                if (max_kh_here - tree_lhs_offset[boot] > orig_diff)
                    info[tid].kh_pvalue += 1.0;
            }
            info[tid].sh_pvalue /= params.topotest_replicates;
            info[tid].kh_pvalue /= params.topotest_replicates;
        }
        
        if (params.do_weighted_test) {
            
            cout << "Computing pairwise logl difference variance ..." << endl;
            /* computing lhdiff_weights as 1/sqrt(lhdiff_variance) */
            for (tid = 0; tid < ntrees; tid++) {
                double *pattern_lh1 = pattern_lhs + (tid * maxnptn);
                lhdiff_weights[tid*ntrees+tid] = 0.0;
                for (tid2 = tid+1; tid2 < ntrees; tid2++) {
                    double lhdiff_variance = tree->computeLogLDiffVariance(pattern_lh1, pattern_lhs + (tid2*maxnptn));
                    lhdiff_weights[tid*ntrees+tid2] = 1.0/sqrt(lhdiff_variance);
                    lhdiff_weights[tid2*ntrees+tid] = lhdiff_weights[tid*ntrees+tid2];
                }
            }
            
            // Weighted KH and SH test
            cout << "Performing WKH and WSH test..." << endl;
            for (tid = 0; tid < ntrees; tid++) {
                double *tree_lhs_offset = tree_lhs + (tid * params.topotest_replicates);
                info[tid].wkh_pvalue = 0.0;
                info[tid].wsh_pvalue = 0.0;
                double worig_diff = -DBL_MAX;
                size_t max_id = -1;
                for (tid2 = 0; tid2 < ntrees; tid2++)
                    if (tid2 != tid) {
                        double wdiff = (orig_tree_lh[tid2] - orig_tree_lh[tid])*lhdiff_weights[tid*ntrees+tid2];
                        if (wdiff > worig_diff) {
                            worig_diff = wdiff;
                            max_id = tid2;
                        }
                    }
                for (size_t boot = 0; boot < params.topotest_replicates; ++boot) {
                    double wmax_diff = -DBL_MAX;
                    for (tid2 = 0; tid2 < ntrees; tid2++)
                        if (tid2 != tid)
                            wmax_diff = max(wmax_diff,
                                            (tree_lhs[tid2*params.topotest_replicates+boot] - avg_lh[tid2] -
                                             tree_lhs_offset[boot] + avg_lh[tid]) * lhdiff_weights[tid*ntrees+tid2]);
                    if (wmax_diff > worig_diff)
                        info[tid].wsh_pvalue += 1.0;
                    wmax_diff = (tree_lhs[max_id*params.topotest_replicates+boot] - avg_lh[max_id] -
                                 tree_lhs_offset[boot] + avg_lh[tid]);
                    if (wmax_diff >  orig_tree_lh[max_id] - orig_tree_lh[tid])
                        info[tid].wkh_pvalue += 1.0;
                }
                info[tid].wsh_pvalue /= params.topotest_replicates;
                info[tid].wkh_pvalue /= params.topotest_replicates;
            }
        }
        
        delete [] avg_lh;
        
        /* now to ELW - Expected Likelihood Weight method */
        cout << "Performing ELW test..." << endl;
        
        for (size_t boot = 0; boot < params.topotest_replicates; ++boot)
            max_lh[boot] = -DBL_MAX;
        for (tid = 0; tid < ntrees; tid++) {
            double *tree_lhs_offset = tree_lhs + (tid * params.topotest_replicates);
            for (size_t boot = 0; boot < params.topotest_replicates; ++boot)
                max_lh[boot] = max(max_lh[boot], tree_lhs_offset[boot]);
        }
        double *sumL = new double[params.topotest_replicates];
        memset(sumL, 0, sizeof(double) * params.topotest_replicates);
        for (tid = 0; tid < ntrees; tid++) {
            double *tree_lhs_offset = tree_lhs + (tid * params.topotest_replicates);
            for (size_t boot = 0; boot < params.topotest_replicates; ++boot) {
                tree_lhs_offset[boot] = exp(tree_lhs_offset[boot] - max_lh[boot]);
                sumL[boot] += tree_lhs_offset[boot];
            }
        }
        for (tid = 0; tid < ntrees; tid++) {
            double *tree_lhs_offset = tree_lhs + (tid * params.topotest_replicates);
            tree_probs[tid] = 0.0;
            for (size_t boot = 0; boot < params.topotest_replicates; ++boot) {
                tree_probs[tid] += (tree_lhs_offset[boot] / sumL[boot]);
            }
            tree_probs[tid] /= params.topotest_replicates;
            info[tid].elw_confident = false;
            info[tid].elw_value = tree_probs[tid];
        }
        
        sort_index(tree_probs, tree_probs + ntrees, tree_ranks);
        prob_sum = 0.0;
        // obtain the confidence set
        for (tid = ntrees-1; tid >= 0; tid--) {
            info[tree_ranks[tid]].elw_confident = true;
            prob_sum += tree_probs[tree_ranks[tid]];
            if (prob_sum > 0.95) break;
        }
        
        // sanity check
        for (tid = 0, prob_sum = 0.0; tid < ntrees; tid++)
            prob_sum += tree_probs[tid];
        if (fabs(prob_sum-1.0) > 0.01)
            outError("Internal error: Wrong ", __func__);
        delete [] sumL;
        
        if (params.do_au_test) {
            cout << "Performing approximately unbiased (AU) test..." << endl;
            performAUTest(params, tree, pattern_lhs, info);
        }
        
        delete [] tree_ranks;
        delete [] tree_probs;
        
    }
    delete [] max_lh;
    delete [] orig_tree_lh;
    aligned_free(pattern_lh);
    aligned_free(pattern_lhs);
    delete [] lhdiff_weights;
    delete [] tree_lhs;
    delete [] boot_samples;
    
    if (params.print_tree_lh) {
        scoreout.close();
    }
    
    treeout.close();
    
    // restore the tree
    tree->readTreeString(saved_tree);
    
    cout << "Time for evaluating all trees: " << getRealTime() - time_start << " sec." << endl;
    
}

void evaluateTrees(string treeset_file, Params &params, IQTree *tree,
                   vector<TreeInfo> &info, IntVector &distinct_ids)
{
    cout << "Reading trees in " << treeset_file << " ..." << endl;
    ifstream in(treeset_file);
    evaluateTrees(in, params, tree, info, distinct_ids);
    in.close();
}

void printTreeTestResults(vector<TreeInfo> &info, IntVector &distinct_trees,
                          IntVector &branch_ids, ostream &out, string out_file)
{

    // print header
    out << "# Test results for rooting positions on every branch" << endl;
    out << "# This file can be read in MS Excel or in R with command:" << endl;
    out << "#    dat=read.csv('" << out_file << "',comment.char='#')" << endl;
    out << "# Columns are comma-separated with following meanings:" << endl;
    out << "#    ID:      Branch ID" << endl;
    out << "#    logL:    Log-likelihood of the tree rooted at this branch" << endl;
    out << "#    deltaL:  logL difference from the maximal logl" << endl;
    if (Params::getInstance().topotest_replicates && info.size() > 1) {
        out
        << "#    bp-RELL: bootstrap proportion using RELL method (Kishino et al. 1990)" << endl
        << "#    p-KH:    p-value of one sided Kishino-Hasegawa test (1989)" << endl
        << "#    p-SH:    p-value of Shimodaira-Hasegawa test (2000)" << endl;
        if (Params::getInstance().do_weighted_test) {
            out
            << "#    p-WKH:   p-value of weighted KH test" << endl
            << "#    p-WSH:   p-value of weighted SH test." << endl;
        }
        out << "#    c-ELW:   Expected Likelihood Weight (Strimmer & Rambaut 2002)" << endl;
        if (Params::getInstance().do_au_test) {
            out
            << "#    p-AU:    p-value of approximately unbiased (AU) test (Shimodaira, 2002)" << endl;
        }
    }

    out << "ID,logL,deltaL";
    if (Params::getInstance().topotest_replicates && info.size() > 1) {
        out << ",bp-RELL,p-KH,p-SH";
        if (Params::getInstance().do_weighted_test)
            out << ",p-WKH,p-WSH";
        out << ",c-ELW";
        if (Params::getInstance().do_au_test)
            out << ",p-AU";
    }
    out << endl;
    
    double maxL = -DBL_MAX;
    int tid, orig_id;
    for (tid = 0; tid < info.size(); tid++)
        if (info[tid].logl > maxL) maxL = info[tid].logl;
    out.precision(10);
    for (orig_id = 0, tid = 0; orig_id < distinct_trees.size(); orig_id++) {
        if (distinct_trees[orig_id] >= 0) {
            continue;
        }
        out << branch_ids[orig_id];
        out << "," << info[tid].logl;
        out << "," << maxL - info[tid].logl;
        if (!Params::getInstance().topotest_replicates || info.size() <= 1) {
            out << endl;
            tid++;
            continue;
        }
        out << "," << info[tid].rell_bp;
        out << "," << info[tid].kh_pvalue;
        out << "," << info[tid].sh_pvalue;

        if (Params::getInstance().do_weighted_test) {
            out << "," << info[tid].wkh_pvalue;
            out << "," << info[tid].wsh_pvalue;
        }
        out << "," << info[tid].elw_value;

        if (Params::getInstance().do_au_test) {
            out << "," << info[tid].au_pvalue;
        }
        out << endl;
        tid++;
    }
    out << endl;

}

void printTreeTestResults(vector<TreeInfo> &info, IntVector &distinct_ids, IntVector &branch_ids, string out_file) {
    try {
        ofstream out(out_file);
        printTreeTestResults(info, distinct_ids, branch_ids, out, out_file);
        out.close();
        cout << "Tree test results printed to " << out_file << endl;
    } catch (...) {
        outError(ERR_WRITE_OUTPUT, out_file);
    }
}
