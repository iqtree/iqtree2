//
//  nnimove.h
//
//  Created by James Barbetti on 6/11/20.
//

#ifndef nnimove_h
#define nnimove_h

#include "phylonode.h"


template<class T=double> class SettingRestorer
{
protected:
    std::vector<T*> variables;
    std::vector<T>  settings;
public:
    void remember(T& rememberMe) {
        variables.emplace_back(&rememberMe);
        settings.emplace_back(rememberMe);
    }
    void restore() {
        size_t i = variables.size();
        while (0 < i) {
            --i;
            (*(variables[i])) = settings[i];
        }
        variables.clear();
        settings.clear();
    }
    ~SettingRestorer() {
        restore();
    }
};

class PhyloTree;
        
class NNIMove {
public:
    // Two nodes representing the central branch
    PhyloNode *node1, *node2;

    // Roots of the two subtree that are swapped
    NeighborVec::iterator node1Nei_it, node2Nei_it;

    // log-likelihood of the tree after applying the NNI
    double newloglh;

    int swap_id;

    // new branch lengths of 5 branches corresponding to the NNI
    DoubleVector newLen[5];

    // pattern likelihoods
    double *ptnlh;

    NNIMove();
    bool operator<(const NNIMove & rhs) const;
    
    void doSwap(PhyloTree* tree);
    void getLengths(bool nni5);
    void optimizeNNIBranches(PhyloTree* tree, bool nni5, int nni5_num_eval);
};


class NNIContext {
protected:
    PhyloTree* tree;
    SettingRestorer<PhyloNeighbor*> s1;
    SettingRestorer<double> s2;
    bool nni5;
    int IT_NUM;
    PhyloNode* node1;
    PhyloNode* node2;
    NeighborVec::iterator it;
    NeighborVec::iterator saved_it[6];
    PhyloNeighbor *saved_nei[6];

public:
    NNIContext(PhyloTree* phylo_tree, PhyloNode* firstNode, PhyloNode* secondNode);
    void resetSubtreeSizes();
    void restore();
};


#endif /* nnimove_hpp */
