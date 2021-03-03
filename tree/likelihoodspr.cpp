//
//  likelihoodspr.cpp
//  Likelihood SPR implementations
//  Created by James Barbetti on 08-Feb-2021.
//  Code in this class was copied here from phylotree.cpp
//  (and struct and class definitions were copied from
//  phylotree.h) and most of it was originally written
//  by Minh Bui.
//

#include <stdio.h>
#include "phylotree.h"

const int SPR_DEPTH = 2;

const int MAX_SPR_MOVES = 20;


/**
        an SPR move.
 */
struct SPRMove {
    PhyloNode* prune_node;
    PhyloNode* prune_dad;
    PhyloNode* regraft_node;
    PhyloNode* regraft_dad;
    double     score;
    inline SPRMove(PhyloNode* pruneNode, PhyloNode* pruneDad,
                   PhyloNode* graftDad,  PhyloNode* graftNode, double sprScore)
        : prune_node(pruneNode),   prune_dad(pruneDad)
        , regraft_node(graftNode), regraft_dad(graftDad)
        , score(sprScore) {}
};

struct SPR_compare {
    bool operator()(const SPRMove& s1, const SPRMove& s2) const {
        return s1.score > s2.score;
    }
};

class SPRMoves : public set<SPRMove, SPR_compare> {
public:
    void add(const SPRMove& move);
};

/****************************************************************************
 SPRMoves class
 ****************************************************************************/

void SPRMoves::add(const SPRMove& spr) {
    if (size() >= MAX_SPR_MOVES && spr.score <= rbegin()->score) {
        return;
    }
    if (size() >= MAX_SPR_MOVES) {
        iterator it = end();
        it--;
        erase(it);
    }
    insert(spr);
}

/****************************************************************************
 Subtree Pruning and Regrafting by maximum likelihood
 ****************************************************************************/

double PhyloTree::optimizeSPR_old(double cur_score,
                                  PhyloNode *node, PhyloNode *dad,
                                  SPRMoves& spr_moves) {
    if (!node) {
        node = getRoot();
    }
    PhyloNeighbor * dad1_nei = NULL;
    PhyloNeighbor * dad2_nei = NULL;
    PhyloNode * sibling1 = NULL;
    PhyloNode * sibling2 = NULL;
    double sibling1_len = 0.0, sibling2_len = 0.0;

    if (dad && !dad->isLeaf()) {

        ASSERT(dad->degree() == 3);
        // assign the sibling of node, with respect to dad

        FOR_EACH_PHYLO_NEIGHBOR(dad, node, it, nei) {
            if (!sibling1) {
                dad1_nei     = nei;
                sibling1     = nei->getNode();
                sibling1_len = nei->length;
            } else {
                dad2_nei     = nei;
                sibling2     = nei->getNode();
                sibling2_len = nei->length;
            }
        }
        // remove the subtree leading to node
        double sum_len = sibling1_len + sibling2_len;
        sibling1->updateNeighbor(dad, sibling2, sum_len);
        sibling2->updateNeighbor(dad, sibling1, sum_len);
        PhyloNeighbor* sibling1_nei = sibling1->findNeighbor(sibling2);
        PhyloNeighbor* sibling2_nei = sibling2->findNeighbor(sibling1);
        sibling1_nei->clearPartialLh();
        sibling2_nei->clearPartialLh();

        // now try to move the subtree to somewhere else
        PhyloNeighborVec spr_path;

        FOR_EACH_PHYLO_NEIGHBOR(sibling1, sibling2, it, nei)
        {
            spr_path.push_back(sibling1_nei);
            double score = swapSPR_old(cur_score, 1, node, dad, sibling1,
                                       sibling2, nei->getNode(), sibling1,
                                       spr_path, spr_moves);
            // if likelihood score improves, return
            if (score > cur_score)

                return score;
            spr_path.pop_back();
        }

        FOR_EACH_PHYLO_NEIGHBOR(sibling2, sibling1, it, nei)
        {
            spr_path.push_back(sibling2_nei);
            double score = swapSPR_old(cur_score, 1, node, dad, sibling1,
                                       sibling2, nei->getNode(), sibling2,
                                       spr_path, spr_moves);
            // if likelihood score improves, return
            if (score > cur_score)

                return score;
            spr_path.pop_back();
        }
        // if likelihood does not improve, swap back
        sibling1->updateNeighbor(sibling2, dad, sibling1_len);
        sibling2->updateNeighbor(sibling1, dad, sibling2_len);
        dad1_nei->node = sibling1;
        dad1_nei->length = sibling1_len;
        dad2_nei->node = sibling2;
        dad2_nei->length = sibling2_len;
        clearAllPartialLH();
    }

    FOR_EACH_ADJACENT_PHYLO_NODE(node, dad, it, child) {
        double score = optimizeSPR_old(cur_score, child, node, spr_moves);
        if (score > cur_score) return score;
    }
    return cur_score;
}

/**
 move the subtree (dad1-node1) to the branch (dad2-node2)
 */
double PhyloTree::swapSPR_old(double cur_score,      int cur_depth,
                              PhyloNode* node1,      PhyloNode* dad1,
                              PhyloNode* orig_node1, PhyloNode* orig_node2,
                              PhyloNode* node2,      PhyloNode* dad2,
                              PhyloNeighborVec &spr_path,
                              SPRMoves& spr_moves) {
    PhyloNeighbor* node1_nei      = node1->findNeighbor(dad1);
    PhyloNeighbor* dad1_nei       = dad1->findNeighbor(node1);
    double         node1_dad1_len = node1_nei->length;
    PhyloNeighbor* node2_nei_out  = node2->findNeighbor(dad2);
    
    if (dad2) {
        // now, connect (node1-dad1) to the branch (node2-dad2)

        bool first = true;
        PhyloNeighbor* node2_nei = node2->findNeighbor(dad2);
        PhyloNeighbor* dad2_nei  = dad2->findNeighbor(node2);
        double         len2      = node2_nei->length;

        FOR_EACH_PHYLO_NEIGHBOR(dad1, node1, it, nei){
            if (first) {
                nei->node = dad2;
                nei->length = len2 * 0.5;
                dad2->updateNeighbor(node2, dad1, len2 / 2);
                first = false;
            } else {
                nei->node = node2;
                nei->length = len2 * 0.5;
                node2->updateNeighbor(dad2, dad1, len2 / 2);
            }
            nei->clearPartialLh();
        }
        node2_nei->clearPartialLh();
        dad2_nei->clearPartialLh();
        node1_nei->clearPartialLh();
        for (auto it2 = spr_path.begin(); it2 != spr_path.end(); ++it2) {
            (*it2)->clearPartialLh();
        }
        
        clearAllPartialLH();
        /* testing different branch optimization */
        optimizeOneBranch(node1, dad1);
        double score = computeLikelihoodFromBuffer();
        if (score > cur_score) {
            return score;
        }
        // else, swap back
        node2->updateNeighbor(dad1, dad2, len2);
        dad2->updateNeighbor(dad1, node2, len2);
        node2_nei->clearPartialLh();
        dad2_nei->clearPartialLh();
        node1_nei->length = node1_dad1_len;
        dad1_nei->length = node1_dad1_len;

        // add to candidate SPR moves
        spr_moves.add(SPRMove(node1, dad1, node2, dad2, score));
    }
    if (cur_depth >= params->spr_radius)
    {
        return cur_score;
    }
    spr_path.push_back(node2_nei_out);

    FOR_EACH_ADJACENT_PHYLO_NODE(node2, dad2, it, child2) {
        double score = swapSPR(cur_score, cur_depth + 1, node1, dad1, orig_node1,
                               orig_node2, child2, node2, spr_path, spr_moves);
        if (score > cur_score) return score;
    }
    spr_path.pop_back();

    return cur_score;

}

double PhyloTree::optimizeSPR(double cur_score, PhyloNode *node,
                              PhyloNode *dad, SPRMoves& spr_moves) {
    if (!node) {
        node = getRoot();
    }
    PhyloNeighbor * dad1_nei = NULL;
    PhyloNeighbor * dad2_nei = NULL;
    PhyloNode * sibling1 = NULL;
    PhyloNode * sibling2 = NULL;
    double sibling1_len = 0.0, sibling2_len = 0.0;

    if (dad && !dad->isLeaf()) {

        ASSERT(dad->degree() == 3);
        // assign the sibling of node, with respect to dad

        FOR_EACH_PHYLO_NEIGHBOR(dad, node, it, nei) {
            if (!sibling1) {
                dad1_nei     = nei;
                sibling1     = nei->getNode();
                sibling1_len = nei->length;
            } else {
                dad2_nei     = nei;
                sibling2     = nei->getNode();
                sibling2_len = nei->length;
            }
        }
        // remove the subtree leading to node
        double sum_len = sibling1_len + sibling2_len;
        sibling1->updateNeighbor(dad, sibling2, sum_len);
        sibling2->updateNeighbor(dad, sibling1, sum_len);
        PhyloNeighbor* sibling1_nei = sibling1->findNeighbor(sibling2);
        PhyloNeighbor* sibling2_nei = sibling2->findNeighbor(sibling1);
        // save partial likelihood
        double* sibling1_partial_lh = sibling1_nei->partial_lh;
        double* sibling2_partial_lh = sibling2_nei->partial_lh;
        sibling1_nei->partial_lh = newPartialLh();
        sibling2_nei->partial_lh = newPartialLh();
        sibling1_nei->clearPartialLh();
        sibling2_nei->clearPartialLh();

        // now try to move the subtree to somewhere else
        PhyloNeighborVec spr_path;

        FOR_EACH_PHYLO_NEIGHBOR(sibling1, sibling2, it, nei)
        {
            spr_path.push_back(sibling1_nei);
            double score = swapSPR(cur_score, 1, node, dad, sibling1,
                                   sibling2, nei->getNode(), sibling1,
                                   spr_path, spr_moves);
            // if likelihood score improves, return
            if (score > cur_score) {
                cout << "cur_score = " << cur_score << endl;
                cout << "Found new BETTER SCORE by SPR: " << score << endl;

                return score;
            }
            spr_path.pop_back();
        }
        FOR_EACH_PHYLO_NEIGHBOR(sibling2, sibling1, it, nei)
        {
            spr_path.push_back(sibling2_nei);
            double score = swapSPR(cur_score, 1, node, dad, sibling1,
                                   sibling2, nei->getNode(), sibling2,
                    spr_path, spr_moves);
            // if likelihood score improves, return
            if (score > cur_score) {
                cout << "cur_score = " << cur_score << endl;
                cout << "Found new BETTER SCORE by SPR: " << score << endl;

                return score;
            }
            spr_path.pop_back();
        }
        // if likelihood does not imporve, swap back
        sibling1->updateNeighbor(sibling2, dad, sibling1_len);
        sibling2->updateNeighbor(sibling1, dad, sibling2_len);
        dad1_nei->node = sibling1;
        dad1_nei->length = sibling1_len;
        dad2_nei->node = sibling2;
        dad2_nei->length = sibling2_len;
        aligned_free(sibling1_nei->partial_lh); //James B. 07-Oct-2020.  Not delete[]s
        aligned_free(sibling2_nei->partial_lh); //as these came from an alligned_alloc.
                                                //see newPartialLh()
        sibling1_nei->partial_lh = sibling1_partial_lh;
        sibling2_nei->partial_lh = sibling2_partial_lh;
    }

    FOR_EACH_ADJACENT_PHYLO_NODE(node, dad, it, child){
        double score = optimizeSPR(cur_score, child, node, spr_moves);
        if (score > cur_score) return score;
    }
    return cur_score;
}

/**
 move the subtree (dad1-node1) to the branch (dad2-node2)
 */
double PhyloTree::swapSPR(double cur_score, int cur_depth,
                          PhyloNode *node1,      PhyloNode *dad1,  PhyloNode *orig_node1,
                          PhyloNode *orig_node2, PhyloNode *node2, PhyloNode *dad2,
                          PhyloNeighborVec &spr_path, SPRMoves& spr_moves) {

    PhyloNeighbor* node1_nei  = node1->findNeighbor(dad1);
    PhyloNeighbor* dad1_nei   = dad1->findNeighbor(node1);
    double node1_dad1_len     = node1_nei->length;
    PhyloNeighbor* node2_nei  = node2->findNeighbor(dad2);
    PhyloNeighbor* dad2_nei   = dad2->findNeighbor(node2);

    double* node2dad2_lh_save = node2_nei->partial_lh;
    double* dad2node2_lh_save = dad2_nei->partial_lh;
    double node2dad2_scale    = node2_nei->lh_scale_factor;
    double dad2node_scale     = dad2_nei->lh_scale_factor;

    double len2    = node2_nei->length;
    double newLen2 = sqrt(len2);

    if (dad2 && cur_depth >= SPR_DEPTH) {
        // now, connect (node1-dad1) to the branch (node2-dad2)

        bool first = true;
        //PhyloNeighbor* node2_nei = node2->findNeighbor(dad2);
        //PhyloNeighbor* dad2_nei  = dad2->findNeighbor(node2);
        //double         len2      = node2_nei->length;

        FOR_EACH_PHYLO_NEIGHBOR(dad1, node1, it, nei){
        // Finding new 2 neighbors for dad1 that are not node1
        if (first) {
            nei->node = dad2;
            //nei->length = len2 / 2;
            nei->length = newLen2;
            dad2->updateNeighbor(node2, dad1, newLen2);
            first = false;
        } else {
            nei->node = node2;
            nei->length = newLen2;
            node2->updateNeighbor(dad2, dad1, newLen2);
        }
        // clear all partial likelihood leading from
        // dad1 to the new neighbors
        nei->clearPartialLh();
    }

    // clear partial likelihood from node2 to dad1
        node2_nei->clearPartialLh();
        // clear partial likelihood from dad2 to dad1
        dad2_nei->clearPartialLh();
        // clear partial likelihood from dad1 to node1
        node1_nei->clearPartialLh();

        // set new legnth as suggested by Alexis
        node1_nei->length = 0.9;
        dad1_nei->length = 0.9;

        //Save the partial likelihood from the removal point to the insertion point
        vector<double*> saved_partial_lhs(spr_path.size());
        for (auto it2 = spr_path.begin(); it2 != spr_path.end(); ++it2) {
            saved_partial_lhs.push_back((*it2)->partial_lh);
            (*it2)->partial_lh = newPartialLh();
            (*it2)->clearPartialLh();
        }

        // optimize relevant branches
        double score;

        /* testing different branch optimization */
        optimizeOneBranch(node1, dad1);
        optimizeOneBranch(dad2, dad1);
        optimizeOneBranch(node2, dad1);
        optimizeOneBranch(orig_node1, orig_node2);
        score = computeLikelihoodFromBuffer();

        /*
         PhyloNode *cur_node = dad2;
         for (int i = spr_path.size()-1; i >= 0; i--) {
         score = optimizeOneBranch(cur_node, spr_path[i]->getNode());
         cur_node = spr_path[i]->getNode();
         }
         */
        //score = optimizeAllBranches(dad1);
        // if score improves, return
        if (score > cur_score) {
            cout << score << endl;
            return score;
        }

        // else, swap back
        node2->updateNeighbor(dad1, dad2, len2);
        dad2->updateNeighbor(dad1, node2, len2);
        //node2_nei->clearPartialLh();
        //dad2_nei->clearPartialLh();
        // restore partial likelihood vectors
        node2_nei->partial_lh = node2dad2_lh_save;
        node2_nei->lh_scale_factor = node2dad2_scale;
        dad2_nei->partial_lh = dad2node2_lh_save;
        dad2_nei->lh_scale_factor = dad2node_scale;
        node2_nei->length = len2;
        dad2_nei->length = len2;
        node1_nei->length = node1_dad1_len;
        dad1_nei->length = node1_dad1_len;
        int index = 0;
        for (auto it2 = spr_path.begin(); it2 != spr_path.end(); ++it2) {
            aligned_free((*it2)->partial_lh); //James B. 07-Oct-2020 (not delete[]).
            (*it2)->partial_lh = saved_partial_lhs.at(index);
            (*it2)->unclearPartialLh();
            ++index;
        }

        // add to candiate SPR moves
        // Tung : why adding negative SPR move ?
        spr_moves.add(SPRMove(node1, dad1, node2, dad2, score));
    }
    if (cur_depth >= params->spr_radius) {
        return cur_score;
    }
    spr_path.push_back(node2_nei);

    FOR_EACH_ADJACENT_PHYLO_NODE(node2, dad2, it, child2) {
        double score = swapSPR(cur_score, cur_depth + 1, node1, dad1,
                               orig_node1, orig_node2, child2, node2,
                               spr_path, spr_moves);
        if (score > cur_score) return score;
    }
    spr_path.pop_back();
    return cur_score;
}

double PhyloTree::assessSPRMove(double cur_score, const SPRMove &spr) {
    PhyloNode *dad = spr.prune_dad;
    PhyloNode *node = spr.prune_node;
    PhyloNode *dad2 = spr.regraft_dad;
    PhyloNode *node2 = spr.regraft_node;

    PhyloNeighbor *dad_nei1 = NULL;
    PhyloNeighbor *dad_nei2 = NULL;
    PhyloNode *sibling1 = NULL;
    PhyloNode *sibling2 = NULL;
    double sibling1_len = 0.0, sibling2_len = 0.0;

    PhyloNeighbor* node1_nei      = node->findNeighbor(dad);
    PhyloNeighbor* dad1_nei       = dad->findNeighbor(node);
    double         node1_dad1_len = node1_nei->length;

    // assign the sibling of node, with respect to dad

    FOR_EACH_PHYLO_NEIGHBOR(dad, node, it, nei) {
        if (!sibling1) {
            dad_nei1     = nei;
            sibling1     = nei->getNode();
            sibling1_len = nei->length;
        } else {

            dad_nei2     = nei;
            sibling2     = nei->getNode();
            sibling2_len = nei->length;
        }
    }
    // remove the subtree leading to node
    double sum_len = sibling1_len + sibling2_len;
    sibling1->updateNeighbor(dad, sibling2, sum_len);
    sibling2->updateNeighbor(dad, sibling1, sum_len);
    // now try to move the subtree to somewhere else

    bool first = true;
    PhyloNeighbor *node2_nei = node2->findNeighbor(dad2);
    //PhyloNeighbor *dad2_nei = dad2->findNeighbor(node2);
    double len2 = node2_nei->length;

    FOR_EACH_PHYLO_NEIGHBOR(dad, node, it, nei)
    {
        if (first) {
            nei->node = dad2;
            nei->length = len2 / 2;
            dad2->updateNeighbor(node2, dad, len2 / 2);
            first = false;
        } else {
            nei->node = node2;
            nei->length = len2 / 2;
            node2->updateNeighbor(dad2, dad, len2 / 2);
        }
        nei->clearPartialLh();
    }

    clearAllPartialLH();
    clearAllPartialParsimony(false);
    // optimize branches
    double score;
    optimizeAllBranches(dad);
    score = computeLikelihoodBranch((PhyloNeighbor*)dad->neighbors.back(), dad,
                                    tree_buffers);

    // if score improves, return
    if (score > cur_score) {
        return score;
    }
    // else, swap back
    node2->updateNeighbor(dad, dad2, len2);
    dad2->updateNeighbor(dad, node2, len2);

    node1_nei->length = node1_dad1_len;
    dad1_nei->length  = node1_dad1_len;

    sibling1->updateNeighbor(sibling2, dad, sibling1_len);
    sibling2->updateNeighbor(sibling1, dad, sibling2_len);
    dad_nei1->node    = sibling1;
    dad_nei1->length  = sibling1_len;
    dad_nei2->node    = sibling2;
    dad_nei2->length  = sibling2_len;
    clearAllPartialLH();
    clearAllPartialParsimony(false);

    return cur_score;
}

double PhyloTree::optimizeSPR() {
    fixNegativeBranch();
    SPRMoves spr_moves;
    double cur_score = computeLikelihood();
    for (int i = 0; i < 100; ++i) {
        cout << "i = " << i << endl;
        spr_moves.clear();
        double score = optimizeSPR_old(cur_score, getRoot()->firstNeighbor()->getNode(),
                                       nullptr, spr_moves);
        clearAllPartialLH();
        clearAllPartialParsimony(false);
        // why this?
        if (score <= cur_score) {
            for (auto it = spr_moves.begin(); it != spr_moves.end(); ++it) {
                //cout << (*it).score << endl;
                score = assessSPRMove(cur_score, *it);
                // if likelihood score improves, apply to SPR
                if (score > cur_score)
                    break;
            }
            if (score <= cur_score) {
                break;
            }
        } else {
            cur_score = optimizeAllBranches();
            cout << "SPR " << i + 1 << " : " << cur_score << endl;
            cur_score = score;
        }
    }
    return cur_score;
}
