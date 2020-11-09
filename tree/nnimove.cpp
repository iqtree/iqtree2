//
//  nnimove.cpp
//  alignment
//
//  Created by James Barbetti on 6/11/20.
//

#include "nnimove.h"
#include "phylotree.h"

NNIMove::NNIMove(): node1(nullptr), node2(nullptr)
    , newloglh(-DBL_MAX), ptnlh(nullptr) {
}
bool NNIMove::operator<(const NNIMove & rhs) const {
    return newloglh > rhs.newloglh;
}

void NNIMove::doSwap(PhyloTree* tree) {
    // do the NNI swap
    NeighborVec::iterator node1_it  = node1Nei_it;
    NeighborVec::iterator node2_it  = node2Nei_it;
    Neighbor*             node1_nei = *node1_it;
    Neighbor*             node2_nei = *node2_it;

    // reorient partial_lh before swap
    tree->reorientPartialLh(node1->findNeighbor(node2), node1);
    tree->reorientPartialLh(node2->findNeighbor(node1), node2);

    node1->updateNeighbor(node1_it, node2_nei);
    node2_nei->node->updateNeighbor(node2, node1);

    node2->updateNeighbor(node2_it, node1_nei);
    node1_nei->node->updateNeighbor(node1, node2);
}

void NNIMove::getLengths(bool nni5) {
    node1->findNeighbor(node2)->getLength(newLen[0]);
    if (nni5) {
        int i = 1;
        FOR_EACH_ADJACENT_PHYLO_NODE(node1, node2, it, node_X) {
            PhyloNeighbor* nei = node_X->findNeighbor(node1);
            nei->getLength(newLen[i]);
            ++i;
        }
        FOR_EACH_ADJACENT_PHYLO_NODE(node2, node1, it, node_Y) {
            PhyloNeighbor* nei = node_Y->findNeighbor(node2);
            nei->getLength(newLen[i]);
            ++i;
        }
    }
}
void NNIMove::optimizeNNIBranches(PhyloTree* tree, bool nni5, int nni5_num_eval) {
    if (nni5) {
        for (int stepsToGo = nni5_num_eval; 0 < stepsToGo; --stepsToGo) {
            FOR_EACH_ADJACENT_PHYLO_NODE(node1, node2, it, node_X) {
                tree->optimizeOneBranch(node1, node_X, false, NNI_MAX_NR_STEP);
            }
            tree->optimizeOneBranch(node1, node2, false, NNI_MAX_NR_STEP);
            FOR_EACH_ADJACENT_PHYLO_NODE(node2, node1, it, node_Y) {
                tree->optimizeOneBranch(node2, node_Y, false, NNI_MAX_NR_STEP);
            }
        }
    }
    tree->optimizeOneBranch(node1, node2, false, NNI_MAX_NR_STEP);
    newloglh = tree->computeLikelihoodFromBuffer();
    TREE_LOG_LINE(*tree, VB_DEBUG, "NNI of nodes with Ids " << node1->id << " - " << node2->id
        << ": scores " << newloglh);
    // compute the pattern likelihoods if wanted
    if (ptnlh != nullptr) {
        tree->computePatternLikelihood(ptnlh, &newloglh);
    }
}

NNIContext::NNIContext(PhyloTree* phylo_tree, PhyloNode* firstNode, PhyloNode* secondNode)
    : tree(phylo_tree), nni5(phylo_tree->params->nni5), IT_NUM(nni5 ? 6 : 2)
    , node1(firstNode), node2(secondNode) {
    s1.remember(tree->current_it);
    s1.remember(tree->current_it_back);

    TREE_LOG_LINE(*tree, VB_DEBUG, "curScore was " << tree->curScore);
    FOR_EACH_PHYLO_NEIGHBOR(node1, nullptr, it3, nei) {
        TREE_LOG_LINE(*tree, VB_DEBUG, "branch [" << nei->id << "] from node1 "
            << pointer_to_hex(nei) << " length was " << nei->length);
        s2.remember(nei->length);
        s2.remember(nei->getNode()->findNeighbor(node1)->length);
    }
    FOR_EACH_PHYLO_NEIGHBOR(node2, nullptr, it4, nei) {
        TREE_LOG_LINE(*tree, VB_DEBUG, "branch [" << nei->id << "] from node2 "
            << pointer_to_hex(nei) << " length was " << nei->length);
        s2.remember(nei->length);
        s2.remember(nei->getNode()->findNeighbor(node2)->length);
    }
    s2.remember(tree->curScore);

    IT_NUM = (tree->params->nni5) ? 6 : 2;
    
    // Upper Bounds ---------------
    /*
    if(params->upper_bound_NNI){
        totalNNIub += 2;
        NNIMove resMove = getBestNNIForBranUB(node1,node2,this);
        // if UB is smaller than the current likelihood, then we don't recompute the likelihood of the swapped topology.
        // Otherwise, follow the normal procedure: evaluate NNIs and compute the likelihood.

        // here, we skip NNI is its UB n times worse than the curLikelihood
        if( resMove.newloglh < (1+params->upper_bound_frac)*this->curScore){
            //cout << "Skipping Likelihood evaluation of NNIs for this branch :) ........................"<<endl;
            return resMove;
        }
    }
    */
    
    //-----------------------------

    int id = 0;

    saved_it[id++] = node1->findNeighborIt(node2);
    saved_it[id++] = node2->findNeighborIt(node1);

    if (nni5) {
        FOR_EACH_ADJACENT_PHYLO_NODE(node1, node2, itx, node0) {
            saved_it[id++] = node0->findNeighborIt(node1);
        }
        FOR_EACH_ADJACENT_PHYLO_NODE(node2, node1, itx, node3) {
            saved_it[id++] = node3->findNeighborIt(node2);
        }
    }
    ASSERT(id == IT_NUM);
        
    if (!nni5) {
        tree->reorientPartialLh(node1->findNeighbor(node2), node1);
        tree->reorientPartialLh(node2->findNeighbor(node1), node2);
    }

    int mem_id = 0;
    // save Neighbor and allocate new Neighbor pointer
    for (id = 0; id < IT_NUM; id++) {
        PhyloNeighbor* oldNei = (PhyloNeighbor*)(*saved_it[id]);
        saved_nei[id]         = oldNei;
        PhyloNeighbor* newNei = saved_nei[id]->newNeighbor();
        *saved_it[id]         = newNei;

        if (oldNei->partial_lh) {
            newNei->partial_lh = tree->nni_partial_lh + mem_id * tree->lh_block_size;
            newNei->scale_num  = tree->nni_scale_num  + mem_id * tree->scale_block_size;
            mem_id++;
            tree->mem_slots.addSpecialNei(newNei);
        }
    }
    if (nni5) {
        ASSERT(mem_id == 2);
    }
}

void NNIContext::resetSubtreeSizes() {
    if (tree->params->lh_mem_save == LM_MEM_SAVE) {
        // reset subtree size to change traversal order
        for (int id = 0; id < IT_NUM; id++) {
            ((PhyloNeighbor*)*saved_it[id])->size = 0;
        }
    }
}

void NNIContext::restore() {
    // restore the Neighbor*
    for (int id = IT_NUM-1; id >= 0; id--) {
        //
        //James B. 02-Nov-2020. Two commented-out lines were not necessarily valid,
        //if current_it and current_it_back were both on one of the new branches
        //after the node swap, they would NOT point at each other at all!
        //It's simplest to just... restore()... what was saved.
        //
        //if (*saved_it[id] == current_it) current_it = saved_nei[id];
        //if (*saved_it[id] == current_it_back) current_it_back = saved_nei[id];
        //
        delete (*saved_it[id]);
        (*saved_it[id]) = saved_nei[id];
    }
    //James B. 02-Nov-2020. Restore current_it and current_it_back.
    s1.restore();

    tree->mem_slots.eraseSpecialNei();

    //James B. 02-Nov-2020 Restore the current score, the lengths of *both*
    //ends of the 5 branches between and around node1 and node2
    s2.restore();
}
