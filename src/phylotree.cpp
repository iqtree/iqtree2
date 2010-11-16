//
// C++ Implementation: phylotree
//
// Description:
//
//
// Author: BUI Quang Minh, Steffen Klaere, Arndt von Haeseler <minh.bui@univie.ac.at>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "phylotree.h"
#include "bionj.h"
#include "rateheterogeneity.h"
#include "alignmentpairwise.h"
#include "acml_mv/acml_mv.h"
#include <algorithm>
#include <limits>

//const static int BINARY_SCALE = floor(log2(1/SCALING_THRESHOLD));
//const static double LOG_BINARY_SCALE = -(log(2) * BINARY_SCALE);


/****************************************************************************
        SPRMoves class
 ****************************************************************************/

void SPRMoves::add(PhyloNode *prune_node, PhyloNode *prune_dad,
        PhyloNode *regraft_node, PhyloNode *regraft_dad, double score) {
    if (size() >= MAX_SPR_MOVES && score <= rbegin()->score)
        return;
    if (size() >= MAX_SPR_MOVES) {
        iterator it = end();
        it--;
        erase(it);
    }
    SPRMove spr;
    spr.prune_node = prune_node;
    spr.prune_dad = prune_dad;
    spr.regraft_node = regraft_node;
    spr.regraft_dad = regraft_dad;
    spr.score = score;
    insert(spr);
}

/****************************************************************************
        PhyloTree class
 ****************************************************************************/

PhyloTree::PhyloTree()
: MTree() {
    aln = NULL;
    model = NULL;
    site_rate = NULL;
    optimize_by_newton = false;
    central_partial_lh = NULL;    
    central_partial_pars = NULL;
    model_factory = NULL;
    tmp_partial_lh1 = NULL;
    tmp_partial_lh2 = NULL;
}

PhyloTree::PhyloTree(Alignment *alignment)
:   ptn_freqs(alignment->size()), p_invar_ptns(alignment->size()),
        lh_ptns (alignment->size()), lh_ptns_log(alignment->size()), MTree() {
    aln = alignment;
    alnSize = aln->size();
    numStates = aln->num_states;
    tranSize = numStates * numStates;
    model = NULL;
    site_rate = NULL;
    optimize_by_newton = false;
    central_partial_lh = NULL;
    central_partial_pars = NULL;
    model_factory = NULL;
    tmp_partial_lh1 = NULL;
    tmp_partial_lh2 = NULL;
    for (int ptn = 0; ptn < alnSize; ++ptn) {
        ptn_freqs[ptn] = (*aln)[ptn].frequency;
    }
    state_freq = NULL;
}

PhyloTree::~PhyloTree() {
    if (central_partial_lh)
        //delete [] central_partial_lh;
        ei_aligned_delete<double>(central_partial_lh, (leafNum - 1) * 4 * block_size);
    central_partial_lh = NULL;
    if (central_partial_pars)
        delete [] central_partial_pars;
    central_partial_pars = NULL;
    if (model && state_freq)
        ei_aligned_delete<double>(state_freq, model->num_states);
    if (model_factory) delete model_factory;
    if (model) delete model;
    if (site_rate) delete site_rate;
//    if (p_invar_ptn)
//        ei_aligned_delete<double>(p_invar_ptn, alnSize);
//    if (ptn_freqs)
//        ei_aligned_delete<double>(ptn_freqs, alnSize);
    if (tmp_partial_lh1)
    	ei_aligned_delete<double>(tmp_partial_lh1, block_size);
        //delete [] tmp_partial_lh1;
    if (tmp_partial_lh2)
    	ei_aligned_delete<double>(tmp_partial_lh2, block_size);
        //delete [] tmp_partial_lh2;
    if (root != NULL)
        freeNode();
//    if (lh_ptns)
//        delete [] lh_ptns;
//    if (lh_ptns_log)
//        delete [] lh_ptns_log;
    root = NULL;
}

void PhyloTree::assignLeafNames(Node *node, Node *dad) {
    if (!node) node = root;
    if (node->isLeaf()) {
        node->id = atoi(node->name.c_str());
        assert(node->id >= 0 && node->id < leafNum);
        node->name = aln->getSeqName(node->id).c_str();
    }
    FOR_NEIGHBOR_IT(node, dad, it)
    assignLeafNames((*it)->node, node);
}

void PhyloTree::copyTree(MTree *tree) {
    MTree::copyTree(tree);
    if (!aln) return;
    // reset the ID with alignment
    setAlignment(aln);
}

void PhyloTree::copyPhyloTree(PhyloTree *tree) {
    MTree::copyTree(tree);
    if (!tree->aln) return;
    setAlignment(tree->aln);
}

void PhyloTree::setAlignment(Alignment *alignment) {
    aln = alignment;
    int nseq = aln->getNSeq();
    for (int seq = 0; seq < nseq; seq++) {
        string seq_name = aln->getSeqName(seq);
        Node *node = findNodeName(seq_name);
        if (!node) {
            string str = "Alignment has a sequence name ";
            str += seq_name;
            str += " which is not in the tree";
            outError(str);
        }
        assert(node->isLeaf());
        node->id = seq;
    }
    alnSize = aln->size();
    ptn_freqs.resize(alnSize);
    p_invar_ptns.resize(alnSize);
    lh_ptns.resize(alnSize);
    lh_ptns_log.resize(alnSize);
    numStates = aln->num_states;
    tranSize = numStates * numStates;
    ptn_freqs.resize(alnSize);
    p_invar_ptns.resize(alnSize);
    lh_ptns.resize(alnSize);
    lh_ptns_log.resize(alnSize);

}

void PhyloTree::rollBack(istream &best_tree_string) {
    best_tree_string.seekg(0);
    freeNode();
    readTree(best_tree_string, rooted);
    assignLeafNames();
    initializeAllPartialLh();
}

void PhyloTree::setModel(SubstModel *amodel) {
    model = amodel;
    calStateFreq();
}

void PhyloTree::setModelFactory(ModelFactory *model_fac) {
    model_factory = model_fac;
}

void PhyloTree::setRate(RateHeterogeneity *rate) {
    site_rate = rate;
    if (!rate) return;
    numCat = site_rate->getNRate();
    p_invar = site_rate->getPInvar();
    p_var_cat = (1.0 - p_invar) / (double) numCat;
    if (aln) {
        block = aln->num_states * numCat;
        lh_size = aln->size() * block;
    }
}

RateHeterogeneity *PhyloTree::getRate() {
    return site_rate;
}

Node* PhyloTree::newNode(int node_id, const char* node_name) {
    return (Node*) (new PhyloNode(node_id, node_name));
}

Node* PhyloTree::newNode(int node_id, int node_name) {
    return (Node*) (new PhyloNode(node_id, node_name));
}

void PhyloTree::clearAllPartialLh() {
    ((PhyloNode*) root->neighbors[0]->node)->clearAllPartialLh((PhyloNode*) root);
}

string PhyloTree::getModelName() {
    return model->name + site_rate->name;
}

/****************************************************************************
        Parsimony function
 ****************************************************************************/



int PhyloTree::getBitsBlockSize() {
    // reserve the last entry for parsimony score
    return (aln->num_states * aln->size() + UINT_BITS - 1) / UINT_BITS + 1;
}

UINT *PhyloTree::newBitsBlock() {
    return new UINT[getBitsBlockSize()];
}

UINT PhyloTree::getBitsBlock(UINT *bit_vec, int index) {
    int nstates = aln->num_states;
    int myindex = (index * nstates);
    int bit_pos_begin = myindex >> BITS_DIV;
    int bit_off_begin = myindex & BITS_MODULO;
    int bit_pos_end;
    int bit_off_end = bit_off_begin + nstates;
    if (bit_off_end > UINT_BITS) {
        bit_off_end -= UINT_BITS;
        bit_pos_end = bit_pos_begin + 1;
    } else
        bit_pos_end = bit_pos_begin;


    if (bit_pos_begin == bit_pos_end)
        return (bit_vec[bit_pos_begin] >> bit_off_begin) & ((1 << nstates) - 1);
    else {
        int part1 = (bit_vec[bit_pos_begin] >> bit_off_begin);
        int part2 = bit_vec[bit_pos_end] & ((1 << bit_off_end) - 1);
        return part1 | (part2 << (UINT_BITS - bit_off_begin));
    }
}

void PhyloTree::setBitsBlock(UINT *bit_vec, int index, UINT value) {
    int nstates = aln->num_states;
    int myindex = (index * nstates);
    int bit_pos_begin = myindex >> BITS_DIV;
    int bit_off_begin = myindex & BITS_MODULO;
    int bit_pos_end;
    int bit_off_end = bit_off_begin + nstates;
    if (bit_off_end > UINT_BITS) {
        bit_off_end -= UINT_BITS;
        bit_pos_end = bit_pos_begin + 1;
    } else
        bit_pos_end = bit_pos_begin;

    int allstates = (1 << nstates) - 1;
    assert(value <= allstates);

    if (bit_pos_begin == bit_pos_end) {
        // first clear the bit between bit_off_begin and bit_off_end
        bit_vec[bit_pos_begin] &= ~(allstates << bit_off_begin);
        // now set the bit
        bit_vec[bit_pos_begin] |= value << bit_off_begin;
    } else {
        int len1 = UINT_BITS - bit_off_begin;
        int allbit1 = (1 << len1) - 1;
        // clear bit from bit_off_begin to UINT_BITS
        bit_vec[bit_pos_begin] &= ~(allbit1 << bit_off_begin);
        // set bit  from bit_off_begin to UINT_BITS
        bit_vec[bit_pos_begin] |= ((value & allbit1) << bit_off_begin);

        // clear bit from 0 to bit_off_end
        bit_vec[bit_pos_end] &= ~((1 << bit_off_end) - 1);
        // now set the bit the value
        bit_vec[bit_pos_end] |= (value >> len1);

    }
}

void PhyloTree::computePartialParsimony(PhyloNeighbor *dad_branch, PhyloNode *dad) {
    // don't recompute the likelihood
    if (dad_branch->partial_lh_computed & 2)
        return;
    Node *node = dad_branch->node;
    assert(node->degree() <= 3);
    int ptn;
    int nstates = aln->num_states;
    int pars_size = getBitsBlockSize();
    assert(dad_branch->partial_pars);

    if (node->isLeaf() && dad) {
        // external node
        memset(dad_branch->partial_pars, 0, pars_size * sizeof (int));
        for (ptn = 0; ptn < aln->size(); ptn++) {
            char state;
            if (node->name == ROOT_NAME) {
                state = STATE_UNKNOWN;
            } else {
                assert(node->id < aln->getNSeq());
                state = (aln->at(ptn))[node->id];
            }
            if (state == STATE_UNKNOWN) {
                // fill all entries with bit 1
                setBitsBlock(dad_branch->partial_pars, ptn, (1 << nstates) - 1);
            } else if (state < nstates) {
                setBitsBlock(dad_branch->partial_pars, ptn, 1 << state);
            } else {
                // ambiguous character, for DNA, RNA
                state = state - (nstates - 1);
                setBitsBlock(dad_branch->partial_pars, ptn, state);
            }
        }
    } else {
        // internal node
        //memset(dad_branch->partial_pars, 127, pars_size * sizeof(int));
        UINT *partial_pars_dad = dad_branch->partial_pars;
        UINT *partial_pars_child1 = NULL, *partial_pars_child2 = NULL;
        // take the intersection of two child states (with &= bit operation)
        FOR_NEIGHBOR_IT(node, dad, it)
        if ((*it)->node->name != ROOT_NAME) {
            computePartialParsimony((PhyloNeighbor*) (*it), (PhyloNode*) node);
            if (!partial_pars_child1)
                partial_pars_child1 = ((PhyloNeighbor*) (*it))->partial_pars;
            else
                partial_pars_child2 = ((PhyloNeighbor*) (*it))->partial_pars;
        }
        assert(partial_pars_child1 && partial_pars_child2);
        for (int i = 0; i < pars_size - 1; i++)
            partial_pars_dad[i] = partial_pars_child1[i] & partial_pars_child2[i];
        int partial_pars = partial_pars_child1[pars_size - 1] + partial_pars_child2[pars_size - 1];
        // now check if some intersection is empty, change to union (Fitch algorithm) and increase the parsimony score
        for (ptn = 0; ptn < aln->size(); ptn++)
            if (getBitsBlock(partial_pars_dad, ptn) == 0) {
                setBitsBlock(partial_pars_dad, ptn, getBitsBlock(partial_pars_child1, ptn) | getBitsBlock(partial_pars_child2, ptn));
                partial_pars += aln->at(ptn).frequency;
            }
        partial_pars_dad[pars_size - 1] = partial_pars;
    }
    dad_branch->partial_lh_computed |= 2;
}

int PhyloTree::computeParsimonyBranch(PhyloNeighbor *dad_branch, PhyloNode *dad) {
    PhyloNode *node = (PhyloNode*) dad_branch->node;
    PhyloNeighbor *node_branch = (PhyloNeighbor*) node->findNeighbor(dad);
    assert(node_branch);
    if (!central_partial_pars)
        initializeAllPartialLh();
    // swap node and dad if dad is a leaf
    if (node->isLeaf()) {
        PhyloNode *tmp_node = dad;
        dad = node;
        node = tmp_node;
        PhyloNeighbor *tmp_nei = dad_branch;
        dad_branch = node_branch;
        node_branch = tmp_nei;
        //cout << "swapped\n";
    }
    if ((dad_branch->partial_lh_computed & 2) == 0)
        computePartialParsimony(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 2) == 0)
        computePartialParsimony(node_branch, node);
    // now combine likelihood at the branch

    int pars_size = getBitsBlockSize();
    int nstates = aln->num_states;
    int i, ptn;
    int tree_pars = node_branch->partial_pars[pars_size - 1] + dad_branch->partial_pars[pars_size - 1];
    UINT *partial_pars = newBitsBlock();
    for (i = 0; i < pars_size - 1; i++)
        partial_pars[i] = (node_branch->partial_pars[i] & dad_branch->partial_pars[i]);

    for (ptn = 0; ptn < aln->size(); ptn++) {
        /*
        for (i = 0; i < aln->getNSeq(); i++)
                cout << aln->convertStateBack(aln->at(ptn)[i]);
        cout << endl;*/

        if (getBitsBlock(partial_pars, ptn) == 0)
            tree_pars += aln->at(ptn).frequency;
    }
    delete [] partial_pars;
    return tree_pars;
}

int PhyloTree::computeParsimony() {
    return computeParsimonyBranch((PhyloNeighbor*) root->neighbors[0], (PhyloNode*) root);
}

int PhyloTree::computeParsimonyScore(int ptn, int &states, PhyloNode *node, PhyloNode *dad) {
    int score = 0;
    states = 0;
    if (!node) node = (PhyloNode*) root;
    if (node->degree() > 3)
        outError("Does not work with multifurcating tree");
    if (verbose_mode == VB_DEBUG)
        cout << ptn << " " << node->id << "  " << node->name << endl;

    if (node->isLeaf()) {
        char state;
        if (node->name == ROOT_NAME) {
            state = STATE_UNKNOWN;
        } else {
            assert(node->id < aln->getNSeq());
            state = (*aln)[ptn][node->id];
        }
        if (state == STATE_UNKNOWN) {
            states = (1 << aln->num_states) - 1;
        } else if (state < aln->num_states)
            states = (1 << state);
        else {
            // ambiguous character, for DNA, RNA
            states = state - 3;
        }
    }
    if (!node->isLeaf() || node == root) {
        int union_states = 0;
        int intersect_states = (1 << aln->num_states) - 1;
        if (states != 0) {
            union_states = states;
            intersect_states = states;
        }

        FOR_NEIGHBOR_IT(node, dad, it) {
            int states_child;
            int score_child = computeParsimonyScore(ptn, states_child, (PhyloNode*) ((*it)->node), node);
            union_states |= states_child;
            intersect_states &= states_child;
            score += score_child;
        }
        if (intersect_states)
            states = intersect_states;
        else {
            states = union_states;
            score++;
        }
    }
    return score;
}

int PhyloTree::computeParsimonyScore() {
    assert(root && root->isLeaf());

    int score = 0;
    for (int ptn = 0; ptn < aln->size(); ptn++)
        if (!aln->at(ptn).is_const) {
            int states;
            score += computeParsimonyScore(ptn, states) * (*aln)[ptn].frequency;
        }
    return score;
}

/****************************************************************************
        Nearest Neighbor Interchange with parsimony
 ****************************************************************************/

double PhyloTree::swapNNI(double cur_score, PhyloNode *node1, PhyloNode *node2) {
    assert(node1->degree() == 3 && node2->degree() == 3);
    FOR_NEIGHBOR_DECLARE(node1, node2, it1)
    break;
    Node *node1_nei = (*it1)->node;

    FOR_NEIGHBOR_IT(node2, node1, it2) {
        // do the NNI swap
        Node *node2_nei = (*it2)->node;
        node1->updateNeighbor(node1_nei, node2_nei);
        node1_nei->updateNeighbor(node1, node2);
        node2->updateNeighbor(node2_nei, node1_nei);
        node2_nei->updateNeighbor(node2, node1);

        // compute the score of the swapped topology
        double score = computeParsimonyScore();
        // if better: return
        if (score < cur_score) return score;
        // else, swap back
        node1->updateNeighbor(node2_nei, node1_nei);
        node1_nei->updateNeighbor(node2, node1);
        node2->updateNeighbor(node1_nei, node2_nei);
        node2_nei->updateNeighbor(node1, node2);
    }
    return cur_score;
}

double PhyloTree::searchNNI(double cur_score, PhyloNode *node, PhyloNode *dad) {
    if (!node) node = (PhyloNode*) root;
    if (!node->isLeaf() && dad && !dad->isLeaf()) {
        double score = swapNNI(cur_score, node, dad);
        if (score < cur_score) return score;
    }

    FOR_NEIGHBOR_IT(node, dad, it) {
        double score = searchNNI(cur_score, (PhyloNode*) (*it)->node, node);
        if (score < cur_score) return score;
    }
    return cur_score;
}

void PhyloTree::searchNNI() {
    cout << "Search with Nearest Neighbor Interchange..." << endl;
    double cur_score = computeParsimonyScore();
    do {
        double score = searchNNI(cur_score);
        if (score >= cur_score) break;
        cout << "Better score found: " << score << endl;
        cur_score = score;
    } while (true);
}

/****************************************************************************
        Stepwise addition (greedy) by maximum parsimony
 ****************************************************************************/

int PhyloTree::addTaxonMP(Node *added_node, Node* &target_node, Node* &target_dad, Node *node, Node *dad) {
    Neighbor *dad_nei = dad->findNeighbor(node);

    // now insert the new node in the middle of the branch node-dad
    double len = dad_nei->length;
    node->updateNeighbor(dad, added_node, len / 2.0);
    dad->updateNeighbor(node, added_node, len / 2.0);
    added_node->updateNeighbor((Node*) 1, node, len / 2.0);
    added_node->updateNeighbor((Node*) 2, dad, len / 2.0);
    // compute the likelihood
    //clearAllPartialLh();
    int best_score = computeParsimonyScore();
    target_node = node;
    target_dad = dad;
    // remove the added node
    node->updateNeighbor(added_node, dad, len);
    dad->updateNeighbor(added_node, node, len);
    added_node->updateNeighbor(node, (Node*) 1, len);
    added_node->updateNeighbor(dad, (Node*) 2, len);

    // now tranverse the tree downwards

    FOR_NEIGHBOR_IT(node, dad, it) {
        Node *target_node2;
        Node *target_dad2;
        double score = addTaxonMP(added_node, target_node2, target_dad2, (*it)->node, node);
        if (score < best_score) {
            best_score = score;
            target_node = target_node2;
            target_dad = target_dad2;
        }
    }
    return best_score;
}

void PhyloTree::growTreeMP(Alignment *alignment) {

    cout << "Stepwise addition using maximum parsimony..." << endl;
    aln = alignment;
    int size = aln->getNSeq();
    if (size < 3)
        outError(ERR_FEW_TAXA);

    root = newNode();
    Node *new_taxon;

    // create initial tree with 3 taxa
    for (leafNum = 0; leafNum < 3; leafNum++) {
        if (verbose_mode >= VB_MAX)
            cout << "Add " << aln->getSeqName(leafNum) << " to the tree" << endl;
        new_taxon = newNode(leafNum, aln->getSeqName(leafNum).c_str());
        root->addNeighbor(new_taxon, 1.0);
        new_taxon->addNeighbor(root, 1.0);
    }
    root = findNodeID(0);
    //optimizeAllBranches();

    // stepwise adding the next taxon
    for (leafNum = 3; leafNum < size; leafNum++) {
        if (verbose_mode >= VB_MAX)
            cout << "Add " << aln->getSeqName(leafNum) << " to the tree";
        // allocate a new taxon and a new ajedcent internal node
        new_taxon = newNode(leafNum, aln->getSeqName(leafNum).c_str());
        Node *added_node = newNode();
        added_node->addNeighbor(new_taxon, 1.0);
        new_taxon->addNeighbor(added_node, 1.0);

        // preserve two neighbors
        added_node->addNeighbor((Node*) 1, 1.0);
        added_node->addNeighbor((Node*) 2, 1.0);


        Node *target_node = NULL;
        Node *target_dad = NULL;
        int score = addTaxonMP(added_node, target_node, target_dad, root->neighbors[0]->node, root);
        if (verbose_mode >= VB_MAX)
            cout << ", score = " << score << endl;
        // now insert the new node in the middle of the branch node-dad
        double len = target_dad->findNeighbor(target_node)->length;
        target_node->updateNeighbor(target_dad, added_node, len / 2.0);
        target_dad->updateNeighbor(target_node, added_node, len / 2.0);
        added_node->updateNeighbor((Node*) 1, target_node, len / 2.0);
        added_node->updateNeighbor((Node*) 2, target_dad, len / 2.0);
        // compute the likelihood
        //clearAllPartialLh();
        //optimizeAllBranches();
        //optimizeNNI();
    }

    nodeNum = 2 * leafNum - 2;
}

/****************************************************************************
        likelihood function
 ****************************************************************************/

void PhyloTree::initializeAllPartialLh() {
    int index;
    block_size = alnSize * numStates * site_rate->getNRate();
    tmp_partial_lh1 = newPartialLh();
    tmp_partial_lh2 = newPartialLh();
    initializeAllPartialLh(index);
    assert(index == (nodeNum - 1)*2);
}

void PhyloTree::initializeAllPartialLh(int &index, PhyloNode *node, PhyloNode *dad) {
    int pars_block_size = getBitsBlockSize();
    if (!node) {
        node = (PhyloNode*) root;
        // allocate the big central partial likelihoods memory
        if (!central_partial_lh) {
            if (verbose_mode >= VB_MED)
                cout << "Allocating " << (leafNum - 1) * 4 * block_size * sizeof (double) << " bytes for partial likelihood vectors" << endl;
            //central_partial_lh = new double[(leafNum-1)*4*block_size];
            central_partial_lh = ei_aligned_new<double>((leafNum - 1) * 4 * block_size);
            if (!central_partial_lh)
                outError("Not enough memory for partial likelihood vectors");

        }
        if (!central_partial_pars) {
            if (verbose_mode >= VB_MED)
                cout << "Allocating " << (leafNum - 1)*4 * pars_block_size * sizeof (UINT) << " bytes for partial parsimony vectors" << endl;
            //central_partial_pars = new UINT[(leafNum-1)*4*pars_block_size];
            central_partial_pars = ei_aligned_new<UINT > ((leafNum - 1)*4 * pars_block_size);
            if (!central_partial_pars)
                outError("Not enough memory for partial parsimony vectors");
        }
        index = 0;
    }
    if (dad) {
        // assign a region in central_partial_lh to both Neihgbors (dad->node, and node->dad)
        PhyloNeighbor *nei = (PhyloNeighbor*) node->findNeighbor(dad);
        //assert(!nei->partial_lh);
        nei->partial_lh = central_partial_lh + (index * block_size);
        nei->partial_pars = central_partial_pars + (index * pars_block_size);
        nei = (PhyloNeighbor*) dad->findNeighbor(node);
        //assert(!nei->partial_lh);
        nei->partial_lh = central_partial_lh + ((index + 1) * block_size);
        nei->partial_pars = central_partial_pars + ((index + 1) * pars_block_size);
        index += 2;
        assert(index < nodeNum * 2 - 1);
    }
    FOR_NEIGHBOR_IT(node, dad, it)
    initializeAllPartialLh(index, (PhyloNode*) (*it)->node, node);
}

double *PhyloTree::newPartialLh() {
    //return new double[aln->size() * aln->num_states * site_rate->getNRate()];
    return ei_aligned_new<double>(block_size);
}

double PhyloTree::computeLikelihood(double *pattern_lh) {
    assert(model);
    assert(site_rate);
    assert(root->isLeaf());
    double *ptn_scale = NULL;
    int nptn = aln->getNPattern();
    if (pattern_lh) {
        double sum_scaling = ((PhyloNeighbor*) root->neighbors[0])->lh_scale_factor;
        if (sum_scaling < 0.0) {
            ptn_scale = new double[nptn];
            memset(ptn_scale, 0, sizeof (double) * nptn);
            //clearAllPartialLh();
            ((PhyloNeighbor*) root->neighbors[0])->clearForwardPartialLh(root);
            computePartialLikelihood((PhyloNeighbor*) root->neighbors[0], (PhyloNode*) root, ptn_scale);
            assert(fabs(sum_scaling - ((PhyloNeighbor*) root->neighbors[0])->lh_scale_factor) < 1e-6);
        }
    }
    double score = computeLikelihoodBranch((PhyloNeighbor*) root->neighbors[0], (PhyloNode*) root, pattern_lh);
    if (ptn_scale) {
        double check_score = 0.0;
        for (int i = 0; i < nptn; i++) {
            pattern_lh[i] += ptn_scale[i];
            check_score += pattern_lh[i] * aln->at(i).frequency;
        }
        delete [] ptn_scale;
        assert(fabs(score - check_score) < 1e-6);
    }
    return score;
}

double PhyloTree::computeLikelihoodBranch(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_lh) {
    if (sse) {
        switch (aln->num_states) {
            case 2: return computeLikelihoodBranchSSE < 2 > (dad_branch, dad, pattern_lh);
            case 4: return computeLikelihoodBranchSSE < 4 > (dad_branch, dad, pattern_lh);
            case 20: return computeLikelihoodBranchSSE < 20 > (dad_branch, dad, pattern_lh);
        }
    } else {
        return computeLikelihoodBranchNaive(dad_branch, dad, pattern_lh);
    }

}

double PhyloTree::computeLikelihoodBranchNaive(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_lh) {
    PhyloNode *node = (PhyloNode*) dad_branch->node;
    PhyloNeighbor *node_branch = (PhyloNeighbor*) node->findNeighbor(dad);
    assert(node_branch);
    if (!central_partial_lh)
        initializeAllPartialLh();
    // swap node and dad if dad is a leaf
    if (node->isLeaf()) {
        PhyloNode *tmp_node = dad;
        dad = node;
        node = tmp_node;
        PhyloNeighbor *tmp_nei = dad_branch;
        dad_branch = node_branch;
        node_branch = tmp_nei;
        //cout << "swapped\n";
    }
    if ((dad_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihood(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihood(node_branch, node);
    // now combine likelihood at the branch

    double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;
    int ncat = site_rate->getNRate();
    double p_invar = site_rate->getPInvar();
    double p_var_cat = (1.0 - p_invar) / (double) ncat;
    int nstates = aln->num_states;
    int block = ncat * nstates;
    int trans_size = nstates * nstates;
    int ptn, cat, state1, state2;

    double trans_mat[ncat * trans_size];
    double state_freq[nstates];
    model->getStateFrequency(state_freq);

    for (cat = 0; cat < ncat; cat++) {
        //trans_mat[cat] = model->newTransMatrix();
        double *trans_cat = trans_mat + (cat * trans_size);
        model_factory->computeTransMatrix(dad_branch->length * site_rate->getRate(cat), trans_cat);
        for (state1 = 0; state1 < nstates; state1++) {
            double *trans_mat_state = trans_cat + (state1 * nstates);
            for (state2 = 0; state2 < nstates; state2++)
                trans_mat_state[state2] *= state_freq[state1];
        }
    }

    for (ptn = 0; ptn < aln->size(); ptn++) {
        double lh_ptn = 0.0; // likelihood of the pattern
        int dad_state = 1000; // just something big enough
        if (dad->isLeaf())
            dad_state = (*aln)[ptn][dad->id];
        for (cat = 0; cat < ncat; cat++) {
            double *partial_lh_site = node_branch->partial_lh + (ptn * block + cat * nstates);
            double *partial_lh_child = dad_branch->partial_lh + (ptn * block + cat * nstates);
            if (dad_state < nstates) { // single state
                // external node
                double *trans_state = trans_mat + (cat * trans_size + dad_state * nstates);
                for (state2 = 0; state2 < nstates; state2++)
                    lh_ptn += partial_lh_child[state2] * trans_state[state2];
            } else {
                // internal node, or external node but ambiguous character
                for (state1 = 0; state1 < nstates; state1++) {
                    double lh_state = 0.0; // likelihood of state1
                    double *trans_state = trans_mat + (cat * trans_size + state1 * nstates);
                    for (state2 = 0; state2 < nstates; state2++)
                        lh_state += partial_lh_child[state2] * trans_state[state2];
                    lh_ptn += lh_state * partial_lh_site[state1];
                }
            }
        }
        lh_ptn *= p_var_cat;
        if ((*aln)[ptn].is_const && (*aln)[ptn][0] < nstates) {
            lh_ptn += p_invar * state_freq[(*aln)[ptn][0]];
        }
        assert(lh_ptn > 0);
        lh_ptn = log(lh_ptn);
        tree_lh += lh_ptn * (*aln)[ptn].frequency;
        if (pattern_lh) pattern_lh[ptn] = lh_ptn;
    }
    //for (cat = ncat-1; cat >= 0; cat--)
    //delete trans_mat[cat];
    //delete state_freq;
    return tree_lh;
}

template<int NSTATES>
inline double PhyloTree::computeLikelihoodBranchSSE(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_lh) {
    PhyloNode *node = (PhyloNode*) dad_branch->node; // Node A
    PhyloNeighbor *node_branch = (PhyloNeighbor*) node->findNeighbor(dad); // Node B
    assert(node_branch);
    if (!central_partial_lh)
        initializeAllPartialLh();
    // swap node and dad if dad is a leaf
    if (node->isLeaf()) {
        PhyloNode *tmp_node = dad;
        dad = node;
        node = tmp_node;
        PhyloNeighbor *tmp_nei = dad_branch;
        dad_branch = node_branch;
        node_branch = tmp_nei;
        //cout << "swapped\n";
    }
    if ((dad_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihoodSSE<NSTATES > (dad_branch, dad);
    if ((node_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihoodSSE<NSTATES > (node_branch, node);
    // now combine likelihood at the branch

    double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;
    int ptn, cat, state1, state2;
    double *partial_lh_site;
    double *partial_lh_child;
    double *trans_state;

    EIGEN_ALIGN16 double trans_mat[numCat * tranSize];
    EIGEN_ALIGN16 double state_freq[NSTATES];
    model->getStateFrequency(state_freq);

    for (cat = 0; cat < numCat; cat++) {
        double *trans_cat = trans_mat + (cat * tranSize);
        model_factory->computeTransMatrix(dad_branch->length * site_rate->getRate(cat), trans_cat);
        for (state1 = 0; state1 < NSTATES; state1++) {
            double *trans_mat_state = trans_cat + (state1 * NSTATES);
            for (state2 = 0; state2 < NSTATES; state2++)
                trans_mat_state[state2] *= state_freq[state1];
        }
    }    
    for (ptn = 0; ptn < alnSize; ++ptn) {
        double lh_ptn = 0.0; // likelihood of the pattern
        for (cat = 0; cat < numCat; cat++) {
            partial_lh_site = node_branch->partial_lh + (ptn * block + cat * NSTATES);
            partial_lh_child = dad_branch->partial_lh + (ptn * block + cat * NSTATES);
            trans_state = trans_mat + cat * tranSize;
            Map<Matrix<double, 1, NSTATES>, Aligned> eigen_partial_lh_child(&partial_lh_child[0]);
            Map<Matrix<double, 1, NSTATES>, Aligned> eigen_partial_lh_site(&partial_lh_site[0]);
            Map<Matrix<double, NSTATES, NSTATES>, Aligned> eigen_trans_state(&trans_state[0]);
            lh_ptn += (eigen_partial_lh_child * eigen_trans_state).dot(eigen_partial_lh_site);
        }

        lh_ptn *= p_var_cat;
        lh_ptn += p_invar_ptns[ptn];
        lh_ptns[ptn] = lh_ptn;
        //tree_lh += log(lh_ptn) * ptn_freqs[ptn];
        if (pattern_lh) pattern_lh[ptn] = lh_ptn;
    }
    vrda_log(alnSize, lh_ptns.data(), lh_ptns_log.data());
    tree_lh += (lh_ptns_log * ptn_freqs).sum();
    return tree_lh;
}

void PhyloTree::computePartialLikelihoodNaive(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_scale) {
    // don't recompute the likelihood
    if (dad_branch->partial_lh_computed & 1)
        return;
    Node *node = dad_branch->node;
    int ptn, cat;
    int ncat = site_rate->getNRate();
    int nstates = aln->num_states;
    int block = nstates * site_rate->getNRate();
    int lh_size = aln->size() * block;
    double *partial_lh_site;

    dad_branch->lh_scale_factor = 0.0;

    assert(dad_branch->partial_lh);
    //if (!dad_branch->partial_lh)
    //	dad_branch->partial_lh = newPartialLh();
    if (node->isLeaf() && dad) {
        /* external node */
        memset(dad_branch->partial_lh, 0, lh_size * sizeof (double));
        for (ptn = 0; ptn < aln->size(); ptn++) {
            char state;
            partial_lh_site = dad_branch->partial_lh + (ptn * block);
            if (node->name == ROOT_NAME) {
                state = STATE_UNKNOWN;
            } else {
                assert(node->id < aln->getNSeq());
                state = (aln->at(ptn))[node->id];
            }
            if (state == STATE_UNKNOWN) {
                // fill all entries (also over rate category) with 1.0
                for (int state2 = 0; state2 < block; state2++) {
                    partial_lh_site[state2] = 1.0;
                }
            } else if (state < nstates) {
                for (cat = 0; cat < ncat; cat++)
                    partial_lh_site[cat * nstates + state] = 1.0;
            } else {
                // ambiguous character, for DNA, RNA
                state = state - (nstates - 1);
                for (int state2 = 0; state2 < nstates; state2++)
                    if (state & (1 << state2)) {
                        for (cat = 0; cat < ncat; cat++)
                            partial_lh_site[cat * nstates + state2] = 1.0;
                    }
            }
        }
    } else {
        /* internal node */
        double *trans_mat[ncat];
        for (cat = 0; cat < ncat; cat++) trans_mat[cat] = model->newTransMatrix();
        for (ptn = 0; ptn < lh_size; ptn++)
            dad_branch->partial_lh[ptn] = 1.0;
        FOR_NEIGHBOR_IT(node, dad, it)
        if ((*it)->node->name != ROOT_NAME) {
            computePartialLikelihoodNaive((PhyloNeighbor*) (*it), (PhyloNode*) node, pattern_scale);

            dad_branch->lh_scale_factor += ((PhyloNeighbor*) (*it))->lh_scale_factor;

            for (cat = 0; cat < ncat; cat++)
                model_factory->computeTransMatrix((*it)->length * site_rate->getRate(cat), trans_mat[cat]);

            for (ptn = 0; ptn < aln->size(); ptn++) {
                for (cat = 0; cat < ncat; cat++) {
                    partial_lh_site = dad_branch->partial_lh + (ptn * block + cat * nstates);
                    double *partial_lh_child = ((PhyloNeighbor*) (*it))->partial_lh + (ptn * block + cat * nstates);
                    for (int state = 0; state < nstates; state++) {
                        double lh_child = 0.0;
                        double *trans_state = trans_mat[cat] + (state * nstates);
                        for (int state2 = 0; state2 < nstates; state2++)
                            lh_child += trans_state[state2] * partial_lh_child[state2];

                        partial_lh_site[state] *= lh_child;
                    }
                }
                // check if one should scale partial likelihoods
                bool do_scale = true;
                partial_lh_site = dad_branch->partial_lh + (ptn * block);
                for (cat = 0; cat < block; cat++)
                    if (partial_lh_site[cat] > SCALING_THRESHOLD) {
                        do_scale = false;
                        break;
                    }
                if (!do_scale) continue;
                // now do the likelihood scaling
                /*
                double lh_max = partial_lh_site[0];
                for (cat = 1; cat < block; cat++)
                        if (lh_max < partial_lh_site[cat]) lh_max = partial_lh_site[cat];
                for (cat = 0; cat < block; cat++)
                        partial_lh_site[cat] /= lh_max;
                dad_branch->lh_scale_factor += log(lh_max) * (*aln)[ptn].frequency;
				
                 */
                for (cat = 0; cat < block; cat++)
                    partial_lh_site[cat] /= SCALING_THRESHOLD;
                dad_branch->lh_scale_factor += LOG_SCALING_THRESHOLD * (*aln)[ptn].frequency;
                if (pattern_scale)
                    pattern_scale[ptn] += LOG_SCALING_THRESHOLD;
            }
        }
        for (cat = ncat - 1; cat >= 0; cat--)
            delete [] trans_mat[cat];
    }
    dad_branch->partial_lh_computed |= 1;

}

void PhyloTree::computePartialLikelihood(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_scale) {
    if (sse) {
        switch (aln->num_states) {
            case 2: return computePartialLikelihoodSSE < 2 > (dad_branch, dad, pattern_scale);
            case 4: return computePartialLikelihoodSSE < 4 > (dad_branch, dad, pattern_scale);
            //case 4: return computePartialLikelihoodNaive(dad_branch, dad, pattern_scale);
            case 20:return computePartialLikelihoodSSE < 20 > (dad_branch, dad, pattern_scale);
        }
    } else {
        return computePartialLikelihoodNaive(dad_branch, dad, pattern_scale);
    }    
}

template<int NSTATES>
inline void PhyloTree::computePartialLikelihoodSSE(PhyloNeighbor *dad_branch, PhyloNode *dad, double *pattern_scale) {
    // don't recompute the likelihood
    if (dad_branch->partial_lh_computed & 1)
        return;
    Node *node = dad_branch->node;
    int ptn, cat;
    double *trans_state;
    double *partial_lh_site;
    double *partial_lh_child;
    double *partial_lh_block;
    bool do_scale = true;
    double freq;
    dad_branch->lh_scale_factor = 0.0;

    if (node->isLeaf() && dad) {
        // external node
        memset(dad_branch->partial_lh, 0, lh_size * sizeof (double));
        //double *partial_lh_site;
        for (ptn = 0; ptn < alnSize; ++ptn) {
            char state;
            partial_lh_site = dad_branch->partial_lh + (ptn * block);

            if (node->name == ROOT_NAME) {
                state = STATE_UNKNOWN;
            } else {                
                state = (aln->at(ptn))[node->id];
            }

            if (state == STATE_UNKNOWN) {
                for (int state2 = 0; state2 < block; state2++) {
                    partial_lh_site[state2] = 1.0;
                }
            } else if (state < NSTATES) {
                cat = 0;
                double *_par_lh_site = partial_lh_site + state;
                while (true) {
                    *_par_lh_site = 1.0;
                    ++cat;
                    if (cat == numCat)
                        break;
                    _par_lh_site += NSTATES;
                }
            } else {
                // ambiguous character, for DNA, RNA
                state = state - (NSTATES - 1);
                for (int state2 = 0; state2 < NSTATES; state2++)
                    if (state & (1 << state2)) {
                        cat = 0;
                        double *_par_lh_site = partial_lh_site + state2;
                        while (true) {
                            *_par_lh_site = 1.0;
                            ++cat;
                            if (cat == numCat)
                                break;
                            _par_lh_site += NSTATES;
                        }
                    }
            }
        }
    } else {
        // internal node
        EIGEN_ALIGN16 double trans_mat[numCat * tranSize];
        for (ptn = 0; ptn < lh_size; ++ptn)
            dad_branch->partial_lh[ptn] = 1.0;

        FOR_NEIGHBOR_IT(node, dad, it)
        if ((*it)->node->name != ROOT_NAME) {
            computePartialLikelihoodSSE<NSTATES > ((PhyloNeighbor*) (*it), (PhyloNode*) node, pattern_scale);
            dad_branch->lh_scale_factor += ((PhyloNeighbor*) (*it))->lh_scale_factor;
            for (cat = 0; cat < numCat; cat++) {
                model_factory->computeTransMatrix((*it)->length * site_rate->getRate(cat), &trans_mat[cat * tranSize]);
            }
            partial_lh_site = dad_branch->partial_lh;
            partial_lh_child = ((PhyloNeighbor*) (*it))->partial_lh;            
            for (ptn = 0; ptn < alnSize; ++ptn) {
                partial_lh_block = partial_lh_site;
                freq = ptn_freqs[ptn];
                trans_state = trans_mat;
                cat = 0;
                do_scale = true;
                while (true) {
                    ++cat;
                    MappedRowVec(NSTATES) ei_partial_lh_child(partial_lh_child);
                    MappedRowVec(NSTATES) ei_partial_lh_site(partial_lh_site);
                    MappedMat(NSTATES) ei_trans_state(trans_state);
                    ei_partial_lh_site.noalias() = (ei_partial_lh_child * ei_trans_state).cwiseProduct(ei_partial_lh_site);
                    partial_lh_site += NSTATES;
                    partial_lh_child += NSTATES;
                    if (cat == numCat)
                        break;
                    else
                        trans_state += tranSize;
                }
                for (cat = 0; cat < block; cat++)
                    if (partial_lh_block[cat] > SCALING_THRESHOLD) {
                        do_scale = false;
                        break;
                    }
                if (do_scale) {
                    Map<ArrayXd, Aligned> ei_lh_block(partial_lh_block, block);
                    ei_lh_block *= SCALING_THRESHOLD_INVER;                    
                    dad_branch->lh_scale_factor += LOG_SCALING_THRESHOLD * freq;
                    if (pattern_scale)
                        pattern_scale[ptn] += LOG_SCALING_THRESHOLD;
                }                
            }
        }
    }
     
    dad_branch->partial_lh_computed |= 1;
}

/****************************************************************************
        computing derivatives of likelihood function
 ****************************************************************************/
template<int NSTATES>
inline double PhyloTree::computeLikelihoodDervSSE(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf) {
    PhyloNode *node = (PhyloNode*) dad_branch->node;
    PhyloNeighbor *node_branch = (PhyloNeighbor*) node->findNeighbor(dad);
    //assert(node_branch);
    // swap node and dad if node is a leaf
    if (node->isLeaf()) {
        PhyloNode *tmp_node = dad;
        dad = node;
        node = tmp_node;
        PhyloNeighbor *tmp_nei = dad_branch;
        dad_branch = node_branch;
        node_branch = tmp_nei;
    }
    if ((dad_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihoodSSE<NSTATES > (dad_branch, dad);
        //computePartialLikelihoodNaive(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihoodSSE<NSTATES > (node_branch, node);
        //computePartialLikelihoodNaive(node_branch, node);

    // now combine likelihood at the branch
    double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;
    df = ddf = 0.0;

    int ptn = 0;
    int cat = 0;
    double *partial_lh_site = node_branch->partial_lh;
    double *partial_lh_child = dad_branch->partial_lh;
    double lh_ptn; // likelihood of the pattern
    double lh_ptn_derv1;
    double lh_ptn_derv2;
    double derv1_frac;
    double derv2_frac;
    double *trans_state;
    double *derv1_state;
    double *derv2_state;

    EIGEN_ALIGN16 double trans_mat[numCat * tranSize];
    EIGEN_ALIGN16 double trans_derv1[numCat * tranSize];
    EIGEN_ALIGN16 double trans_derv2[numCat * tranSize];

    Map<Array<double, 1, NSTATES>, Aligned> ei_state_freq(state_freq);
    Array<double, NSTATES, NSTATES> ei_state_freq_mat = ei_state_freq.colwise().replicate(NSTATES);

    trans_state = trans_mat;
    derv1_state = trans_derv1;
    derv2_state = trans_derv2;    
    while (true) {        
        double rate_val = site_rate->getRate(cat);
        double rate_sqr = rate_val * rate_val;
        model_factory->computeTransDerv(dad_branch->length * rate_val, trans_state, derv1_state, derv2_state);
        MappedArr2D(NSTATES) ei_trans_cat(trans_state);
        MappedArr2D(NSTATES) ei_derv1_cat(derv1_state);
        MappedArr2D(NSTATES) ei_derv2_cat(derv2_state);
        ei_trans_cat *= ei_state_freq_mat;
        ei_derv1_cat *= (ei_state_freq_mat * rate_val);
        ei_derv2_cat *= (ei_state_freq_mat * rate_sqr);
        ++cat;
        if ( cat == numCat)
            break;
        trans_state += tranSize;
        derv1_state += tranSize;
        derv2_state += tranSize;
    }
    int dad_state = STATE_UNKNOWN;
    for ( ; ptn < alnSize; ++ptn) {
        lh_ptn = 0.0;
        lh_ptn_derv1 = 0.0;
        lh_ptn_derv2 = 0.0;
        double freq = ptn_freqs[ptn];
        int padding;
        if (dad->isLeaf()) {
            dad_state = (*aln)[ptn][dad->id];
            padding = dad_state * NSTATES;
        }
        if (dad_state < NSTATES) {
            //external node
            trans_state = trans_mat + padding;
            derv1_state = trans_derv1 + padding;
            derv2_state = trans_derv2 + padding;
            cat = 0;
            while (true) {
                ++cat;
                MappedVec(NSTATES) ei_partial_lh_child(partial_lh_child);
                MappedVec(NSTATES) ei_trans_state(trans_state);
                MappedVec(NSTATES) ei_derv1_state(derv1_state);
                MappedVec(NSTATES) ei_derv2_state(derv2_state);
                lh_ptn += ei_partial_lh_child.dot(ei_trans_state);
                lh_ptn_derv1 += ei_partial_lh_child.dot(ei_derv1_state);
                lh_ptn_derv2 += ei_partial_lh_child.dot(ei_derv2_state);
                partial_lh_child += NSTATES;
                partial_lh_site += NSTATES;
                if ( cat == numCat )
                    break;
                trans_state += tranSize;
                derv1_state += tranSize;
                derv2_state += tranSize;                
            }
        } else {
            // internal node, or external node but ambiguous character
            cat = 0;
            trans_state = trans_mat;
            derv1_state = trans_derv1;
            derv2_state = trans_derv2;            
            while (true) {
                ++cat;
                MappedRowVec(NSTATES) ei_partial_lh_site(partial_lh_site);
                MappedRowVec(NSTATES) ei_partial_lh_child(partial_lh_child);
                MappedMat(NSTATES) ei_trans_state(trans_state);
                MappedMat(NSTATES) ei_derv1_state(derv1_state);
                MappedMat(NSTATES) ei_derv2_state(derv2_state);
                lh_ptn += (ei_partial_lh_child * ei_trans_state).dot(ei_partial_lh_site);
                lh_ptn_derv1 += (ei_partial_lh_child * ei_derv1_state).dot(ei_partial_lh_site);
                lh_ptn_derv2 += (ei_partial_lh_child * ei_derv2_state).dot(ei_partial_lh_site);
                partial_lh_site += NSTATES;
                partial_lh_child += NSTATES;
                if (cat == numCat)
                    break;
                trans_state += tranSize;
                derv1_state += tranSize;
                derv2_state += tranSize;
            }
        }

        lh_ptn = lh_ptn * p_var_cat + p_invar_ptns[ptn];
        //lh_ptn += p_invar_ptns[ptn];        
        double tmp = p_var_cat / lh_ptn;
        derv1_frac = lh_ptn_derv1 * tmp;
        derv2_frac = lh_ptn_derv2 * tmp;
        df += derv1_frac * freq;
        ddf += (derv2_frac - derv1_frac * derv1_frac) * freq;
        lh_ptns[ptn] = lh_ptn;
        //tree_lh += log(lh_ptn) * freq;
    }
    vrda_log(alnSize, lh_ptns.data(), lh_ptns_log.data());
    tree_lh += (lh_ptns_log * ptn_freqs).sum();
    return tree_lh;
}

double PhyloTree::computeLikelihoodDervNaive(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf) {
    PhyloNode *node = (PhyloNode*) dad_branch->node;
    PhyloNeighbor *node_branch = (PhyloNeighbor*) node->findNeighbor(dad);
    assert(node_branch);
    // swap node and dad if dad is a leaf
    if (node->isLeaf()) {
        PhyloNode *tmp_node = dad;
        dad = node;
        node = tmp_node;
        PhyloNeighbor *tmp_nei = dad_branch;
        dad_branch = node_branch;
        node_branch = tmp_nei;
        //cout << "swapped\n";
    }
    if ((dad_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihoodNaive(dad_branch, dad);
    if ((node_branch->partial_lh_computed & 1) == 0)
        computePartialLikelihoodNaive(node_branch, node);
    // now combine likelihood at the branch

    double tree_lh = node_branch->lh_scale_factor + dad_branch->lh_scale_factor;
    df = ddf = 0.0;
    int ncat = site_rate->getNRate();
    double p_invar = site_rate->getPInvar();
    double p_var_cat = (1.0 - p_invar) / (double) ncat;
    int nstates = aln->num_states;
    int block = ncat * nstates;
    int trans_size = nstates * nstates;
    int ptn, cat, state1, state2;

    double trans_mat[ncat * trans_size];
    double trans_derv1[ncat * trans_size];
    double trans_derv2[ncat * trans_size];
    double state_freq[nstates];
    model->getStateFrequency(state_freq);

    for (cat = 0; cat < ncat; cat++) {
        //trans_mat[cat] = model->newTransMatrix();
        double *trans_cat = trans_mat + (cat * trans_size);
        double *derv1_cat = trans_derv1 + (cat * trans_size);
        double *derv2_cat = trans_derv2 + (cat * trans_size);
        double rate_val = site_rate->getRate(cat);
        double rate_sqr = rate_val * rate_val;
        model_factory->computeTransDerv(dad_branch->length * rate_val, trans_cat, derv1_cat, derv2_cat);
        for (state1 = 0; state1 < nstates; state1++) {
            double *trans_mat_state = trans_cat + (state1 * nstates);
            double *trans_derv1_state = derv1_cat + (state1 * nstates);
            double *trans_derv2_state = derv2_cat + (state1 * nstates);

            for (state2 = 0; state2 < nstates; state2++) {
                trans_mat_state[state2] *= state_freq[state1];
                trans_derv1_state[state2] *= (state_freq[state1] * rate_val);
                trans_derv2_state[state2] *= (state_freq[state1] * rate_sqr);
            }
        }
    }

    for (ptn = 0; ptn < aln->size(); ptn++) {
        double lh_ptn = 0.0; // likelihood of the pattern
        double lh_ptn_derv1 = 0.0;
        double lh_ptn_derv2 = 0.0;
        int dad_state = STATE_UNKNOWN;
        if (dad->isLeaf())
            dad_state = (*aln)[ptn][dad->id];
        for (cat = 0; cat < ncat; cat++) {
            double *partial_lh_site = node_branch->partial_lh + (ptn * block + cat * nstates);
            double *partial_lh_child = dad_branch->partial_lh + (ptn * block + cat * nstates);
            if (dad_state < nstates) {
                // external node
                double *trans_state = trans_mat + (cat * trans_size + dad_state * nstates);
                double *derv1_state = trans_derv1 + (cat * trans_size + dad_state * nstates);
                double *derv2_state = trans_derv2 + (cat * trans_size + dad_state * nstates);
                for (state2 = 0; state2 < nstates; state2++) {
                    lh_ptn += partial_lh_child[state2] * trans_state[state2];
                    lh_ptn_derv1 += partial_lh_child[state2] * derv1_state[state2];
                    lh_ptn_derv2 += partial_lh_child[state2] * derv2_state[state2];
                }
            } else {
                // internal node, or external node but ambiguous character
                for (state1 = 0; state1 < nstates; state1++) {
                    double lh_state = 0.0; // likelihood of state1
                    double lh_state_derv1 = 0.0;
                    double lh_state_derv2 = 0.0;
                    double *trans_state = trans_mat + (cat * trans_size + state1 * nstates);
                    double *derv1_state = trans_derv1 + (cat * trans_size + state1 * nstates);
                    double *derv2_state = trans_derv2 + (cat * trans_size + state1 * nstates);
                    for (state2 = 0; state2 < nstates; state2++) {
                        lh_state += partial_lh_child[state2] * trans_state[state2];
                        lh_state_derv1 += partial_lh_child[state2] * derv1_state[state2];
                        lh_state_derv2 += partial_lh_child[state2] * derv2_state[state2];
                    }
                    lh_ptn += lh_state * partial_lh_site[state1];
                    lh_ptn_derv1 += lh_state_derv1 * partial_lh_site[state1];
                    lh_ptn_derv2 += lh_state_derv2 * partial_lh_site[state1];
                }
            }
        }                
        lh_ptn *= p_var_cat;
        lh_ptn_derv1 *= p_var_cat;
        lh_ptn_derv2 *= p_var_cat;
        if ((*aln)[ptn].is_const && (*aln)[ptn][0] < nstates) {
            lh_ptn += p_invar * state_freq[(*aln)[ptn][0]];
        }
        assert(lh_ptn > 0);
        double derv1_frac = lh_ptn_derv1 / lh_ptn;
        double derv2_frac = lh_ptn_derv2 / lh_ptn;
        tree_lh += log(lh_ptn) * (*aln)[ptn].frequency;
        df += derv1_frac * (*aln)[ptn].frequency;
        ddf += (derv2_frac - derv1_frac * derv1_frac) * (*aln)[ptn].frequency;
    }
    //for (cat = ncat-1; cat >= 0; cat--)
    //delete trans_mat[cat];
    //delete state_freq;
    return tree_lh;
}

double PhyloTree::computeLikelihoodDerv(PhyloNeighbor *dad_branch, PhyloNode *dad, double &df, double &ddf) {
    if (sse) {
        switch (aln->num_states) {
            case 2: return computeLikelihoodDervSSE < 2 > (dad_branch, dad, df, ddf);
            case 4: return computeLikelihoodDervSSE < 4 > (dad_branch, dad, df, ddf);
            case 20:return computeLikelihoodDervSSE < 20 > (dad_branch, dad, df, ddf);
        }
    } else {
        return computeLikelihoodDervNaive(dad_branch, dad, df, ddf);
    }
}

/****************************************************************************
        Branch length optimization by maximum likelihood
 ****************************************************************************/

double PhyloTree::computeFunction(double value) {
    current_it->length = value;
    current_it_back->length = value;
    return -computeLikelihoodBranch(current_it, (PhyloNode*) current_it_back->node);
}

double PhyloTree::computeFuncDerv(double value, double &df, double &ddf) {
    current_it->length = value;
    current_it_back->length = value;
    double lh = -computeLikelihoodDerv(current_it, (PhyloNode*) current_it_back->node, df, ddf);
    df = -df;
    ddf = -ddf;
    return lh;
}

double PhyloTree::optimizeOneBranch(PhyloNode *node1, PhyloNode *node2, bool clearLH) {
    double negative_lh;
    current_it = (PhyloNeighbor*) node1->findNeighbor(node2);
    assert(current_it);
    current_it_back = (PhyloNeighbor*) node2->findNeighbor(node1);
    assert(current_it_back);
    double current_len = current_it->length;
    double ferror, optx;
    /*if (verbose_mode == VB_DEBUG) {
            cout << "For branch " << node1->name << "," << node2->name << endl;
    }*/
    if (optimize_by_newton) // Newton-Raphson method
        optx = minimizeNewton(MIN_BRANCH_LEN, current_len, MAX_BRANCH_LEN, TOL_BRANCH_LEN, negative_lh);
    else // Brent methodcase 20:return computeLikelihoodDervSSE<20> (dad_branch, dad, df, ddf);
        optx = minimizeOneDimen(MIN_BRANCH_LEN, current_len, MAX_BRANCH_LEN, TOL_BRANCH_LEN, &negative_lh, &ferror);
    if (fabs(current_len - optx) <= TOL_BRANCH_LEN) // if nothing changes, return
        return -negative_lh;
    current_it->length = optx;
    current_it_back->length = optx;
    if (clearLH) {
        node1->clearReversePartialLh(node2);
        node2->clearReversePartialLh(node1);
    }
    return -negative_lh;
}

double PhyloTree::optimizeChildBranches(PhyloNode *node, PhyloNode *dad) {
    double tree_lh = 0.0;

    FOR_NEIGHBOR_IT(node, dad, it) {
        tree_lh = optimizeOneBranch((PhyloNode*) node, (PhyloNode*) (*it)->node);
    }
    return tree_lh;
}

double PhyloTree::optimizeAllBranches(PhyloNode *node, PhyloNode *dad) {
    //double tree_lh = optimizeChildBranches(node, dad);
    double tree_lh;

    FOR_NEIGHBOR_IT(node, dad, it) {
        //if (!(*it)->node->isLeaf())
        tree_lh = optimizeAllBranches((PhyloNode*) (*it)->node, node);
    }
    if (dad) tree_lh = optimizeOneBranch(node, dad);
    return tree_lh;
}

double PhyloTree::optimizeAllBranches(int iterations) {
    if (verbose_mode > VB_MED)
        cout << "Optimizing branch lenths (max " << iterations << " loops)..." << endl;
    double tree_lh = computeLikelihood();
    //cout << tree_lh << endl;
    for (int i = 0; i < iterations; i++) {
        double new_tree_lh = optimizeAllBranches((PhyloNode*) root);
        //clearAllPartialLh();
        //new_tree_lh = computeLikelihood();
        if (verbose_mode > VB_MAX) {
            cout << "BRANCH LEN " << i + 1 << " : ";
            cout.precision(10);
            cout << new_tree_lh << endl;
        }
        if (new_tree_lh <= tree_lh + TOL_LIKELIHOOD)
            return (new_tree_lh > tree_lh) ? new_tree_lh : tree_lh;
        tree_lh = new_tree_lh;
    }
    return tree_lh;
}

/****************************************************************************
        Stepwise addition (greedy) by maximum likelihood
 ****************************************************************************/

double PhyloTree::addTaxonML(Node *added_node, Node* &target_node, Node* &target_dad, Node *node, Node *dad) {
    Neighbor *dad_nei = dad->findNeighbor(node);

    // now insert the new node in the middle of the branch node-dad
    double len = dad_nei->length;
    node->updateNeighbor(dad, added_node, len / 2.0);
    dad->updateNeighbor(node, added_node, len / 2.0);
    added_node->updateNeighbor((Node*) 1, node, len / 2.0);
    added_node->updateNeighbor((Node*) 2, dad, len / 2.0);
    // compute the likelihood
    clearAllPartialLh();
    double best_score = optimizeChildBranches((PhyloNode*) added_node);
    target_node = node;
    target_dad = dad;
    // remove the added node
    node->updateNeighbor(added_node, dad, len);
    dad->updateNeighbor(added_node, node, len);
    added_node->updateNeighbor(node, (Node*) 1, len);
    added_node->updateNeighbor(dad, (Node*) 2, len);

    // now tranverse the tree downwards

    FOR_NEIGHBOR_IT(node, dad, it) {
        Node *target_node2;
        Node *target_dad2;
        double score = addTaxonML(added_node, target_node2, target_dad2, (*it)->node, node);
        if (score > best_score) {
            best_score = score;
            target_node = target_node2;
            target_dad = target_dad2;
        }
    }
    return best_score;
}

void PhyloTree::growTreeML(Alignment *alignment) {

    cout << "Stepwise addition using ML..." << endl;
    aln = alignment;
    int size = aln->getNSeq();
    if (size < 3)
        outError(ERR_FEW_TAXA);

    root = newNode();
    Node *new_taxon;

    // create initial tree with 3 taxa
    for (leafNum = 0; leafNum < 3; leafNum++) {
        cout << "Add " << aln->getSeqName(leafNum) << " to the tree" << endl;
        new_taxon = newNode(leafNum, aln->getSeqName(leafNum).c_str());
        root->addNeighbor(new_taxon, 1.0);
        new_taxon->addNeighbor(root, 1.0);
    }
    root = findNodeID(0);
    optimizeAllBranches();

    // stepwise adding the next taxon
    for (leafNum = 3; leafNum < size; leafNum++) {
        cout << "Add " << aln->getSeqName(leafNum) << " to the tree" << endl;
        // allocate a new taxon and a new ajedcent internal node
        new_taxon = newNode(leafNum, aln->getSeqName(leafNum).c_str());
        Node *added_node = newNode();
        added_node->addNeighbor(new_taxon, 1.0);
        new_taxon->addNeighbor(added_node, 1.0);

        // preserve two neighbors
        added_node->addNeighbor((Node*) 1, 1.0);
        added_node->addNeighbor((Node*) 2, 1.0);


        Node *target_node = NULL;
        Node *target_dad = NULL;
        addTaxonML(added_node, target_node, target_dad, root->neighbors[0]->node, root);
        // now insert the new node in the middle of the branch node-dad
        double len = target_dad->findNeighbor(target_node)->length;
        target_node->updateNeighbor(target_dad, added_node, len / 2.0);
        target_dad->updateNeighbor(target_node, added_node, len / 2.0);
        added_node->updateNeighbor((Node*) 1, target_node, len / 2.0);
        added_node->updateNeighbor((Node*) 2, target_dad, len / 2.0);
        // compute the likelihood
        clearAllPartialLh();
        optimizeAllBranches();
        optimizeNNI();
    }

    nodeNum = 2 * leafNum - 2;
}

/****************************************************************************
        Distance function
 ****************************************************************************/

double PhyloTree::computeDist(int seq1, int seq2, double initial_dist) {
    // if no model or site rate is specified, return JC distance
    if (initial_dist == 0.0)
        initial_dist = aln->computeDist(seq1, seq2);
    if (!model_factory || !site_rate) return initial_dist;

    // now optimize the distance based on the model and site rate
    AlignmentPairwise aln_pair(this, seq1, seq2);
    return aln_pair.optimizeDist(initial_dist);
}

void PhyloTree::computeDist(double *dist_mat) {
    time_t begin_time = clock();
    int nseqs = aln->getNSeq();
    int pos = 0;
    double longest_dist = 0.0;
    for (int seq1 = 0; seq1 < nseqs; seq1++)
        for (int seq2 = 0; seq2 < nseqs; seq2++, pos++) {
            if (seq1 == seq2)
                dist_mat[pos] = 0.0;
            else if (seq2 > seq1) {
                dist_mat[pos] = computeDist(seq1, seq2, dist_mat[pos]);
            } else dist_mat[pos] = dist_mat[seq2 * nseqs + seq1];
            if (dist_mat[pos] > longest_dist)
                longest_dist = dist_mat[pos];
        }
    if (longest_dist > MAX_GENETIC_DIST * 0.99)
        outWarning("Some distances are saturated. Please check your alignment again");
    cout << "Time: " << (double) (clock() - begin_time) / CLOCKS_PER_SEC << " seconds" << endl;
}

void PhyloTree::computeDist(Params &params, Alignment *alignment, double* &dist_mat, string &dist_file) {

    aln = alignment;
    dist_file = params.out_prefix;
    if (!model_factory)
        dist_file += ".jcdist";
    else
        dist_file += ".mldist";

    if (!dist_mat) {
        dist_mat = new double[alignment->getNSeq() * alignment->getNSeq()];
        memset(dist_mat, 0, sizeof (double) * alignment->getNSeq() * alignment->getNSeq());
    }
    if (!params.dist_file) {
        computeDist(dist_mat);
        alignment->printDist(dist_file.c_str(), dist_mat);
    } else {
        alignment->readDist(params.dist_file, dist_mat);
        dist_file = params.dist_file;
    }
}

/****************************************************************************
        compute BioNJ tree, a more accurate extension of Neighbor-Joining
 ****************************************************************************/

void PhyloTree::computeBioNJ(Params &params, Alignment *alignment, string &dist_file) {
    string bionj_file = params.out_prefix;
    bionj_file += ".bionj";

    cout << "Computing BIONJ tree..." << endl;
    BioNj bionj;
    bionj.create(dist_file.c_str(), bionj_file.c_str());
    bool my_rooted = false;
    bool non_empty_tree = (root != NULL);
    if (root) freeNode();
    readTree(bionj_file.c_str(), my_rooted);
    if (non_empty_tree) initializeAllPartialLh();
    setAlignment(alignment);
}

int PhyloTree::fixNegativeBranch(double fixed_length, Node *node, Node *dad) {
    if (!node) node = root;
    int fixed = 0;

    FOR_NEIGHBOR_IT(node, dad, it) {
        if ((*it)->length < 0.0) { // negative branch length detected
            if (verbose_mode == VB_DEBUG)
                cout << "Negative branch length " << (*it)->length << " was set to ";
            //(*it)->length = ((double)(rand())/RAND_MAX)*0.1+TOL_BRANCH_LEN;
            (*it)->length = fixed_length;
            if (verbose_mode == VB_DEBUG)
                cout << (*it)->length << endl;
            // set the backward branch length
            (*it)->node->findNeighbor(node)->length = (*it)->length;
            fixed++;
        }
        fixed += fixNegativeBranch(fixed_length, (*it)->node, node);
    }
    return fixed;
}

/****************************************************************************
        Nearest Neighbor Interchange by maximum likelihood
 ****************************************************************************/

double PhyloTree::swapNNIBranch(double cur_score, PhyloNode *node1, PhyloNode *node2) {
    assert(node1->degree() == 3 && node2->degree() == 3);

    PhyloNeighbor *node12_it = (PhyloNeighbor*) node1->findNeighbor(node2);
    PhyloNeighbor *node21_it = (PhyloNeighbor*) node2->findNeighbor(node1);
    double node12_len = node12_it->length;

    // save the likelihood vector at the two ends of node1-node2
    double *node1_lh_save = node12_it->partial_lh;
    double *node2_lh_save = node21_it->partial_lh;
    node12_it->partial_lh = newPartialLh();
    node21_it->partial_lh = newPartialLh();

    // TUNG save the first found neighbor (2 Neighbor total) of node 1 (excluding node2) in node1_it
    FOR_NEIGHBOR_DECLARE(node1, node2, node1_it)
    break;
    Neighbor *node1_nei = *node1_it;
    double node1_len = node1_nei->length;
    NeighborVec::iterator node1_nei_it = node1_nei->node->findNeighborIt(node1);


    // Neighbors of node2 which are not node1

    FOR_NEIGHBOR_IT(node2, node1, node2_it) {

        // do the NNI swap
        Neighbor *node2_nei = *node2_it;
        // TUNG unused variable ?
        NeighborVec::iterator node2_nei_it = node2_nei->node->findNeighborIt(node2);

        double node2_len = node2_nei->length;

        node1->updateNeighbor(node1_it, node2_nei);
        node2_nei->node->updateNeighbor(node2, node1);

        node2->updateNeighbor(node2_it, node1_nei);
        node1_nei->node->updateNeighbor(node1, node2);

        // partial_lhclear partial likelihood vector
        node12_it->clearPartialLh();
        node21_it->clearPartialLh();

        // compute the score of the swapped topology
        double score;
        score = optimizeOneBranch(node1, node2);
        // if better: return
        if (score > cur_score) {
            node2->clearReversePartialLh(node1);
            node1->clearReversePartialLh(node2);
            cur_score = score;
            cout << "Swapped neighbors :" << node1_nei->node->id << " and " << node2_nei->node->id << endl;
            break;

        }
        //cout << "Could not find better score for NNI with Node " << node1->id << "->" << "Node " << node2->id << endl;

        // else, swap back, also recover the branch lengths
        node1->updateNeighbor(node1_it, node1_nei, node1_len);
        node1_nei->node->updateNeighbor(node2, node1, node1_len);
        node2->updateNeighbor(node2_it, node2_nei, node2_len);
        node2_nei->node->updateNeighbor(node1, node2, node2_len);
        node12_it->length = node12_len;
        node21_it->length = node12_len;
    }

    // restore the partial likelihood vector
    delete [] node21_it->partial_lh;
    delete [] node12_it->partial_lh;
    node12_it->partial_lh = node1_lh_save;
    node21_it->partial_lh = node2_lh_save;
    return cur_score;
}

double PhyloTree::optimizeNNI(double cur_score, PhyloNode *node, PhyloNode *dad) {
    if (!node) node = (PhyloNode*) root;
    if (!node->isLeaf() && dad && !dad->isLeaf()) {
        double score = swapNNIBranch(cur_score, node, dad);
        if (score > cur_score) return score;
    }

    FOR_NEIGHBOR_IT(node, dad, it) {
        double score = optimizeNNI(cur_score, (PhyloNode*) (*it)->node, node);
        if (score > cur_score) return score;
    }
    return cur_score;
}

double PhyloTree::optimizeNNI() {
    double cur_score = computeLikelihood();
    for (int i = 0; i < 100; i++) {
        double score = optimizeNNI(cur_score);
        if (score <= cur_score) break;
        if (verbose_mode > VB_MED)
            cout << "NNI " << i + 1 << " : " << score << endl;
        //cur_score = score;
        cur_score = optimizeAllBranches((PhyloNode*) root);
    }
    return optimizeAllBranches();
}

double PhyloTree::optimizeNNIBranches() {
    if (verbose_mode > VB_MIN)
        cout << "Search with Nearest Neighbor Interchange (NNI) using ML..." << endl;
    double cur_score = computeLikelihood();
    for (int i = 0; i < 100; i++) {
        double score = optimizeNNI();
        if (score <= cur_score + TOL_LIKELIHOOD) break;
        cur_score = score;
    }
    return cur_score;
}

/****************************************************************************
        Subtree Pruning and Regrafting by maximum likelihood
 ****************************************************************************/

double PhyloTree::optimizeSPR(double cur_score, PhyloNode *node, PhyloNode *dad) {
    if (!node) node = (PhyloNode*) root;
    PhyloNeighbor *dad1_nei = NULL;
    PhyloNeighbor *dad2_nei = NULL;
    PhyloNode *sibling1 = NULL;
    PhyloNode *sibling2 = NULL;
    double sibling1_len = 0.0, sibling2_len = 0.0;

    if (dad && !dad->isLeaf()) {
        assert(dad->degree() == 3);
        // assign the sibling of node, with respect to dad

        FOR_NEIGHBOR_DECLARE(dad, node, it) {
            if (!sibling1) {
                dad1_nei = (PhyloNeighbor*) (*it);
                sibling1 = (PhyloNode*) (*it)->node;
                sibling1_len = (*it)->length;
            } else {
                dad2_nei = (PhyloNeighbor*) (*it);
                sibling2 = (PhyloNode*) (*it)->node;
                sibling2_len = (*it)->length;
            }
        }
        // remove the subtree leading to node
        double sum_len = sibling1_len + sibling2_len;
        sibling1->updateNeighbor(dad, sibling2, sum_len);
        sibling2->updateNeighbor(dad, sibling1, sum_len);
        PhyloNeighbor* sibling1_nei = (PhyloNeighbor*) sibling1->findNeighbor(sibling2);
        PhyloNeighbor* sibling2_nei = (PhyloNeighbor*) sibling2->findNeighbor(sibling1);
        sibling1_nei->clearPartialLh();
        sibling2_nei->clearPartialLh();

        // now try to move the subtree to somewhere else
        vector<PhyloNeighbor*> spr_path;

        FOR_NEIGHBOR(sibling1, sibling2, it) {
            spr_path.push_back(sibling1_nei);
            double score = swapSPR(cur_score, 1, node, dad, sibling1, sibling2, (PhyloNode*) (*it)->node, sibling1, spr_path);
            // if likelihood score improves, return
            if (score > cur_score) return score;
            spr_path.pop_back();
        }

        FOR_NEIGHBOR(sibling2, sibling1, it) {
            spr_path.push_back(sibling2_nei);
            double score = swapSPR(cur_score, 1, node, dad, sibling1, sibling2, (PhyloNode*) (*it)->node, sibling2, spr_path);
            // if likelihood score improves, return
            if (score > cur_score) return score;
            spr_path.pop_back();
        }
        // if likelihood does not imporve, swap back
        sibling1->updateNeighbor(sibling2, dad, sibling1_len);
        sibling2->updateNeighbor(sibling1, dad, sibling2_len);
        dad1_nei->node = sibling1;
        dad1_nei->length = sibling1_len;
        dad2_nei->node = sibling2;
        dad2_nei->length = sibling2_len;
        clearAllPartialLh();
    }

    FOR_NEIGHBOR_IT(node, dad, it) {
        double score = optimizeSPR(cur_score, (PhyloNode*) (*it)->node, node);
        if (score > cur_score) return score;
    }
    return cur_score;
}

double PhyloTree::swapSPR(double cur_score, int cur_depth, PhyloNode *node1, PhyloNode *dad1,
        PhyloNode *orig_node1, PhyloNode *orig_node2,
        PhyloNode *node2, PhyloNode *dad2, vector<PhyloNeighbor*> &spr_path) {

    PhyloNeighbor *node1_nei = (PhyloNeighbor*) node1->findNeighbor(dad1);
    PhyloNeighbor *dad1_nei = (PhyloNeighbor*) dad1->findNeighbor(node1);
    double node1_dad1_len = node1_nei->length;
    PhyloNeighbor *node2_nei = (PhyloNeighbor*) node2->findNeighbor(dad2);

    if (dad2) {
        // now, connect (node1-dad1) to the branch (node2-dad2)
        bool first = true;
        PhyloNeighbor *node2_nei = (PhyloNeighbor*) node2->findNeighbor(dad2);
        PhyloNeighbor *dad2_nei = (PhyloNeighbor*) dad2->findNeighbor(node2);
        double len2 = node2_nei->length;

        FOR_NEIGHBOR_IT(dad1, node1, it) {
            if (first) {
                (*it)->node = dad2;
                (*it)->length = len2 / 2;
                dad2->updateNeighbor(node2, dad1, len2 / 2);
                first = false;
            } else {
                (*it)->node = node2;
                (*it)->length = len2 / 2;
                node2->updateNeighbor(dad2, dad1, len2 / 2);
            }
            ((PhyloNeighbor*) (*it))->clearPartialLh();
        }
        node2_nei->clearPartialLh();
        dad2_nei->clearPartialLh();
        node1_nei->clearPartialLh();
        vector<PhyloNeighbor*>::iterator it2;
        for (it2 = spr_path.begin(); it2 != spr_path.end(); it2++)
            (*it2)->clearPartialLh();
        clearAllPartialLh();
        // optimize relevant branches
        double score;

        /* testing different branch optimization */
        score = optimizeOneBranch(node1, dad1);
        //score = optimizeOneBranch(dad2, dad1);
        //score = optimizeOneBranch(node2, dad1);
        /*
        PhyloNode *cur_node = dad2;
        for (int i = spr_path.size()-1; i >= 0; i--) {
                score = optimizeOneBranch(cur_node, (PhyloNode*)spr_path[i]->node);
                cur_node = (PhyloNode*)spr_path[i]->node;
        }
         */
        //score = optimizeAllBranches(dad1);

        // if score improves, return
        if (score > cur_score) return score;
        // else, swap back
        node2->updateNeighbor(dad1, dad2, len2);
        dad2->updateNeighbor(dad1, node2, len2);
        node2_nei->clearPartialLh();
        dad2_nei->clearPartialLh();
        node1_nei->length = node1_dad1_len;
        dad1_nei->length = node1_dad1_len;

        // add to candiate SPR moves
        spr_moves.add(node1, dad1, node2, dad2, score);
    }
    if (cur_depth >= spr_radius) return cur_score;
    spr_path.push_back(node2_nei);

    FOR_NEIGHBOR_IT(node2, dad2, it) {
        double score = swapSPR(cur_score, cur_depth + 1, node1, dad1, orig_node1, orig_node2, (PhyloNode*) (*it)->node, node2, spr_path);
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

    PhyloNeighbor *node1_nei = (PhyloNeighbor*) node->findNeighbor(dad);
    PhyloNeighbor *dad1_nei = (PhyloNeighbor*) dad->findNeighbor(node);
    double node1_dad1_len = node1_nei->length;

    // assign the sibling of node, with respect to dad

    FOR_NEIGHBOR_DECLARE(dad, node, it) {
        if (!sibling1) {
            dad_nei1 = (PhyloNeighbor*) (*it);
            sibling1 = (PhyloNode*) (*it)->node;
            sibling1_len = (*it)->length;
        } else {
            dad_nei2 = (PhyloNeighbor*) (*it);
            sibling2 = (PhyloNode*) (*it)->node;
            sibling2_len = (*it)->length;
        }
    }
    // remove the subtree leading to node
    double sum_len = sibling1_len + sibling2_len;
    sibling1->updateNeighbor(dad, sibling2, sum_len);
    sibling2->updateNeighbor(dad, sibling1, sum_len);
    // now try to move the subtree to somewhere else

    bool first = true;
    PhyloNeighbor *node2_nei = (PhyloNeighbor*) node2->findNeighbor(dad2);
    PhyloNeighbor *dad2_nei = (PhyloNeighbor*) dad2->findNeighbor(node2);
    double len2 = node2_nei->length;

    FOR_NEIGHBOR(dad, node, it) {
        if (first) {
            (*it)->node = dad2;
            (*it)->length = len2 / 2;
            dad2->updateNeighbor(node2, dad, len2 / 2);
            first = false;
        } else {
            (*it)->node = node2;
            (*it)->length = len2 / 2;
            node2->updateNeighbor(dad2, dad, len2 / 2);
        }
        ((PhyloNeighbor*) (*it))->clearPartialLh();
    }

    clearAllPartialLh();
    // optimize branches
    double score;
    score = optimizeAllBranches(dad);

    // if score improves, return
    if (score > cur_score) return score;
    // else, swap back
    node2->updateNeighbor(dad, dad2, len2);
    dad2->updateNeighbor(dad, node2, len2);

    node1_nei->length = node1_dad1_len;
    dad1_nei->length = node1_dad1_len;

    sibling1->updateNeighbor(sibling2, dad, sibling1_len);
    sibling2->updateNeighbor(sibling1, dad, sibling2_len);
    dad_nei1->node = sibling1;
    dad_nei1->length = sibling1_len;
    dad_nei2->node = sibling2;
    dad_nei2->length = sibling2_len;
    clearAllPartialLh();

    return cur_score;

}

double PhyloTree::optimizeSPR() {
    double cur_score = computeLikelihood();
    //spr_radius = leafNum / 5;
    spr_radius = 2;
    for (int i = 0; i < 100; i++) {
        spr_moves.clear();
        double score = optimizeSPR(cur_score, (PhyloNode*) root->neighbors[0]->node);
        clearAllPartialLh();
        if (score <= cur_score) {
            for (SPRMoves::iterator it = spr_moves.begin(); it != spr_moves.end(); it++) {
                //cout << (*it).score << endl;
                score = assessSPRMove(cur_score, *it);
                // if likelihood score improves, apply to SPR
                if (score > cur_score) break;
            }
            if (score <= cur_score) break;
        }
        cout.precision(10);
        cout << "SPR " << i + 1 << " : " << score << endl;
        //cur_score = score;
        cur_score = optimizeAllBranches();
    }
    return optimizeAllBranches();
}

double PhyloTree::optimizeSPRBranches() {
    cout << "Search with Subtree Pruning and Regrafting (SPR) using ML..." << endl;
    double cur_score = computeLikelihood();
    for (int i = 0; i < 100; i++) {
        double score = optimizeSPR();
        if (score <= cur_score + TOL_LIKELIHOOD) break;
        cur_score = score;
    }
    return cur_score;
}

void PhyloTree::pruneSubtree(PhyloNode *node, PhyloNode *dad, PruningInfo &info) {
    bool first = true;
    info.node = node;
    info.dad = dad;

    FOR_NEIGHBOR_IT(dad, node, it) {
        if (first) {
            info.dad_it_left = it;
            info.dad_nei_left = (*it);
            info.dad_lh_left = ((PhyloNeighbor*) (*it))->partial_lh;
            info.left_node = (*it)->node;
            info.left_len = (*it)->length;
            first = false;
        } else {
            info.dad_it_right = it;
            info.dad_nei_right = (*it);
            info.dad_lh_right = ((PhyloNeighbor*) (*it))->partial_lh;
            info.right_node = (*it)->node;
            info.right_len = (*it)->length;
        }
    }
    info.left_it = info.left_node->findNeighborIt(dad);
    info.right_it = info.right_node->findNeighborIt(dad);
    info.left_nei = (*info.left_it);
    info.right_nei = (*info.right_it);

    info.left_node->updateNeighbor(info.left_it, info.dad_nei_right);
    info.right_node->updateNeighbor(info.right_it, info.dad_nei_left);
    ((PhyloNeighbor*) info.dad_nei_right)->partial_lh = newPartialLh();
    ((PhyloNeighbor*) info.dad_nei_left)->partial_lh = newPartialLh();
}

void PhyloTree::regraftSubtree(PruningInfo &info,
        PhyloNode *in_node, PhyloNode *in_dad) {
    NeighborVec::iterator in_node_it = in_node->findNeighborIt(in_dad);
    NeighborVec::iterator in_dad_it = in_dad->findNeighborIt(in_node);
    Neighbor *in_dad_nei = (*in_dad_it);
    Neighbor *in_node_nei = (*in_node_it);
    double in_len = in_dad_nei->length;
    info.dad->updateNeighbor(info.dad_it_right, in_dad_nei);
    info.dad->updateNeighbor(info.dad_it_left, in_node_nei);
    // SOMETHING NEED TO BE DONE
    //in_dad->updateNeighbor(in_dad_it,

}

/****************************************************************************
        Approximate Likelihood Ratio Test with SH-like interpretation
 ****************************************************************************/


void PhyloTree::computeNNIPatternLh(
        double cur_lh,
        double &lh2, double *pattern_lh2,
        double &lh3, double *pattern_lh3,
        PhyloNode *node1, PhyloNode *node2) {
    assert(node1->degree() == 3 && node2->degree() == 3);

    // recompute pattern scaling factors if necessary
    PhyloNeighbor *node12_it = (PhyloNeighbor*) node1->findNeighbor(node2);
    PhyloNeighbor *node21_it = (PhyloNeighbor*) node2->findNeighbor(node1);
    double *ptn_scale = NULL;
    double new_scale[aln->getNPattern()];
    double sum_scaling = node12_it->lh_scale_factor + node21_it->lh_scale_factor;
    NeighborVec::iterator it;
    ptn_scale = new double[aln->getNPattern()];
    memset(ptn_scale, 0, sizeof (double) * aln->getNPattern());

    FOR_NEIGHBOR(node1, node2, it)
    if (((PhyloNeighbor*) * it)->lh_scale_factor < 0.0) {
        ((PhyloNeighbor*) * it)->clearForwardPartialLh(node1);
        computePartialLikelihood((PhyloNeighbor*) (*it), node1, ptn_scale);
    }
    FOR_NEIGHBOR(node2, node1, it)
    if (((PhyloNeighbor*) * it)->lh_scale_factor < 0.0) {
        ((PhyloNeighbor*) * it)->clearForwardPartialLh(node2);
        computePartialLikelihood((PhyloNeighbor*) (*it), node2, ptn_scale);
    }

    const int IT_NUM = 6;

    NeighborVec::iterator saved_it[IT_NUM];
    int id = 0;

    FOR_NEIGHBOR(node1, node2, it) {
        saved_it[id++] = (*it)->node->findNeighborIt(node1);
    } else {
        saved_it[id++] = it;
    }

    FOR_NEIGHBOR(node2, node1, it) {
        saved_it[id++] = (*it)->node->findNeighborIt(node2);
    } else {
        saved_it[id++] = it;
    }
    assert(id == IT_NUM);


    Neighbor * saved_nei[IT_NUM];
    // save Neighbor and allocate new Neighbor pointer
    for (id = 0; id < IT_NUM; id++) {
        saved_nei[id] = (*saved_it[id]);
        *saved_it[id] = new PhyloNeighbor(saved_nei[id]->node, saved_nei[id]->length);
        ((PhyloNeighbor*) (*saved_it[id]))->partial_lh = newPartialLh();
    }

    // get the Neighbor again since it is replaced for saving purpose
    node12_it = (PhyloNeighbor*) node1->findNeighbor(node2);
    node21_it = (PhyloNeighbor*) node2->findNeighbor(node1);

    PhyloNeighbor *node2_lastnei;

    // save the first found neighbor of node 1 (excluding node2) in node1_it
    FOR_NEIGHBOR_DECLARE(node1, node2, node1_it) break;
    Neighbor *node1_nei = *node1_it;

    bool first = true;

    FOR_NEIGHBOR_IT(node2, node1, node2_it) {
        /* do the NNI swap */
        Neighbor *node2_nei = *node2_it;
        node1->updateNeighbor(node1_it, node2_nei);
        node2_nei->node->updateNeighbor(node2, node1);

        node2->updateNeighbor(node2_it, node1_nei);
        node1_nei->node->updateNeighbor(node1, node2);

        // re-optimize five adjacent branches
        double old_score = -INFINITY, new_score = old_score;

        // clear partial likelihood vector
        node12_it->clearPartialLh();
        node21_it->clearPartialLh();
        int i;
        for (i = 0; i < 2; i++) {

            new_score = optimizeOneBranch(node1, node2, false);

            FOR_NEIGHBOR(node1, node2, it) {
                //for (id = 0; id < IT_NUM; id++)
                //((PhyloNeighbor*)(*saved_it[id]))->clearPartialLh();
                ((PhyloNeighbor*) (*it)->node->findNeighbor(node1))->clearPartialLh();
                new_score = optimizeOneBranch(node1, (PhyloNode*) (*it)->node, false);
            }

            node21_it->clearPartialLh();

            FOR_NEIGHBOR(node2, node1, it) {
                //for (id = 0; id < IT_NUM; id++)
                //((PhyloNeighbor*)(*saved_it[id]))->clearPartialLh();
                ((PhyloNeighbor*) (*it)->node->findNeighbor(node2))->clearPartialLh();
                new_score = optimizeOneBranch(node2, (PhyloNode*) (*it)->node, false);
                node2_lastnei = (PhyloNeighbor*) (*it);
            }
            node12_it->clearPartialLh();
            if (new_score < old_score + TOL_LIKELIHOOD) break;
            old_score = new_score;
        }
        //new_score = optimizeOneBranch(node1, node2, false);
        if (new_score > cur_lh + TOL_LIKELIHOOD)
            cout << "Alternative NNI shows better likelihood " << new_score << " > " << cur_lh << endl;
        double *result_lh;
        if (first) {
            result_lh = pattern_lh2;
            lh2 = new_score;
        } else {
            result_lh = pattern_lh3;
            lh3 = new_score;
        }

        memcpy(new_scale, ptn_scale, sizeof (double) * aln->getNPattern());

        node12_it->clearPartialLh();
        node21_it->clearPartialLh();

        computePartialLikelihood(node12_it, node1, new_scale);
        computePartialLikelihood(node21_it, node2, new_scale);

        old_score = new_score;
        // compute the score of the swapped topology
        //new_score = computeLikelihoodBranch(node2_lastnei, node2, result_lh);
        //double x[aln->getNPattern()];
        new_score = computeLikelihoodBranch(node12_it, node1, result_lh);

        if (ptn_scale) {
            double check_score = 0.0;
            for (i = 0; i < aln->getNPattern(); i++) {
                result_lh[i] += new_scale[i];
                check_score += result_lh[i] * aln->at(i).frequency;
            }
            // make sure that the pattern likelihoods were just computed correctly
            assert(fabs(new_score - check_score) < 1e-4);
        }
        /*
                        if (fabs(new_score - old_score) > TOL_LIKELIHOOD)
                                cout << "something wrong" << endl;
                        for (i = 0; i < aln->getNPattern(); i++)
                                if (fabs(result_lh[i] - x[i]) > TOL_LIKELIHOOD)
                                        cout << "wrong " << i << endl;*/
        // swap back and recover the branch lengths
        node1->updateNeighbor(node1_it, node1_nei);
        node1_nei->node->updateNeighbor(node2, node1);
        node2->updateNeighbor(node2_it, node2_nei);
        node2_nei->node->updateNeighbor(node1, node2);
        first = false;
    }

    // restore the Neighbor*
    for (id = 0; id < IT_NUM; id++) {
        delete [] ((PhyloNeighbor*) * saved_it[id])->partial_lh;
        delete (*saved_it[id]);
        (*saved_it[id]) = saved_nei[id];
    }

    // restore the length of 4 branches around node1, node2
    FOR_NEIGHBOR(node1, node2, it)
            (*it)->length = (*it)->node->findNeighbor(node1)->length;
    FOR_NEIGHBOR(node2, node1, it)
            (*it)->length = (*it)->node->findNeighbor(node2)->length;
    if (ptn_scale) delete [] ptn_scale;
}

void PhyloTree::resampleLh(double **pat_lh, double *lh_new) {
    int nsite = aln->getNSite();
    int nptn = aln->getNPattern();
    memset(lh_new, 0, sizeof (double) * 3);
    int freq[nptn], i;
    memset(freq, 0, nptn * sizeof (int));
    for (i = 0; i < nsite; i++) {
        int site_id = floor(((double) (rand()) / RAND_MAX) * nsite);
        int ptn_id = aln->getPatternID(site_id);
        freq[ptn_id]++;
    }
    for (i = 0; i < nptn; i++) {
        lh_new[0] += freq[i] * pat_lh[0][i];
        lh_new[1] += freq[i] * pat_lh[1][i];
        lh_new[2] += freq[i] * pat_lh[2][i];
    }
}

// Implementation of testBranch follows Guindon et al. (2010)

double PhyloTree::testOneBranch(
        double best_score, double *pattern_lh,
        int times, PhyloNode *node1, PhyloNode *node2) {
    int NUM_NNI = 3;
    double lh[NUM_NNI];
    double *pat_lh[NUM_NNI];
    lh[0] = best_score;
    pat_lh[0] = pattern_lh;
    pat_lh[1] = new double[aln->getNPattern()];
    pat_lh[2] = new double[aln->getNPattern()];
    computeNNIPatternLh(best_score, lh[1], pat_lh[1], lh[2], pat_lh[2], node1, node2);
    double aLRT;
    if (lh[1] > lh[2])
        aLRT = (lh[0] - lh[1]);
    else
        aLRT = (lh[0] - lh[2]);

    int support = 0;

    for (int i = 0; i < times; i++) {
        double lh_new[NUM_NNI];
        // resampling estimated log-likelihood (RELL)
        resampleLh(pat_lh, lh_new);
        double cs[NUM_NNI], cs_best, cs_2nd_best;
        cs[0] = lh_new[0] - lh[0];
        cs[1] = lh_new[1] - lh[1];
        cs[2] = lh_new[2] - lh[2];
        if (cs[0] >= cs[1] && cs[0] >= cs[2]) {
            cs_best = cs[0];
            if (cs[1] > cs[2]) cs_2nd_best = cs[1];
            else cs_2nd_best = cs[2];
        } else if (cs[1] >= cs[2]) {
            cs_best = cs[1];
            if (cs[0] > cs[2]) cs_2nd_best = cs[0];
            else cs_2nd_best = cs[2];
        } else {
            cs_best = cs[2];
            if (cs[0] > cs[1]) cs_2nd_best = cs[0];
            else cs_2nd_best = cs[1];
        }
        if (aLRT > (cs_best - cs_2nd_best) + 0.05) support++;
    }
    delete [] pat_lh[2];
    delete [] pat_lh[1];
    return ((double) support) / times;
}

int PhyloTree::testAllBranches(int threshold, double best_score, double *pattern_lh,
        int times, PhyloNode *node, PhyloNode *dad) {
    int num_low_support = 0;
    if (!node) {
        node = (PhyloNode*) root;
        root->neighbors[0]->node->name = "";
    }
    if (dad && !node->isLeaf() && !dad->isLeaf()) {
        int support = round(testOneBranch(best_score, pattern_lh, times, node, dad)*100);
        node->name = convertIntToString(support);
        if (support < threshold) num_low_support = 1;
        ((PhyloNeighbor*) node->findNeighbor(dad))->partial_pars[0] = support;
        ((PhyloNeighbor*) dad->findNeighbor(node))->partial_pars[0] = support;
    }
    FOR_NEIGHBOR_IT(node, dad, it)
    num_low_support += testAllBranches(threshold, best_score, pattern_lh, times, (PhyloNode*) (*it)->node, node);
    return num_low_support;
}

/****************************************************************************
        Collapse stable (highly supported) clades by one representative
 ****************************************************************************/

void PhyloTree::deleteLeaf(Node *leaf) {
    Node *near_node = leaf->neighbors[0]->node;
    assert(leaf->isLeaf() && near_node->degree() == 3);
    Node *node1 = NULL;
    Node *node2 = NULL;
    double sum_len = 0.0;

    FOR_NEIGHBOR_IT(near_node, leaf, it) {
        sum_len += (*it)->length;
        if (!node1)
            node1 = (*it)->node;
        else
            node2 = (*it)->node;
    }
    // make sure that the returned node1 and node2 are correct
    assert(node1 && node2);
    // update the neighbor
    node1->updateNeighbor(near_node, node2, sum_len);
    node2->updateNeighbor(near_node, node1, sum_len);
}

void PhyloTree::reinsertLeaf(Node *leaf, Node *node,
        Node *dad) {
    bool first = true;
    Node *adjacent_node = leaf->neighbors[0]->node;
    Neighbor *nei = node->findNeighbor(dad);
    double len = nei->length;

    FOR_NEIGHBOR_IT(adjacent_node, leaf, it) {
        if (first) {
            (*it)->node = node;
            (*it)->length = len / 2;
            node->updateNeighbor(dad, adjacent_node, len / 2);
        } else {
            (*it)->node = dad;
            (*it)->length = len / 2;
            dad->updateNeighbor(node, adjacent_node, len / 2);
        }
        first = false;
    }
}

bool PhyloTree::isSupportedNode(PhyloNode* node, int min_support) {
    FOR_NEIGHBOR_IT(node, NULL, it)
    if (!(*it)->node->isLeaf())
        if (((PhyloNeighbor*) * it)->partial_pars[0] < min_support) {
            return false;
        }
    return true;
}

int PhyloTree::collapseStableClade(int min_support, NodeVector &pruned_taxa, StrVector &linked_name, double* &dist_mat) {
    NodeVector taxa;
    NodeVector::iterator tax_it;
    StrVector::iterator linked_it;
    getTaxa(taxa);
    IntVector linked_taxid;
    linked_taxid.resize(leafNum, -1);
    int num_pruned_taxa; // global num of pruned taxa
    int ntaxa = leafNum;
    do {
        num_pruned_taxa = 0;
        for (tax_it = taxa.begin(); tax_it != taxa.end(); tax_it++)
            if (linked_taxid[(*tax_it)->id] < 0) {
                Node *taxon = (*tax_it);
                PhyloNode *near_node = (PhyloNode*) taxon->neighbors[0]->node;
                Node *adj_taxon = NULL;
                FOR_NEIGHBOR_DECLARE(near_node, taxon, it)
                if ((*it)->node->isLeaf()) {
                    adj_taxon = (*it)->node;
                    break;
                }
                // if it is not a cherry
                if (!adj_taxon) continue;
                assert(linked_taxid[adj_taxon->id] < 0);
                PhyloNeighbor *near_nei = NULL;
                FOR_NEIGHBOR(near_node, taxon, it)
                if ((*it)->node != adj_taxon) {
                    near_nei = (PhyloNeighbor*) (*it);
                    break;
                }
                assert(near_nei);
                // continue if the cherry is not stable, or distance between two taxa is near ZERO
                if (!isSupportedNode((PhyloNode*) near_nei->node, min_support) && dist_mat[taxon->id * ntaxa + adj_taxon->id] > 2e-6) continue;
                // now do the taxon pruning
                Node *pruned_taxon = taxon, *stayed_taxon = adj_taxon;
                // prune the taxon that is far away
                if (adj_taxon->neighbors[0]->length > taxon->neighbors[0]->length) {
                    pruned_taxon = adj_taxon;
                    stayed_taxon = taxon;
                }
                deleteLeaf(pruned_taxon);
                linked_taxid[pruned_taxon->id] = stayed_taxon->id;
                pruned_taxa.push_back(pruned_taxon);
                linked_name.push_back(stayed_taxon->name);
                num_pruned_taxa++;
                // do not prune more than n-4 taxa
                if (pruned_taxa.size() >= ntaxa - 4) break;
            }
    } while (num_pruned_taxa && pruned_taxa.size() < ntaxa - 4);

    if (pruned_taxa.empty()) return 0;


    if (verbose_mode >= VB_MED)
        for (tax_it = pruned_taxa.begin(), linked_it = linked_name.begin(); tax_it != pruned_taxa.end(); tax_it++, linked_it++)
            cout << "Delete " << (*tax_it)->name << " from " << (*linked_it) << endl;

    // set root to the first taxon which was not deleted
    for (tax_it = taxa.begin(); tax_it != taxa.end(); tax_it++)
        if (linked_taxid[(*tax_it)->id] < 0) {
            root = (*tax_it);
            break;
        }
    // extract the sub alignment
    IntVector stayed_id;
    int i, j;
    for (i = 0; i < taxa.size(); i++)
        if (linked_taxid[i] < 0) stayed_id.push_back(i);
    assert(stayed_id.size() + pruned_taxa.size() == leafNum);
    Alignment *pruned_aln = new Alignment();
    pruned_aln->extractSubAlignment(aln, stayed_id, 2); // at least 2 informative characters
    nodeNum = leafNum = stayed_id.size();
    initializeTree();
    setAlignment(pruned_aln);

    double *pruned_dist = new double [leafNum * leafNum];
    for (i = 0; i < leafNum; i++)
        for (j = 0; j < leafNum; j++)
            pruned_dist[i * leafNum + j] = dist_mat[stayed_id[i] * ntaxa + stayed_id[j]];
    dist_mat = pruned_dist;
    return pruned_taxa.size();
}

int PhyloTree::restoreStableClade(Alignment *original_aln, NodeVector &pruned_taxa, StrVector &linked_name) {
    int num_inserted_taxa;
    NodeVector::reverse_iterator tax_it;
    StrVector::reverse_iterator linked_it;
    tax_it = pruned_taxa.rbegin();
    linked_it = linked_name.rbegin();
    for (; tax_it != pruned_taxa.rend(); tax_it++, linked_it++) {
        //cout << "Reinsert " << (*tax_it)->name << " to " << (*linked_it) << endl;
        Node *linked_taxon = findNodeName((*linked_it));
        assert(linked_taxon);
        assert(linked_taxon->isLeaf());
        leafNum++;
        reinsertLeaf((*tax_it), linked_taxon, linked_taxon->neighbors[0]->node);
    }
    assert(leafNum == original_aln->getNSeq());
    nodeNum = leafNum;
    initializeTree();
    setAlignment(original_aln);
    root = findNodeName(aln->getSeqName(0));
    //if (verbose_mode >= VB_MED) drawTree(cout);
}

bool PhyloTree::checkEqualScalingFactor(double &sum_scaling, PhyloNode *node, PhyloNode *dad) {
    if (!node) node = (PhyloNode*) root;
    if (dad) {
        double scaling = ((PhyloNeighbor*) node->findNeighbor(dad))->lh_scale_factor +
                ((PhyloNeighbor*) dad->findNeighbor(node))->lh_scale_factor;
        if (sum_scaling > 0) sum_scaling = scaling;
        if (fabs(sum_scaling - scaling) > 1e-6) {
            cout << sum_scaling << " " << scaling << endl;
            return false;
        }
    }
    FOR_NEIGHBOR_IT(node, dad, it)
    if (!checkEqualScalingFactor(sum_scaling, (PhyloNode*) (*it)->node, node)) return false;
    return true;
}
