/***************************************************************************
 *   Copyright (C) 2009 by BUI Quang Minh   *
 *   minh.bui@univie.ac.at   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "tinatree.h"

TinaTree::TinaTree()
 : PhyloTree()
{
}
TinaTree::TinaTree(Alignment *alignment) : PhyloTree(alignment) {
}


TinaTree::~TinaTree()
{
}


int TinaTree::computeParsimonyScore(int ptn, int &states, PhyloNode *node, PhyloNode *dad) {
    int score = 0;
    states = 0;
    if (!node) {
        node = getRoot();
    }
    if (node->degree() > 3)
        outError("Does not work with multifurcating tree");
    if (verbose_mode == VB_DEBUG)
        cout << ptn << " " << node->id << "  " << node->name << endl;

    if (node->isLeaf()) {
        char state;
        if (node->name == ROOT_NAME) {
            state = aln->STATE_UNKNOWN;
        } else {
            ASSERT(node->id < aln->getNSeq());
            state = (*aln)[ptn][node->id];
        }
        if (state == aln->STATE_UNKNOWN) {
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

        FOR_EACH_ADJACENT_PHYLO_NODE(node, dad, it, child) {
            int states_child;
            int score_child = computeParsimonyScore(ptn, states_child, child, node);
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


int TinaTree::computeParsimonyScore() {
    ASSERT(root && root->isLeaf());

    int score = 0;
    for (size_t ptn = 0; ptn < aln->size(); ++ptn)
        if (!aln->at(ptn).isConst()) {
            int states;
            int ptn_score = computeParsimonyScore(ptn, states);
            score += ptn_score * (*aln)[ptn].frequency;
            if (verbose_mode >= VB_MAX) {
            	for (size_t seq=0; seq < aln->getNSeq(); ++seq)
            		cout << aln->convertStateBackStr(aln->at(ptn)[seq]);
            	cout << " " << ptn_score << endl;
            }
        }
    if (verbose_mode >= VB_MAX)
    	cout << endl;
    return score;
}

void TinaTree::initializeAllPartialLh() {
    int index, indexlh;
    initializeAllPartialLh(index, indexlh);
    ASSERT(index == (nodeNum - 1)*2);
}

void TinaTree::initializeAllPartialLh(int &index, int &indexlh,
                                      PhyloNode *node, PhyloNode *dad) {
    int pars_block_size = getBitsBlockSize();
    if (!node) {
        node = getRoot();
        // allocate the big central partial likelihoods memory

        if (!central_partial_pars) {
            if (verbose_mode >= VB_MAX)
                cout << "Allocating " << (leafNum - 1)*4 * pars_block_size * sizeof (UINT) << " bytes for partial parsimony vectors" << endl;
            central_partial_pars = new UINT[(leafNum-1)*4*pars_block_size];
            if (!central_partial_pars)
                outError("Not enough memory for partial parsimony vectors");
        }
        index = 0;
    }
    if (dad) {
        // assign a region in central_partial_lh to both Neihgbors (dad->node, and node->dad)
        PhyloNeighbor *nei = node->findNeighbor(dad);
        //assert(!nei->partial_lh);
        nei->partial_pars = central_partial_pars + (index * pars_block_size);
        nei = dad->findNeighbor(node);
        //assert(!nei->partial_lh);
        nei->partial_pars = central_partial_pars + ((index + 1) * pars_block_size);
        index += 2;
        ASSERT(index < nodeNum * 2 - 1);
    }
    FOR_EACH_ADJACENT_PHYLO_NODE(node, dad, it, child) {
        initializeAllPartialLh(index, indexlh, child, node);
    }
}
