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


int TinaTree::computeParsimonyScore() {
    assert(root && root->isLeaf());

    int score = 0;
    for (int ptn = 0; ptn < aln->size(); ptn++)
        if (!aln->at(ptn).is_const) {
            int states;
            score += computeParsimonyScore(ptn, states) * (*aln)[ptn].frequency;
        }
    return score;
}

