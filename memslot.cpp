/***************************************************************************
 *   Copyright (C) 2009-2016 by                                            *
 *   BUI Quang Minh <minh.bui@univie.ac.at>                                *
 *                                                                         *
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

#include "phylotree.h"
#include "memslot.h"

const int MEM_LOCKED = 1;

void MemSlotVector::init(PhyloTree *tree, int num_slot) {
    if (Params::getInstance().lh_mem_save != LM_MEM_SAVE)
        return;
    resize(num_slot);
    size_t lh_size = tree->getPartialLhSize();
    size_t scale_size = tree->getScaleNumSize();
    reset();
    for (iterator it = begin(); it != end(); it++) {
        it->partial_lh = tree->central_partial_lh + lh_size*(it-begin());
        it->scale_num = tree->central_scale_num + scale_size*(it-begin());
    }
}

void MemSlotVector::reset() {
    if (Params::getInstance().lh_mem_save != LM_MEM_SAVE)
        return;
    for (iterator it = begin(); it != end(); it++) {
        it->status = 0;
        it->nei = NULL;
    }
    nei_id_map.clear();
    free_count = 0;
}


MemSlotVector::iterator MemSlotVector::findNei(PhyloNeighbor *nei) {
    auto it = nei_id_map.find(nei);
    assert(it != nei_id_map.end());
//    assert(at(it->second).nei == nei);
    return it->second;
}

void MemSlotVector::assignNei(iterator it, PhyloNeighbor* nei) {
    it->nei = nei;
}

void MemSlotVector::addNei(PhyloNeighbor *nei, iterator it) {
    nei->partial_lh = it->partial_lh;
    nei->scale_num = it->scale_num;
    assignNei(it, nei);
    nei_id_map[nei] = it;
}

bool MemSlotVector::lock(PhyloNeighbor *nei) {
    if (Params::getInstance().lh_mem_save != LM_MEM_SAVE)
        return false;
    if (nei->node->isLeaf())
        return false;
    iterator id = findNei(nei);
    assert((id->status & MEM_LOCKED) == 0);
    id->status |= MEM_LOCKED;
    return true;
}

void MemSlotVector::unlock(PhyloNeighbor *nei) {
    if (Params::getInstance().lh_mem_save != LM_MEM_SAVE)
        return;
    if (nei->node->isLeaf())
        return;
    iterator id = findNei(nei);
    assert((id->status & MEM_LOCKED) != 0);
    id->status &= ~MEM_LOCKED;
}

bool MemSlotVector::locked(PhyloNeighbor *nei) {
    if (Params::getInstance().lh_mem_save != LM_MEM_SAVE)
        return false;
    if (nei->node->isLeaf())
        return false;
    iterator id = findNei(nei);

    if ((id->status & MEM_LOCKED) == 0)
        return false;
    else
        return true;
}

int MemSlotVector::allocate(PhyloNeighbor *nei) {
    if (Params::getInstance().lh_mem_save != LM_MEM_SAVE)
        return -1;

    // first find a free slot
    if (free_count < size()) {
        iterator it = begin() + free_count;
        assert(it->nei == NULL);
        addNei(nei, it);
        free_count++;
        return it-begin();
    }

    int min_size = INT_MAX;
    iterator best = end();


    // no free slot found, find an unlocked slot with minimal size
    for (iterator it = begin(); it != end(); it++)
        if ((it->status & MEM_LOCKED) == 0 && min_size > it->nei->size) {
            best = it;
            min_size = it->nei->size;
            // 2 is the minimum size
            if (min_size == 2)
                break;
        }

    assert(best != end() && "No free/unlocked mem slot found!");

    // clear mem assigned to it->nei
    best->nei->clearPartialLh();

    // assign mem to nei
    addNei(nei, best);
    return best-begin();

}

void MemSlotVector::update(PhyloNeighbor *nei) {
    if (Params::getInstance().lh_mem_save != LM_MEM_SAVE)
        return;

    iterator it = findNei(nei);
    if (it->nei != nei) {
        // clear mem assigned to it->nei
        it->nei->clearPartialLh();

        // assign mem to nei
        addNei(nei, it);
    }
}


void MemSlotVector::cleanup() {
    if (Params::getInstance().lh_mem_save != LM_MEM_SAVE)
        return;
    unordered_map<PhyloNeighbor*, iterator> new_map;
    for (auto it = nei_id_map.begin(); it != nei_id_map.end(); it++)
        if (it->first != it->second->nei) {
            it->first->partial_lh_computed &= ~1; // clear bit
            it->first->partial_lh = NULL;
            it->first->scale_num = NULL;
        } else {
            new_map[it->first] = it->second;
        }
    nei_id_map = new_map;
    assert(nei_id_map.size() == size());
}

void MemSlotVector::takeover(PhyloNeighbor *nei, PhyloNeighbor *taken_nei) {
    assert(taken_nei->partial_lh);
    nei->partial_lh = taken_nei->partial_lh;
    nei->scale_num = taken_nei->scale_num;
    taken_nei->partial_lh = NULL;
    taken_nei->scale_num = NULL;
    taken_nei->partial_lh_computed &= ~1; // clear bit
    if (Params::getInstance().lh_mem_save != LM_MEM_SAVE)
        return;
    iterator id = findNei(taken_nei);
    nei_id_map.erase(nei_id_map.find(taken_nei));
    nei_id_map[nei] = id;
    if (id->nei == taken_nei) {
        assignNei(id, nei);
    }
}
