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

#ifndef MEMSLOT_H
#define MEMSLOT_H

#ifndef PHYLOTREE_H
#error "Please #include phylotree.h before including this header file" 
#endif

/**
    one memory slot, used for memory saving technique
*/
struct MemSlot {
    int status; // status of this slot
    PhyloNeighbor *nei; // neighbor assigned to this slot
    double *partial_lh; // partial_lh assigned to this slot
    UBYTE *scale_num; // scale_num assigned to this slot

    PhyloNeighbor *saved_nei;
};

/**
    all memory slots, used for memory saving technique
*/
class MemSlotVector : public vector<MemSlot> {
public:

    /** initialize with a specified number of slots */
    void init(PhyloTree *tree, int num_slot);

    /** 
        lock the memory assigned to nei
        @param nei neighbor to lock
        @return TRUE if successfully locked, FALSE otherwise 
    */
    bool lock(PhyloNeighbor *nei);

    /** unlock the memory assigned to nei */
    void unlock(PhyloNeighbor *nei);

    /** test if the memory assigned to nei is locked or not */
    bool locked(PhyloNeighbor *nei);

    /** allocate free or unlocked memory to nei */
    int allocate(PhyloNeighbor *nei);

    /** update neighbor */
    void update(PhyloNeighbor *nei);

    /** find ID the a neighbor */
    iterator findNei(PhyloNeighbor *nei);

    /** add neighbor into a specified iterator */
    void addNei(PhyloNeighbor *nei, iterator it);

    /** reset everything */
    void reset();

    /** clean up all neighbors where partial_lh_computed = 0 */
    void cleanup();

    /** take over neighbor from another one */
    void takeover(PhyloNeighbor *nei, PhyloNeighbor *taken_nei);

    /** add special neihbor e.g. for NNI */
    void addSpecialNei(PhyloNeighbor *nei);

    /** erase special neihbor e.g. for NNI */
    void eraseSpecialNei();

    /** replace a neighbor, used for NNI */
    void replace(PhyloNeighbor *new_nei, PhyloNeighbor *old_nei);

    /** restore neighbor, after calling replace */
    void restore(PhyloNeighbor *new_nei, PhyloNeighbor *old_nei);

protected:


    /** 
        map from neighbor to slot ID for fast lookup
        IMPORTANT: mapping to ID instead of (unsafe) iterator
    */
    unordered_map<PhyloNeighbor*, int> nei_id_map;

    /** counter of free slot ID */
    int free_count;

};


#endif // MEMSLOT_H
