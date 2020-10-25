//
// blockallocator.cpp
// Does allocation (and tracking of) parsimony, likelihood, and
// scale vectors.
//
// Created by James Barbetti on 08-Oct-2020.
//

#include "blockallocator.h"

BlockAllocator::BlockAllocator(PhyloTree& tree,
               int parsimonyIndex)
    : phylo_tree(tree), index_parsimony(parsimonyIndex) {
        tree.determineBlockSizes();
}
BlockAllocator::~BlockAllocator() {
    //C++17 or later, merely: for ( block : uint_blocks ) aligned_free(block);
    for (auto it = uint_blocks.rbegin(); it!=uint_blocks.rend(); ++it) {
        aligned_free(*it);
    }
}
void BlockAllocator::allocateParsimonyBlock(UINT*& partial_pars) {
    partial_pars = phylo_tree.central_partial_pars
                   + (index_parsimony * phylo_tree.pars_block_size);
    ASSERT( partial_pars < phylo_tree.tip_partial_pars );
    ++index_parsimony;
}
void BlockAllocator::allocateMemoryFor(PhyloNeighbor* nei) {
    if (nei->partial_pars==nullptr) {
        allocateParsimonyBlock(nei->partial_pars);
    }
}
PhyloTree& BlockAllocator::getTree() {
    return phylo_tree;
}
int BlockAllocator::getParsimonyBlockCount() const {
    return index_parsimony;
}
void BlockAllocator::handOverComputedState(PhyloNeighbor* from_nei, PhyloNeighbor* to_nei) {
    std::swap(to_nei->partial_lh         , from_nei->partial_lh);
    std::swap(to_nei->partial_pars       , from_nei->partial_pars);
    std::swap(to_nei->scale_num          , from_nei->scale_num);
    std::swap(to_nei->partial_lh_computed, from_nei->partial_lh_computed);
    allocateMemoryFor(from_nei);
    to_nei->partial_lh_computed = from_nei->partial_lh_computed;
    from_nei->clearComputedFlags();
}
int BlockAllocator::getLikelihoodBlockCount() const {
    return 0;
}
void BlockAllocator::allocateLikelihoodBlocks(double*& partial_lh, UBYTE*& scale_num) {
    //Don't!
}
bool BlockAllocator::usesLikelihood() {
    return false;
}
void BlockAllocator::makeTreeReady(PhyloNode* first, PhyloNode* second) {
}

LikelihoodBlockPair::LikelihoodBlockPair()
    : partial_lh(nullptr), scale_num(nullptr), blocks_are_owned(false) {
}

LikelihoodBlockPair::LikelihoodBlockPair(const LikelihoodBlockPair& rhs)
    : blocks_are_owned(false) {
    //If the other is borrowing, this one can borrow too
    if (!rhs.blocks_are_owned) {
        partial_lh = rhs.partial_lh;
        scale_num  = rhs.scale_num;
    } else {
        partial_lh = nullptr;
        scale_num  = nullptr;
    }
}

LikelihoodBlockPair& LikelihoodBlockPair::operator = (const LikelihoodBlockPair& rhs) {
    if (&rhs==this) {
        return *this;   //self-referential assignment isn't a clear!
    }
    clear();
    //If the other is borrowing, this one can borrow too
    if (!rhs.blocks_are_owned) {
        partial_lh = rhs.partial_lh;
        scale_num  = rhs.scale_num;
    } else {
        partial_lh = nullptr;
        scale_num  = nullptr;
    }
    return *this;
}

LikelihoodBlockPair::LikelihoodBlockPair(double* lh, UBYTE* scale, bool takeOwnership)
    : partial_lh(lh), scale_num(scale), blocks_are_owned(takeOwnership) {}

void LikelihoodBlockPair::clear() {
    if (blocks_are_owned) {
        aligned_free(partial_lh);
        aligned_free(scale_num);
    }
    else {
        partial_lh = nullptr;
        scale_num  = nullptr;
    }
    blocks_are_owned = false;
}

LikelihoodBlockPair::~LikelihoodBlockPair() {
    clear();
}

void    LikelihoodBlockPair::copyFrom(PhyloTree& tree, PhyloNeighbor* source) {
    allocate(tree);
    if (source->partial_lh != nullptr) {
        memcpy(partial_lh, source->partial_lh, tree.lh_block_size * sizeof(partial_lh[0]));
    } else {
        memset(partial_lh, 0, tree.lh_block_size * sizeof(partial_lh[0]) );
    }
    if (source->scale_num != nullptr) {
        memcpy(scale_num, source->scale_num, tree.scale_block_size * sizeof(scale_num[0]));
    } else {
        memset(scale_num, 0, tree.scale_block_size * sizeof(scale_num[0]) );
    }
}

void LikelihoodBlockPair::allocate(PhyloTree& tree)
{
    //
    //Assumed: If partial_lh and scale_num are already allocated, and
    //         owned by this instance, they're the right size for tree,
    //         (as per tree.pars_block_size and tree.scale_block_size).
    //
    if (!blocks_are_owned) {
        clear();
    }
    ensure_aligned_allocated(partial_lh, tree.lh_block_size);
    ensure_aligned_allocated(scale_num,  tree.scale_block_size);
    blocks_are_owned = true;
}

void    LikelihoodBlockPair::lendTo(PhyloNeighbor* borrower) {
    borrower->partial_lh = this->partial_lh;
    borrower->scale_num  = this->scale_num;
}


LikelihoodBlockAllocator::LikelihoodBlockAllocator(PhyloTree& tree,
                                                   int parsimonyIndex,
                                                   int likelihoodIndex)
    : super(tree, parsimonyIndex), index_lh(likelihoodIndex) {}
LikelihoodBlockAllocator::~LikelihoodBlockAllocator() {
    //C++17 or later, merely: for ( block : double_blocks ) aligned_free(block);
    for (auto it = double_blocks.rbegin(); it!=double_blocks.rend(); ++it) {
        aligned_free(*it);
    }
    for (auto it = ubyte_blocks.rbegin(); it!=ubyte_blocks.rend(); ++it) {
        aligned_free(*it);
    }
}

int LikelihoodBlockAllocator::getLikelihoodBlockCount() const {
    return index_lh;
}

void LikelihoodBlockAllocator::allocateLikelihoodBlocks(double*& partial_lh,
                                                        UBYTE*& scale_num) {
    partial_lh = phylo_tree.central_partial_lh
               + (index_lh * phylo_tree.lh_block_size);
    scale_num  = phylo_tree.central_scale_num
               + (index_lh * phylo_tree.scale_block_size);
    ++index_lh;
}
void LikelihoodBlockAllocator::allocateMemoryFor(PhyloNeighbor* nei) {
    super::allocateMemoryFor(nei);
    if (nei->partial_lh==nullptr) {
        double* ph = nullptr;
        allocateLikelihoodBlocks(ph, nei->scale_num);
        nei->partial_lh = ph;
    }
}
void LikelihoodBlockAllocator::handOverComputedState(PhyloNeighbor* from_nei,
                                                     PhyloNeighbor* to_nei) {
    std::swap(to_nei->partial_lh         , from_nei->partial_lh);
    std::swap(to_nei->partial_pars       , from_nei->partial_pars);
    std::swap(to_nei->scale_num          , from_nei->scale_num);
    std::swap(to_nei->partial_lh_computed, from_nei->partial_lh_computed);
    allocateMemoryFor(from_nei);
    to_nei->partial_lh_computed = from_nei->partial_lh_computed;
    from_nei->clearComputedFlags();
}
bool LikelihoodBlockAllocator::usesLikelihood() {
    return true;
}
void LikelihoodBlockAllocator::findMissingInputs(PhyloNode* first, PhyloNode* second,
                                                 PhyloNeighborVec& list) {
    PhyloNeighbor* wanted = first->findNeighbor(second);
    if (second->isLeaf()) {
        return;
    }
    FOR_EACH_ADJACENT_PHYLO_NODE(second, first, it, third) {
        findMissingInputs(second, third, list);
    }
    if (wanted->partial_lh == nullptr) {
        list.emplace_back(wanted);
        //TREE_LOG_LINE(phylo_tree, VB_MIN, "Input " << pointer_to_hex(wanted)
        //              << " has partial_lh " << pointer_to_hex(wanted->partial_lh));
    }
}
void LikelihoodBlockAllocator::findUselessOutputs(PhyloNode* first, PhyloNode* second,
                                                  PhyloNeighborVec& list) {
    FOR_EACH_ADJACENT_PHYLO_NODE(second, first, it, third) {
        findUselessOutputs(second, third, list);
        PhyloNeighbor* unwanted = third->findNeighbor(second);
        if (unwanted->partial_lh != nullptr) {
            list.emplace_back(unwanted);
            //TREE_LOG_LINE(phylo_tree, VB_MIN, "Non-Input " << pointer_to_hex(unwanted)
            //              << " has partial_lh " << pointer_to_hex(unwanted->partial_lh));
        }
    }
}
void LikelihoodBlockAllocator::makeTreeReady(PhyloNode* first,
                                             PhyloNode* second) {
    if (phylo_tree.params->lh_mem_save != LM_PER_NODE) {
        return;
    }
    PhyloNeighborVec missing_inputs;
    findMissingInputs(first, second, missing_inputs);
    findMissingInputs(second, first, missing_inputs);
    size_t countMissing = missing_inputs.size();

    PhyloNeighborVec useless_outputs;
    findUselessOutputs(first, second, useless_outputs);
    findUselessOutputs(second, first, useless_outputs);
    size_t countUseless  = useless_outputs.size();
    
    //TREE_LOG_LINE(phylo_tree, VB_MIN, "missing " << countMissing
    //              << ", useless " << countUseless);
    ASSERT( countMissing <= countUseless + spare_block_pairs.size());

    for (size_t i=0; i<countMissing; ++i) {
        PhyloNeighbor* recipient = missing_inputs[i];
        if ( i < countUseless ) {
            //recycle blocks from newly useless PhyloNeighbor instances
            PhyloNeighbor* donor = useless_outputs[i];
            std::swap(recipient->partial_lh, donor->partial_lh);
            std::swap(recipient->scale_num,  donor->scale_num);
            donor->setLikelihoodComputed(false);
        } else {
            //recycle older unused blocks
            auto block_pair = spare_block_pairs.back();
            spare_block_pairs.pop_back();
            block_pair.lendTo(recipient);
        }
        recipient->setLikelihoodComputed(false);
    }
    for (size_t i=countMissing; i<countUseless; ++i) {
        //remove unused blocks from newly useless
        //PhyloNeighbor instances (that weren't recycled above)
        //so they will be available for recycling, if when needed
        PhyloNeighbor* donor     = useless_outputs[i];
        spare_block_pairs.emplace_back(donor->partial_lh, donor->scale_num);
        donor->partial_lh = nullptr;
        donor->scale_num  = nullptr;
        donor->setLikelihoodComputed(false);
    }
}

