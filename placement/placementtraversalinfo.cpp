//
// placementtraversalinfo.cpp
// Implementation of the PlacementTraversalInfo class.
//
// Created by James Barbetti on 08-Oct-2020.
//

#include "placementtraversalinfo.h"

PlacementTraversalInfo::PlacementTraversalInfo(
        PhyloTree& tree,LikelihoodBufferSet& buffersToUse,
        PhyloNeighbor* dad_branch, PhyloNode* dad)
        : super(dad_branch,dad), phylo_tree(tree), buffers(buffersToUse) {
    size_t nstates         = phylo_tree.aln->num_states;
    size_t ncat            = phylo_tree.site_rate->getNRate();
    size_t ncat_mix        = (phylo_tree.model_factory->fused_mix_rate) ? ncat : ncat*phylo_tree.model->getNMixtures();
    size_t block           = nstates * ncat_mix;
    size_t children_size   = get_safe_upper_limit(block*nstates*2);
    size_t leaf_size_floor = (phylo_tree.aln->STATE_UNKNOWN+1)*block*2;
    size_t lh_leaf_size    = get_safe_upper_limit(leaf_size_floor);
    echildren              = aligned_alloc<double>(children_size);
    partial_lh_leaves      = aligned_alloc<double>(lh_leaf_size);
    buffer_tmp             = aligned_alloc<double>(phylo_tree.aln->num_states);
}

void PlacementTraversalInfo::computePartialLikelihood(PhyloNeighbor* nei, PhyloNode* node) {
    if (nei != nullptr) {
        dad_branch = nei;
    }
    if (node != nullptr) {
        dad = node;
    }
    size_t orig_nptn = roundUpToMultiple(phylo_tree.aln->size(), 8);
    size_t nptn      = roundUpToMultiple(orig_nptn+phylo_tree.model_factory->unobserved_ptns.size(), 8);

    phylo_tree.computePartialInfoDouble(*this, buffer_tmp);
    
    //
    //Note: Using 8 here is dodgy.  It would be better if we knew the
    //      size (in doubles) of the vector class the tree is using.
    //      Because we'd prefer to use that.
    //
    std::vector<size_t> limits;
    phylo_tree.computePatternPacketBounds(8, phylo_tree.num_threads,
                                          phylo_tree.num_packets, nptn, limits);
    #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,1) num_threads(phylo_tree.num_threads)
    #endif
    for (int packet_id = 0; packet_id < phylo_tree.num_packets; ++packet_id) {
        phylo_tree.computePartialLikelihood(*this, limits[packet_id],
                                            limits[packet_id+1], packet_id,
                                            buffers);
    }
}

double PlacementTraversalInfo::getBranchScore() {
    return phylo_tree.computeLikelihoodBranch(dad_branch, dad,
                                              buffers);
}

PlacementTraversalInfo::~PlacementTraversalInfo() {
    aligned_free(buffer_tmp);
    aligned_free(echildren);
    aligned_free(partial_lh_leaves);
}


