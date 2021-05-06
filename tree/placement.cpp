//
//  placement.cpp
//  Evolutionary placement: adding taxa to a phylogenetic tree
//  This file created by James Barbetti on 25/8/20, but:
//  1. addNewTaxaToTree was formerly in phylotree.cpp;
//  2. addTaxonML likewise (and was the work of BUI Quang Minh)
//
#include <tree/phylotree.h>
#include <placement/placement.h> //for PlacementParameters
#include <placement/blockallocator.h>
#include <placement/placement.h>
#include <placement/searchheuristic.h>
#include <placement/targetbranch.h>
#include <placement/taxontoplace.h>
#include <placement/placementcostcalculator.h>
#include <placement/placementoptimizer.h>
#include <placement/placementrun.h>
#include <utils/timekeeper.h>
#include <utils/timeutil.h>    //for getRealTime
#include "parsimonyspr.h"      //for ParsimonySPRMove
#include "parsimonymove.h"             //for ParsimonyPathVector
#include "phylotreethreadingcontext.h" //for PhyloTreeThreadingContext

/****************************************************************************
 Stepwise addition (greedy) by maximum likelihood
 ****************************************************************************/

double PhyloTree::recomputeParsimonyBranchLength(PhyloNode* fromNode, PhyloNode* toNode) {
    PhyloNeighbor* nei     = fromNode->findNeighbor(toNode);
    PhyloNeighbor* backnei = toNode->findNeighbor(fromNode);
    int       branch_subst = 0;
    computeParsimonyBranchFast(nei, fromNode, &branch_subst);
    double uncorrected_length = (branch_subst > 0)
                    ? ((double) branch_subst / getAlnNSite())
                    : (1.0 / getAlnNSite());
    double alpha    = (site_rate) ? site_rate->getGammaShape() : 1.0;
    nei->length     = correctBranchLengthF81(uncorrected_length, alpha);
    backnei->length = nei->length;
    return nei->length;
}

double PhyloTree::addTaxonML(PhyloNode* added_taxon,     PhyloNode *added_node,
                             PhyloNode* node,            PhyloNode* dad,
                             bool isAddedAtMidpoint,
                             PhyloNode* &target_node,    PhyloNode* &target_dad,
                             double& len_to_new_taxon,
                             double& len_to_target_node, double& len_to_target_dad) {

    Neighbor *dad_nei = dad->findNeighbor(node);

    //link the new interior node into the middle of the branch node-dad:
    //
    //   dad <---*---> added_node <---*---> node
    //                      ^
    //                      |
    //                      V
    //                 added_taxon
    //
    double len     = dad_nei->length;
    double halfLen = 0.5 * len;
    node->updateNeighbor(dad, added_node, halfLen);
    dad->updateNeighbor(node, added_node, halfLen);
    added_node->updateNeighbor(DUMMY_NODE_1, node, halfLen);
    added_node->updateNeighbor(DUMMY_NODE_2, dad,  halfLen);
    added_node->updateNeighbor(added_taxon, added_taxon, -1);
    added_taxon->updateNeighbor(added_node, added_node, -1);
    
    LOG_LINE(VerboseMode::VB_DEBUG, "  Placement branch length " << len);

    FOR_EACH_PHYLO_NEIGHBOR(added_node, nullptr, it, nei) {
        nei->clearComputedFlags();
        nei->getNode()->findNeighbor(added_node)->clearComputedFlags();
    }
    
    //compute the likelihood
    PhyloNeighbor* nei;
    double best_score = 0;
    if (isAddedAtMidpoint) {
        len_to_new_taxon   = recomputeParsimonyBranchLength(added_taxon, added_node);
        LOG_LINE(VerboseMode::VB_DEBUG, 
                 "  Parsimony taxon->interior length " << len_to_new_taxon );
        nei                = added_taxon->findNeighbor(added_node);
        best_score         = computeLikelihoodBranch(nei, added_taxon,
                                                     tree_buffers);
        LOG_LINE(VerboseMode::VB_DEBUG, 
                 "  Traversal info size is " << traversal_info.size());
        LOG_LINE(VerboseMode::VB_DEBUG, 
                 "  Likelihood before optimization " << best_score);
        optimizeOneBranch(added_taxon, added_node, false, 20);
        len_to_target_dad  = halfLen;
        len_to_target_node = halfLen;
        len_to_new_taxon   = nei->length;
        best_score         = computeLikelihoodFromBuffer();
        LOG_LINE(VerboseMode::VB_DEBUG, "  Likelihood after optimization " << best_score
                           << " (len = " << len_to_new_taxon << ")");
    }
    else {
        len_to_new_taxon = recomputeParsimonyBranchLength(added_taxon, added_node);
        optimizeOneBranch(added_node,  dad, false, 20);
        nei                = added_node->findNeighbor(dad);
        len_to_target_dad  = nei->length;

        optimizeOneBranch(added_node,  node, false, 20);
        nei                = added_node->findNeighbor(node);
        len_to_target_node = nei->length;
        
        optimizeOneBranch(added_taxon,  added_node, false, 20);
        nei                = added_node->findNeighbor(added_taxon);
        best_score         = computeLikelihoodFromBuffer();
        len_to_new_taxon   = nei->length;
    }
    target_node        = node;
    target_dad         = dad;
    LOG_LINE(VerboseMode::VB_DEBUG, "  ML Lengths " << len_to_target_dad
             << ", " << len_to_target_node << ", " << len_to_new_taxon << std::endl);

    //unlink the added node
    node->updateNeighbor(added_node, dad, len);
    dad->updateNeighbor(added_node, node, len);
    added_node->updateNeighbor(node, DUMMY_NODE_1, halfLen);
    added_node->updateNeighbor(dad,  DUMMY_NODE_2, halfLen);
    node->findNeighbor(dad)->clearComputedFlags();
    dad->findNeighbor(node)->clearComputedFlags();
    trackProgress(1.0);

    //now traverse the tree downwards
    FOR_EACH_ADJACENT_PHYLO_NODE(node, dad, it, child) {
        PhyloNode* target_node2 = nullptr;
        PhyloNode* target_dad2  = nullptr;
        double     len_child    = 0;
        double     len_node     = 0;
        double     len_dad      = 0;
        double score = addTaxonML(added_taxon, added_node, child, node,
                                  isAddedAtMidpoint,
                                  target_node2, target_dad2,
                                  len_child, len_node, len_dad);
        if (score > best_score) {
            best_score         = score;
            target_node        = target_node2;
            target_dad         = target_dad2;
            len_to_new_taxon   = len_child;
            len_to_target_node = len_node;
            len_to_target_dad  = len_dad;
        }
    }
    return best_score;
}

void PhyloTree::removeSampleTaxaIfRequested() {
    int nseq = static_cast<int>(aln->getNSeq());
    size_t countOfTaxaToRemove = PlacementParameters().getNumberOfTaxaToRemoveAndReinsert(nseq);
    if (0<countOfTaxaToRemove) {
        LOG_LINE(VerboseMode::VB_DEBUG, 
                 "Will mark " << countOfTaxaToRemove << " (out of " << nseq
                 << ") sequences, to be removed.");
        map<string, Node*> mapNameToNode;
        getMapOfTaxonNameToNode(nullptr, nullptr, mapNameToNode);
        size_t r = 0;
        LOG_LINE(VerboseMode::VB_MAX, "Before removing sequences:");
        for (auto it=mapNameToNode.begin(); it!=mapNameToNode.end(); ++it) {
            LOG_LINE(VerboseMode::VB_MAX, "sequence " << it->first 
                     << " has id " << it->second->id );
        }

        for (int seq = 0; seq < nseq; ++seq) {
            r += countOfTaxaToRemove;
            if ( nseq <= r ) {
                r -= nseq;
                string seq_name = aln->getSeqName(seq);
                auto it = mapNameToNode.find(seq_name);
                if (it!=mapNameToNode.end()) {
                    Node* node = it->second;
                    auto newName = node->name + "_Removed";
                    if (mapNameToNode.find(newName) == mapNameToNode.end()) {
                        node->name = newName;
                        mapNameToNode[newName] = node;
                    }
                    LOG_LINE(VerboseMode::VB_MED, 
                             "Marking sequence " << seq << " (" << node->name << ")"
                             " as to be removed.");
                }
            }
        }
    }
}

bool PhyloTree::shouldPlacementUseSankoffParsimony() const {
    return PlacementParameters().doesPlacementUseSankoffParsimony();
}

bool PhyloTree::shouldPlacementUseLikelihood() const {
    return PlacementParameters().doesPlacementUseLikelihood();
}

void PhyloTree::reinsertTaxaViaStepwiseParsimony(const IntVector& taxaIdsToAdd) {
    constraintTree.readConstraint(*this);
    //clearing all the nodes...
    freeNode();
    root = nullptr;
    cout << "Creating fast initial parsimony tree by random order stepwise addition..." << endl;
    double  start  = getRealTime();
    Params& params = Params::getInstance();
    const char* doing_what;
    double  score  = computeParsimonyTree(aln, randstream, params.out_prefix.c_str(), doing_what);
    cout << getRealTime() - start << " seconds, parsimony score: " << score
        << " (based on " << aln->num_parsimony_sites << " sites)"<< endl;

    //Note that this score tends to disagree.
    double parsimonyStart = getRealTime();
    clearAllPartialParsimony(false);
    double parsimonyScore = computeParsimony("Recalculating parsimony score");
    LOG_LINE(VerboseMode::VB_MED, 
             "Recalculated parsimony score " << parsimonyScore
             << " (recalculation cost " << (getRealTime() - parsimonyStart) << " sec)" );
}

typedef TaxonToPlace TaxonTypeInUse;
//typedef LessFussyTaxon TaxonTypeInUse;

void PhyloTree::optimizePlacementRegion(ParsimonySearchParameters& s,
                                        TargetBranchRange& targets,
                                        size_t region_target_index,
                                        ParsimonyPathVector& per_thread_path_parsimony,
                                        PhyloTreeThreadingContext& context,
                                        LikelihoodBlockPairs &blocks) {
    //TimeKeeper copyingIn("Copying in");
    //copyingIn.start();
    s.initializing.start();
    //Clone all the nodes that correspond to the placement region,
    //and build a map from the node (ids) back to those nodes
    //(that's map_to_real_node).
    //Clone all the target branches too (but with new branch numbers,
    //between the cloned nodes and cloned branches)
    //and construct a mapping to the IDs of targets in
    //TargetBranchRange (that's fake_to_real_branch_id).
    
    //1. Set up local_targets (subset copy of targets)
    //   This takes time proportional to B (the number of
    //   branches that will be in local_targets).
    //
    std::vector<size_t> fake_to_real_branch_id; 
    TargetBranch&  tb = targets[region_target_index];
    targets.getFinalReplacementBranchIndexes(region_target_index, fake_to_real_branch_id);
    TargetBranchRange local_targets(targets, fake_to_real_branch_id);

    //2. Find out which nodes are referenced by local_targets.
    //   (and bring the boundary nodes to the front! of the array!)
    //   This ought to be *very* quick indeed.
    //
    NodeVector real_nodes;
    local_targets.getNodes(real_nodes);
    std::partition(real_nodes.begin(), real_nodes.end(),
                   [tb](Node* n) { return n==tb.first || n==tb.second; });
#if (0)
    LOG_LINE(VerboseMode::VB_MIN, "targ bran " << region_target_index
             << " and first->id " << tb.first->id
             << " and second->id " << tb.second->id
             << " node count " << real_nodes.size());
#endif
    
    //3. Set up mappings (id #s to nodes), set up
    //   local nodes (they point to entries in a sequential
    //   array, fake_nodes).  This does 2N inserts into maps.
    std::map<int, PhyloNode*> map_to_real_node;
    std::map<int, PhyloNode*> map_to_fake_node;
    std::vector<PhyloNode> fake_nodes;
    intptr_t node_count = real_nodes.size();
    fake_nodes.resize(node_count);
    for (intptr_t i=0; i<node_count; ++i) {
        int node_id               = real_nodes[i]->id;
        map_to_real_node[node_id] = (PhyloNode*)real_nodes[i];
        map_to_fake_node[node_id] = &fake_nodes[i];
        fake_nodes[i].id          = node_id;
    }
    
    //4. Copy subtree structure (neighbor relationships)
    //   from the nodes in the target region in the
    //   original tree (EXCEPT those referring to nodes
    //   that are *outside* the region of interest).
    //   (the "EXCEPT" leaves any local copy, of a boundary
    //   node that was an interior node (in the tree),
    //   into an exterior node (in the sub-region copy).
    //
    //   This does N node map lookups (ouch!) but they can
    //   be done in parallel, because each each "i" only adds
    //   neighbors to the ith fake node (that match the ith
    //   real one).  And read-only map lookups may execute
    //   in parallel.
    //
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (intptr_t i=0; i<node_count; ++i) {
        FOR_EACH_PHYLO_NEIGHBOR(real_nodes[i], nullptr, it, real_nei) {
            PhyloNode* real_adjacent = real_nei->getNode();
            auto fake_it = map_to_fake_node.find(real_adjacent->id);
            if (fake_it != map_to_fake_node.end()) {
                fake_nodes[i].addNeighbor(fake_it->second, real_nei->length);
                PhyloNeighbor* fake_nei = fake_nodes[i].lastNeighbor();
                fake_nei->copyComputedState(real_nei);
            } else {
                //This is a neighbor, referring to a node that is outside
                //of the region of interest, that we don't want to copy.
            }
        }
    }
    
    //5. Remap node references from local_targets to fake nodes
    //   And, remap all branch numbers, between fake nodes,
    //   so that they are now mapped back to branch ids in
    //   local_targets.
    //
    //   This does 2B node map lookups (ouch!),
    //   but it can do them in parallel, because the
    //   *writes* are to the ids of the neighbors
    //   between the bth branch's nodes (no overlaps!).
    //
    intptr_t branch_count = local_targets.size();
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (intptr_t b=0; b<branch_count; ++b) {
        TargetBranch& branch = local_targets.at(b);
        //BEGIN: workaround for "outward" parsimony not computed
        bool left_nei_computed = false;
        if (branch.second==tb.first || branch.second==tb.second) {
            PhyloNeighbor* checkNei = branch.first->findNeighbor(branch.second);
            if (!checkNei->isParsimonyComputed()) {
                computeParsimonyBranch(checkNei, branch.first);
                left_nei_computed = true;
            }
        }
        bool right_nei_computed = false;
        if (branch.first==tb.first || branch.first==tb.second ) {
            PhyloNeighbor* checkNei = branch.second->findNeighbor(branch.first);
            if (!checkNei->isParsimonyComputed()) {
                computeParsimonyBranch(checkNei, branch.second);
                right_nei_computed = true;
            }
        }
        //FINISH: workaround
        branch.first  = map_to_fake_node[branch.first->id];
        branch.second = map_to_fake_node[branch.second->id];
        PhyloNeighbor* left_nei  = branch.first->findNeighbor(branch.second);
        PhyloNeighbor* right_nei = branch.second->findNeighbor(branch.first);
        left_nei->id = right_nei->id = static_cast<int>(b);
        //BEGIN: workaround for "outward" parsimony not computed
        if (left_nei_computed) {
            left_nei->setParsimonyComputed(true);
        }
        if (right_nei_computed) {
            right_nei->setParsimonyComputed(true);
        }
        //FINISH: workaround
    }
    //copyingIn.stop();
    s.initializing.stop();

    for (int step=0; step<2; ++step) {
        s.rescoring.start();
        PhyloNode*     firstNode  = local_targets[0].first;
        PhyloNode*     secondNode = local_targets[0].second;
        PhyloNeighbor* firstNeigh = firstNode->findNeighbor(secondNode);
        computeParsimony("Rescoring parsimony for subtree", true, false, firstNeigh, firstNode);
        s.rescoring.stop();
    
        if (0<step) continue;
        TimeKeeper optimizing("optimizing");
        optimizing.start();
        optimizeSubtreeParsimony<ParsimonySPRMove>(s, local_targets,
                                                   per_thread_path_parsimony,
                                                   context, false);
        optimizing.stop();
    }

    //TimeKeeper copyingOut("copying subtree back");
    //copyingOut.start();
    //Copy subtree state back!  Most of this can be done in parallel (!)
    //1. Copy for all (N-2) nodes that aren't boundary nodes.
    //This requires (N-2)+2B look-ups of real nodes (via map_to_real_node)
    //((N-2) for the nodes, 2B for the neighbors at each end of each branch),
    //but the look-ups are done in parallel.
    //
    intptr_t fake_node_count = map_to_fake_node.size();
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (intptr_t fake_node_index = 2; fake_node_index < fake_node_count;
         ++fake_node_index) {
        //Copy state for nodes in the *interior* of the subtree
        PhyloNode* fake = &fake_nodes[fake_node_index];
        PhyloNode* real = map_to_real_node[fake->id];
        int   nei_count = static_cast<int>(real->neighbors.size());
        ASSERT(nei_count == fake->neighbors.size());
        for (int nei_no = 0; nei_no < nei_count; ++nei_no) {
            PhyloNeighbor* fake_nei = fake->getNeighborByIndex(nei_no);
            PhyloNeighbor* real_nei = real->getNeighborByIndex(nei_no);
            real_nei->copyComputedState(fake_nei);
            real_nei->node = map_to_real_node[fake_nei->node->id];
            real_nei->id   = static_cast<int>(fake_to_real_branch_id[fake_nei->id]);
        }
    }
    //2. Copy for the 2 boundary nodes.
    //   This requires (at most) 8 node lookups (it requires only 6
    //   if one of the boundary nodes is, in the real tree, a leaf node)
    for (intptr_t fake_node_index =0; fake_node_index < 2; ++fake_node_index) {
        //Copy state for subtree *boundary* nodes
        PhyloNode*     fake          = &fake_nodes[fake_node_index];
        ASSERT(fake->neighbors.size()==1);
        PhyloNode*     real          = map_to_real_node[fake->id];
        PhyloNeighbor* fake_nei      = fake->firstNeighbor();
        PhyloNode*     fake_adjacent = fake_nei->getNode();
        for (Neighbor* real_nei : real->neighbors) {
            PhyloNeighbor* real_phylo_nei = (PhyloNeighbor*)real_nei;
            PhyloNode*     real_adjacent  = real_phylo_nei->getNode();
            if ( real_adjacent->id == fake_adjacent->id ) {
                //Each boundary node might have 3 neighbors
                //but only one will correspond to a node in the
                //region of interest.  And this one must be it!
                real_phylo_nei->copyComputedState(fake_nei);
            }
        }
    }
    for (intptr_t fake_branch_index = 0;
         fake_branch_index < branch_count; ++fake_branch_index) {
        TargetBranch& fake_branch = local_targets[fake_branch_index];
        auto real_branch_id = fake_to_real_branch_id[fake_branch_index];
        TargetBranch& real_branch = targets[real_branch_id];
        real_branch.first  = map_to_real_node[fake_branch.first->id];
        real_branch.second = map_to_real_node[fake_branch.second->id];
        real_branch.copyComputedState(fake_branch);
    }
    //copyingOut.stop();
    /*
    LOG_LINE(VerboseMode::VB_MIN, "Opt " << branch_count << " branches "
             << "In "  << copyingIn.elapsed_wallclock_time << " " << copyingIn.elapsed_cpu_time
             << "Opt " << optimizing.elapsed_wallclock_time << " " << optimizing.elapsed_cpu_time
             << "Out " << copyingOut.elapsed_wallclock_time << " " << copyingOut.elapsed_cpu_time
             );
    */
}

int  PhyloTree::renumberInternalNodes() {
    //Reassign ids of internal nodes
    NodeVector pnv;
    this->getInternalNodes(pnv);
    int first_old_interior_id = aln->getNSeq32();
    int node_count = static_cast<int>(pnv.size());
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (int i=0; i<node_count; ++i) {
        pnv[i]->id = first_old_interior_id + i;
    }
    return first_old_interior_id + static_cast<int>(pnv.size());
}



void PhyloTree::addNewTaxaToTree(const IntVector& taxaIdsToAdd,
                                 const char* description,
                                 const char* placement_parameter_string,
                                 bool be_quiet) {
    //
    //Assumes: The tree is rooted.
    //
    if (0==strlen(placement_parameter_string)) {
        placement_parameter_string = Params::getInstance().incremental_method.c_str();
    }
    PlacementParameters place_params(placement_parameter_string);
    PlacementRun pr(*this, place_params, taxaIdsToAdd, be_quiet);
    deleteAllPartialLhAndParsimony();
    
    bool trackLikelihood      = shouldPlacementUseLikelihood();
    int  additional_sequences = static_cast<int>(taxaIdsToAdd.size());

    //Todo: Change how the caller selects step-wise parsimony!
    if ( pr.taxa_per_batch == 1 && pr.heuristic->isGlobalSearch()  &&
         !pr.calculator->usesLikelihood() && !trackLikelihood &&
         3 < leafNum ) {
        //The minimum for leafNum is because ConstraintTree::initFromTree
        //(called from ConstraintTree::readConstraint, called from
        //reinsertTaxViaStepwiseParsimony(), refuses to run if leafNum
        //is less than 4).
        
        pr.setUpAllocator(additional_sequences*4, false, 0);
        reinsertTaxaViaStepwiseParsimony(taxaIdsToAdd);
        pr.donePlacement();
        return;
    }

    //Todo: What about # of likelihood vectors when
    //      likelihood vectors aren't being allocated per node.
    //      extra_lh_blocks should be even less, in that case.
    int sequences              = static_cast<int>(aln->getNSeq());
    int target_branch_count    = sequences * 2 - 3;
    int extra_parsimony_blocks = target_branch_count + additional_sequences * 4;
        //each target branch (including those for the sequences that
        //aren't yet in the tree) needs a parsimony block, and
        //for each additional taxon to place, there are 4 additional
        //                        parsimony blocks needed.
        //              A         (Two on the branch between the added
        //              |          taxon and its interior, I, two on each
        //              | (2)      of the branches that connect the
        //              |          taxon's interior, I, to existing interior
        //              I          nodes, L and R, minus the two on the branch
        //             / \         that formerly connected L and R directly
        //        (2) /   \ (2)    to each other.
        //           /     \
        //          L       R
        //           x(-2)-x
    
    intptr_t fixed_spr_radius = 3;
    intptr_t blocksPerThread  = ParsimonySPRMove::getParsimonyVectorSize(fixed_spr_radius);
    intptr_t threadsNeeded    = ParsimonySPRMove::getMinimumPathVectorCount();
    intptr_t threadCount      = PhyloTreeThreadingContext::getMaximumThreadCount();
    ParsimonyPathVector pv(blocksPerThread, threadsNeeded, threadCount);
    intptr_t spr_blocks       = pv.getTotalNumberOfBlocksRequired();
    extra_parsimony_blocks   += static_cast<int>(spr_blocks);
    
    int extra_lh_blocks = trackLikelihood
                          ? (target_branch_count + additional_sequences * 4)
                          : 0;
        //each target branch needs its own lh block, and
        //although each TaxonToPlace just needs one lh block
        //each addition of a taxon later creates 3 *new* target branches
        //(each of which needs its own lh block)
    
    intptr_t newTaxaCount           = taxaIdsToAdd.size();
    double   estimate_per_placement = 1.0 + floor(newTaxaCount / leafNum / 2.0);
    double   estimate               = newTaxaCount * (2.0 + estimate_per_placement);
    initProgress(estimate, description, "", "");

    pr.overall.start();
    pr.initializing.start();
    pr.prepareForPlacementRun();
    pr.setUpAllocator(extra_parsimony_blocks, trackLikelihood, extra_lh_blocks);
    
    //TimeKeeper optimizing("optimizing");
    //Allocate per-thread parsimony vector work areas used to calculate
    //modified parsimony scores along the path between the
    //pruning and regrafting points.
    auto num_path_vectors = pv.getNumberOfPathsRequired();
    pv.resize(num_path_vectors);
    for (int vector=0; vector<num_path_vectors; ++vector) {
        pr.block_allocator->allocateVectorOfParsimonyBlocks
            (blocksPerThread, pv[vector]);
    }
    pr.initializing.stop();

    double setUpStartTime = getRealTime();
    pr.searchTime.start();
    LOG_LINE ( VerboseMode::VB_DEBUG, "Before allocating TaxonToPlace array"
              << ", index_lh was "
              << pr.block_allocator->getLikelihoodBlockCount() );
    
    TypedTaxaToPlace<TaxonTypeInUse> candidates(newTaxaCount);
    int first_new_interior_id = renumberInternalNodes();
    for (int i=0; i<newTaxaCount; ++i) {
        int         taxonId   = taxaIdsToAdd[i];
        std::string taxonName = aln->getSeqName(taxonId);        
        candidates.emplace_back(pr.block_allocator, first_new_interior_id+i,
                                taxonId, taxonName, true);
    }
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (intptr_t i=0; i<newTaxaCount; ++i) {
        candidates[i].computeParsimony(this);
        if ((i%1000) == 999) {
            trackProgress(1000.0);
        }
    }
    trackProgress((double)(newTaxaCount % 1000));
    pr.searchTime.stop();

    LOG_LINE ( VerboseMode::VB_DEBUG, "After allocating TaxonToPlace"
               << ", index_lh was " << pr.block_allocator->getLikelihoodBlockCount()
               << ", index_pars was " << pr.block_allocator->getParsimonyBlockCount());

    TargetBranchRange targets(*this, pr.block_allocator, pr.calculator, false);
    LOG_LINE ( VerboseMode::VB_DEBUG, "After allocating TargetBranchRange"
               << ", index_lh was " << pr.block_allocator->getLikelihoodBlockCount()
               << ", index_pars was " << pr.block_allocator->getParsimonyBlockCount());
    
    if (!be_quiet) {
        LOG_LINE ( VerboseMode::VB_MIN, "Placement set-up time was "
                   << (getRealTime() - setUpStartTime) << " sec");
    }
    
    ParsimonySearchParameters s("SPR");
    s.iterations                 = 3;
    s.lazy_mode                  = false;
    s.radius                     = static_cast<int>(fixed_spr_radius);
    s.calculate_connection_costs = true;
    s.be_quiet                   = true;

    LikelihoodBlockPairs spare_blocks(2);
    for (; 0<newTaxaCount; newTaxaCount = candidates.size() ) {
        if (newTaxaCount<static_cast<intptr_t>(pr.taxa_per_batch)) {
            pr.taxa_per_batch = newTaxaCount;
        }
        size_t batchStart=0;
        for (; static_cast<intptr_t>(batchStart+pr.taxa_per_batch) <= newTaxaCount
             ; batchStart+=pr.taxa_per_batch) {
            pr.prepareForBatch();

            size_t batchStop  = batchStart + pr.taxa_per_batch;
            pr.doBatchPlacementCosting(candidates, batchStart, batchStop, targets);

            size_t insertStop = batchStart;
            pr.selectPlacementsForInsertion(candidates, batchStart, batchStop,
                                            insertStop);

            pr.startBatchInsert();
            pr.doBatchInsert(candidates, batchStart, insertStop, spare_blocks,
                             estimate_per_placement, targets,
                             s, pv);
            pr.doneBatch(candidates, batchStart, batchStop, targets);
        } //batches of items
        
        pr.donePass(candidates, batchStart, targets);
    }
    doneProgress();
    
    if (!be_quiet) {
        LOG_LINE ( VerboseMode::VB_MED, "Total number of blocked inserts was "  << pr.taxa_inserted_nearby );
        LOG_LINE ( VerboseMode::VB_MED, "At the end of addNewTaxaToTree, index_lhs was "
                  << pr.block_allocator->getLikelihoodBlockCount() << ", index_pars was "
                  << pr.block_allocator->getParsimonyBlockCount() << ".");
    }
    
    pr.refreshTime.start();
    clearAllPartialParsimony(false);
    clearAllPartialLH();
    pr.refreshTime.stop();
        
    pr.donePlacement();
    
    if (VerboseMode::VB_MIN <= verbose_mode && !be_quiet) {
        pr.reportActivity();
    }
}
