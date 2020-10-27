//
//  placement.cpp
//  Evolutionary placement: adding taxa to a phylogenetic tree
//  This file created by James Barbetti on 25/8/20, but:
//  1. addNewTaxaToTree was formerly in phylotree.cpp;
//  2. addTaxonML likewise (and was the work of BUI Quang Minh)
//
#include <tree/phylotree.h>
#include <placement/blockallocator.h>
#include <placement/placement.h>
#include <placement/searchheuristic.h>
#include <placement/targetbranch.h>
#include <placement/taxontoplace.h>
#include <placement/placementcostcalculator.h>
#include <placement/placementoptimizer.h>
#include <placement/placementrun.h>

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
    
    LOG_LINE(VB_DEBUG, "  Placement branch length " << len);

    FOR_EACH_PHYLO_NEIGHBOR(added_node, nullptr, it, nei) {
        nei->clearComputedFlags();
        nei->getNode()->findNeighbor(added_node)->clearComputedFlags();
    }
    
    //compute the likelihood
    PhyloNeighbor* nei;
    double best_score = 0;
    if (isAddedAtMidpoint) {
        len_to_new_taxon   = recomputeParsimonyBranchLength(added_taxon, added_node);
        LOG_LINE(VB_DEBUG, "  Parsimony taxon->interior length " << len_to_new_taxon );
        nei                = added_taxon->findNeighbor(added_node);
        best_score         = computeLikelihoodBranch(nei, added_taxon,
                                                     tree_buffers);
        LOG_LINE(VB_DEBUG, "  Traversal info size is " << traversal_info.size());
        LOG_LINE(VB_DEBUG, "  Likelihood before optimization " << best_score);
        optimizeOneBranch(added_taxon, added_node, false, 20);
        len_to_target_dad  = halfLen;
        len_to_target_node = halfLen;
        len_to_new_taxon   = nei->length;
        best_score         = computeLikelihoodFromBuffer();
        LOG_LINE(VB_DEBUG, "  Likelihood after optimization " << best_score
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
    LOG_LINE(VB_DEBUG, "  ML Lengths " << len_to_target_dad
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
    size_t nseq = aln->getNSeq();
    size_t countOfTaxaToRemove = Placement::getNumberOfTaxaToRemoveAndReinsert(nseq);
    if (0<countOfTaxaToRemove) {
        LOG_LINE(VB_DEBUG, "Will mark " << countOfTaxaToRemove << " (out of " << nseq
                    << ") sequences, to be removed.");
        map<string, Node*> mapNameToNode;
        getMapOfTaxonNameToNode(nullptr, nullptr, mapNameToNode);
        size_t r = 0;
        LOG_LINE(VB_MAX, "Before removing sequences:");
        for (auto it=mapNameToNode.begin(); it!=mapNameToNode.end(); ++it) {
            LOG_LINE( VB_MAX, "sequence " << it->first << " has id " << it->second->id );
        }

        for (size_t seq = 0; seq < nseq; ++seq) {
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
                    LOG_LINE(VB_MED, "Marking sequence " << seq << " (" << node->name << ") to be removed.");
                }
            }
        }
    }
}

double PhyloTree::taxaAdditionWorkEstimate(size_t newTaxaCount,
                                           size_t taxaPerBatch,
                                           size_t insertsPerBatch) {
    if ( newTaxaCount <= taxaPerBatch || taxaPerBatch == 0 ) {
        if ( newTaxaCount <= insertsPerBatch || insertsPerBatch == 0 ) {
            return 3.0 * newTaxaCount * leafNum;
        }
        return 3.0 * newTaxaCount * leafNum
                   * newTaxaCount / insertsPerBatch;
    }
    size_t batchesThisPass  = newTaxaCount / taxaPerBatch;
    double workThisPass     = batchesThisPass * taxaPerBatch * leafNum;
    double progressThisPass = batchesThisPass * insertsPerBatch;
    //Optimistic if inserts = 100% and batches are large.
    return (3.0 * workThisPass / progressThisPass) * newTaxaCount;
}

bool PhyloTree::shouldPlacementUseSankoffParsimony() const {
    return Placement::doesPlacementUseParsimony();
}

bool PhyloTree::shouldPlacementUseLikelihood() const {
    return Placement::doesPlacementUseLikelihood();
}

void PhyloTree::reinsertTaxaViaStepwiseParsimony(const IntVector& taxaIdsToAdd) {
    constraintTree.readConstraint(*this);
    //clearing all the nodes...
    freeNode();
    root = nullptr;
    cout << "Creating fast initial parsimony tree by random order stepwise addition..." << endl;
    double start = getRealTime();
    double score = computeParsimonyTree(Params::getInstance().out_prefix, aln, randstream);
    cout << getRealTime() - start << " seconds, parsimony score: " << score
        << " (based on " << aln->num_parsimony_sites << " sites)"<< endl;

    //Note that this score tends to disagree.
    double parsimonyStart = getRealTime();
    clearAllPartialParsimony(false);
    double parsimonyScore = computeParsimony("Recalculating parsimony score");
    LOG_LINE( VB_MED, "Recalculated parsimony score " << parsimonyScore
                << " (recalculation cost " << (getRealTime() - parsimonyStart) << " sec)" );
}

typedef TaxonToPlace TaxonTypeInUse;
//typedef LessFussyTaxon TaxonTypeInUse;


void PhyloTree::addNewTaxaToTree(const IntVector& taxaIdsToAdd) {
    //
    //Assumes: The tree is rooted.
    //
    PlacementRun pr(*this, taxaIdsToAdd);
    deleteAllPartialLh();
    
    if ( pr.taxa_per_batch == 1 && pr.heuristic->isGlobalSearch()  &&
         !pr.calculator->usesLikelihood() ) {
        reinsertTaxaViaStepwiseParsimony(taxaIdsToAdd);
        pr.donePlacement();
        return;
    }
    
    bool trackLikelihood  = shouldPlacementUseLikelihood();

    //Todo: What about # of likelihood vectors when
    //      likelihood vectors aren't being allocated per node.
    //      extra_lh_blocks should be even less, in that case.
    size_t   sequences              = aln->getNSeq();
    size_t   target_branch_count    = sequences * 2 - 1;
    size_t   additional_sequences   = taxaIdsToAdd.size();
    size_t   extra_parsimony_blocks = target_branch_count + additional_sequences * 4;
        //each target branch needs a parsimony block, and
        //although each taxon to place just needs one parsimony block
        //(looking down from the new interior to the new exterior), 
        //each addition of a taxon later creates 3 *new* target branches, 
        //each of which needs a new "upward" looking parsimony block.  
        //(hence the: multiplication, of additional_sequences, by 4).

    size_t   extra_lh_blocks        = trackLikelihood
                                    ? (target_branch_count + additional_sequences * 4)
                                    : 0;
        //each target branch needs its own lh block, and
        //although each TaxonToPlace just needs one lh block
        //each addition of a taxon later creates 3 *new* target branches
        //(each of which needs its own lh block)
    
    pr.setUpAllocator(extra_parsimony_blocks, trackLikelihood, extra_lh_blocks);
    pr.prepareForPlacementRun();
    
    double setUpStartTime = getRealTime();
    size_t newTaxaCount = taxaIdsToAdd.size();
    
    TypedTaxaToPlace<TaxonTypeInUse> candidates(newTaxaCount);
    LOG_LINE ( VB_DEBUG, "Before allocating TaxonToPlace array"
              << ", index_lh was " << pr.block_allocator->getLikelihoodBlockCount() );
    for (size_t i=0; i<newTaxaCount; ++i) {
        int         taxonId   = taxaIdsToAdd[i];
        std::string taxonName = aln->getSeqName(taxonId);
        candidates.emplace_back(pr.block_allocator, taxonId, taxonName);
    }
    LOG_LINE ( VB_DEBUG, "After allocating TaxonToPlace"
              << ", index_lh was " << pr.block_allocator->getLikelihoodBlockCount()
              << ", index_pars was " << pr.block_allocator->getParsimonyBlockCount());

    TargetBranchRange targets(*this, pr.block_allocator, pr.calculator);
    LOG_LINE ( VB_DEBUG, "After allocating TargetBranchRange"
              << ", index_lh was " << pr.block_allocator->getLikelihoodBlockCount()
              << ", index_pars was " << pr.block_allocator->getParsimonyBlockCount());
    LOG_LINE ( VB_MIN, "Set up time was " << (getRealTime() - setUpStartTime) << " sec");
        
    double estimate = taxaAdditionWorkEstimate
    ( newTaxaCount, pr.taxa_per_batch, pr.inserts_per_batch );
    double timeSpentOnRefreshes = 0.0; //Time spent recalculating parsimony &/or likelihood
                                       //for the entire tree
    double timeSpentOnSearches  = 0.0;
    double timeSpentOnRanking   = 0.0;
    double timeSpentOnInserts   = 0.0;
    double timeSpentOnCleanUp   = 0.0;
    
    LikelihoodBlockPairs spare_blocks(2);
    initProgress(estimate, "Adding new taxa to tree", "", "");
    for (; 0<newTaxaCount; newTaxaCount = candidates.size() ) {
        if (newTaxaCount<pr.taxa_per_batch) {
            pr.taxa_per_batch = newTaxaCount;
        }
        size_t batchStart=0;
        for (; batchStart+pr.taxa_per_batch <= newTaxaCount
             ; batchStart+=pr.taxa_per_batch) {
            timeSpentOnRefreshes -= getRealTime();
            pr.prepareForBatch();
            timeSpentOnRefreshes += getRealTime();

            timeSpentOnSearches  -= getRealTime();
            size_t batchStop      = batchStart + pr.taxa_per_batch;
            pr.doBatchPlacementCosting(candidates, batchStart, batchStop, targets);
            timeSpentOnSearches  += getRealTime();
            
            timeSpentOnRanking  -= getRealTime();
            size_t insertStop    = batchStart;
            pr.selectPlacementsForInsertion( candidates, batchStart, batchStop, insertStop);
            timeSpentOnRanking  += getRealTime();
            
            timeSpentOnInserts -= getRealTime();
            pr.startBatchInsert();
            for ( size_t i = batchStart; i<insertStop; ++i) {
                pr.insertTaxon(candidates, i, targets, spare_blocks);
            }
            timeSpentOnInserts += getRealTime();
            
            timeSpentOnCleanUp -= getRealTime();
            pr.doneBatch(candidates, batchStart, batchStop, targets);
            timeSpentOnCleanUp += getRealTime();

        } //batches of items
        
        timeSpentOnCleanUp -= getRealTime();
        pr.donePass(candidates, batchStart, targets);
        timeSpentOnCleanUp += getRealTime();

        auto workLeft   = taxaAdditionWorkEstimate
                          ( newTaxaCount, pr.taxa_per_batch, pr.inserts_per_batch );
        this->progress->setWorkRemaining(workLeft);
    }
    doneProgress();
    
    LOG_LINE ( VB_MIN, "Time spent on refreshes was "          << timeSpentOnRefreshes << " sec");
    LOG_LINE ( VB_MIN, "Time spent on searches was "           << timeSpentOnSearches << " sec");
    LOG_LINE ( VB_MIN, "Time spent on actual inserts was "     << timeSpentOnInserts << " sec");
    LOG_LINE ( VB_MIN, "Time for post-batch optimization was " << timeSpentOnCleanUp << " sec");
    LOG_LINE ( VB_MIN, "Total number of blocked inserts was "  << pr.taxa_inserted_nearby );
    LOG_LINE ( VB_MED, "At the end of addNewTaxaToTree, index_lhs was "
              << pr.block_allocator->getLikelihoodBlockCount() << ", index_pars was "
              << pr.block_allocator->getParsimonyBlockCount() << ".");

    pr.donePlacement();
}

