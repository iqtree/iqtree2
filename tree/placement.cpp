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
    size_t countOfTaxaToRemove = Placement::getNumberOfTaxaToRemove(nseq);
    if (0<countOfTaxaToRemove) {
        map<string, Node*> mapNameToNode;
        getMapOfTaxonNameToNode(nullptr, nullptr, mapNameToNode);
        size_t r = 0;
        for (size_t seq = 0; seq < nseq; ++seq) {
            r += countOfTaxaToRemove;
            if ( r >= nseq ) {
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
    return Placement::getCostFunction() == Placement::SANKOFF_PARSIMONY;
}

bool PhyloTree::shouldPlacementUseLikelihood() const {
    Placement::CostFunction       costFunction    = Placement::getCostFunction();
    return ( costFunction != Placement::MAXIMUM_PARSIMONY
             && costFunction != Placement::SANKOFF_PARSIMONY);
}

namespace {
    void logInsert(PhyloTree* tree, Params& params,
                   Placement::CostFunction costFunction,
                   size_t totalInsertCount, const char* verb,
                   TaxonToPlace & c, const char * where) {
        if (( verbose_mode >= VB_MIN && !params.suppress_list_of_sequences)
            || verbose_mode >= VB_MED ) {
            stringstream s;
            s << totalInsertCount << ". " << verb << " "
                << c.taxonName << " " << where << ". It had ";
            const PossiblePlacement& p = c.getBestPlacement();
            if (costFunction==Placement::MAXIMUM_PARSIMONY ||
                costFunction==Placement::SANKOFF_PARSIMONY) {
                s << "parsimony score " << (int)(p.score);
            } else {
                s << "likelihood score " << p.score;
            }
            s << " (and path lengths " << p.lenToNode1
                << ", " << p.lenToNode2 << ", " << p.lenToNewTaxon << ")";
            tree->logLine(s.str());
        }
    }
}

typedef TaxonToPlace TaxonTypeInUse;
//typedef LessFussyTaxon TaxonTypeInUse;

#define NEW_TAXON_MAJOR (0)

void PhyloTree::addNewTaxaToTree(const IntVector& taxaIdsToAdd) {
    //
    //Assumes: The tree is rooted.
    //
    Params&            params          = Params::getInstance();
    size_t             taxaPerBatch    = Placement::getTaxaPerBatch(taxaIdsToAdd.size());
                                         //Must be 1 or more
    size_t             insertsPerBatch = Placement::getInsertsPerBatch(taxaIdsToAdd.size(), taxaPerBatch);
                                         //Must be 1 or more
    PlacementRun pr(*this);
    deleteAllPartialLh();
    
    if ( taxaPerBatch == 1 && pr.heuristic->isGlobalSearch()  &&
        ( pr.costFunction == Placement::MAXIMUM_PARSIMONY
         || pr.costFunction == Placement::SANKOFF_PARSIMONY ) ) {
        //For now, we might as well use the existing step-wise
        //parsimony stuff for adding to a constraint tree, eh?
        //Since, for now, it is a lot faster.
        constraintTree.readConstraint(*this);
        //clearing all the nodes...
        freeNode();
        root = nullptr;
        cout << "Creating fast initial parsimony tree by random order stepwise addition..." << endl;
        double start = getRealTime();
        double score = computeParsimonyTree(params.out_prefix, aln, randstream);
        cout << getRealTime() - start << " seconds, parsimony score: " << score
            << " (based on " << aln->num_parsimony_sites << " sites)"<< endl;
        
        //Note that this score tends to disagree.
        double parsimonyStart = getRealTime();
        clearAllPartialParsimony(false);
        double parsimonyScore = computeParsimony("Recalculating parsimony score");
        LOG_LINE( VB_MED, "Recalculated parsimony score " << parsimonyScore
                    << " (recalculation cost " << (getRealTime() - parsimonyStart) << " sec)" );

        finishUpAfterTaxaAddition();
        return;
    }
    
    bool trackLikelihood  = shouldPlacementUseLikelihood();

    //Todo: Recalculate what both of these should be.
    //      Their formulas were once spot-on, but that was
    //      a while back, and the pacement code now uses
    //      fewer likelihood vectors (in particular).
    //      The correct # of extras required is probably
    //      more like nodeNum - 1.
    //Todo: What about # of likelihood vectors when
    //      likelihood vectors aren't being allocated per node.
    //      extra_lh_blocks should be even less, in that case.
    //
    size_t   extra_parsimony_blocks = leafNum * 2 - 4;
    size_t   extra_lh_blocks        = trackLikelihood
                                    ? (leafNum * 4 - 4 + taxaIdsToAdd.size())
                                    : 0;
    pr.setUpAllocator(extra_parsimony_blocks, trackLikelihood, extra_lh_blocks);
    if (pr.costFunction == Placement::SANKOFF_PARSIMONY) {
        computeTipPartialParsimony();
    }
    
    LOG_LINE ( VB_MED, "After overallocating lh blocks, index_lh was "
              << pr.block_allocator->getLikelihoodBlockCount() );
    if (VB_MED <= verbose_mode) {
        curScore = computeLikelihood();
        LOG_LINE ( VB_MED, "Likelihood score before insertions was " << curScore );
        #if (0)
            curScore = optimizeAllBranches(2);
            LOG_LINE ( VB_MED, "Optimized likelihood score before insertions was " << curScore);
        #endif
    }
    LOG_LINE ( VB_MED, "Batch size is " << taxaPerBatch
              << " and the number of inserts per batch is " << insertsPerBatch);
    
    double setUpStartTime = getRealTime();
    size_t newTaxaCount = taxaIdsToAdd.size();
    
    TaxaToPlace<TaxonTypeInUse> candidates(newTaxaCount);
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
                      ( newTaxaCount, taxaPerBatch, insertsPerBatch );
    size_t totalInsertCount     = 0;
    size_t blockedInsertCount   = 0;
    double timeSpentOnRefreshes = 0.0; //Time spent recalculating parsimony &/or likelihood
                                       //for the entire tree
    double timeSpentOnSearches  = 0.0;
    double timeSpentOnInserts   = 0.0;
    initProgress(estimate, "Adding new taxa to tree", "", "");
    while (0<newTaxaCount) {
        if (newTaxaCount<taxaPerBatch) {
            taxaPerBatch = newTaxaCount;
        }
        size_t batchStart=0;
        for (; batchStart+taxaPerBatch <= newTaxaCount; batchStart+=taxaPerBatch) {
            timeSpentOnRefreshes -= getRealTime();
            if (trackLikelihood) {
                clearAllPartialLH(false);
                clearAllScaleNum(false);
                double likelihoodScore = computeLikelihood();
                LOG_LINE( VB_MIN, "Log-likelihood is currently " << likelihoodScore);
            }
            size_t batchStop = batchStart + taxaPerBatch;
            TargetBranch* pointStart = targets.data();
            TargetBranch* pointStop  = pointStart + targets.size();
            clearAllPartialParsimony(false);
            timeSpentOnRefreshes += getRealTime();
            timeSpentOnSearches -= getRealTime();
            TaxonTypeInUse* candidateStart = candidates.data() + batchStart;
            TaxonTypeInUse* candidateStop  = candidates.data() + batchStop;
#if (NEW_TAXON_MAJOR)
            for (TargetBranch* point = pointStart; point<pointStop; ++point) {
                point->computeState(*this);
            }
            for (TaxonTypeInUse* c = candidateStart; c<candidateStop; ++c) {
                LOG_LINE(VB_DEBUG, "Scoring ... " << c->taxonName);
                c->findPlacement(*this, targets,
                                 pr.heuristic, pr.calculator);
                const PossiblePlacement& p = c->getBestPlacement();
                LOG_LINE(VB_DEBUG, "Scored " << p.score << " for placement"
                         << " of " << c->taxonName << " with lengths "
                         << p.lenToNode1 << ", " << p.lenToNode2 << ", " << p.lenToNewTaxon);
            }
#else //INSERTION_POINT_MAJOR
            for (TargetBranch* point = pointStart; point<pointStop; ++point) {
                LOG_LINE(VB_DEBUG, "Scoring target branch " << (point-pointStart) << " of " << (pointStop-pointStart));
                point->computeState(*this);
                point->costPlacementOfTaxa(*this, &targets, point-pointStart,
                                           candidateStart, candidateStop,
                                           pr.heuristic, pr.calculator,
                                           point==pointStart );
                point->forgetState();
            }
#endif
            timeSpentOnSearches += getRealTime();
            insertsPerBatch      = Placement::getInsertsPerBatch(taxaIdsToAdd.size(), batchStop-batchStart);
            size_t insertStop    = batchStart + insertsPerBatch;
            std::sort( candidates.begin() + batchStart, candidates.begin() + batchStop);
            if (batchStop <= insertStop) {
                insertStop = batchStop; //Want them all
            }
            timeSpentOnInserts -= getRealTime();
            size_t insertCount = 0;
            for ( size_t i = batchStart; i<insertStop; ++i) {
                TaxonToPlace& c = candidates[i];
                if (c.canInsert()) {
                    ++insertCount;
                    ++totalInsertCount;
                    c.insertIntoTree(*this, pr.block_allocator, targets, *pr.calculator);
                    logInsert(this, params, pr.costFunction, totalInsertCount,
                              "Inserted", c, "at its preferred branch");
                } else {
                    //Another candidate taxon has gotten there first
                    ++blockedInsertCount;
                    ++insertCount;
                    ++totalInsertCount;
                    c.insertNearby(*this, pr.block_allocator, targets, *pr.calculator);
                    logInsert(this, params, pr.costFunction, totalInsertCount,
                              "Inserted", c, "near its preferred branch");
                }
                pr.taxon_placement_optimizer->cleanUpAfterTaxonPlacement(c, this);
            }
            timeSpentOnInserts += getRealTime();
            if ( 1 < batchStop - batchStart ) {
                LOG_LINE ( VB_MED,  "Inserted " << (insertCount)
                          << " out of a batch of " << (batchStop - batchStart) << "." );
            }
            pr.batch_placement_optimizer->cleanUpAfterBatch(candidates, batchStart, batchStop, this);
            if (trackLikelihood) {
                fixNegativeBranch();
            }
            if (insertCount == 0) {
                outError("No taxa inserted in batch");
                break;
            }
        } //batches of items
        
        targets.removeUsed();
        //Remove all the candidates that we were able to place
        std::vector<TaxonTypeInUse> oldCandidates;
        std::swap(oldCandidates, candidates);
        //1. Any candidates not considered this time go to the
        //   first batch to consider in the next pass.
        for (size_t r=batchStart; r<newTaxaCount; ++r) {
            candidates.emplace_back(oldCandidates[r]);
        }
        //2. Any candidates that were considered, but were not
        //   inserted, are to be considered in the next pass.
        for (size_t r=0; r<batchStart; ++r) {
            if (!oldCandidates[r].inserted) {
                //Keep this one to be considered next time
                candidates.emplace_back(oldCandidates[r]);
            }
        }
        newTaxaCount = candidates.size();
        insertsPerBatch = Placement::getInsertsPerBatch(taxaIdsToAdd.size(), taxaPerBatch);
        auto workLeft   = taxaAdditionWorkEstimate
                          ( newTaxaCount, taxaPerBatch, insertsPerBatch );
        this->progress->setWorkRemaining(workLeft);
        LOG_LINE ( VB_MAX, "At the end of this pass, index_lhs was "
                  << pr.block_allocator->getLikelihoodBlockCount() << ", index_pars was "
                  << pr.block_allocator->getParsimonyBlockCount());
    }
    doneProgress();
    
    LOG_LINE ( VB_MED, "Tidying up tree after inserting taxa.");
    pr.global_placement_optimizer->cleanUpAfterPlacement(this);
    
    LOG_LINE ( VB_MIN, "Time spent on refreshes was "      << timeSpentOnRefreshes << " sec");
    LOG_LINE ( VB_MIN, "Time spent on searches was "       << timeSpentOnSearches << " sec");
    LOG_LINE ( VB_MIN, "Time spent on actual inserts was " << timeSpentOnInserts << " sec");
    LOG_LINE ( VB_MIN, "Total number of blocked inserts was " << blockedInsertCount );
    LOG_LINE ( VB_MED, "At the end of addNewTaxaToTree, index_lhs was "
              << pr.block_allocator->getLikelihoodBlockCount() << ", index_pars was "
              << pr.block_allocator->getParsimonyBlockCount() << ".");
    if (!trackLikelihood) {
        fixNegativeBranch();
    }
    finishUpAfterTaxaAddition();
    
    if (VB_MED <= verbose_mode) {
        PhyloNodeVector idToNode;
        getArrayOfTaxaNodesById(nullptr, nullptr, idToNode);
        for (size_t i=0; i<taxaIdsToAdd.size(); ++i) {
            PhyloNode*     leaf      = idToNode[taxaIdsToAdd[i]];
            if ( leaf != nullptr ) {
                PhyloNeighbor* upLink    = leaf->firstNeighbor();
                PhyloNode*     upNode    = upLink->getNode();
                PhyloNeighbor* leftLink  = nullptr;
                PhyloNeighbor* rightLink = nullptr;
                double leftLength  = -1;
                double rightLength = -1;
                FOR_EACH_PHYLO_NEIGHBOR(upNode, leaf, itNei, nei) {
                    if (leftLink==nullptr) {
                        leftLink   = nei;
                        leftLength = nei->length;
                    } else {
                        rightLink   = nei;
                        rightLength = nei->length;
                    }
                }
                auto length = leaf->firstNeighbor()->length;
                std::cout << (i+1) << "." << "Node [" << leaf->id << "]=" << leaf->name
                    << " now has branch length " << length
                    << " (interior left branch " << leftLength
                    << " , and right branch " << rightLength << ")"
                    << std::endl;
            }
        }
    }
}

void PhyloTree::finishUpAfterTaxaAddition() {
    initializeTree();
    deleteAllPartialLh();
    initializeAllPartialLh();
    LOG_LINE ( VB_MED, "Number of leaves " << this->leafNum
              << ", of nodes " << this->nodeNum );
    Placement::CostFunction costFunction = Placement::getCostFunction();
    if (costFunction==Placement::MAXIMUM_LIKELIHOOD_ANYWHERE ||
        costFunction==Placement::MAXIMUM_LIKELIHOOD_MIDPOINT) {
        double score = optimizeAllBranches();
        LOG_LINE ( VB_MIN, "After optimizing, likelihood score was " << score );
    }
}
