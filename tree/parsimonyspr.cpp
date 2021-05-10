//
//  parsimonyspr.cpp
//  Parsimony SPR implementation
//  (ParsimonyLazySPRMove scores possible SPR moves by considering
//   only the length of the moving branch before and after a  move,
//   ParsimonySPRMove scores them by properly calculating the
//   parsimony score that will obtain if a move is done)
//  (In practice ParsimonySPRMove is a lot more useful, because
//   ParsimonyLazySPRMove suggests too many moves that don't really
//   result in a benefit) (but with a better heuristic
//   ParsimonyLazySPRMove might be competitive!)
//  Created by James Barbetti on 8-Dec-2020.
//

#include "parsimonyspr.h"
#include "parsimonysearch.h"
#include <placement/targetbranch.h>            //for TargetBranchRange
#include <placement/placementcostcalculator.h> //for ParsimonyCostCalculator
#include <utils/timeutil.h>                    //for getRealTime

PhyloNodeVector PhyloTree::getTaxaNodesInIDOrder() const {
    PhyloNodeVector taxa;
    getTaxa(taxa);
    PhyloNodeVector result;
    result.resize(leafNum, nullptr);
    for (auto it=taxa.begin(); it!=taxa.end(); ++it) {
        auto node = *it;
        ASSERT(0<=node->id && node->id<static_cast<int>(leafNum));
        ASSERT(result[node->id] == nullptr);
        result[node->id] = node;
    }
    for (auto itResult=result.begin(); itResult!=result.end(); ++itResult) {
        ASSERT((*itResult)!=nullptr);
    }
    return result;
}

ParsimonyLazySPRMove::LazySPRSearch::LazySPRSearch
    (const PhyloTree& phylo_tree, const TargetBranchRange& target_branches,
     double disconnection_benefit, ParsimonyLazySPRMove& output)
    : tree(phylo_tree), branches(target_branches)
    , source(target_branches[output.source_branch_id])
    , discon(disconnection_benefit)
    , put_answer_here(output) {
}

void ParsimonyLazySPRMove::LazySPRSearch::searchForForwardsSPR
    (PhyloNode* current, PhyloNode* prev, int radius) {
    --radius;
    FOR_EACH_ADJACENT_PHYLO_NODE(current, prev, it, next) {
        if (0<radius) {
            searchForForwardsSPR(next, current, radius);
        }
        int target_branch_id = (*it)->id;
        const TargetBranch& target = branches[target_branch_id];
        double cost    = source.getForwardConnectionCost(tree, target);
        double benefit = discon - cost;
        if (put_answer_here.benefit<benefit) {
            put_answer_here.benefit          = benefit;
            put_answer_here.target_branch_id = target_branch_id;
            put_answer_here.isForward        = true;
        }
        ++put_answer_here.positions_considered;
    }
}

void ParsimonyLazySPRMove::LazySPRSearch::searchForBackwardsSPR
    ( PhyloNode* current, PhyloNode* prev, int radius) {
    --radius;
    FOR_EACH_ADJACENT_PHYLO_NODE(current, prev, it, next) {
        if (0<radius) {
            searchForBackwardsSPR(next, current, radius);
        }
        int target_branch_id = (*it)->id;
        const TargetBranch& target = branches[target_branch_id];
        double cost    = source.getBackwardConnectionCost(tree, target);
        double benefit = discon - cost;
        if (put_answer_here.benefit<benefit) {
            put_answer_here.benefit          = benefit;
            put_answer_here.target_branch_id = target_branch_id;
            put_answer_here.isForward        = false;
        }
        ++put_answer_here.positions_considered;
    }
}

ParsimonyLazySPRMove::ParsimonyLazySPRMove(): super() {
    initialize(0, true);
}

/*static*/ intptr_t ParsimonyLazySPRMove::getParsimonyVectorSize(intptr_t radius) {
    return 0;
}

void ParsimonyLazySPRMove::initialize(intptr_t source_branch, bool beLazy) {
    source_branch_id = source_branch;
    lazy             = beLazy;
    benefit          = -1.0;
    target_branch_id = -1;
    isForward        = false;
    source_first     = nullptr;
    source_second    = nullptr;
    target_first     = nullptr;
    target_second    = nullptr;
    positions_considered = 0;
}

std::string ParsimonyLazySPRMove::getDescription() const {
    std::stringstream s;
    s << " linking branch " << source_branch_id
      << ((isForward) ? " forward" : " backward" )
      << " to branch " << target_branch_id;
    return s.str();
}

void ParsimonyLazySPRMove::finalize(PhyloTree& tree,
              const TargetBranchRange& branches) {
    if (target_branch_id < 0) {
        return;
    }
    if (0 < benefit) {
        TREE_LOG_LINE(tree, VerboseMode::VB_DEBUG, 
            "move s=" << source_branch_id
            << ",d=" << target_branch_id
            << ", f=" << isForward
            << ", b=" << benefit);
    }
    auto source   = branches[source_branch_id];
    source_first  = source.first;
    source_second = source.second;
    auto target   = branches[target_branch_id];
    target_first  = (0<benefit) ? target.first  : nullptr;
    target_second = (0<benefit) ? target.second : nullptr;
}

void ParsimonyLazySPRMove::findForwardLazySPR(const PhyloTree& tree, const TargetBranchRange& branches,
                    int radius, double disconnection_benefit) {
    LazySPRSearch s(tree, branches, disconnection_benefit, *this);
    const TargetBranch& tb = branches[source_branch_id];
    PhyloNode* left;
    PhyloNode* right;
    GET_OTHER_ADJACENT_PHYLO_NODES(tb.first, tb.second,
                                   left, right);
    s.searchForForwardsSPR(left,  tb.first, radius);
    s.searchForForwardsSPR(right, tb.first, radius);
}

void ParsimonyLazySPRMove::findBackwardLazySPR(const PhyloTree& tree, const TargetBranchRange& branches,
                     int radius, double disconnection_benefit) {
    LazySPRSearch s(tree, branches, disconnection_benefit, *this);
    const TargetBranch& tb = branches[source_branch_id];
    PhyloNode* left;
    PhyloNode* right;
    GET_OTHER_ADJACENT_PHYLO_NODES(tb.second, tb.first,
                                   left, right);
    s.searchForBackwardsSPR(left,  tb.second, radius);
    s.searchForBackwardsSPR(right, tb.second, radius);
}

bool ParsimonyLazySPRMove::isStillPossible(const TargetBranchRange& branches,
                             PhyloBranchVector& path) const {
    path.clear();
    const TargetBranch& source = branches[source_branch_id];
    if (source.first  != source_first)  return false;
    if (source.second != source_second) return false;
    const TargetBranch& target = branches[target_branch_id];
    if (target.first  != target_first)  return false;
    if (target.second != target_second) return false;
    if (isForward) {
        return isAConnectedThroughBToC(source.second, source.first, target.first, path);
    } else {
        return isAConnectedThroughBToC(source.first, source.second, target.first, path);
    }
}

double ParsimonyLazySPRMove::recalculateBenefit
    ( PhyloTree& tree, double parsimony_score,
      TargetBranchRange& branches, LikelihoodBlockPairs &blocks,
      ParsimonyPathVector& parsimony_path_vectors) const {
    TargetBranch& source = branches[source_branch_id];
    TargetBranch& target = branches[target_branch_id];
    source.computeState(tree, parsimony_score, source_branch_id, blocks);
    target.computeState(tree, parsimony_score, target_branch_id, blocks);
    BenefitPair benefitPair = source.getPartialDisconnectionBenefit(tree, branches);
    double updated_benefit  = isForward ? benefitPair.forwardBenefit : benefitPair.backwardBenefit;
    double updated_cost     = isForward
                            ? source.getForwardConnectionCost(tree, target)
                            : source.getBackwardConnectionCost(tree, target);
    return updated_benefit - updated_cost;
}
//Note: callers may assume that ParsimonyLazySPRMove::apply
//      is its own inverse (that calling it a second time,
//      immediately after it is called the first time, will
//      reverse the changes that calling it the first time
//      made).
//
double ParsimonyLazySPRMove::apply(PhyloTree& tree,  double parsimony_score,
                                   TargetBranchRange& branches, LikelihoodBlockPairs blocks,
                                   ParsimonyPathVector& parsimony_path_vectors) {
    TargetBranch& source     = branches[source_branch_id];
    TargetBranch& target     = branches[target_branch_id];
    PhyloNode*    moved_node = isForward ? source.first : source.second;
    PhyloNode*    other_node = isForward ? source.second : source.first;
    PhyloNode*    new_left   = target.first;
    PhyloNode*    new_right  = target.second;
    PhyloNode*    snip_left;
    PhyloNode*    snip_right;
    GET_OTHER_ADJACENT_PHYLO_NODES(moved_node, other_node, snip_left, snip_right);
    int           snip_left_id  = snip_left->findNeighbor(moved_node)->id;
    int           snip_right_id = snip_right->findNeighbor(moved_node)->id;
    
    //Update node linkage in the tree (note, the neighbors
    //of moved_node have to be done in two phases, steps 3+4
    //and 5+6, just in case there's overlap, with one of
    //new_left, new_right being the same as one of snip_left,
    //and snip_right).
    //
    snip_left ->updateNeighbor(moved_node,   snip_right);
    snip_right->updateNeighbor(moved_node,   snip_left);
    moved_node->updateNeighbor(snip_left,    DUMMY_NODE_1);
    moved_node->updateNeighbor(snip_right,   DUMMY_NODE_2);
    moved_node->updateNeighbor(DUMMY_NODE_1, new_left);
    moved_node->updateNeighbor(DUMMY_NODE_2, new_right);
    new_left  ->updateNeighbor(new_right,    moved_node);
    new_right ->updateNeighbor(new_left,     moved_node);
    
    TargetBranch& left_branch  = branches[snip_left_id];
    left_branch.updateMapping(snip_left_id, new_left, moved_node, true);
            
    TargetBranch& right_branch = branches[snip_right_id];
    right_branch.updateMapping(snip_right_id, new_right, moved_node, true);
    
    target.updateMapping(target_branch_id, snip_left, snip_right, true);
    //Note: The branch that target_branch_id now refers to,
    //      is the branch that, were we reversing the SPR,
    //      would be the "new" target branch (it's the branch
    //      that came to be when we "snipped out" moved_node).
    //      It has to be this way, if apply() is to be its own
    //      inverse.
    
    double score=-1;
    score = left_branch .computeState (tree, score, snip_left_id,     blocks);
    score = right_branch.computeState (tree, score, snip_right_id,    blocks);
    score = target      .computeState (tree, score, target_branch_id, blocks);
    score = source      .computeState (tree, score, source_branch_id, blocks);
    
    left_branch.setParsimonyLength(tree);
    right_branch.setParsimonyLength(tree);
    target.setParsimonyLength(tree);
    source.setParsimonyLength(tree);

    TREE_LOG_LINE(tree, VerboseMode::VB_MAX,
                  "Updated parsimony score"
                  << " after applying SPR move was " << score);
    return score;
} //ParsimonyLazySPRMove::apply

/*static*/ intptr_t ParsimonySPRMove::getParsimonyVectorSize(intptr_t radius) {
    return radius + 1;
}

ParsimonySPRMove::ProperSPRSearch::ProperSPRSearch
    ( const PhyloTree& phylo_tree,
      const TargetBranchRange& target_branches,
      double disconnection_benefit,
      std::vector<UINT*>& path_parsimony_to_use,
      ParsimonyLazySPRMove& output)
    : super(phylo_tree, target_branches, disconnection_benefit, output),
      path_parsimony(path_parsimony_to_use) /*proper*/ {
}

PhyloNode* ParsimonySPRMove::ProperSPRSearch::other_adjacent_node(PhyloNode*a, PhyloNode*b, PhyloNode*c) {
    //Return the other node, adjacent to a, that isn't b or c
    FOR_EACH_ADJACENT_PHYLO_NODE(a, b, it, x) {
        if (x!=c) return x;
    }
    ASSERT(0 && "could not find other adjacent node");
    return nullptr;
}

void ParsimonySPRMove::ProperSPRSearch::prepareToSearch(PhyloNode* left, PhyloNode* right,
                     PhyloNode* snipped, int radius) {
    //Sets up path_parsimony[radius+1] to the view from
    //snipped, of the subtree containing right.
    path_parsimony.resize(radius+2, nullptr);
    path_parsimony[radius+1] = snipped->findNeighbor(right)->get_partial_pars();
    discon = tree.getSubTreeParsimony(left->findNeighbor(snipped))
           - tree.getSubTreeParsimony(snipped->findNeighbor(right));
}

void ParsimonySPRMove::ProperSPRSearch::searchForForwardsSPR(PhyloNode* current, PhyloNode* prev,
                          int radius, double parsimony) {
    FOR_EACH_ADJACENT_PHYLO_NODE(current, prev, it, next) {
        //path_parsimony[radius+1] has already been set "above"
        //(in the callstack), or in prepareToSearch.
        UINT*      on_path_vector  = path_parsimony[radius+1];
        PhyloNode* off_path_node   = other_adjacent_node(current, next, prev);
        UINT*      off_path_vector = current->findNeighbor(off_path_node)->get_partial_pars();
        tree.computePartialParsimonyOutOfTree
            ( on_path_vector, off_path_vector
            , path_parsimony[radius] );
        if (1<radius) {
            searchForForwardsSPR(next, current, radius-1, parsimony);
        }
        int target_branch_id = (*it)->id;
        
        double pruned_tree_score = tree.computePartialParsimonyOutOfTree
                                   ( current->findNeighbor(next)->get_partial_pars()
                                   , path_parsimony[radius], path_parsimony[0] );
        int new_branch_cost = 0;
        tree.computeParsimonyOutOfTree
            ( path_parsimony[0]
            , source.first->findNeighbor(source.second)->get_partial_pars()
            , &new_branch_cost );
        auto   subtree_root_nei = source.first->findNeighbor(source.second);
        double subtree_cost = tree.getSubTreeParsimony(subtree_root_nei);
        double benefit = parsimony       - pruned_tree_score
                       - new_branch_cost - subtree_cost;
        
        if (put_answer_here.benefit<benefit) {
            put_answer_here.benefit          = benefit;
            put_answer_here.target_branch_id = target_branch_id;
            put_answer_here.isForward        = true;
        }
        ++put_answer_here.positions_considered;
    }
}

void ParsimonySPRMove::ProperSPRSearch::searchForBackwardsSPR
    ( PhyloNode* current, PhyloNode* prev, int radius, double parsimony) {
    FOR_EACH_ADJACENT_PHYLO_NODE(current, prev, it, next) {
        UINT*      on_path_vector  = path_parsimony[radius+1];
        PhyloNode* off_path_node   = other_adjacent_node(current, next, prev);
        UINT*      off_path_vector = current->findNeighbor(off_path_node)->get_partial_pars();
        tree.computePartialParsimonyOutOfTree
            ( on_path_vector, off_path_vector
            , path_parsimony[radius] );
        if (1<radius) {
            searchForBackwardsSPR(next, current, radius-1, parsimony);
        }
        int target_branch_id = (*it)->id;
        
        double pruned_tree_score = tree.computePartialParsimonyOutOfTree
                                   ( current->findNeighbor(next)->get_partial_pars()
                                   , path_parsimony[radius], path_parsimony[0] );
        int new_branch_cost = 0;
        tree.computeParsimonyOutOfTree
            ( path_parsimony[0]
            , source.second->findNeighbor(source.first)->get_partial_pars()
            , &new_branch_cost );
        auto   subtree_root_nei = source.second->findNeighbor(source.first);
        double subtree_cost = tree.getSubTreeParsimony(subtree_root_nei);
        double benefit = parsimony       - pruned_tree_score
                       - new_branch_cost - subtree_cost;
        
        if (put_answer_here.benefit<benefit) {
            put_answer_here.benefit          = benefit;
            put_answer_here.target_branch_id = target_branch_id;
            put_answer_here.isForward        = false;
        }
        ++put_answer_here.positions_considered;
    }
}

void ParsimonySPRMove::findForwardSPR(const PhyloTree& tree, const TargetBranchRange& branches,
                                int radius, double disconnection_benefit,
                                std::vector<UINT*> &path_parsimony,
                                double parsimony_score) {
    if (lazy) {
        super::findForwardLazySPR(tree, branches, radius, disconnection_benefit);
        return;
    }
    ProperSPRSearch s(tree, branches, disconnection_benefit,
                      path_parsimony, *this);
    const TargetBranch& tb = branches[source_branch_id];
    PhyloNode* left;
    PhyloNode* right;
    GET_OTHER_ADJACENT_PHYLO_NODES(tb.first, tb.second,
                                   left, right);
    s.prepareToSearch     (left,  right,    tb.first, radius);
    s.searchForForwardsSPR(left,  tb.first, radius,   parsimony_score);
    s.prepareToSearch     (right, left,     tb.first, radius);
    s.searchForForwardsSPR(right, tb.first, radius,   parsimony_score);
}

void ParsimonySPRMove::findBackwardSPR(const PhyloTree& tree, const TargetBranchRange& branches,
                     int radius, double disconnection_benefit,
                     std::vector<UINT*> &path_parsimony,
                     double parsimony_score) {
    if (lazy) {
        super::findBackwardLazySPR(tree, branches, radius, disconnection_benefit);
        return;
    }
    ProperSPRSearch s(tree, branches, disconnection_benefit,
                      path_parsimony, *this);
    const TargetBranch& tb = branches[source_branch_id];
    PhyloNode* left;
    PhyloNode* right;
    GET_OTHER_ADJACENT_PHYLO_NODES(tb.second, tb.first,
                                   left, right);
    s.prepareToSearch      (left,  right,     tb.second, radius);
    s.searchForBackwardsSPR(left,  tb.second, radius,    parsimony_score);
    s.prepareToSearch      (right, left,      tb.second, radius);
    s.searchForBackwardsSPR(right, tb.second, radius,    parsimony_score);
}

void ParsimonySPRMove::findMove(const PhyloTree& tree,
                                const TargetBranchRange& branches,
                                int radius,
                                std::vector<UINT*> &path_parsimony,
                                double parsimony_score) {
    auto source = branches[source_branch_id];
    BenefitPair   benefit = source.getPartialDisconnectionBenefit(tree, branches);

    if (source.first->isInterior()) {
        findForwardSPR(tree, branches, radius, benefit.forwardBenefit,
                       path_parsimony, parsimony_score);
    }
    if (source.second->isInterior()) {
        findBackwardSPR(tree, branches, radius, benefit.backwardBenefit,
                        path_parsimony, parsimony_score);
    }
}

int PhyloTree::doParsimonySPR(VerboseMode how_loud) {
    return doParsimonySPR(params->parsimony_spr_iterations,
                          params->use_lazy_parsimony_spr,
                          params->spr_radius, verbose_mode < how_loud );
}


int PhyloTree::doParsimonySPR(intptr_t iterations, bool lazy,
                               int radius, bool quiet) {

    ParsimonySearchParameters s("SPR");
        
    s.iterations                 = iterations;
    s.lazy_mode                  = lazy;
    s.radius                     = radius;
    s.calculate_connection_costs = lazy;
    s.be_quiet                   = quiet;

    return doParsimonySearch<ParsimonySPRMove>(s);
}

void IQTree::doPLLParsimonySPR(VerboseMode how_loud) {
    double    init_start = getRealTime();
    StrVector oldNames;
    bool      areNamesDummied = false;

    deleteAllPartialParsimony(); 
    if (!isInitializedPLL()) {
        //Names containing '/' characters give PLL trouble.
        //Simplest way to avoid the risk... is to dummy all
        //the names for the (PLL) duration.
        areNamesDummied = true;
        oldNames        = aln->getSeqNames();
        PhyloNodeVector taxa ( getTaxaNodesInIDOrder() );
        ASSERT(taxa.size() == aln->getSeqNames().size());
        intptr_t taxa_count = taxa.size();
        for (intptr_t i=0; i<taxa_count; ++i) {
            //Note: the... if (areNamesDummied)... block below
            //depends on the exact format used here.
            std::stringstream dummy;
            dummy << "S" << i;
            std::string dummyName = dummy.str();
            taxa[i]->name = dummyName;
            aln->setSeqName(static_cast<int>(i), dummyName);
        }
        initializePLL(*params);
    }

    stringstream tree_stream;
    setRootNode(params->root);
    printTree(tree_stream, WT_SORT_TAXA);
    string constructedTreeString = tree_stream.str();;
    
    pllReadNewick(constructedTreeString);
    double opt_start = getRealTime();
    int iterations = pllOptimizeWithParsimonySPR(pllInst, pllPartitions,
        params->parsimony_spr_iterations,  
        params->spr_radius);
    double read_start = getRealTime();
    pllTreeToNewick(pllInst->tree_string, pllInst, pllPartitions,
        pllInst->start->back, PLL_FALSE, PLL_TRUE, PLL_FALSE,
        PLL_FALSE, PLL_FALSE, PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);
    string pllTreeString = string(pllInst->tree_string);

    if (areNamesDummied) {
        //Destroying PLL instance NOW avoids PLL being told about
        //the updated version of the tree when we read it back(!).
        //That would be pointless...
        pllDestroyInstance(pllInst);
        pllInst = nullptr;
    }

    PhyloTree::readTreeString(pllTreeString, areNamesDummied);
    
    if (areNamesDummied) {
        PhyloNodeVector taxa;
        getTaxa(taxa);
        intptr_t seq_name_count = aln->getSeqNames().size();
        intptr_t taxa_count     = taxa.size();
        ASSERT(taxa_count == seq_name_count);
        for (intptr_t i=0; i<taxa_count; ++i) {
            aln->setSeqName(static_cast<int>(i), oldNames[i]);
        }
        for (auto taxon_node: taxa) {
            taxon_node->id = atoi(taxon_node->name.c_str()+1); 
            //Names are of the form Sn, where n is decimal number,
            //so need to skip over the S at the front of the name.
            //See the... if (!isInitializedPLL())... block above.
            taxon_node->name = oldNames[taxon_node->id];
        }
    }

    double pars_start = getRealTime();
    initializeAllPartialPars();
    double optimized_parsimony = computeParsimony("Computing post optimization parsimony", true);
    //Note: bidirectional=true, because setAllBranchLengthsFromParsimony() needs 
    //      every branch's parsimony views calculated in *both* directions.
    LOG_LINE(how_loud, "After " << iterations << " rounds of Parsimony SPR,"
             << " parsimony score was " << optimized_parsimony);
    double pars_time             = (getRealTime() - pars_start);
    double pll_setup_overhead    = (opt_start - init_start);
    double pll_teardown_overhead = (pars_start - read_start);
    double spr_time              = (read_start - opt_start);
    LOG_LINE(how_loud, "SPR optimization  took " << spr_time              << " wall-clock seconds");
    LOG_LINE(how_loud, "PLL set-up        took " << pll_setup_overhead    << " wall-clock seconds");
    LOG_LINE(how_loud, "PLL tear-down     took " << pll_teardown_overhead << " wall-clock seconds");   
    LOG_LINE(how_loud, "Parsimony scoring took " << pars_time             << " wall-clock seconds");
    
    double fix_start = getRealTime();
    setAllBranchLengthsFromParsimony(false, optimized_parsimony);
    double fix_time  = getRealTime()-fix_start;
    LOG_LINE(how_loud, "Setting branch lengths took " << fix_time        << " wall-clock seconds");
}
