#include "phylotree.h"
#include <utils/gzstream.h>
#include <model/modeldivergent.h>

void PhyloTree::setUpSubtreesForDivergentModels(ModelSubst* top_model) {

    #if (0)
        std::cout << "Tree (in setUpSubtreesForDivergentModels):" << std::endl;
        PhyloNodeVector all_nodes (getAllNodesInTree());
        for (PhyloNode* node : all_nodes ) {
            std::cout << "  " << node->id 
                    << " (" << node->name << ")" << std::endl;
        }
    #endif

    DivergentModels div_models;
    top_model->getDivergentModels(div_models);
    if (div_models.empty()) {
        //No divergent models?  Nothing to do!
        return;
    }
    std::string div_file_path = params->divergence_graph_file_path;
    if (div_file_path.empty()) {
        outError("Cannot use divergent models without a divergence graph");
    }

    //Load up divergence graph
    MTree div_tree;
    bool  is_div_tree_rooted = false;
    std::cout << "Loading divergence graph"
              << " from file: " << div_file_path 
              << "." << std::endl;
    div_tree.readTree(div_file_path.c_str(), is_div_tree_rooted);
    Node* div_root = div_tree.getRoot();

    //Set up a mapping from taxon name to taxon id.
    NameToIDMap name_to_id;
    intptr_t taxa_in_alignment = aln->getMapFromNameToID(name_to_id);
    intptr_t taxa_in_div_graph = div_tree.setNodeIdsFromMap(name_to_id);
    if (taxa_in_alignment!=taxa_in_div_graph) {
        std::stringstream complaint;
        complaint << "There were " << taxa_in_alignment << " taxa"
                  << " in the alignment, but only " 
                  << taxa_in_div_graph
                  << " of those taxa were mentioned,"
                  << " in the divergence graph.";
        outWarning(complaint.str());
    }

    //For each divergent model, identify the way that it
    //partitions the taxa in the tree
    std::vector<IntVector> subsets_for_all_models;
    ASSERT(div_models.size()==1);
    for (ModelDivergent* div_model: div_models) {
        std::vector<IntVector> subsets_for_model;        
        div_model->identifyTaxonSubsets(div_root,
                                        subsets_for_model);
        //Todo: splitSubsets(subsets_for_model, subsets_for_all_models);
        std::swap(subsets_for_all_models, subsets_for_model);
    }

    //Set up an array that maps taxon id to subset,
    //and an array that maps taxon id to taxon node
    //in the tree
    PhyloNodeVector taxa(getTaxaNodesInIDOrder());
    IntVector       taxon_to_subset(aln->getNSeq());
    int subset_number=0;
    for (IntVector& subset : subsets_for_all_models) {
        for (int taxon_id : subset) {
            taxon_to_subset[taxon_id] = subset_number;
            if (taxa.size()<=taxon_id) {
                continue;
            }
            if (taxa[taxon_id]==nullptr) {
                //If taxon is not currently in the tree, never
                //mind.  Perhaps it will be added later.
                continue;
            }
            taxa[taxon_id]->setSubsetNumber(subset_number);
            aln->setSequenceSubset(taxon_id, subset_number);
        }
        ++subset_number;
    }

    setSubsetNumbersForAllNodes();

    hideProgress();
    for (ModelDivergent* div_model: div_models) {
        div_model->mapTaxonSubsetsToModels
            (div_root, subset_number, taxon_to_subset);
        div_model->calculateSubtreeFrequencyEstimates
            (aln, this);
    }

    showProgress();
}

PhyloNodeVector PhyloTree::getAllNodesInTree() const {
    auto root = getRoot();
    PhyloNodeVector answer;
    getAllNodesInSubtree(root, nullptr, answer);
    return answer;
}

void PhyloTree::setSubsetNumbersForLeafNodes() {
    for (PhyloNode* leaf: getTaxaNodesInIDOrder()) {
        int subset_number = aln->getSequenceSubset(leaf->id);
        leaf->setSubsetNumber(subset_number);
    }
}

void PhyloTree::setSubsetNumbersForAllNodes() {
    setSubsetNumbersForLeafNodes();
    PhyloNodeVector all_nodes (getAllNodesInTree());
    PhyloNodeVector layer; //starts out with leaf nodes
    for (PhyloNode* visit: all_nodes) {
        //Mark all interior node subsets as uncalculated
        //Add all exterior nodes to the across-the-breadth 
        //vector of nodes to look at next.
        if (!visit->isLeaf()) {
            visit->setSubsetNumber(SUBSET_UNKNOWN);
        } else {     
            LOG_LINE(VerboseMode::VB_MAX, 
                     "Leaf " << visit->id  
                     << " (" << visit->name << ")"
                     << " in subset " << visit->getSubsetNumber());
            layer.push_back(visit);
        }
    }
    do {
        PhyloNodeVector next_layer; //interior nodes with subset number
                                    //set, in this do-loop iteration.
        for (PhyloNode* visited : layer ) {
            int subset = visited->getSubsetNumber();
            FOR_EACH_ADJACENT_PHYLO_NODE
                (visited, nullptr, it, interior) {
                if (interior->getSubsetNumber()==SUBSET_UNKNOWN) {
                    FOR_EACH_ADJACENT_PHYLO_NODE
                        (interior, visited, it2, next_door ) {
                        if (next_door->getSubsetNumber() == subset) {
                            //A second node, adjacent to interior, is
                            //in the same subset that visited was.
                            interior->setSubsetNumber(subset);
                            next_layer.push_back(interior);
                            LOG_LINE(VerboseMode::VB_MAX,
                                      "Interior " << interior->id 
                                      << " (" << interior->name << ")"
                                      << " in subset "
                                      << interior->getSubsetNumber()
                                      << " (because adjacent node with ids " 
                                      << visited->id 
                                      << " and "
                                      << next_door->id
                                      << " were)." );
                            break;
                        }
                    }
                }
            }
        }
        std::swap(layer, next_layer);
        //Next layer of interior nodes (that just had their
        //subset numbers calculated)... will be the layer of
        //most recently visited nodes, in the next iteration 
        //of the do-loop.
    }
    while (!layer.empty());

    for (PhyloNode* visit: all_nodes) {
        int subset_number = visit->getSubsetNumber();
        if (subset_number != SUBSET_UNKNOWN) {
            continue;
        }
        LOG_LINE(VerboseMode::VB_QUIET, 
                 "Node " << visit->id << "(" << visit->name << ")"
                 " could not be assigned to a subset.");
        //Dump the tree structure
        dumpSubsetStructure(getRoot(),nullptr,0);

        ASSERT(subset_number != SUBSET_UNKNOWN);
    }
}

void PhyloTree::dumpSubsetStructure
        (PhyloNode* node, PhyloNode* prev, int indent) {
    LOG_LINE(VerboseMode::VB_MIN, std::string(indent, ' ')
        << node->id << " " << node->name 
        << " (subset " << node->getSubsetNumber() << ")" );
    FOR_EACH_ADJACENT_PHYLO_NODE(node, prev, it, adj) {
        dumpSubsetStructure(adj, node, indent+1);
    }
}

int PhyloTree::getSubTreeNumberForBranch(PhyloNode* dad, 
                                         PhyloNode* node) const {
    if (!model->isDivergentModel()) {
        return 0;
    }
    ModelDivergent* model_div = dynamic_cast<ModelDivergent*>(model);
    int dad_subset            = dad->getSubsetNumber();
    ModelMarkov* dad_model    = model_div->getSubsetModel(dad_subset);
    int child_subset          = node->getSubsetNumber();
    ModelMarkov* child_model  = model_div->getSubsetModel(child_subset);
    if (dad_model == child_model) {
        return model_div->getSubtreeNumberOfSubset(child_subset);
    } else {
        return 0;
    }
}

ModelSubst* PhyloTree::getModelForBranch
                (PhyloNode* dad, PhyloNode* node,
                 ModelSubst*& other_model) const {
    if (!model->isDivergentModel()) {
        return model;
    }
    ModelDivergent* model_div = dynamic_cast<ModelDivergent*>(model);
    int dad_subset            = dad->getSubsetNumber();
    ModelMarkov* dad_model    = model_div->getSubsetModel(dad_subset);
    int child_subset          = node->getSubsetNumber();
    ModelMarkov* child_model  = model_div->getSubsetModel(child_subset);

    if (child_model == dad_model) {
        return child_model;
    }
    return model_div->getBranchJoiningModel(dad_subset, child_subset);
}

RateHeterogeneity* PhyloTree::getRateModelForBranch
    (PhyloNode* dad, PhyloNode* node, 
     RateHeterogeneity*& other_rate_model) const {
    if (!model->isDivergentModel()) {
        return site_rate;
    }
    auto model_div        = dynamic_cast<ModelDivergent*>(model);
    int  dad_subset       = dad->getSubsetNumber();
    auto rate_models      = model_div->getSubtreeRateModels();
    other_rate_model      = rate_models[dad_subset];
    other_rate_model      = (other_rate_model!=nullptr)
                          ? other_rate_model : site_rate;
    int  child_subset     = node->getSubsetNumber();
    auto child_rate_model = rate_models[child_subset];
    child_rate_model      = (child_rate_model!=nullptr)
                          ? child_rate_model : site_rate;

    return child_rate_model;
}

ModelSubst* PhyloTree::getModelForCurrentBranch() const {
    PhyloNode*  dad  = current_it_back->getNode();
    PhyloNode*  node = current_it->getNode();
    ASSERT(dad!=nullptr);
    ASSERT(node!=nullptr);
    ModelSubst* other_model = nullptr;
    return getModelForBranch(dad, node, other_model);
}

RateHeterogeneity* PhyloTree::getRateModelForCurrentBranch() const {
    PhyloNode*  dad  = current_it_back->getNode();
    PhyloNode*  node = current_it->getNode();
    ASSERT(dad!=nullptr);
    ASSERT(node!=nullptr);
    RateHeterogeneity* other_rate_model = nullptr;
    return getRateModelForBranch(dad, node, other_rate_model);
}

void PhyloTree::getModelAndTipLikelihood
        (PhyloNode*   dad, PhyloNode*   node, 
         ModelSubst*& model_to_use, ModelSubst*& other_model, 
         RateHeterogeneity*& rate_model, RateHeterogeneity*& other_rate,
         double*&     tip_lh) const {
    model_to_use = model;
    other_model  = model_to_use;
    rate_model   = site_rate;
    other_rate   = rate_model;
    tip_lh       = tip_partial_lh;
    if (model_to_use->isDivergentModel()) {
        ModelDivergent* div_model = 
            dynamic_cast<ModelDivergent*>(model_to_use);
        getSubtreeModelAndTipLikelihood(dad, node, div_model, 
                                        model_to_use, other_model,
                                        rate_model, other_rate, tip_lh);
    }
}

void PhyloTree::getSubtreeModelAndTipLikelihood
        (PhyloNode* dad, PhyloNode* node, ModelDivergent* div_model,
         ModelSubst*& model_to_use, ModelSubst*& other_model, 
         RateHeterogeneity*& rate_model, RateHeterogeneity*& other_rate,
         double*&     tip_lh) const {
    int subtree_number = getSubTreeNumberForBranch(dad, node);
    model_to_use = div_model->getNthSubtreeModel(subtree_number);
    rate_model   = div_model->getNthSubtreeRateModel(subtree_number);
    rate_model   = (rate_model!=nullptr) ? rate_model : site_rate;
    tip_lh       = tip_partial_lh 
                    + subtree_number * tip_partial_lh_size_per_model;
    int other_subtree_number = getSubTreeNumberForBranch(node, dad);
    other_model  = div_model->getNthSubtreeModel(other_subtree_number);
    other_rate   = div_model->getNthSubtreeRateModel(other_subtree_number);
    other_rate   = (other_rate!=nullptr) ? other_rate : site_rate;

    if (other_model != model_to_use) {
        LOG_LINE(VerboseMode::VB_DEBUG,
                    "Model for dad (" << dad->id << ")"
                    << " different from model for node (" << node->id << ")" );
    }
}


void PhyloTree::handleDivergentModelBoundary
        (TraversalInfo&     info, 
         ModelSubst*        model_to_use, ModelSubst*        other_model, 
         RateHeterogeneity* rate_model,   RateHeterogeneity* other_rate,
         intptr_t ptn_lower, intptr_t ptn_upper, intptr_t packet_id,
         int start_cat_no,   int stop_cat_no,
         const LikelihoodBufferSet& buffers) {

    //Note:
    //This assumes that organization of partial_lh 
    //(major to minor) is as follows:
    //
    //   Rank             Number of choices         Indexing Multiplier
    //-- -------------    -----------------         -------------------
    //1. Pattern-Major    (PatternNo / VectorSize)  block
    //2. Mixture Category (model->getNMixtures())   
    //3. Rate Category    (rate_model->getNRate())  v_by_s
    //4. State            (aln->num_states)         vector_size
    //5. Pattern-Minor    (PatternNo % VectorSize)  1
    //

    PhyloNeighbor* dad_branch     = info.dad_branch;
    int            nstates        = aln->num_states;
    const int      ncat           = rate_model->getNRate();
    const bool     fused          = model_factory->fused_mix_rate;
    const int      n_mix          = model->getNMixtures();
    const int      ncat_mix       = fused ? ncat : ncat * n_mix;
    const int      v_by_s         = nstates * vector_size;
    const int      block          = v_by_s  * ncat_mix;
    double*        start_lh       = dad_branch->partial_lh 
                                  + ptn_lower 
                                    * static_cast<intptr_t>(nstates)
                                    * static_cast<intptr_t>(ncat_mix);

    ASSERT((ptn_lower % vector_size) == 0);
    //Otherwise, formula for start_lh is invalid.
    ASSERT((ptn_upper % vector_size) == 0);

    LOG_LINE(VerboseMode::VB_MIN, "nstates=" << nstates
             << ", ncat=" << ncat << ", fused=" << fused
             << ", n_mix=" << n_mix << ", ncat_mix=" << ncat_mix
             << ", v_by_s=" << v_by_s << ", block=" << v_by_s * ncat_mix
             << ", ptn_lower=" << ptn_lower
             << ", ptn_upper=" << ptn_upper);

    //
    //Question: for mixture models, is information stored 
    //(in each block) in mixture category major order, 
    //or in rate category major order? The code below 
    //assumes mixture category major order (as it's simpler!)
    //I think that's right, but I'm not 100% sure.
    //
    if (!fused && 1<ncat && 1<n_mix) {
        start_cat_no *= ncat;
        stop_cat_no  *= ncat;
    }

    std::vector<double> log_lh_total(nstates, 0.0);

    //Get category proportions
    auto cat_count = stop_cat_no - start_cat_no;
    std::vector<double> cat_prop(cat_count, 1.0);
    for (auto c=start_cat_no; c<stop_cat_no; ++c) {
        int     denom          = fused ? 1 : ncat;
        int     mymix          = c/denom;
        int     mycat          = c%ncat;
        cat_prop[c - start_cat_no ] 
            = rate_model->getProp(mycat) 
            * model_to_use->getMixtureWeight(mymix);
        LOG_LINE(VerboseMode::VB_MIN, "cat_prop for cat " << c
             << " is " << cat_prop[c - start_cat_no ]);
    }

    double* partial_lh   = start_lh;
    double  total_log_lh = 0;

    for (intptr_t ptn = ptn_lower; ptn < ptn_upper
         ; partial_lh += block, ptn += vector_size) {
        int    state_offset = 0;
        double log_lh_here = 0.0; //log of likelihood for the 
                                  //sites in this pattern block
        int v_size_here = (ptn + vector_size < ptn_upper) 
                        ? vector_size : (ptn_upper-ptn);

        for (int vector_index = 0; vector_index < v_size_here; ++vector_index) { 
            double lh_for_vector_posn = 0.0; //likelihood (over all states), 
                                             //this pattern            
            for (int state = 0; state < nstates; ++state, state_offset+=vector_size) {
                double lh_for_state = 0.0; //likelihood over all rate categories
                                           //and mixtures, for this pattern
                int    cat          = start_cat_no;                
                for (int cat_offset = start_cat_no * v_by_s
                     ; cat < stop_cat_no
                     ; cat_offset += v_by_s, ++cat) {
                    auto   ix      = cat_offset + state_offset + vector_index;
                    double lh_here = abs(partial_lh [ ix ])
                                   * cat_prop   [ cat - start_cat_no ];
                    //Todo: what about scaling?!  Doesn't that need to be
                    //taken into account?  Somehow?
                    //LOG_LINE(VerboseMode::VB_MIN, "ptn+i=" << (ptn+vector_index)
                    //         << ", ix=" << ix << ", lh=" << lh_here);
                    lh_for_state += lh_here;
                }
                lh_for_vector_posn += lh_for_state;
            }
            //LOG_LINE(VerboseMode::VB_MIN, "ptn[" << (ptn+vector_index) <<"]"
            //         << " has lh " << lh_for_vector_posn);
            ASSERT(0<lh_for_vector_posn);
            log_lh_here += log(lh_for_vector_posn) * ptn_freq[ptn+vector_index];        
        }
        total_log_lh += log_lh_here;
    }

    double total_ptn_freq = 0.0;
    for (auto ptn=ptn_lower; ptn<ptn_upper; ++ptn) {
        total_ptn_freq += static_cast<double>(ptn_freq[ptn]);
    }

    double ptn_count      = ptn_upper - ptn_lower;
    double multiplier     = 1.0 / total_ptn_freq;;
    double log_lh_per_ptn = total_log_lh * multiplier;
    double lh_per_ptn     = exp(log_lh_per_ptn);    
    LOG_LINE(VerboseMode::VB_MIN, "ptn_count is " << ptn_count
             << ", total_ptn_freq is " << total_ptn_freq
             << ", multiplier is "     << multiplier
             << ", log_lh_per_ptn is " << log_lh_per_ptn 
             << ", and lh_per_ptn is " << lh_per_ptn);
    
    partial_lh = start_lh;
    int start_cat_offset = start_cat_no * v_by_s;
    int stop_cat_offset  = stop_cat_no  * v_by_s;
    LOG_LINE(VerboseMode::VB_MIN, 
             "start_cat_offset is " << start_cat_offset
             << ", stop_cat_offset is " << stop_cat_offset);
    for (intptr_t ptn = ptn_lower; ptn < ptn_upper
         ; partial_lh += block, ptn += vector_size) {
        for (int cat_offset = start_cat_offset
             ; cat_offset < stop_cat_offset
             ; cat_offset += v_by_s) {
            int state_offset = 0;
            for (int state = 0; state < nstates; ++state, state_offset+=vector_size) {
                for (int vector_index = 0; vector_index < vector_size; ++vector_index) {
                    intptr_t ix = cat_offset + state_offset + vector_index;                   
                    partial_lh [  ix ] = lh_per_ptn;
                }
            }
        }
    }
}
