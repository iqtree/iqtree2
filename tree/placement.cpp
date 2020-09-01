//
//  placement.cpp
//  Evolutionary placement: adding taxa to a phylogenetic tree
//  This file created by James Barbetti on 25/8/20, but:
//  1. addNewTaxaToTree was formerly in phylotree.cpp;
//  2. addTaxonML likewise (and was the work of BUI Quang Minh)
//
#include "phylotree.h"

namespace {
    enum CostFunction {
        MAXIMUM_PARSIMONY,
        MAXIMUM_LIKELIHOOD
    };
    std::string getIncrementalParameter(const char letter, const char* defaultValue) {
        const std::string& inc = Params::getInstance().incremental_method;
        std::string answer = defaultValue;
        int braceLevel = 0;
        int i;
        for (i=0; i<inc.length(); ++i) {
            if (inc[i]==letter && braceLevel==0) {
                break;
            } else if (inc[i]=='{') {
                ++braceLevel;
            } else if (inc[i]=='}') {
                --braceLevel;
            }
        }
        if (i==inc.length()) {
            return answer;  //Didn't find it
        }
        ++i;
        defaultValue = "";
        int j;
        for (j=i; j<inc.length(); ++j) {
            if (inc[j]=='+' && braceLevel==0) {
                break;
            } else if (inc[j]=='-' && braceLevel==0) {
                break;
            } else if (inc[j]=='{') {
                ++braceLevel;
            } else if (inc[j]=='}') {
                --braceLevel;
            }
        }
        answer = inc.substr(i, j-i);
        if (!answer.empty() && answer[0]=='{'
            && answer[answer.length()-1]=='}' ) {
            answer = answer.substr(1, answer.length()-2);
        }
        return answer;
    }
    size_t getIncrementalParameter(const char letter, size_t defaultValue) {
        auto s = getIncrementalParameter(letter, "");
        if (s.empty()) {
            return defaultValue;
        }
        int i = convert_int_nothrow(s.c_str(), defaultValue);
        if (i<0) {
            return defaultValue;
        }
        return static_cast<size_t>(i);
    }
    CostFunction getCostFunction() {
        auto cf = getIncrementalParameter('C', "MP");
        if (cf=="ML") {
            return MAXIMUM_LIKELIHOOD;
        }
        return MAXIMUM_PARSIMONY;
    }
    enum SearchHeuristic {
        GLOBAL_SEARCH
    };
    SearchHeuristic getSearchHeuristic() {
        return GLOBAL_SEARCH;
    }
    enum LocalCleanup {
        NO_LOCAL_CLEANUP
    };
    LocalCleanup getLocalCleanupAlgorithm() {
        auto f = getIncrementalParameter('L', "");
        return NO_LOCAL_CLEANUP;
    }
    size_t getTaxaPerBatch(size_t totalTaxa) {
        size_t taxaPerBatch = getIncrementalParameter('B', 1);
        if (taxaPerBatch==0) {
            taxaPerBatch = totalTaxa;
        }
        return taxaPerBatch;
    }
    size_t getInsertsPerBatch(size_t totalTaxa, size_t taxaPerBatch) {
        size_t insertsPerBatch = getIncrementalParameter('I', taxaPerBatch);
        if (insertsPerBatch ==0 ) {
            insertsPerBatch = taxaPerBatch;
        }
        return insertsPerBatch;
    }
    enum BatchCleanup {
        NO_BATCH_CLEANUP
    };
    BatchCleanup getBatchCleanupAlgorithm() {
        auto f = getIncrementalParameter('A', "");
        return NO_BATCH_CLEANUP;
    }
    enum GlobalCleanup {
        NO_GLOBAL_CLEANUP
    };
    GlobalCleanup getGlobalCleanupAlgorithm() {
        auto f = getIncrementalParameter('T', "");
        return NO_GLOBAL_CLEANUP;
    }
};

/****************************************************************************
 Stepwise addition (greedy) by maximum likelihood
 ****************************************************************************/

double PhyloTree::addTaxonML(PhyloNode* added_taxon,     PhyloNode *added_node,
                             PhyloNode* node,            PhyloNode* dad,
                             PhyloNode* &target_node,    PhyloNode* &target_dad,
                             double& len_to_new_taxon,
                             double& len_to_target_node, double& len_to_target_dad) {
    Neighbor *dad_nei = dad->findNeighbor(node);

    // now insert the new node in the middle of the branch node-dad
    double len = dad_nei->length;
    double halfLen = 0.5 * len;
    node->updateNeighbor(dad, added_node, halfLen);
    dad->updateNeighbor(node, added_node, halfLen);
    added_node->updateNeighbor((Node*) 1, node, halfLen);
    added_node->updateNeighbor((Node*) 2, dad, halfLen);
    // compute the likelihood

    double best_score = optimizeChildBranches((PhyloNode*) added_node);
    target_node = node;
    target_dad = dad;
    // remove the added node
    node->updateNeighbor(added_node, dad, len);
    dad->updateNeighbor(added_node, node, len);
    added_node->updateNeighbor(node, (Node*) 1, len);
    added_node->updateNeighbor(dad, (Node*) 2, len);

    // now traverse the tree downwards
    FOR_NEIGHBOR_IT(node, dad, it){
        PhyloNode *target_node2 = nullptr;
        PhyloNode *target_dad2  = nullptr;
        double len_child = 0;
        double len_node  = 0;
        double len_dad   = 0;
        double score = addTaxonML(added_taxon,  added_node,
                                  (PhyloNode*)(*it)->node, (PhyloNode*)node,
                                  target_node2, target_dad2,
                                  len_child, len_node, len_dad);
        if (score > best_score) {
            best_score = score;
            target_node = target_node2;
            target_dad  = target_dad2;
            len_to_new_taxon   = len_child;
            len_to_target_node = len_node;
            len_to_target_dad  = len_dad;
        }
    }
    return best_score;
}

void PhyloTree::addTaxonMP(PhyloNode* added_taxon, PhyloNode *added_node,
                           PhyloNode *node, PhyloNode *dad,
                           PhyloNode* &target_node, PhyloNode* &target_dad,
                           int& bestScore) {
    Neighbor *dad_nei = dad->findNeighbor(node);
    double len = dad_nei->length;
    double halfLen = 0.5 * len;
    int    score;
    {
        // insert the new node in the middle of the branch node-dad
        node->updateNeighbor(dad, added_node, halfLen);
        dad->updateNeighbor(node, added_node, halfLen);
        added_node->updateNeighbor((Node*) 1, node, halfLen);
        added_node->updateNeighbor((Node*) 2, dad, halfLen);
        
        ((PhyloNeighbor*)added_taxon->findNeighbor(added_node))->partial_lh_computed = 0;
        ((PhyloNeighbor*)added_node->findNeighbor(added_taxon))->partial_lh_computed = 0;
        ((PhyloNeighbor*)added_node->findNeighbor(node))->partial_lh_computed = 0;
        ((PhyloNeighbor*)added_node->findNeighbor(dad))->partial_lh_computed = 0;

        // score it
        score = computeParsimonyBranch((PhyloNeighbor*)added_taxon->findNeighbor(added_node),
                                       (PhyloNode*)added_taxon);

        // remove the added node
        node->updateNeighbor ( added_node, dad, len);
        dad->updateNeighbor  ( added_node, node, len);
        added_node->updateNeighbor ( node, (Node*) 1, 0);
        added_node->updateNeighbor ( dad,  (Node*) 2, 0);
    }

    if (score < bestScore) {
        best_pars_score = score; //So shortcutting will work
        target_node = node;
        target_dad  = dad;
        bestScore = score;
    }

    // now traverse the tree downwards
    FOR_NEIGHBOR_IT(node, dad, it){
        addTaxonMP( added_taxon,  added_node,
                    (PhyloNode*)((*it)->node),  node,
                    target_node, target_dad, bestScore);
    }
}

class CandidateTaxon {
public:
    PhyloTree* phylo_tree;
    int taxonId;
    std::string taxonName;
    PhyloNode* new_taxon;
    PhyloNode* added_node;
    PhyloNode* target_node;
    PhyloNode* target_dad;
    double  score;         //score (the likelihood, or minus
                           //the parsimony score)
    double  lenToNewTaxon; //(best scoring) length of the edge between
                           //new_taxon and added_node
    double  lenToDad;      //(best scoring) length of edge between
                           //target_dad and added_node
    double  lenToNode;     //(best scoring) length of edge between
                           //target_child and added_node
    bool    inserted;
    
    CandidateTaxon(const CandidateTaxon& rhs) = default; //copy everything!
    
    CandidateTaxon(PhyloTree* tree, int id, std::string name,
                   int& index_parsimony, int& index_lh,
                   uint64_t lh_block_size, uint64_t scale_block_size)
        : phylo_tree(tree), taxonId(id), taxonName(name) {
        new_taxon  = phylo_tree->newNode(taxonId, taxonName.c_str()); //leaf
        added_node = phylo_tree->newNode(); //interior

        //link them together
        added_node->addNeighbor(new_taxon, 1.0);
        new_taxon->addNeighbor(added_node, 1.0);

        //add dummy neighbors (need this, because of how addTaxonML works)
        added_node->addNeighbor((Node*) 1, 1.0);
        added_node->addNeighbor((Node*) 2, 1.0);
            
        target_node   = nullptr;
        target_dad    = nullptr;
        score         = 0;
        lenToDad      = 0;
        lenToNode     = 0;
        lenToNewTaxon = 1.0;
        inserted      = false;
            
        // assign partial_pars and partial_lh storage
        phylo_tree->initializeAllPartialLh(index_parsimony, index_lh, added_node, new_taxon);
        phylo_tree->initializeAllPartialLh(index_parsimony, index_lh, new_taxon,  added_node);
    }
    void findPlacement(CostFunction costFunction,
                       PhyloNode *search_start, PhyloNode *search_backstop) {
        switch (costFunction) {
            case MAXIMUM_LIKELIHOOD:
                score = phylo_tree->addTaxonML(new_taxon,    added_node,
                                               search_start, search_backstop,
                                               target_node,  target_dad,
                                               lenToNewTaxon, lenToNode, lenToDad);
                break;
            case MAXIMUM_PARSIMONY:
                //parsimony scores are change counts: the lower the better
                //(hence the unary minus - here).
                phylo_tree->best_pars_score = INT_MAX;
                int bestParsimonyScore = INT_MAX;
                phylo_tree->addTaxonMP(new_taxon, added_node,
                                       search_start, search_backstop,
                                       target_node,  target_dad,
                                       bestParsimonyScore);
                score         = -bestParsimonyScore;
                lenToDad      = 0.5 * target_dad->findNeighbor(target_node)->length;
                lenToNode     = lenToDad;
                break;
        }
    }
    bool canInsert() {
        //Returns true if and only if this candidate can be
        //inserted into the tree in the spot that was chosen
        //for it (it can't, if another candidate has already
        //been inserted there!).
        return target_dad!=nullptr && target_node!=nullptr
            && target_dad->isNeighbor(target_node);
    }
    void insertIntoTree() {
        target_node->updateNeighbor(target_dad,  added_node, lenToNode);
        target_dad->updateNeighbor (target_node, added_node, lenToDad);
        
        //replace dummy neighbours (also needed, due to how (both)
        //                          addTaxonML and addTaxonMP work)
        added_node->updateNeighbor((Node*) 1, target_node, lenToNode);
        added_node->updateNeighbor((Node*) 2, target_dad,  lenToDad);
        
        //Set length of link to new taxon too
        added_node->findNeighbor(new_taxon)->length = lenToNewTaxon;
        new_taxon->findNeighbor(added_node)->length = lenToNewTaxon;
        
        inserted = true;
    }
    bool operator < (const CandidateTaxon& rhs) const {
        return rhs.score < score; //sort best (highest) score to front (not lowest!)
    }
    bool operator <= (const CandidateTaxon& rhs) const {
        return rhs.score <= score; //sort best (highest) score to front (not lowest!)
    }
};

class CandidateTaxa: public std::vector<CandidateTaxon>
{
public:
    explicit CandidateTaxa(size_t reservation) {
        reserve(reservation);
    }
};

void PhyloTree::addNewTaxaToTree(const IntVector& taxaIdsToAdd,
                                 int index_parsimony, int index_lh) {
    //
    //Assumes: The tree is rooted.
    //
    Params             params          = Params::getInstance();
    SearchHeuristic    heuristic       = getSearchHeuristic();    //global? or localized?
    CostFunction       placement       = getCostFunction(); //parsimony or likelihood
    LocalCleanup       localCleanup    = getLocalCleanupAlgorithm();
    size_t             taxaPerBatch    = getTaxaPerBatch(taxaIdsToAdd.size());
                                         //Must be 1 or more
    size_t             insertsPerBatch = getInsertsPerBatch(taxaIdsToAdd.size(), taxaPerBatch);
                                         //Must be 1 or more
    BatchCleanup       batchCleanup    = getBatchCleanupAlgorithm();
    GlobalCleanup      globalCleanup   = getGlobalCleanupAlgorithm();
    double             score           = 0;
    
    size_t   nptn;
    uint64_t pars_block_size, lh_block_size, scale_block_size;
    getBlockSizes( nptn, pars_block_size, lh_block_size, scale_block_size );
    
    switch (placement) {
        case MAXIMUM_PARSIMONY:
            //Todo: still need this?  I think not: initializeAllPartialPars();
            break;
        case MAXIMUM_LIKELIHOOD:
            break;
    }
    if ( verbose_mode >= VB_MED) {
        std::cout << "Batch size is " << taxaPerBatch
            << " and the number of inserts per batch is "
            << insertsPerBatch << std::endl;
    }
    size_t newTaxaCount = taxaIdsToAdd.size();
    CandidateTaxa candidates(newTaxaCount);
    for (size_t i=0; i<newTaxaCount; ++i) {
        int         taxonId   = taxaIdsToAdd[i];
        std::string taxonName = aln->getSeqName(taxonId);
        candidates.emplace_back(this, taxonId, taxonName,
                                index_parsimony, index_lh,
                                lh_block_size, scale_block_size);
    }
    initProgress(newTaxaCount,
                 "Adding new taxa to tree",
                 "added", "taxon");
    while (0<newTaxaCount) {
        if (newTaxaCount<taxaPerBatch) {
            taxaPerBatch = newTaxaCount;
        }
        size_t batchStart=0;
        for (; batchStart+taxaPerBatch <= newTaxaCount; batchStart+=taxaPerBatch) {
            switch (placement) {
                case MAXIMUM_PARSIMONY:
                    clearAllPartialLH();
                    score = - computeParsimony();
                    if (true) { //if ( verbose_mode >= VB_DEBUG) {
                        hideProgress();
                        std::cout << "Parsimony score is currently "
                            << -score << std::endl;
                        showProgress();
                    }
                    break;
                case MAXIMUM_LIKELIHOOD:
                    score = computeLikelihood();
                    if (true) { //if ( verbose_mode >= VB_DEBUG) {
                        hideProgress();
                        std::cout << "Log-likelihood is currently "
                            << score << std::endl;
                        showProgress();
                    }
                    break;
            }
            size_t batchStop = batchStart + taxaPerBatch;
            //Todo: see note 3. Search heuristics
            //      (for restricting the subtree to search)
            PhyloNode* search_start;
            PhyloNode* search_backstop;
            switch (heuristic) {
                case GLOBAL_SEARCH:
                    search_start    = (PhyloNode*) root->neighbors[0]->node;
                    search_backstop = (PhyloNode*) root;
                    break;
            }
            for (size_t i=batchStart; i<batchStop; ++i) {
                CandidateTaxon& c = candidates[i];
                c.findPlacement(placement, search_start, search_backstop);
            }
            size_t insertStop = batchStart + insertsPerBatch;
            std::sort( candidates.begin() + batchStart, candidates.begin() + batchStop);
            if (batchStop <= insertStop) {
                insertStop = batchStop; //Want them all
            }
            size_t insertCount = 0;
            for (size_t i = batchStart; i<insertStop; ++i) {
                CandidateTaxon& c = candidates[i];
                if (c.canInsert()) {
                    //if ( verbose_mode >= VB_DEBUG) {
                        hideProgress();
                        std::cout << "Inserting " << c.taxonName
                            << " which had (tree) score " << c.score << std::endl;
                        showProgress();
                    //}
                    c.insertIntoTree();
                    switch (localCleanup) {
                        //Todo: implement local cleanup options
                        default:
                            break;
                    }
                    trackProgress(1.0);
                    ++insertCount;
                }
            }
            if ( 1 < batchStop - batchStart && verbose_mode >= VB_MED ) {
                hideProgress();
                std::cout << "Inserted " << (insertCount)
                    << " out of a batch of " << (batchStop - batchStart) << "." << std::endl;
                showProgress();
            }
            switch (batchCleanup) {
                //Todo: implement per-batch cleanup options
                default:
                    break;
            }
        } //batches of items
        
        //Remove all the candidates that we were able to place
        std::vector<CandidateTaxon> oldCandidates;
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
                candidates.emplace_back(oldCandidates[r]);
            }
        }
        newTaxaCount = candidates.size();
        /*
        hideProgress();
        std::cout << "There are now " << newTaxaCount << " taxa left to insert." << std::endl;
        showProgress();
        */
    }
    doneProgress();
    if ( verbose_mode >= VB_MED ) {
        std::cout << "Tidying up tree after inserting taxa." << std::endl;
    }
    switch (globalCleanup) {
        //Todo: implement global cleanup options
        default:
            break;
    }
    initializeTree();
    deleteAllPartialLh();
    initializeAllPartialLh();
    optimizeAllBranches();
}
