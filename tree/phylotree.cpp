/***************************************************************************
 *   Copyright (C) 2009-2015 by                                            *
 *   BUI Quang Minh <minh.bui@univie.ac.at>                                *
 *   Lam-Tung Nguyen <nltung@gmail.com>                                    *
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
#include "phylotreethreadingcontext.h"
#include <utils/progress.h>  //for progress_display
#include <utils/starttree.h>
#include <utils/stringfunctions.h> //for convert_string_vec
//#include "rateheterogeneity.h"
#include "alignment/alignmentpairwise.h"
#include "alignment/alignmentsummary.h"
#include <algorithm>
#include <limits>
#include "utils/timeutil.h"
#include "utils/pllnni.h"
#include "phylosupertree.h"
#include "phylosupertreeplen.h"
#include "upperbounds.h"
#include "utils/MPIHelper.h"
#include "utils/hammingdistance.h"
#include "model/modelmixture.h"
#include "phylonodemixlen.h"
#include "phylotreemixlen.h"
#include "model/modelfactory.h" //for readModelsDefinition
#include <placement/parallelparsimonycalculator.h> //for ParallelParsimonyCalculator

#ifdef _MSC_VER
#include <boost/scoped_array.hpp>
#endif

const int LH_MIN_CONST = 1;

//const static int BINARY_SCALE = floor(log2(1/SCALING_THRESHOLD));
//const static double LOG_BINARY_SCALE = -(log(2) * BINARY_SCALE);



/****************************************************************************
 PhyloTree class
 ****************************************************************************/

PhyloTree::PhyloTree() : MTree(), CheckpointFactory() {
    init();
}

void PhyloTree::init() {
    aln = NULL;
    model = NULL;
    site_rate = NULL;
    optimize_by_newton              = true;
    central_partial_lh              = nullptr;
    lh_block_size                   = 0; //will be set, later, by determineBlockSizes()
    nni_partial_lh                  = nullptr;
    tip_partial_lh                  = nullptr;
    tip_partial_lh_computed         = 0;
    ptn_freq_computed               = false;
    central_scale_num               = nullptr;
    scale_block_size                = 0;
    central_scale_num_size_in_bytes = 0; //will be set, later, by determineBlockSizes()
    nni_scale_num                   = nullptr;
    central_partial_pars            = nullptr; //points to per-node partial parsimony vectors
                                               //followed by tip partial parsimony vectors.
    total_parsimony_mem_size        = 0;       //will be set, when central_partial_pars is allocated
    pars_block_size                 = 0;       //will be set, later, by determineBlockSizes()
    tip_partial_pars                = nullptr; //points to the last part of central_partial_pars
                                               //(and is set when central_partial_pars is).
    
    cost_matrix = NULL;
    model_factory = NULL;
    discard_saturated_site = true;
    _pattern_lh_cat_state = NULL;
    _site_lh = NULL;
    //root_state = STATE_UNKNOWN;
    root_state = 126;
    ptn_freq = NULL;
    ptn_freq_pars = NULL;
    ptn_invar = NULL;
    subTreeDistComputed = false;
    dist_matrix = nullptr;
    dist_matrix_rank = 0;
    is_dist_file_read = false;
    params = NULL;
    setLikelihoodKernel(LK_SSE2);  // FOR TUNG: you forgot to initialize this variable!
    setNumThreads(1);
    num_threads = 0;
    num_packets = 0;
    max_lh_slots = 0;
    save_all_trees = 0;
    nodeBranchDists = NULL;
    // FOR: upper bounds
    mlCheck = 0;
    skippedNNIub = 0;
    totalNNIub = 0;
    minStateFreq = 0.0;
    //minUB = 0.0;
    //meanUB = 0.0;
    //maxUB = 0.0;
    pllInst = NULL;
    pllAlignment = NULL;
    pllPartitions = NULL;
//    lhComputed = false;
    curScore = -DBL_MAX;
    root = NULL;
    current_scaling = 1.0;
    is_opt_scaling = false;
    num_partial_lh_computations = 0;
    vector_size = 0;
    safe_numeric = false;
    summary = nullptr;
    isSummaryBorrowed = false;
    progress = nullptr;
    progressStackDepth = 0;
    isShowingProgressDisabled = false;
    warnedAboutThreadCount = false;
    warnedAboutNumericalUnderflow = false;
}

PhyloTree::PhyloTree(Alignment *aln) : MTree(), CheckpointFactory() {
    init();
    this->aln = aln;
}

PhyloTree::PhyloTree(string& treeString, Alignment* aln, bool isRooted) : MTree() {
    init();
    stringstream str;
    str << treeString;
    str.seekg(0, ios::beg);
    freeNode();
    readTree(str, isRooted);
    setAlignment(aln);
}

void PhyloTree::startCheckpoint() {
    checkpoint->startStruct("PhyloTree");
}


void PhyloTree::saveCheckpoint() {
    startCheckpoint();
//    StrVector leafNames;
//    getTaxaName(leafNames);
//    CKP_VECTOR_SAVE(leafNames);
    string newick = PhyloTree::getTreeString();
    CKP_SAVE(newick);
//    CKP_SAVE(curScore);
    endCheckpoint();
    CheckpointFactory::saveCheckpoint();
}

void PhyloTree::restoreCheckpoint() {
    CheckpointFactory::restoreCheckpoint();
    startCheckpoint();
//    StrVector leafNames;
//    if (CKP_VECTOR_RESTORE(leafNames)) {
//        if (leafNames.size() +(int)rooted != leafNum)
//            outError("Alignment mismatched from checkpoint!");
//
//        StrVector taxname;
//        getTaxaName(taxname);
//        for (int i = 0; i < leafNames.size(); i++)
//            if (taxname[i] != leafNames[i])
//                outError("Sequence name " + taxname[i] + " mismatched from checkpoint");
//    }    
    string newick;
    if (CKP_RESTORE(newick)) {
        PhyloTree::readTreeString(newick);
    }
    endCheckpoint();
}

void PhyloTree::discardSaturatedSite(bool val) {
    discard_saturated_site = val;
}

void myPartitionsDestroy(partitionList *pl) {
    int i;
    for (i = 0; i < pl->numberOfPartitions; i++) {
        rax_free(pl->partitionData[i]->partitionName);
        rax_free(pl->partitionData[i]);
    }
    rax_free(pl->partitionData);
    rax_free(pl);
}

PhyloTree::~PhyloTree() {
    doneComputingDistances();
    aligned_free(nni_scale_num);
    aligned_free(nni_partial_lh);
    aligned_free(central_partial_lh);
    aligned_free(central_scale_num);
    aligned_free(central_partial_pars);
    aligned_free(cost_matrix);

    delete model_factory;
    model_factory = NULL;
    delete model;
    model = NULL;
    delete site_rate;
    site_rate = NULL;
    aligned_free(_site_lh);
    aligned_free(ptn_freq);
    aligned_free(ptn_freq_pars);
    ptn_freq_computed = false;
    aligned_free(ptn_invar);
    delete[] dist_matrix;
    dist_matrix = nullptr;
    dist_matrix_rank = 0;

    if (pllPartitions!=nullptr) {
        myPartitionsDestroy(pllPartitions);
        pllPartitions = nullptr;
    }
    if (pllAlignment!=nullptr) {
        pllAlignmentDataDestroy(pllAlignment);
        pllAlignment = nullptr;
    }
    if (pllInst) {
        pllDestroyInstance(pllInst);
        pllInst = nullptr;
    }
    if (!isSummaryBorrowed) {
        delete summary;
    }
    summary = nullptr;
    delete progress;
    progress = nullptr;
    progressStackDepth = 0;
}

void PhyloTree::readTree(const char *infile, bool &is_rooted) {
    MTree::readTree(infile, is_rooted);
    // 2015-10-14: has to reset this pointer when read in
    current_it = current_it_back = NULL;
    if (rooted && root)
    {
        computeBranchDirection();
    }
}

void PhyloTree::readTree(istream &in, bool &is_rooted) {
    MTree::readTree(in, is_rooted);
    // 2015-10-14: has to reset this pointer when read in
    current_it = current_it_back = nullptr;
    // remove taxa if necessary
    if (removed_seqs.size() > 0) {
        removeTaxa(removed_seqs, true, "");
    }
    // collapse any internal node of degree 2
    NodeVector nodes;
    getInternalNodes(nodes);
    int num_collapsed = 0;
    for (NodeVector::iterator it = nodes.begin(); it != nodes.end(); it++) {
        if ((*it)->degree() == 2) {
            Node*  left  = (*it)->neighbors[0]->node;
            Node*  right = (*it)->neighbors[1]->node;
            double len   = (*it)->neighbors[0]->length+(*it)->neighbors[1]->length;
            left->updateNeighbor((*it), right, len);
            right->updateNeighbor((*it), left, len);
            delete (*it);
            ++num_collapsed;
            if (verbose_mode >= VerboseMode::VB_MED) {
                cout << "Node of degree 2 collapsed" << endl;
            }
        }
    }
    if (num_collapsed) {
        initializeTree();
    }
    if (rooted) {
        computeBranchDirection();
    }
}

void PhyloTree::assignLeafNames(Node *node, Node *dad) {
    if (!node)
        node = root;
    if (node->isLeaf()) {
        if (rooted && node == root) {
            ASSERT(node->id == leafNum-1);
            root->name = ROOT_NAME;
        } else {
            node->id = atoi(node->name.c_str());
            node->name = aln->getSeqName(node->id);
        }
        ASSERT(node->id >= 0 && node->id < static_cast<int>(leafNum));
    }
    FOR_NEIGHBOR_IT(node, dad, it)assignLeafNames((*it)->node, node);
}

void PhyloTree::copyTree(MTree *tree) {
    MTree::copyTree(tree);
    if (aln==nullptr) {
        return;
    }
    // reset the ID with alignment
    setAlignment(aln);
}

void PhyloTree::copyTree(MTree *tree, string &taxa_set) {
    MTree::copyTree(tree, taxa_set);
    if (rooted)
        computeBranchDirection();
    if (aln==nullptr) {
        return;
    }
    // reset the ID with alignment
    setAlignment(aln);
}

void PhyloTree::copyPhyloTree(PhyloTree *tree, bool borrowSummary) {
    MTree::copyTree(tree);
    if (!tree->aln)
        return;
    setAlignment(tree->aln);
    if (borrowSummary && summary!=tree->summary && tree->summary!=nullptr) {
        if (!isSummaryBorrowed) {
            delete summary;
        }
        summary           = tree->summary;
        isSummaryBorrowed = (summary!=nullptr);
    }
}

void PhyloTree::copyPhyloTreeMixlen(PhyloTree *tree, int mix, bool borrowSummary) {
    if (tree->isMixlen()) {
        ((PhyloTreeMixlen*)tree)->cur_mixture = mix;
    }
    copyPhyloTree(tree, borrowSummary);
    if (tree->isMixlen()) {
        ((PhyloTreeMixlen*)tree)->cur_mixture = -1;
    }
}

#define FAST_NAME_CHECK 1
void PhyloTree::setAlignment(Alignment *alignment) {
    aln = alignment;
    //double checkStart = getRealTime();
    int nseq = aln->getNSeq32();
    bool err = false;
#if FAST_NAME_CHECK
    map<string, Node*> mapNameToNode;
    getMapOfTaxonNameToNode(nullptr, nullptr, mapNameToNode);
    for (int seq = 0; seq < nseq; seq++) {
        string seq_name = aln->getSeqName(seq);
        auto it = mapNameToNode.find(seq_name);
        if (it==mapNameToNode.end()) {
            string str = "Alignment sequence ";
            str += seq_name;
            str += " does not appear in the tree";
            err = true;
            outError(str, false);
        } else {
            (*it).second->id = seq;
            mapNameToNode.erase(it);
        }
    }
#else
    //For each sequence, find the correponding node, via
    //the map (if there is one), tag it, and remove it from
    //the map.
    for (size_t seq = 0; seq < nseq; seq++) {
        string seq_name = aln->getSeqName(seq);
        Node *node = findLeafName(seq_name);
        if (!node) {
            string str = "Alignment sequence ";
            str += seq_name;
            str += " does not appear in the tree";
            err = true;
            outError(str, false);
        } else {
            ASSERT(node->isLeaf());
            node->id = seq;
        }
    }
#endif
    if (err) {
        printTree(cout, WT_NEWLINE);
        outError("Tree taxa and alignment sequence do not match (see above)");
    }
#if FAST_NAME_CHECK
    //Since every sequence in the alignment has had its
    //corresponding node tagged with an id, and that node's
    //entry has been removed from the map, it follows:
    //any node still referenced from the map does NOT
    //correspond to a sequence in the alignment.
    for (auto it = mapNameToNode.begin(); it != mapNameToNode.end(); ++it) {
        if (it->first != ROOT_NAME) {
            outError((string)"Tree taxon " + it->first + " does not appear in the alignment", false);
            err = true;
        }
    }
#else
    StrVector taxname;
    getTaxaName(taxname);
    for (StrVector::iterator it = taxname.begin(); it != taxname.end(); it++)
        if ((*it) != ROOT_NAME && alignment->getSeqID(*it) < 0) {
            outError((string)"Tree taxon " + (*it) + " does not appear in the alignment", false);
            err = true;
        }
#endif
    if (err) {
        outError("Tree taxa and alignment sequence do not match (see above)");
    }
    if (rooted) {
        ASSERT(root->name == ROOT_NAME);
        root->id = aln->getNSeq32();
    }
    /*
    if (verbose_mode >= VerboseMode::VB_MED) {
        cout << "Alignment namecheck took " << (getRealTime()-checkStart)
        << " sec (of wall-clock time)" << endl;
    }
    */
}

void PhyloTree::configureLikelihoodKernel(const Params& params, bool force) {
    if (computeLikelihoodBranchPointer==nullptr || force)
    {
        setLikelihoodKernel(params.SSE);
        optimize_by_newton = params.optimize_by_newton;
        setNumThreads(params.num_threads);
    } else if ( num_threads <=0 ) {
        setNumThreads(params.num_threads);
    }
}

void PhyloTree::configureModel(Params& params, PhyloTree* report_to_tree) {
    if (model==nullptr) {
        ModelsBlock *models_block = readModelsDefinition(params);
        if (model_factory==nullptr && aln!=nullptr) {
            initializeModel(params, aln->model_name,
                            models_block, report_to_tree);
        }
        if (getRate()->isHeterotachy() && !isMixlen()) {
            ASSERT(0 && "Heterotachy tree not properly created");
        }
        delete models_block;
    }
}

void PhyloTree::initializeModel(Params &params, string model_name,
                                ModelsBlock *models_block,
                                PhyloTree* report_to_tree) {
    //Must be overridden in IQTree etc.
    throw "Uh-oh!  Called initializeModel"
          " on PhyloTree rather than a subclass";
}

void PhyloTree::prepareForDeletes() {
    LOG_LINE ( VerboseMode::VB_DEBUG, 
               "Before delete, number of leaves " << this->leafNum
               << ", of nodes "       << this->nodeNum
               << ", of branches "    << this->branchNum);
    deleteAllPartialLh(); //Todo: Only clear partial Lh
    clearAllPartialParsimony(false);
    if (shouldPlacementUseSankoffParsimony()) {
        computeTipPartialParsimony();
    }
    int initial_parsimony = computeParsimony("Computing initial parsimony");
    LOG_LINE ( VerboseMode::VB_MED, 
               "Parsimony score before deletions was " << initial_parsimony );
    if (shouldPlacementUseLikelihood()) {
        initializeAllPartialLh();
        fixNegativeBranch();
        double likelihoodStart    = getRealTime();
        double likelihoodCPUStart = getCPUTime();
        double likelihood         = computeLikelihood();
        LOG_LINE ( VerboseMode::VB_MIN, 
                   "Computed initial likelihood in " << (getRealTime() - likelihoodStart) << " wall-clock secs"
                  << " and " << (getCPUTime() - likelihoodCPUStart) << " cpu secs");
        LOG_LINE ( VerboseMode::VB_MIN, "Likelihood score before deletions was " << likelihood );
    }
}

void PhyloTree::doneDeletes(size_t countRemoved, bool are_placements_to_follow) {
    //Todo: Carry out any requested "after-batch-of-deletes" global tidy-up
    //Todo: Carry out any requested "after-all-deletes" global tidy-up of the tree.
    leafNum   -= static_cast<unsigned int>(countRemoved);
    branchNum -= static_cast<int>(countRemoved*2);
    nodeNum   -= static_cast<int>(countRemoved*2);
    LOG_LINE ( VerboseMode::VB_DEBUG, "After delete, number of leaves " << this->leafNum
                    << ", of nodes "       << this->nodeNum
                    << ", of branches "    << this->branchNum);

    deleteAllPartialLhAndParsimony();
    if (!are_placements_to_follow) {
        if (isUsingSankoffParsimony() && !params->sankoff_cost_file) {
            stopUsingSankoffParsimony();
        }
        std::stringstream sx;
        sx << "Computing parsimony (after deleting taxa)";
        int parsimony_score = computeParsimony(sx.str().c_str());
        LOG_LINE(VerboseMode::VB_MIN, 
                 "Parsimony score after deleting taxa was " << parsimony_score);
        fixNegativeBranch();
        double likelihood = computeLikelihood();
        LOG_LINE(VerboseMode::VB_MIN, 
                 "Likelihood score after deleting taxa was " << likelihood);
        deleteAllPartialLhAndParsimony();
    }
    else if (VerboseMode::VB_MED <= verbose_mode ) {
        int post_delete_parsimony = computeParsimony("Computing post-delete parsimony");
        LOG_LINE (VerboseMode::VB_MIN, 
                  "Parsimony score after deletions was " << post_delete_parsimony );
    }

}

void PhyloTree::prepareForPlacement() {
    bool using_sankoff = shouldPlacementUseSankoffParsimony();
    if (using_sankoff) {
        initCostMatrix(CM_UNIFORM);
        //Force the use of Sankoff parsimony kernel
    }
    if (aln!=nullptr) {
        if ( aln->ordered_pattern.empty() ) {
            aln->orderPatternByNumChars(PAT_VARIANT);
        }
    }
    configureLikelihoodKernel(*params, true);
    configureModel(*params, this);
    setParsimonyKernel(params->SSE);
}

void PhyloTree::logTaxaToBeRemoved(const map<string, Node*>& mapNameToNode) {
    for (auto it=mapNameToNode.begin(); it!=mapNameToNode.end(); ++it) {
        PhyloNode*     leaf      = (PhyloNode*)(it->second);
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
                rightLength = rightLink->length;
            }
        }
        LOG_LINE (VerboseMode::VB_MED, "Before deletion, " << it->second->name
                  << " had branch length " << upLink->length
                  << " (left branch " << leftLength
                  << ", and right branch " << rightLength << ")");
    }
}

bool PhyloTree::updateToMatchAlignment(Alignment* alignment) {
    PhyloTreeThreadingContext context(*this, false);
    aln      = alignment;
    int nseq = aln->getNSeq32();

    removeSampleTaxaIfRequested();

    //For sequence names that are found in both the
    //tree and the alignment, set the node id.
    //(Also: identify sequences NOT found in the tree)
    map<string, Node*> mapNameToNode;
    double mapStart = getRealTime();
    getMapOfTaxonNameToNode(nullptr, nullptr, mapNameToNode);
    LOG_LINE (VerboseMode::VB_DEBUG, "Mapping taxa names to nodes took "
              << (getRealTime()-mapStart) << " wall-clock secs" );
    if (VerboseMode::VB_MAX <= verbose_mode) {
        for (auto it=mapNameToNode.begin(); it!=mapNameToNode.end(); ++it) {
            LOG_LINE(VerboseMode::VB_MAX, "sequence " << it->first << 
                     " has id " << it->second->id );
        }
    }

    IntVector taxaIdsToAdd; //not found in tree
    for (int seq = 0; seq < nseq; seq++) {
        string seq_name = aln->getSeqName(seq);
        auto   it       = mapNameToNode.find(seq_name);
        if (it==mapNameToNode.end()) {
            if ( !params->suppress_list_of_sequences
                 || VerboseMode::VB_MAX<=verbose_mode ) {
            LOG_LINE(VerboseMode::VB_MED, "Could not find sequence " << seq
                     << " " << seq_name << " in the tree."
                     << " It will be added.");
            }
            taxaIdsToAdd.push_back((int)seq);
        } else {
            (*it).second->id = seq;
            mapNameToNode.erase(it);
        }
    }
    auto rootInMap = mapNameToNode.find(ROOT_NAME);
    if (rootInMap != mapNameToNode.end()) {
        root     = (*rootInMap).second;
        root->id = nseq;
        mapNameToNode.erase(rootInMap);
    } else if (root==nullptr && !taxaIdsToAdd.empty()) {
        //Todo: make sure we HAVE a root, since later steps will need one.
        throw "Cannot add new taxa to unrooted tree";
    }
    bool will_delete = !mapNameToNode.empty();
    bool will_add    = !taxaIdsToAdd.empty();
    bool modified    = will_delete || will_add;
    
    if (!modified) {
        return false;
    }
    prepareForPlacement();
    
    if (will_delete) {
        //Necessary, if we read in a tree that came out of a
        //distance matrix algorithm (which might have
        //negative branch lengths (which would bust the stuff below).
        //Todo: figure out what happens for PhyloSuperTree here
        //Todo: What if -fixbr is supplied.  Isn't that problematic?
        
        //mapNameToNode now lists leaf nodes to be removed
        //from the tree.  Remove them.
        StrVector taxaToRemove;
        for (auto it=mapNameToNode.begin(); it!=mapNameToNode.end(); ++it) {
            taxaToRemove.push_back(it->first);
        }
        if ( !params->suppress_list_of_sequences ||
             VerboseMode::VB_MAX <= verbose_mode ) {
            logTaxaToBeRemoved(mapNameToNode);
        }
        mapNameToNode.clear();
        
        prepareForDeletes();
        auto countRemoved = removeTaxa(taxaToRemove, !will_add,
                                       "Removing deleted taxa from tree");
        //Todo: Carry out any requested "after-each-delete" local tidy-up of the tree
        //      (via an extra parameter to removeTaxa, perhaps?)
        doneDeletes(countRemoved, will_add);
    }
    if (will_add) {
        addNewTaxaToTree(taxaIdsToAdd, "Adding new taxa to tree", "");
        //Note: If the Sankoff cost matrix isn't wanted after placement, it
        //      will be deallocated (switching back to regular parsimony)
        //      in GlobalPlacementOptimizer::optimizeAfterPlacement.
    }
    return true;
}

namespace {
    //Mapping names to ids
    class NameToIDMap: public std::map<std::string, int> {
    public:
        explicit NameToIDMap(const vector<std::string>& src) {
            int count = static_cast<int>(src.size());
            for (int id = 0 ; id < count; ++id) {
                operator[](src[id]) = id;
            }
        }
    };

    typedef std::set<std::string>          NameSet;
    typedef std::pair<size_t, NameSet>     CountAndNames;
    typedef std::pair<size_t, std::string> HashedSequence;

    class UniqueSequenceMap: public std::map<HashedSequence,CountAndNames> {
        public:
        size_t trackHashAndSequence(size_t hash,
                                    std::string sequence,
                                    std::string name) {
            HashedSequence hashed(hash, sequence);
            auto it = find(hashed);
            if (it==end()) {
                NameSet name_set;
                name_set.insert(name);
                (*this)[hashed] = CountAndNames(1, name_set);
                return 1;
            } else {
                it->second.second.insert(name); //record name
                return ++it->second.first;   //bump up count
            }
        }
        size_t untrackHashAndSequence(size_t      hash,
                                      std::string sequence,
                                      std::string name) {
            HashedSequence hashed(hash, sequence);
            auto it = find(hashed);
            ASSERT( it!=end() );
            if (it==end()) {
                return 0;
            }
            it->second.second.erase(name);
            size_t count =  --it->second.first ;
            if ( count == 0) {
                erase(it);
            }
            return count;
        }
        bool getHashAndSequenceFirstName(size_t hash,
                                         std::string sequence,
                                         std::string& put_name_here) {
            HashedSequence hashed(hash, sequence);
            auto it = find(hashed);
            if (it==end() || it->second.second.empty()) {
                return false;
            }
            put_name_here = *(it->second.second.begin());
            return true;
        }
        void loadAlignment(const StrVector& state_strings,
                           const std::vector<size_t>& old_hashes,
                           Alignment* aln,
                           const NameToIDMap& name_to_id,
                           const StrVector&   removed_seqs,
                           const StrVector&   twin_seqs) {
            int old_count = aln->getNSeq32();
            //Count the number of copies,
            //of each sequence in the old alignment
            for (int old_seq_id = 0 ;
                 old_seq_id < old_count; ++old_seq_id) {
                std::string old_sequence;
                aln->getOneSequence (state_strings, old_seq_id,
                                     old_sequence);
                trackHashAndSequence(old_hashes[old_seq_id],
                                     old_sequence,
                                     aln->getSeqName(old_seq_id));
            }
            //Also count the duplicate sequences,
            //that have been dropped from the alignment.
            size_t old_dupe_count = removed_seqs.size();
            for ( size_t i=0; i < old_dupe_count; ++i ) {
                auto f = name_to_id.find(twin_seqs[i]);
                if ( f != name_to_id.end() ) {
                    size_t      old_seq_id = f->second;
                    std::string old_sequence;
                    aln->getOneSequence (state_strings, old_seq_id,
                                         old_sequence);
                    trackHashAndSequence(old_hashes[old_seq_id],
                                         old_sequence,
                                         removed_seqs[i]);
                }
            }
        }
    };

    //Mapping names of duplicate sequences to the names of
    //the old sequences, they are copies of
    class NameToNameMap: public std::map<std::string, std::string> {
    public:
        NameToNameMap(const StrVector& keys, const StrVector& values) {
            ASSERT(keys.size() == values.size());
            for (size_t i = 0; i < keys.size(); ++i) {
                operator[](keys[i]) = values[i];
            }
        }
        void getKeysAndValues(StrVector& keys, StrVector& values) const {
            keys.clear();
            values.clear();
            for (auto it = begin(); it!=end(); ++it ) {
                keys.push_back(it->first);
                values.push_back(it->second);
            }
        }
    };
}

void PhyloTree::mergeAlignment(const Alignment* new_aln) {

    const std::string& path     = new_aln->aln_file;
    int            old_count    = aln->getNSeq32();
    int            new_count    = new_aln->getNSeq32();
    size_t         site_count   = aln->getNSite();
    std::string    task_name    = "Merging alignment " + path;
    StrVector      state_strings;
    
    aln->getStateStrings(state_strings);
    double work_estimate = static_cast<double>(old_count + new_count * 2 + 2 * site_count);
    initProgress(work_estimate,
                 task_name.c_str(), "", "");
    
    vector<size_t>   old_hashes     = aln->getPatternIndependentSequenceHashes(progress);
    vector<size_t>   new_hashes     = new_aln->getPatternIndependentSequenceHashes(progress);
    const StrVector& old_names      = aln->getSeqNames();
    auto             new_names      = new_aln->getSeqNames();
    size_t           count_the_same = 0;
    
    NameToIDMap   name_to_id   (old_names);
    NameToNameMap dupe_to_twin (removed_seqs, twin_seqs);

    UniqueSequenceMap unique_seq_map;
    if (params->ignore_identical_seqs) {
        unique_seq_map.loadAlignment(state_strings, old_hashes, aln,
                                     name_to_id, removed_seqs, twin_seqs);
    }
    
    IntVector added_sequences;   //ids of sequences in alignment just added
    std::vector<std::pair<int,int>> updated_sequences;
        //first is the sequence number in the original alignment
        //second is the sequence number in the additional (second) alignment
    StrVector names_of_updated_sequences;
    IntVector duplicated_sequence_ids;
        //ids of sequences, in alignment getting
        //merged, that were duplicated in the original sequence
    
    for (int new_seq_id = 0; new_seq_id < new_count; ++new_seq_id) {
        trackProgress(1.0);
        std::string new_name = new_names[new_seq_id];
        std::string new_sequence;
        new_aln->getOneSequence(state_strings, new_seq_id, new_sequence);
        auto find = name_to_id.find(new_name);
        bool was_old_dupe = false;
        if (find==name_to_id.end()) {
            auto find_dupe = dupe_to_twin.find(new_name);
            if (find_dupe!=dupe_to_twin.end()) {
                find = name_to_id.find(find_dupe->second);
                was_old_dupe = true;
            }
        }
        if (find==name_to_id.end()) {
            bool suppress = false;
            if (params->ignore_identical_seqs) {
                size_t copies = unique_seq_map.trackHashAndSequence
                                ( new_hashes[new_seq_id], new_sequence,
                                  new_name );
                if (2<copies) {
                    duplicated_sequence_ids.push_back(new_seq_id);
                    if (!params->suppress_list_of_sequences) {
                        LOG_LINE(VerboseMode::VB_MIN, "Sequence " << new_name
                                 << " duplicated in " << path);
                    }
                    suppress = true;
                }
            }
            if (!suppress) {
                added_sequences.push_back(new_seq_id);
                if (!params->suppress_list_of_sequences) {
                    LOG_LINE(VerboseMode::VB_MIN, "Sequence " << new_name
                             << " will be added from " << path);
                }
                int merged_seq_id = old_count + static_cast<int>(added_sequences.size());
                name_to_id[new_name] = merged_seq_id;
            }
        } else {
            int old_seq_id = find->second;
            if (old_count <= old_seq_id) {
                //This sequence name duplicates that of a sequence found
                //earlier in the same file.
                LOG_LINE(VerboseMode::VB_MIN, "Ignoring sequence with duplicated name "
                         << new_name << " found in " << path);
                continue;
            }
            //Name already seen in the alignment. Has the sequence changed?
            bool same = (old_hashes[old_seq_id] == new_hashes[new_seq_id]);
            std::string old_sequence;
            aln->getOneSequence(state_strings, old_seq_id, old_sequence);
            ASSERT(old_sequence.length() == new_sequence.length());
            if (same) {
                //
                //Todo: It would be better to ask Alignment to do this comparison
                //      (wouldn't have to construct the two sequence strings, then).
                //
                size_t len = old_sequence.length();
                for (size_t i = 0 ; i < len; ++i ) {
                    if (old_sequence[i]!=new_sequence[i]) {
                        LOG_LINE(VerboseMode::VB_DEBUG,
                                 "Sequence " << new_name << ": first changed "
                                 << " site was " << i << " (of " << len << ")");
                        same = false;
                        break;
                    }
                }
                if (same) {
                    ++count_the_same;
                    if (!params->suppress_list_of_sequences) {
                        LOG_LINE(VerboseMode::VB_MIN, "Sequence " << new_name
                                 << " unchanged in " << path);
                    }
                    continue;
                }
            } else {
                LOG_LINE(VerboseMode::VB_DEBUG, "Sequence " << new_name
                         << " had different hash (" << new_hashes[new_seq_id] << ")"
                         << " from before (" << old_hashes[old_seq_id] << ")");
            }
            bool still_a_dupe = false;
            if (params->ignore_identical_seqs) {
                unique_seq_map.untrackHashAndSequence(old_hashes[old_seq_id],
                                                      old_sequence, new_name);
                if (2<unique_seq_map.trackHashAndSequence  (new_hashes[new_seq_id],
                                                            new_sequence, new_name)) {
                    still_a_dupe = true;
                    std::string old_name;
                    unique_seq_map.getHashAndSequenceFirstName(new_hashes[new_seq_id],
                                                               new_sequence, old_name);
                    dupe_to_twin[new_name] = old_name;
                }
            }
            if (!params->suppress_list_of_sequences) {
                const char* preamble = was_old_dupe ? "Duplicated sequence " : "Sequence ";
                LOG_LINE(VerboseMode::VB_MIN, preamble << new_name << " updated in " << path);
            }
            if (was_old_dupe) {
                //this was an update to a sequence that was a duplicate
                if (!still_a_dupe) {
                    dupe_to_twin.erase(new_name);
                    added_sequences.push_back(new_seq_id);
                    if (!params->suppress_list_of_sequences) {
                        LOG_LINE(VerboseMode::VB_MIN, "Sequence " << new_name
                                 << " will be added from " << path);
                    }
                    int merged_seq_id = old_count + static_cast<int>(added_sequences.size());
                    name_to_id[new_name] = merged_seq_id;
                }
            } else {
                updated_sequences.emplace_back(old_seq_id, new_seq_id);
                names_of_updated_sequences.push_back(new_name);
                name_to_id[new_name] = old_count;
            }
        }
    }
    size_t count_updated          = updated_sequences.size();
    size_t count_added            = added_sequences.size();
    size_t count_updated_or_added = count_updated + count_added;
    size_t count_duplicated       = duplicated_sequence_ids.size();
    if (count_updated_or_added) {
        LOG_LINE(VerboseMode::VB_MIN, "Adding " << count_added << ","
                << " and Updating " << count_updated << " sequences");
        if (0<count_the_same) {
            LOG_LINE(VerboseMode::VB_MIN, "Leaving " << count_the_same
                     << " sequences unchanged" );
        }
        if (0<count_duplicated) {
            for (auto it=duplicated_sequence_ids.begin();
                 it!=duplicated_sequence_ids.end(); ++it) {
                size_t      new_seq_id = *it;
                std::string new_name   = new_names[new_seq_id];
                size_t      new_hash   = new_hashes[new_seq_id];
                std::string new_sequence;
                std::string old_name;
                new_aln->getOneSequence(state_strings,
                                        new_seq_id, new_sequence);
                if (unique_seq_map.getHashAndSequenceFirstName(new_hash,
                                                               new_sequence,
                                                               old_name ) ) {
                    dupe_to_twin[new_name] = old_name;
                } else {
                    //Nasty corner case.  This was a duplicate, of an existing sequence
                    //but that sequence has been replaced by an updated version (!).
                    //Ee-uw!  So we *can* add this sequence, after all!
                    --count_duplicated;
                    unique_seq_map.trackHashAndSequence(new_hash, new_sequence,
                                                        new_name);
                    added_sequences.push_back(static_cast<int>(new_seq_id));
                }
            }
            LOG_LINE(VerboseMode::VB_MIN, "Dropped " << count_duplicated
                     << " duplicated sequences from consideration ");
            LOG_LINE(VerboseMode::VB_MIN, "These will be added back after analysis is complete.");
        }
        aln->updateFrom(new_aln, updated_sequences,
                        added_sequences, progress);
        doneProgress();
        LOG_LINE(VerboseMode::VB_MED, "Incremental Method: "
                 << params->incremental_method);
        prepareForPlacement();
        IntVector ids_of_sequences_to_place;
        if (0<count_updated) {
            prepareForDeletes();
            removeTaxa(names_of_updated_sequences, false, "");
            doneDeletes(count_updated, true);
            for (auto it=updated_sequences.begin();
                 it!=updated_sequences.end(); ++it) {
                ids_of_sequences_to_place.push_back(it->first);
                //Each taxon that was removed (because it was updated)
                //has to be put back.
            }
        }
        for (auto it=added_sequences.begin();
             it!=added_sequences.end(); ++it) {
            //Each taxon that was added to the alignment,
            //needs to be added to the tree.
            //added_sequences gives sequence id
            //in the *new* alignment, not the *old*
            ids_of_sequences_to_place.push_back(old_count);
            ++old_count;
        }
        dupe_to_twin.getKeysAndValues(removed_seqs, twin_seqs);
        addNewTaxaToTree(ids_of_sequences_to_place, "Adding new taxa to tree", "");
    }
}

void PhyloTree::mergeAlignments(const StrVector& paths) {
    for (auto it=paths.begin(); it!=paths.end(); ++it ) {
        auto path = *it;
        Alignment *new_aln = createAlignment(path, aln->sequence_type.c_str(),
                                              params->intype, aln->model_name);
        //new_aln->isShowingProgressDisabled = true;
        if (params->remove_empty_seq) {
            Alignment* dead_aln = new_aln;
            new_aln = new_aln->removeGappySeq();
            if (dead_aln!=new_aln) {
                delete dead_aln;
            }
        }
        try {
            if (new_aln->getNSite()==aln->getNSite()) {
                new_aln->name             = path; //Just want the file
                                                  //name part of the path
                new_aln->model_name       = aln->model_name;
                new_aln->aln_file         = path;
                new_aln->position_spec    = "";
                new_aln->showNoProgress();
                ASSERT(aln->seq_type == new_aln->seq_type);
                if (aln->seq_type == new_aln->seq_type) {
                    mergeAlignment(new_aln);
                }
            }
        }
        catch (...) {
            delete new_aln;
            throw;
        }
        delete new_aln;
    }
    auto mergedTree = getTreeString();
    CKP_SAVE(mergedTree);
    saveCheckpoint();
    checkpoint->dump();
}

void PhyloTree::setRootNode(const char *my_root, bool multi_taxa) {
    if (rooted) {
        computeBranchDirection();
        return;
    }
    
    if (!my_root) {
        root = findNodeName(aln->getSeqName(0));
        ASSERT(root);
        return;
    }

    if (strchr(my_root, ',') == NULL) {
        string root_name = my_root;
        root = findNodeName(root_name);
        ASSERT(root);
        return;
    }

    // my_root is a list of taxa
    StrVector taxa;
    convert_string_vec(my_root, taxa);
    root = findNodeName(taxa[0]);
    ASSERT(root);
    if (!multi_taxa) {
        return;
    }
    unordered_set<string> taxa_set;
    for (auto it = taxa.begin(); it != taxa.end(); it++)
        taxa_set.insert(*it);
    pair<Node*,Neighbor*> res = {NULL, NULL};
    findNodeNames(taxa_set, res, root->neighbors[0]->node, root);
    if (res.first)
        root = res.first;
    else
        outWarning("Branch separating outgroup is not found");
}

void PhyloTree::readTreeString(const string &tree_string, bool keep_node_names) {
    stringstream str(tree_string);
    freeNode();
    
    MTree::readTree(str, rooted);
    if (!keep_node_names) {
        // bug fix 2016-04-14: in case taxon name happens to be ID
        assignLeafNames();
    } else {
        // bug fix 2021-05-10: in case taxon name *isn't* ID.
    }

    setRootNode(Params::getInstance().root);
    if (isSuperTree()) {
        ((PhyloSuperTree*) this)->mapTrees();
    }
    if (Params::getInstance().pll) {
        pllReadNewick(getTreeString());
    }
    resetCurScore();
    if (Params::getInstance().fixStableSplits ||
        Params::getInstance().adaptPertubation) {
        buildNodeSplit();
    }
    current_it = current_it_back = NULL;
}

void PhyloTree::readTreeStringSeqName(const string &tree_string) {
    stringstream str(tree_string);
    freeNode();
    if (rooted) {
        rooted = false;
        readTree(str, rooted);
        if (!rooted)
            convertToRooted();
    } else {
        readTree(str, rooted);
    }
    setAlignment(aln);
    setRootNode(params->root);

    if (isSuperTree()) {
        ((PhyloSuperTree*) this)->mapTrees();
    }
    if (params->pll) {
        pllReadNewick(getTreeString());
    }
    resetCurScore();
//    lhComputed = false;
    if (params->fixStableSplits) {
        buildNodeSplit();
    }
    current_it = current_it_back = NULL;
}

int PhyloTree::wrapperFixNegativeBranch(bool force_change) {
    // Initialize branch lengths for the parsimony tree
    initializeAllPartialPars();
    clearAllPartialLH();
    int numFixed = fixNegativeBranch(force_change);
    if (params->pll) {
        pllReadNewick(getTreeString());
    }
    resetCurScore();
    if (verbose_mode >= VerboseMode::VB_MAX) {
        printTree(cout);
    }
    return numFixed;
}

void PhyloTree::pllReadNewick(string newickTree) {
    pllNewickTree *newick = pllNewickParseString(newickTree.c_str());
    pllTreeInitTopologyNewick(pllInst, newick, PLL_FALSE);
    pllNewickParseDestroy(&newick);
}

void PhyloTree::readTreeFile(const string &file_name) {
    ifstream str;
    str.open(file_name.c_str());
    freeNode();
    if (rooted) {
        rooted = false;
        readTree(str, rooted);
        if (!rooted) {
            convertToRooted();
        }
    } else {
        readTree(str, rooted);
    }
    setAlignment(aln);
    if (isSuperTree()) {
        ((PhyloSuperTree*) this)->mapTrees();
    } else {
        clearAllPartialLH(false);
        clearAllScaleNum(false);
        clearAllPartialParsimony(false);
    }
    str.close();
    current_it = current_it_back = NULL;
}

string PhyloTree::getTreeString() {
    stringstream tree_stream;
    setRootNode(params->root);
    printTree(tree_stream, WT_TAXON_ID + WT_BR_LEN + WT_SORT_TAXA);
    return tree_stream.str();
}

string PhyloTree::getTopologyString(bool printBranchLength) {
    stringstream tree_stream;
    // important: to make topology string unique
    setRootNode(params->root);
    //printTree(tree_stream, WT_TAXON_ID + WT_SORT_TAXA);
    if (printBranchLength) {
        printTree(tree_stream, WT_SORT_TAXA + WT_BR_LEN + WT_TAXON_ID);
    } else {
        printTree(tree_stream, WT_SORT_TAXA);
    }
    return tree_stream.str();
}

void PhyloTree::rollBack(istream &best_tree_string) {
    best_tree_string.seekg(0, ios::beg);
    freeNode();
    readTree(best_tree_string, rooted);
    assignLeafNames();
    if (params->compute_likelihood) {
        initializeAllPartialLh();
        clearAllPartialLH();
        clearAllScaleNum(false);
    }
    clearAllPartialParsimony(false);
}

void PhyloTree::setModel(ModelSubst *amodel) {
    model = amodel;
    //state_freqs = new double[numStates];
    //model->getStateFrequency(state_freqs);
}

void PhyloTree::setModelFactory(ModelFactory *model_fac) {
    model_factory = model_fac;
    if (model_fac) {
        model = model_factory->model;
        site_rate = model_factory->site_rate;
        if (!isSuperTree()) {
            setLikelihoodKernel(sse);
            setNumThreads(num_threads);
        }
    } else {
        model = NULL;
        site_rate = NULL;
    }
}

bool PhyloTree::hasModelFactory() const {
    return model_factory != nullptr;
}

bool PhyloTree::hasModel() const {
    return model != nullptr;
}

void PhyloTree::setRate(RateHeterogeneity *rate) {
    site_rate = rate;
}

bool PhyloTree::hasRateHeterogeneity() const {
    return site_rate != nullptr;
}

RateHeterogeneity *PhyloTree::getRate() {
    return site_rate;
}

PhyloNode* PhyloTree::getRoot() {
    return (PhyloNode*)root;
}

PhyloNode* PhyloTree::newNode(int node_id, const char* node_name) {
    return new PhyloNode(node_id, node_name);
}

PhyloNode* PhyloTree::newNode(int node_id, int node_name) {
    return new PhyloNode(node_id, node_name);
}

bool PhyloTree::isDummyNode(PhyloNode* node) const {
    return node == DUMMY_NODE_1 || node == DUMMY_NODE_2;
}

void PhyloTree::clearAllPartialLH(bool set_to_null) {
    if (!root) {
        return;
    }
    PhyloNode*     atRoot     = getRoot();
    PhyloNeighbor* nei        = atRoot->firstNeighbor();
    PhyloNode*     nextToRoot = nei->getNode();
    nextToRoot->clearAllPartialLh(set_to_null, atRoot);
    atRoot->clearAllPartialLh(set_to_null, nextToRoot);
    tip_partial_lh_computed = 0;
    current_it      = nullptr;
    current_it_back = nullptr;
}

void PhyloTree::clearAllPartialParsimony(bool set_to_null) {
    if (!root) {
        return;
    }
    PhyloNode*     r          = getRoot();
    PhyloNeighbor* nei        = r->firstNeighbor();
    PhyloNode*     nextToRoot = nei->getNode();
    nextToRoot->clearAllPartialParsimony(set_to_null, r);
    r->clearAllPartialParsimony(set_to_null, nextToRoot);
    nei->setParsimonyComputed(false);
    if (set_to_null) {
        nei->partial_pars = nullptr;
    }
}

void PhyloTree::clearAllScaleNum(bool set_to_null) {
    if (!root) {
        return;
    }
    PhyloNode*     r          = getRoot();
    PhyloNeighbor* nei        = r->firstNeighbor();
    PhyloNode*     nextToRoot = nei->getNode();
    if (set_to_null) {
        nei->scale_num        = nullptr;
    }
    tip_partial_lh_computed   = 0;
    current_it                = nullptr;
    current_it_back           = nullptr;
    nextToRoot->clearAllScaleNum(set_to_null, r);
    r->clearAllScaleNum(set_to_null, nextToRoot);
    if (!set_to_null && central_scale_num != nullptr) {
        memset(central_scale_num, 0, central_scale_num_size_in_bytes);
    }
}


string getASCName(ASCType ASC_type) {
    switch (ASC_type) {
        case ASC_NONE:
            return "";
        case ASC_VARIANT:
            return "+ASC";
        case ASC_VARIANT_MISSING:
            return "+ASC_MIS";
        case ASC_INFORMATIVE:
            return "+ASC_INF";
        case ASC_INFORMATIVE_MISSING:
            return "+ASC_INF_MIS";
        default:
            ASSERT(0 && "ASC_type not recognized");
            return "";
    }
}

string PhyloTree::getSubstName() {
    return model->getName() + getASCName(model_factory->getASC());
}

string PhyloTree::getRateName() {
    if (model_factory->fused_mix_rate) {
        return "*" + site_rate->name.substr(1);
    } else {
        return site_rate->name;
    }
}

string PhyloTree::getModelName() {
    return getSubstName() + getRateName();
}

string PhyloTree::getModelNameParams() {
    string name = model->getNameParams();
    name += getASCName(model_factory->getASC());
    string rate_name = site_rate->getNameParams();

    if (model_factory->fused_mix_rate) {
        name += "*" + rate_name.substr(1);
    } else {
        name += rate_name;
    }

    return name;
}

void PhyloTree::saveBranchLengths(DoubleVector &lenvec, int startid,
                                  PhyloNode *node, PhyloNode *dad) {
    if (!node) {
        node = getRoot();
        ASSERT(branchNum == nodeNum-1);
        if (lenvec.empty()) {
            lenvec.resize(branchNum*getMixlen() + startid);
        }
    }
    FOR_EACH_PHYLO_NEIGHBOR(node, dad, it, nei){
        nei->getLength(lenvec, nei->id*getMixlen() + startid);
        PhyloTree::saveBranchLengths(lenvec, startid, nei->getNode(), node);
    }
}

void PhyloTree::restoreBranchLengths(DoubleVector &lenvec, int startid,
                                     PhyloNode *node, PhyloNode *dad) {
    if (!node) {
        node = getRoot();
        ASSERT(!lenvec.empty());
    }
    size_t mixLen = getMixlen();
    FOR_EACH_PHYLO_NEIGHBOR(node, dad, it, nei){
        PhyloNode* child = nei->getNode();
        intptr_t indexIntoVector = nei->id * mixLen;
        PhyloNeighbor* back_nei =  child->findNeighbor(node);
        ASSERT(nei->id == back_nei->id);
        nei->setLength(lenvec, static_cast<int>(indexIntoVector) + startid, static_cast<int>(mixLen));
        back_nei->setLength(lenvec, static_cast<int>(indexIntoVector) + startid, static_cast<int>(mixLen));
        PhyloTree::restoreBranchLengths(lenvec, startid, child, node);
    }
}


/****************************************************************************
 Parsimony function
 ****************************************************************************/

/*
 double PhyloTree::computeCorrectedParsimonyBranch(PhyloNeighbor *dad_branch, PhyloNode *dad) {
 //    double corrected_bran = 0;
 //    int parbran;
 //    int parscore = computeParsimonyBranch(node21_it, node2, &parbran);
 //    if (site_rate->getGammaShape() != 0) {
 //        corrected_bran = (aln->num_states - 1.0) / aln->num_states
 //                * site_rate->getGammaShape()
 //                * (pow( 1.0 - aln->num_states / (aln->num_states - 1.0) * ((double) parbran / aln->getNSite()),
 //                        -1.0 / site_rate->getGammaShape()) - 1.0);
 //    } else {
 //        corrected_bran = -((aln->num_states - 1.0) / aln->num_states)
 //                * log(1.0 - (aln->num_states / (aln->num_states - 1.0)) * ((double) parbran / aln->getNSite()));
 //    }
 //    return corrected_bran;
 }
 */
int PhyloTree::initializeAllPartialPars() {
    if (!ptn_freq_pars) {
        ptn_freq_pars = aligned_alloc<UINT>(get_safe_upper_limit_float(getAlnNPattern()));
    }
    int index = 0;
    initializeAllPartialPars(index);
    clearAllPartialParsimony(false);
    return index;
}

void PhyloTree::ensureCentralPartialParsimonyIsAllocated(size_t extra_vector_count) {
    if (central_partial_pars != nullptr) {
        return;
    }
    determineBlockSizes();
    uint64_t tip_partial_pars_size = get_safe_upper_limit_float(aln->num_states * (aln->STATE_UNKNOWN+1));
    uint64_t vector_count          = aln->getNSeq() * 4 - 2 + extra_vector_count;
    //2N-3 branches in an unrooted tree, 2N-1 in a rooted tree, and each branch needs
    //two vectors, so allocate 4N-2, just in case the tree is rooted.
    total_parsimony_mem_size       = vector_count * pars_block_size + tip_partial_pars_size;

    LOG_LINE(VerboseMode::VB_DEBUG, "Allocating " << total_parsimony_mem_size * sizeof(UINT)
             << " bytes for " << vector_count << " partial parsimony vectors"
             << " of " << (pars_block_size * sizeof(UINT)) << " bytes each, and "
             << (tip_partial_pars_size * sizeof(UINT)) << " additional bytes for tip vectors");
    
    try {
        central_partial_pars = aligned_alloc<UINT>(total_parsimony_mem_size);
    } catch (std::bad_alloc &) {
        outError("Not enough memory for partial parsimony vectors (bad_alloc)");
    }
    if (!central_partial_pars) {
        outError("Not enough memory for partial parsimony vectors");
    }
    tip_partial_pars = central_partial_pars + total_parsimony_mem_size - tip_partial_pars_size;
    
    LOG_LINE(VerboseMode::VB_MAX, "central_partial_pars is " << pointer_to_hex(central_partial_pars)
             << ", tip_partial_pars is " << pointer_to_hex(tip_partial_pars)
             << ", end of allocation is " << pointer_to_hex(tip_partial_pars+tip_partial_pars_size));
    LOG_LINE(VerboseMode::VB_MAX, "parsimony_index upper bound is "
             << ((tip_partial_pars - central_partial_pars) / pars_block_size)
             << " block size is "  << pars_block_size);
    return;
}

void PhyloTree::initializeAllPartialPars(int &index, PhyloNode *node, PhyloNode *dad) {
    if (!node) {
        ensureCentralPartialParsimonyIsAllocated(0);
        node  = getRoot();
        index = 0;
    }
    if (dad) {
        // assign blocks in central_partial_lh to both Neighbors (dad->node, and node->dad)
        PhyloNeighbor* backNei = node->findNeighbor(dad);
        backNei->partial_pars  = central_partial_pars + (index * pars_block_size);
        PhyloNeighbor* nei     = dad->findNeighbor(node);
        nei->partial_pars      = central_partial_pars + ((index + 1) * pars_block_size);
        index                 += 2;
        ASSERT( tip_partial_pars==nullptr || nei->partial_pars < tip_partial_pars );
    }
    FOR_EACH_ADJACENT_PHYLO_NODE(node, dad, it, child) {
        initializeAllPartialPars(index, child, node);
    }
}

size_t PhyloTree::getBitsBlockSize() {
    return pars_block_size;
}

UINT *PhyloTree::newBitsBlock() {
    return aligned_alloc<UINT>(pars_block_size);
}

void PhyloTree::computePartialParsimony(PhyloNeighbor *dad_branch, PhyloNode *dad) {
    (this->*computePartialParsimonyPointer)(dad_branch, dad);
}

int PhyloTree::getSubTreeParsimony(PhyloNeighbor* dad_branch) const {
    return (this->*getSubTreeParsimonyPointer)(dad_branch);
}

int PhyloTree::computeMarginalParsimony(PhyloNeighbor* dad_branch, PhyloNode* dad) {
    PhyloNode* node = dad_branch->getNode();
    computePartialParsimony(dad_branch, dad);
    int cost = this->getSubTreeParsimony(dad_branch);
    FOR_EACH_PHYLO_NEIGHBOR(node, dad, it, nei) {
        cost -= this->getSubTreeParsimony(nei);
    }
    return cost;
}

double PhyloTree::computePartialParsimonyOutOfTree(const UINT* left_partial_pars,
                                      const UINT* right_partial_pars,
                                                 UINT* dad_partial_pars) const {
    return (this->*computePartialParsimonyOutOfTreePointer)
        ( left_partial_pars, right_partial_pars, dad_partial_pars );
}

void PhyloTree::computePartialInfoDouble(TraversalInfo &info, double* buffer) {
    (this->*computePartialInfoPointer)(info, buffer);
}

void PhyloTree::computeReversePartialParsimony(PhyloNode *node, PhyloNode *dad) {
    PhyloNeighbor* node_nei = node->findNeighbor(dad);
    ASSERT(node_nei);
    computePartialParsimony(node_nei, node);
    FOR_EACH_ADJACENT_PHYLO_NODE(node, dad, it, child) {
        computeReversePartialParsimony(child, node);
    }
}

int PhyloTree::computeParsimonyBranch(PhyloNeighbor* dad_branch,
                                      PhyloNode* dad, int* branch_subst) {
    return (this->*computeParsimonyBranchPointer)(dad_branch, dad, branch_subst);
}

int PhyloTree::computeParsimonyOutOfTree(const UINT* dad_partial_pars,
                                         const UINT* node_partial_pars,
                                         int* branch_subst) const {
    return (this->*computeParsimonyOutOfTreePointer)
           (dad_partial_pars, node_partial_pars, branch_subst);
}

int PhyloTree::computeParsimony(const char* taskDescription,
                                bool bidirectional, bool countProgress,
                                PhyloNeighbor* neighbor,
                                PhyloNode* starting_node) {
    if (central_partial_pars == nullptr) {
        initializeAllPartialPars();
    }
    PhyloNode*     r   = (starting_node!=nullptr) ? starting_node : getRoot();
    PhyloNeighbor* nei = (neighbor!=nullptr)      ? neighbor : r->firstNeighbor();
    if (taskDescription==nullptr || taskDescription[0]=='\0') {
        return computeParsimonyBranch(r->firstNeighbor(), r);
    }
    ParallelParsimonyCalculator calculator(*this, countProgress);
    if (bidirectional) {
        return calculator.computeAllParsimony(nei, r);
    }
    return calculator.computeParsimonyBranch
           ( nei, r, taskDescription );
}

int PhyloTree::computeParsimony(const std::string& taskDescription,
                                bool bidirectional, bool countProgress,
                                PhyloNeighbor* neighbor,
                                PhyloNode* starting_node) {
    return computeParsimony(taskDescription.c_str(), bidirectional,
                            countProgress, neighbor, starting_node);
}

/****************************************************************************
 likelihood function
 ****************************************************************************/

size_t PhyloTree::getBufferPartialLhSize() {
    const size_t VECTOR_SIZE = 8; // TODO, adjusted
    // 2017-12-13: make sure that num_threads was already set
    ASSERT(num_threads > 0);
    size_t ncat_mix    = site_rate->getNRate() * ((model_factory->fused_mix_rate)? 1 : model->getNMixtures());
    size_t block       = model->num_states * ncat_mix;
    size_t buffer_size = 0;

    // buffer for traversal_info.echildren and partial_lh_leaves
    if (!Params::getInstance().buffer_mem_save) {
        buffer_size += get_safe_upper_limit(block * model->num_states * 2) * aln->getNSeq();
        buffer_size += get_safe_upper_limit(block *(aln->STATE_UNKNOWN+1)) * aln->getNSeq();
    }

    buffer_size += get_safe_upper_limit(block *(aln->STATE_UNKNOWN+1));
    buffer_size += (block*2+model->num_states)*VECTOR_SIZE*num_packets;

    // always more buffer for non-rev kernel, in case switching between kernels
    buffer_size += get_safe_upper_limit(block)*(aln->STATE_UNKNOWN+1)*2;
    buffer_size += block*2*VECTOR_SIZE*num_packets;
    buffer_size += get_safe_upper_limit(3*block*model->num_states);

    if (isMixlen()) {
        size_t nmix = max(getMixlen(), getRate()->getNRate());
        buffer_size += nmix*(nmix+1)*VECTOR_SIZE + (nmix+3)*nmix*VECTOR_SIZE*num_packets;
    }
    return buffer_size;
}

void PhyloTree::ensurePartialLHIsAllocated(size_t count_of_extra_parsimony_blocks,
                                           size_t count_of_extra_lh_blocks) {
    determineBlockSizes();

    int numStates = model->num_states;
    // Minh's question: why getAlnNSite() but not getAlnNPattern() ?
    //size_t mem_size = ((getAlnNSite() % 2) == 0) ? getAlnNSite() : (getAlnNSite() + 1);
    // extra #numStates for ascertainment bias correction
    size_t mem_size = get_safe_upper_limit(getAlnNPattern()) + max(get_safe_upper_limit(numStates),
        get_safe_upper_limit(model_factory->unobserved_ptns.size()));

    tree_buffers.theta_block_size = mem_size * numStates * site_rate->getNRate()
                                    * ((model_factory->fused_mix_rate)? 1 : model->getNMixtures());
    // make sure _pattern_lh size is divisible by 4 (e.g., 9->12, 14->16)
    size_t lhcat_size = mem_size * site_rate->getNDiscreteRate()
                      * ((model_factory->fused_mix_rate)? 1 : model->getNMixtures());
    
    tree_buffers.ensurePatternLhAllocated(mem_size);
    tree_buffers.ensurePatternLhCatAllocated(lhcat_size);
    tree_buffers.ensurePartialLhAllocated(getBufferPartialLhSize());

    if (!_site_lh && (params->robust_phy_keep < 1.0 || params->robust_median)) {
        _site_lh = aligned_alloc<double>(getAlnNSite());
    }
    ptn_freq_computed &= (ptn_freq!=nullptr);
    
    tree_buffers.ensureThetaAllocated(tree_buffers.theta_block_size);
    tree_buffers.ensureScaleAllAllocated(mem_size);
    
    ensure_aligned_allocated(ptn_freq,          mem_size);
    ensure_aligned_allocated(ptn_freq_pars,     mem_size);
    ensure_aligned_allocated(ptn_invar,         mem_size);

    allocateCentralBlocks(count_of_extra_parsimony_blocks,
                          count_of_extra_lh_blocks);
}

void PhyloTree::initializeAllPartialLh() {
    ASSERT(params->compute_likelihood);
    clearAllPartialLH(true);
    clearAllScaleNum(true);
    deleteAllPartialParsimony();
    getMemoryRequired();
    ensurePartialLHIsAllocated(0,0);

    int index_parsimony = 0;
    int index_lh = 0;
    initializeAllPartialLh(index_parsimony, index_lh);
    if (params->lh_mem_save == LM_MEM_SAVE) {
        mem_slots.init(this, static_cast<int>(max_lh_slots));
    }
    ASSERT(index_parsimony == (nodeNum - 1) * 2);
    if (params->lh_mem_save == LM_PER_NODE) {
        ASSERT(index_lh == nodeNum-leafNum);
    }
    clearAllPartialLH();
}

void PhyloTree::deleteAllPartialParsimony() {
    aligned_free(central_partial_pars);
    tip_partial_pars        = nullptr;
    clearAllPartialParsimony(true);
    tip_partial_lh_computed &= ~2;
}

void PhyloTree::deleteAllPartialLh() {
    //Note: aligned_free now sets the pointer to nullptr
    //      (so there's no need to do that explicitly any more)
    aligned_free(central_partial_lh);
    aligned_free(central_scale_num);
    aligned_free(nni_scale_num);
    aligned_free(nni_partial_lh);
    aligned_free(ptn_invar);
    aligned_free(ptn_freq);
    aligned_free(ptn_freq_pars);
    tree_buffers.freeBuffers();
    aligned_free(_site_lh);

    ptn_freq_computed        = false;
    tip_partial_lh           = nullptr;
    tip_partial_lh_computed &= ~1; //was missing!

    clearAllPartialLH(true);
    clearAllScaleNum(true);
}

void PhyloTree::deleteAllPartialLhAndParsimony() {
    deleteAllPartialLh();
    deleteAllPartialParsimony();
}
 
uint64_t PhyloTree::getMemoryRequired(size_t ncategory, bool full_mem) {
    // +num_states for ascertainment bias correction
    int64_t nptn = get_safe_upper_limit(aln->getNPattern()) + get_safe_upper_limit(aln->num_states);
    if (model_factory) {
        nptn = get_safe_upper_limit(aln->getNPattern()) + max(get_safe_upper_limit(aln->num_states), get_safe_upper_limit(model_factory->unobserved_ptns.size()));
    }
    int64_t scale_block_size = nptn;
    if (site_rate) {
        scale_block_size *= site_rate->getNRate();
    }
    else {
        scale_block_size *= ncategory;
    }
    if (model && !model_factory->fused_mix_rate) {
        scale_block_size *= model->getNMixtures();
    }

    int64_t block_size = scale_block_size * aln->num_states;
    int64_t mem_size;
    // memory to tip_partial_lh
    if (model) {
        mem_size = aln->num_states * (aln->STATE_UNKNOWN+1) * model->getNMixtures() * sizeof(double);
    }
    else {
        mem_size = aln->num_states * (aln->STATE_UNKNOWN+1) * sizeof(double);
    }

    // memory for UFBoot
    if (params->gbo_replicates) {
        mem_size += params->gbo_replicates*nptn*sizeof(BootValType);
    }

    // memory for model
    if (model) {
        mem_size += model->getMemoryRequired();
    }

    int64_t lh_scale_size = block_size * sizeof(double) + scale_block_size * sizeof(UBYTE);

    int64_t max_leaves = leafNum;
    if (aln!=nullptr && max_leaves<static_cast<int64_t>(aln->getNSeq())+1) {
        //If -mlnj-only has been supplied, then it is possible
        //this function might have been called before a tree
        //has been constructed, in which case leafNum will still
        //be zero.  But we know what leafNum might become...
        //(number of sequences, plus one for a root).
        max_leaves = aln->getNSeq()+1;
    }
    max_lh_slots = max_leaves - 2;

    if (!full_mem && params->lh_mem_save == LM_MEM_SAVE) {
        int64_t min_lh_slots = static_cast<int64_t>(floor(log2(max_leaves)))+LH_MIN_CONST;
        double  fraction = (params->max_mem_is_in_bytes)
                         ? ((double)params->max_mem_size / (double)mem_size)
                         : params->max_mem_size;
        int64_t max_bytes = (params->max_mem_is_in_bytes)
                          ? static_cast<int64_t>(params->max_mem_size) : getMemorySize();
         
        if (params->max_mem_size == 0.0) {
            max_lh_slots = min_lh_slots;
        } else if (fraction<= 1) {
            max_lh_slots = static_cast<int64_t>(floor(fraction*(max_leaves-2)));
        } else {
            int64_t rest_mem = max_bytes - mem_size;
            
            // include 2 blocks for nni_partial_lh
            max_lh_slots = rest_mem / lh_scale_size - 2;

            // RAM over requirement, reset to LM_PER_NODE
            if (max_lh_slots + 2 > max_leaves)
            {
                max_lh_slots = max_leaves-2;
                if (max_lh_slots < min_lh_slots) {
                    max_lh_slots = min_lh_slots;
                }
            }
        }
        if (max_lh_slots < min_lh_slots) {
            double newSizeInMegabytes = (mem_size + (min_lh_slots+2)*lh_scale_size)/1048576.0 ;
            if (1<newSizeInMegabytes) {
                hideProgress();
                cout << "WARNING: Too low -mem, automatically increased to "
                    << newSizeInMegabytes << " MB" << endl;
                showProgress();
            }
            max_lh_slots = min_lh_slots;
        }
    }

    // also count MEM for nni_partial_lh
    mem_size += (max_lh_slots+2) * lh_scale_size;
    return mem_size;
}

uint64_t PhyloTree::getMemoryRequiredThreaded(size_t ncategory, bool full_mem) {
    return getMemoryRequired(ncategory, full_mem);
}

void PhyloTree::getMemoryRequired(uint64_t &partial_lh_entries,
                                  uint64_t &scale_num_entries,
                                  uint64_t &partial_pars_entries) {
    // +num_states for ascertainment bias correction
    uint64_t block_size = get_safe_upper_limit(aln->getNPattern()) + get_safe_upper_limit(aln->num_states);
    if (model_factory)
        block_size = get_safe_upper_limit(aln->getNPattern()) + max(get_safe_upper_limit(aln->num_states), get_safe_upper_limit(model_factory->unobserved_ptns.size()));
    size_t scale_size = block_size;
    block_size = block_size * aln->num_states;
    if (site_rate) {
        block_size *= site_rate->getNRate();
        scale_size *= site_rate->getNRate();
    }
    if (model!=nullptr && model_factory!=nullptr
        && !model_factory->fused_mix_rate) {
        block_size *= model->getNMixtures();
        scale_size *= model->getNMixtures();
    }

    uint64_t tip_partial_lh_size   = aln->num_states * (aln->STATE_UNKNOWN+1) * model->getNMixtures();
    uint64_t tip_partial_pars_size = aln->num_states * (aln->STATE_UNKNOWN+1);

    int64_t max_leaves = leafNum;
    if (aln!=nullptr && max_leaves<static_cast<int64_t>(aln->getNSeq())+1) {
        //If -mlnj-only has been supplied, then it is possible
        //this function might have been called before a tree
        //has been constructed, in which case leafNum will still
        //be zero.  But we know what leafNum might become...
        //(number of sequences, plus one for a root).
        max_leaves = aln->getNSeq()+1;
    }

    // TODO mem save
    partial_lh_entries = ((uint64_t)max_leaves - 2) * (uint64_t) block_size + 4 + tip_partial_lh_size;
    scale_num_entries  = (max_leaves - 2) * scale_size;

    partial_pars_entries = (max_leaves - 1) * 4 * pars_block_size + tip_partial_pars_size;
}

void PhyloTree::determineBlockSizes() {

    size_t old_pars_block_size  = pars_block_size;
    size_t old_scale_block_size = scale_block_size;
    size_t old_lh_block_size    = lh_block_size;

    // reserve the last entry for parsimony score
    // no longer: pars_block_size = (aln->num_states * aln->size() + UINT_BITS - 1) / UINT_BITS + 1;
    if (isUsingSankoffParsimony()) {
        //Sankoff parsimony tracks a number (a UINT) for each state for each site
        //(and should have 1 additional UINT for a total score)
        intptr_t ptnCount = aln->ordered_pattern.size();
        size_t nstates  = aln->num_states;
        pars_block_size = ptnCount * nstates + 1;
    } else {
        //Otherwise, parsimony tracks a bit per site per state
        //(and should have 1 additional UINTs for a total score)
        size_t bits_per_state  = max(aln->size(), (size_t)aln->num_variant_sites);
        size_t uints_per_state = (bits_per_state + SIMD_BITS - 1) / UINT_BITS;
        pars_block_size = aln->getMaxNumStates() * uints_per_state + 1;
    }
    pars_block_size = get_safe_upper_limit_float(pars_block_size);

    // +num_states for ascertainment bias correction
    size_t states     = aln->num_states;
    size_t unobserved = (model_factory == nullptr)
                        ? 0 : model_factory->unobserved_ptns.size();
    intptr_t nptn     = get_safe_upper_limit(aln->size())
                        + get_safe_upper_limit(max(states,unobserved));
    size_t mixtures   = (model_factory == nullptr || model_factory->fused_mix_rate || model==nullptr)
                        ? 1 : model->getNMixtures();
    size_t rates      = (site_rate!=nullptr) ? site_rate->getNRate() : 1;
    scale_block_size  = nptn * rates * mixtures;
    size_t mod_states = (model==nullptr) ? aln->num_states : model->num_states;
    lh_block_size     = scale_block_size * mod_states;

    if (0 < old_scale_block_size && old_scale_block_size < scale_block_size
        && this->central_scale_num != nullptr) {
        deleteAllPartialLh();
    }
    else if (0 < old_lh_block_size && old_lh_block_size < lh_block_size
             && this->central_partial_lh != nullptr) {
        deleteAllPartialLh();
    }
    if (0 < old_pars_block_size && old_pars_block_size < pars_block_size
        && this->central_partial_pars != nullptr) {
        deleteAllPartialParsimony();
    }
}

size_t PhyloTree::getLhBlockSize() {
    return lh_block_size;
}

void PhyloTree::allocateCentralBlocks(size_t extra_parsimony_block_count,
                                      size_t extra_lh_block_count) {
    determineBlockSizes();

    // allocate the big central partial likelihoods memory
    size_t IT_NUM = 2;
    if (!nni_partial_lh) {
        // allocate memory only once!
        nni_partial_lh = aligned_alloc<double>(IT_NUM*lh_block_size);
        nni_scale_num = aligned_alloc<UBYTE>(IT_NUM*scale_block_size);
    }
    
    int mixtures = (model == nullptr) ? 1 : model->getNMixtures();
    uint64_t tip_partial_lh_size = get_safe_upper_limit(aln->num_states * (aln->STATE_UNKNOWN+1) * mixtures );
    if (!central_partial_lh) {
        if (model!=nullptr && model->isSiteSpecificModel())
        {
            tip_partial_lh_size = get_safe_upper_limit(aln->size()) * model->num_states * leafNum;
        }
        if (max_lh_slots == 0) {
            getMemoryRequired();
        }
        uint64_t mem_size = (uint64_t)(max_lh_slots + extra_lh_block_count) * lh_block_size
            + 4 + tip_partial_lh_size;

        if (0<extra_lh_block_count) {
            LOG_LINE(VerboseMode::VB_MED, "max_lh_slots was " << max_lh_slots
                << " and extra_lh_block_count was " << extra_lh_block_count);
        }
        LOG_LINE (VerboseMode::VB_DEBUG, 
                  "Allocating " << mem_size * sizeof(double) << " bytes for "
                  << (max_lh_slots + extra_lh_block_count) << " partial likelihood vectors "
                  << "( " << (tip_partial_lh_size+4) * sizeof(double) << " bytes for tips and scores)");
        try {
            central_partial_lh = aligned_alloc<double>(mem_size);
        } catch (std::bad_alloc &) {
            outError("Not enough memory for partial likelihood vectors (bad_alloc)");
        }
        if (!central_partial_lh) {
            outError("Not enough memory for partial likelihood vectors");
        }
        LOG_LINE (VerboseMode::VB_DEBUG, 
                  "central_partial_lh is " << pointer_to_hex(central_partial_lh)
                  << " through " << pointer_to_hex(central_partial_lh + mem_size));
        memset(central_partial_lh, 0, mem_size*sizeof(double));
                
        // We need not treat params->lh_mem_save == LM_PER_NODE as a special case.
        tip_partial_lh = central_partial_lh + ((max_lh_slots + extra_lh_block_count) * lh_block_size);
        LOG_LINE (VerboseMode::VB_DEBUG, 
                  "tip_partial_lh is " << pointer_to_hex(tip_partial_lh)
                  << " through " << pointer_to_hex(tip_partial_lh + tip_partial_lh_size ));
    }


    if (central_scale_num == nullptr) {
        auto slots_wanted = max_lh_slots + extra_lh_block_count;
        uint64_t mem_size = slots_wanted * scale_block_size;
        central_scale_num_size_in_bytes = mem_size * sizeof(UBYTE);
        LOG_LINE(VerboseMode::VB_MAX, 
                 "Allocating " << central_scale_num_size_in_bytes 
                 << " bytes for scale num vectors");
        try {
            central_scale_num = aligned_alloc<UBYTE>(mem_size);
        } catch (std::bad_alloc &) {
            outError("Not enough memory for scale num vectors (bad_alloc)");
        }
        if (!central_scale_num) {
            outError("Not enough memory for scale num vectors");
        }
        LOG_LINE (VerboseMode::VB_DEBUG, 
                  "Allocated " << central_scale_num_size_in_bytes << " bytes"
                  << " for " << slots_wanted << " scale blocks");
        LOG_LINE (VerboseMode::VB_DEBUG, 
                  "Address range for scale blocks is " << pointer_to_hex(central_scale_num)
                  << " to " << pointer_to_hex(central_scale_num + mem_size) );
    }
    ensureCentralPartialParsimonyIsAllocated(extra_parsimony_block_count);
}
    
void PhyloTree::initializeAllPartialLh(int &index_pars, int &index_lh,
                                       PhyloNode *node, PhyloNode *dad) {
    ASSERT(params->compute_likelihood);
    if (!node) {
        node = getRoot();
        allocateCentralBlocks(0, 0);
        index_pars = 0;
        index_lh   = 0;
    }
    if (dad) {
        //assign blocks of memory in central_partial_pars and
        //central_partial_lh to both Neighbors (dad->node, and node->dad)
        PhyloNeighbor *nei  = isDummyNode(node) ? nullptr : node->findNeighbor(dad);
        PhyloNeighbor *nei2 = isDummyNode(node) ? nullptr : dad->findNeighbor(node);
        
        // first initialize partial_pars
        if (nei!=nullptr && nei->partial_pars==nullptr) {
            nei->partial_pars = central_partial_pars + (index_pars * pars_block_size);
            ++index_pars;
        }
        if (nei2!=nullptr && nei2->partial_pars==nullptr) {
            nei2->partial_pars = central_partial_pars + (index_pars * pars_block_size);
            ++index_pars;
        }
        ASSERT(index_pars < nodeNum * 2 - 1);
        
        // now initialize partial_lh and scale_num
        if (params->lh_mem_save == LM_PER_NODE) {
            if (nei2!=nullptr && !node->isLeaf()) {
                // only allocate memory to internal node
                if (nei2->partial_lh==nullptr) {
                    nei2->scale_num  = central_scale_num  + (index_lh * scale_block_size);
                    nei2->partial_lh = central_partial_lh + (index_lh * lh_block_size);
#if (0)
                    LOG_LINE ( VerboseMode::VB_MAX, "allocating partial_lh block " << index_lh
                        << " " << pointer_to_hex(nei2->partial_lh)
                        << " to front-neighbour " << pointer_to_hex(nei2)
                              << " of node " << pointer_to_hex(dad));
#endif
                    ++index_lh;
                }
            }
        }
    }
    FOR_EACH_PHYLO_NEIGHBOR(node, dad, it, nei) {
        PhyloNode*   child = nei->getNode();
        if (!isDummyNode(child)) {
            initializeAllPartialLh(index_pars, index_lh, child, node);
        } else {
            if (nei->partial_pars==nullptr) {
                nei->partial_pars = central_partial_pars + (index_pars * pars_block_size);
                ++index_pars;
            }
            if (nei->partial_lh==nullptr) {
                nei->scale_num  = central_scale_num  + (index_lh * scale_block_size);
                nei->partial_lh = central_partial_lh + (index_lh * lh_block_size);
#if (0)
                LOG_LINE( VB_MAX, "allocating partial_lh block " << index_lh
                    << " " << pointer_to_hex(nei->partial_lh)
                    << " to root-neighbour " << pointer_to_hex(nei)
                    << " of node " << pointer_to_hex(node));
#endif
                ++index_lh;
            }
        }
    }
}

double *PhyloTree::newPartialLh() {
    return aligned_alloc<double>(getPartialLhSize());
}

size_t PhyloTree::getPartialLhSize() {
    // +num_states for ascertainment bias correction
    size_t block_size = get_safe_upper_limit(aln->size())+max(get_safe_upper_limit(aln->num_states),
        get_safe_upper_limit(model_factory->unobserved_ptns.size()));
    block_size *= model->num_states * site_rate->getNRate() * ((model_factory->fused_mix_rate)? 1 : model->getNMixtures());
    return block_size;
}

size_t PhyloTree::getPartialLhBytes() {
    // +num_states for ascertainment bias correction
    return getPartialLhSize() * sizeof(double);
}

size_t PhyloTree::getScaleNumSize() {
    size_t block_size = get_safe_upper_limit(aln->size())+max(get_safe_upper_limit(aln->num_states),
        get_safe_upper_limit(model_factory->unobserved_ptns.size()));
    return (block_size) * site_rate->getNRate() * ((model_factory->fused_mix_rate)? 1 : model->getNMixtures());
}

size_t PhyloTree::getScaleNumBytes() {
    return getScaleNumSize()*sizeof(UBYTE);
}

UBYTE *PhyloTree::newScaleNum() {
    return aligned_alloc<UBYTE>(getScaleNumSize());
}

double PhyloTree::computeLikelihood(double *pattern_lh) {
    ASSERT(model);
    ASSERT(site_rate);
    ASSERT(root->isLeaf());
    if (current_it==nullptr) {
        PhyloNode* leaf = findFarthestLeaf();
        current_it      = leaf->firstNeighbor();
        current_it_back = current_it->getNode()->findNeighbor(leaf);
    }
    PhyloNode* node = current_it_back->getNode();
    curScore        = computeLikelihoodBranch(current_it, node, tree_buffers);
    if (!pattern_lh) {
        return curScore;
    }
    memmove(pattern_lh, tree_buffers._pattern_lh, aln->size() * sizeof(double));
    if (current_it->lh_scale_factor < 0.0) {
        intptr_t nptn = aln->getNPattern();
        for (intptr_t i = 0; i < nptn; i++) {
            pattern_lh[i] += max(current_it->scale_num[i], UBYTE(0)) * LOG_SCALING_THRESHOLD;
        }
    }
    return curScore;
}

//double PhyloTree::computeLikelihoodRooted(PhyloNeighbor *dad_branch, PhyloNode *dad) {
//    double score = computeLikelihoodBranchNaive(dad_branch, dad);
//    if (verbose_mode >= VerboseMode::VB_DEBUG) {
//        printTransMatrices(dad_branch->node, dad);
//        /*
//         FOR_EACH_PHYLO_NEIGHBOR(dad_branch->node, dad, it, pit) {
//              cout << pit->node->name << "\t" << pit->partial_lh[0] << endl;
//         }*/
//    }
//    double* state_freq = new double[aln->num_states];
//    model->getStateFrequency(state_freq);
//    score -= log(state_freq[(int) root_state]);
//    delete[] state_freq;
//    return score;
//}

int PhyloTree::getNumLhCat(SiteLoglType wsl) {
    int ncat = 0;
    switch (wsl) {
    case WSL_NONE: ASSERT(0 && "is not WSL_NONE"); return 0;
    case WSL_SITE: ASSERT(0 && "is not WSL_SITE"); return 0;
    case WSL_MIXTURE_RATECAT: 
        ncat = getRate()->getNDiscreteRate();
        if (getModel()->isMixture() && !getModelFactory()->fused_mix_rate)
            ncat *= getModel()->getNMixtures();
        return ncat;
    case WSL_RATECAT:
        return getRate()->getNDiscreteRate();
    case WSL_MIXTURE:
        return getModel()->getNMixtures();
    default:
        ASSERT(0 && "SiteLoglType recognized");
        return 0;
    }
}

void PhyloTree::transformPatternLhCat() {
    if (vector_size == 1) {
        return;
    }

    intptr_t nptn = ((aln->size()+vector_size-1)/vector_size)*vector_size;
//    size_t nstates = aln->num_states;
    size_t ncat = site_rate->getNRate();
    if (!model_factory->fused_mix_rate) ncat *= model->getNMixtures();

    double *mem = aligned_alloc<double>(nptn*ncat);
    memcpy(mem, tree_buffers._pattern_lh_cat, nptn*ncat*sizeof(double));
    double *memptr = mem;

    for (intptr_t ptn = 0; ptn < nptn; ptn+=vector_size) {
        double *lh_cat_ptr = &tree_buffers._pattern_lh_cat[ptn*ncat];
        for (size_t cat = 0; cat < ncat; cat++) {
            for (size_t i = 0; i < vector_size; i++) {
                lh_cat_ptr[i*ncat+cat] = memptr[i];
            }
            memptr += vector_size;
        }
    }
    aligned_free(mem);
}

PhyloNode* PhyloTree::findFirstFarLeaf(PhyloNode* node) const {
    PhyloNode* dad = nullptr;
    do {
        FOR_EACH_PHYLO_NEIGHBOR(node, dad, it, nei) {
            dad = node;
            node = nei->getNode();
            break;
        }
    } while (!node->isLeaf());
    return node;
}

double PhyloTree::computePatternLhCat(SiteLoglType wsl) {
    if (!current_it) {
        PhyloNode *leaf = findFirstFarLeaf(getRoot());
        current_it      = leaf->firstNeighbor();
        current_it_back = current_it->getNode()->findNeighbor(leaf);
    }
    double score  = computeLikelihoodBranch(current_it,
                                            current_it_back->getNode(),
                                            tree_buffers);
    transformPatternLhCat();

    if (!getModel()->isSiteSpecificModel() && getModel()->isMixture() 
        && !getModelFactory()->fused_mix_rate) {
        if (wsl == WSL_MIXTURE || wsl == WSL_RATECAT) {
            double *lh_cat = tree_buffers._pattern_lh_cat;
            double *lh_res = tree_buffers._pattern_lh_cat;
            intptr_t ptn, nptn = aln->getNPattern();
            size_t m, nmixture = getModel()->getNMixtures();
            size_t c, ncat = getRate()->getNRate();
            if (wsl == WSL_MIXTURE && ncat > 1) {
                // transform to lh per mixture class
                for (ptn = 0; ptn < nptn; ptn++) {
                    for (m = 0; m < nmixture; m++) {
                        double lh = lh_cat[0];
                        for (c = 1; c < ncat; c++)
                            lh += lh_cat[c];
                        lh_res[m] = lh;
                        lh_cat += ncat;
                    }
                    lh_res += nmixture;
                }
            } else if (wsl == WSL_RATECAT && nmixture > 1) {
                // transform to lh per rate category
                for (ptn = 0; ptn < nptn; ptn++) {
                    if (lh_res != lh_cat)
                        memcpy(lh_res, lh_cat, ncat*sizeof(double));
                    lh_cat += ncat;
                    for (m = 1; m < nmixture; m++) {
                        for (c = 0; c < ncat; c++)
                            lh_res[c] += lh_cat[c];
                        lh_cat += ncat;
                    }
                    lh_res += ncat;
                }
            }
        }
    }
    return score;
}

void PhyloTree::computePatternStateFreq(double *ptn_state_freq) {
    ASSERT(getModel()->isMixture());
    computePatternLhCat(WSL_MIXTURE);
    double*  lh_cat   = tree_buffers._pattern_lh_cat;
    intptr_t nptn     = getAlnNPattern();
    int      nmixture = getModel()->getNMixtures();
    double*  freq_ptr = ptn_state_freq;
    int      nstates  = aln->num_states;
//    ModelMixture *models = (ModelMixture*)model;
    
    if (params->print_site_state_freq == WSF_POSTERIOR_MEAN) {
        cout << "Computing posterior mean site frequencies...." << endl;
        // loop over all site-patterns
        for (intptr_t ptn = 0; ptn < nptn; ++ptn) {
        
            // first compute posterior for each mixture component
            double sum_lh = 0.0;
            for (size_t m = 0; m < nmixture; ++m) {
                sum_lh += lh_cat[m];
            }
            sum_lh = 1.0/sum_lh;
            for (size_t m = 0; m < nmixture; ++m) {
                lh_cat[m] *= sum_lh;
            }
            
            // now compute state frequencies
            for (int state = 0; state < nstates; ++state) {
                double freq = 0;
                for (int m = 0; m < nmixture; ++m)
                    freq += model->getMixtureClass(m)->state_freq[state] * lh_cat[m];
                freq_ptr[state] = freq;
            }
            
            // increase the pointers
            lh_cat   += nmixture;
            freq_ptr += nstates;
        }
    } else if (params->print_site_state_freq == WSF_POSTERIOR_MAX) {
        cout << "Computing posterior max site frequencies...." << endl;
        // loop over all site-patterns
        for (intptr_t ptn = 0; ptn < nptn; ++ptn) {
            // first compute posterior for each mixture component
            int max_comp = 0;
            for (int m = 1; m < nmixture; ++m)
                if (lh_cat[m] > lh_cat[max_comp]) {
                    max_comp = m;
                }
            
            // now compute state frequencies
            memcpy(freq_ptr, model->getMixtureClass(max_comp)->state_freq, static_cast<int>(nstates*sizeof(double)));
            
            // increase the pointers
            lh_cat += nmixture;
            freq_ptr += nstates;
        }
    }
}



void PhyloTree::computePatternLikelihood(double *ptn_lh, double *cur_logl, double *ptn_lh_cat, SiteLoglType wsl) {
    /*    if (!dad_branch) {
     dad = getRoot();
     dad_branch = dad->firstNeighbor();
     }*/
    intptr_t nptn = aln->getNPattern();
    int      ncat = getNumLhCat(wsl);
    if (ptn_lh_cat) {
        // Right now only Naive version store _pattern_lh_cat!
        computePatternLhCat(wsl);
    }
    
    double sum_scaling = current_it->lh_scale_factor + current_it_back->lh_scale_factor;
    if (sum_scaling < 0.0) {
        if (current_it->lh_scale_factor == 0.0) {
            for (intptr_t i = 0; i < nptn; i++) {
                UBYTE scale_up = max(UBYTE(0), current_it_back->scale_num[i]);
                ptn_lh[i] = tree_buffers._pattern_lh[i] + scale_up * LOG_SCALING_THRESHOLD;
            }
        } else if (current_it_back->lh_scale_factor == 0.0){
            for (intptr_t i = 0; i < nptn; i++) {
                UBYTE scale_up = max(UBYTE(0), current_it->scale_num[i]);
                ptn_lh[i] = tree_buffers._pattern_lh[i] + scale_up * LOG_SCALING_THRESHOLD;
            }
        } else {
            for (intptr_t i = 0; i < nptn; i++) {
                UBYTE scale_up_fwd  = max(UBYTE(0), current_it->scale_num[i]);
                UBYTE scale_up_back = max(UBYTE(0), current_it_back->scale_num[i]);
                int   scale_up      = scale_up_fwd + scale_up_back;
                ptn_lh[i] = tree_buffers._pattern_lh[i] + scale_up * LOG_SCALING_THRESHOLD;
            }
        }
    } else {
        memmove(ptn_lh, tree_buffers._pattern_lh, nptn * sizeof(double));
    }

    if (!ptn_lh_cat)
        return;

    /*
    if (ptn_lh_cat && model->isSiteSpecificModel()) {
        int offset = 0;
        if (sum_scaling == 0.0) {
            int nptncat = nptn * ncat;
            for (i = 0; i < nptncat; i++) {
                ptn_lh_cat[i] = log(_pattern_lh_cat[i]);
            }
        } else if (current_it->lh_scale_factor == 0.0) {
            for (i = 0; i < nptn; i++) {
                double scale = (max(UBYTE(0), current_it_back->scale_num[i])) * LOG_SCALING_THRESHOLD;
                for (int j = 0; j < ncat; j++, offset++)
                    ptn_lh_cat[offset] = log(_pattern_lh_cat[offset]) + scale;
            }
        } else if (current_it_back->lh_scale_factor == 0.0) {
            for (i = 0; i < nptn; i++) {
                double scale = (max(UBYTE(0), current_it->scale_num[i])) * LOG_SCALING_THRESHOLD;
                for (int j = 0; j < ncat; j++, offset++)
                    ptn_lh_cat[offset] = log(_pattern_lh_cat[offset]) + scale;
            }
        } else {
            for (i = 0; i < nptn; i++) {
                double scale = (max(UBYTE(0), current_it->scale_num[i]) +
                        max(UBYTE(0), current_it_back->scale_num[i])) * LOG_SCALING_THRESHOLD;
                for (int j = 0; j < ncat; j++, offset++)
                    ptn_lh_cat[offset] = log(_pattern_lh_cat[offset]) + scale;
            }
        }
        return;
    }
    */
    
    // New kernel
    PhyloNeighbor *nei1 = current_it;
    PhyloNeighbor *nei2 = current_it_back;
    if (!nei1->node->isLeaf() && nei2->node->isLeaf()) {
        std::swap(nei1, nei2);
    }
    if (nei1->node->isLeaf()) {
        // external branch
        double* lh_cat     = tree_buffers._pattern_lh_cat;
        double* out_lh_cat = ptn_lh_cat;
        UBYTE*  nei2_scale = nei2->scale_num;
        if (params->lk_safe_scaling || static_cast<int>(leafNum) >= params->numseq_safe_scaling) {
            // per-category scaling
            for (intptr_t ptn = 0; ptn < nptn; ptn++) {
                for (intptr_t i = 0; i < ncat; i++) {
                    out_lh_cat[i] = log(lh_cat[i]) + nei2_scale[i] * LOG_SCALING_THRESHOLD;
                }
                lh_cat += ncat;
                out_lh_cat += ncat;
                nei2_scale += ncat;
            }
        } else {
            // normal scaling
            for (intptr_t ptn = 0; ptn < nptn; ptn++) {
                double scale = nei2_scale[ptn] * LOG_SCALING_THRESHOLD;
                for (intptr_t i = 0; i < ncat; i++)
                    out_lh_cat[i] = log(lh_cat[i]) + scale;
                lh_cat += ncat;
                out_lh_cat += ncat;
            }
        }
    } else {
        // internal branch
        double* lh_cat     = tree_buffers._pattern_lh_cat;
        double* out_lh_cat = ptn_lh_cat;
        UBYTE*  nei1_scale = nei1->scale_num;
        UBYTE*  nei2_scale = nei2->scale_num;
        if (params->lk_safe_scaling || static_cast<int>(leafNum) >= params->numseq_safe_scaling) {
            // per-category scaling
            for (intptr_t ptn = 0; ptn < nptn; ptn++) {
                for (intptr_t i = 0; i < ncat; i++) {
                    out_lh_cat[i] = log(lh_cat[i]) + (nei1_scale[i]+nei2_scale[i]) * LOG_SCALING_THRESHOLD;
                }
                lh_cat += ncat;
                out_lh_cat += ncat;
                nei1_scale += ncat;
                nei2_scale += ncat;
            }
        } else {
            // normal scaling
            for (intptr_t ptn = 0; ptn < nptn; ptn++) {
                double scale = (nei1_scale[ptn] + nei2_scale[ptn]) * LOG_SCALING_THRESHOLD;
                for (intptr_t i = 0; i < ncat; i++)
                    out_lh_cat[i] = log(lh_cat[i]) + scale;
                lh_cat += ncat;
                out_lh_cat += ncat;
            }
        }
    }

//    if (cur_logl) {
//        double check_score = 0.0;
//        for (int i = 0; i < nptn; i++) {
//            check_score += (ptn_lh[i] * (aln->at(i).frequency));
//        }
//        if (fabs(check_score - *cur_logl) > 0.01) {
//            cout << *cur_logl << " " << check_score << endl;
//            assert(0);
//        }
//    }
    //double score = computeLikelihoodBranch(dad_branch, dad, pattern_lh);
    //return score;
}

void PhyloTree::computePatternProbabilityCategory(double *ptn_prob_cat, SiteLoglType wsl) {
    /*    if (!dad_branch) {
     dad = getRoot();
     dad_branch = dad->firstNeighbor();
     }*/
    intptr_t ptn, nptn = aln->getNPattern();
    size_t cat, ncat = getNumLhCat(wsl);
    // Right now only Naive version store _pattern_lh_cat!
    computePatternLhCat(wsl);

    memcpy(ptn_prob_cat, tree_buffers._pattern_lh_cat, sizeof(double)*nptn*ncat);

    for (ptn = 0; ptn < nptn; ptn++) {
        double *lh_cat = ptn_prob_cat + ptn*ncat;
        double sum = lh_cat[0];
        for (cat = 1; cat < ncat; cat++)
            sum += lh_cat[cat];
        sum = 1.0/sum;
        for (cat = 0; cat < ncat; cat++)
            lh_cat[cat] *= sum;
    }
}

int PhyloTree::computePatternCategories(IntVector *pattern_ncat) {
    if (sse != LK_386) {
        // compute _pattern_lh_cat
        computePatternLhCat(WSL_MIXTURE_RATECAT);
    }    
    intptr_t npattern = aln->getNPattern();
    int      ncat     = getRate()->getNRate();
    int      nmixture;
    if (getModel()->isMixture() && !getModelFactory()->fused_mix_rate)
        nmixture = getModel()->getNMixtures();
    else
        nmixture = ncat;
    if (pattern_ncat)
        pattern_ncat->resize(npattern);
    if (ptn_cat_mask.empty())
        ptn_cat_mask.resize(npattern, 0);
    
    size_t num_best_mixture = 0;
    ASSERT(ncat < sizeof(uint64_t)*8 && nmixture < sizeof(uint64_t)*8);

    double* lh_cat            = tree_buffers._pattern_lh_cat;
    double* lh_mixture        = new double[nmixture];
    double* sorted_lh_mixture = new double[nmixture];
    int*    id_mixture        = new int[nmixture];
    
//    for (c = 0; c < ncat; c++)
//        cat_prob[c] = getRate()->getProp(c);
//    cout << "Ptn\tFreq\tNumMix\tBestMix" << endl;
    
    size_t sum_nmix = 0;
    for (intptr_t ptn = 0; ptn < npattern; ptn++) {
        double sum_prob = 0.0, acc_prob = 0.0;
        memset(lh_mixture, 0, nmixture*sizeof(double));
        if (getModel()->isMixture() && !getModelFactory()->fused_mix_rate) {
            for (int m = 0; m < nmixture; m++) {
                for (int c = 0; c < ncat; c++) {
                    lh_mixture[m] += lh_cat[c];
                }
                sum_prob += lh_mixture[m];
                lh_cat += ncat;
                id_mixture[m] = m;
            }
        } else {
            for (int  m = 0; m < nmixture; m++) {
                lh_mixture[m] = lh_cat[m];
                sum_prob += lh_mixture[m];
                id_mixture[m] = m;
            }
            lh_cat += nmixture;
        }
        sum_prob = 1.0 / sum_prob;
        for (size_t m = 0; m < nmixture; m++) {
            lh_mixture[m] *= sum_prob;
            sorted_lh_mixture[m] = -lh_mixture[m];
        }
        quicksort(sorted_lh_mixture, 0, nmixture-1, id_mixture);

        int m;
        for ( m = 0; m < nmixture && acc_prob <= 0.99; m++) {
            acc_prob -= sorted_lh_mixture[m];
            ptn_cat_mask[ptn] |= (uint64_t)1 << id_mixture[m];
        }
        if (m > num_best_mixture) {
            num_best_mixture = m;
        }
        sum_nmix += m;
        if (pattern_ncat) {
            (*pattern_ncat)[ptn] = m;
        }
        if (verbose_mode >= VerboseMode::VB_MED) {
            cout << ptn << "\t" << (int)ptn_freq[ptn] << "\t" << m << "\t" << id_mixture[0];
            for (size_t c = 0; c < m; c++) {
                cout  << "\t" << id_mixture[c] << "\t" << -sorted_lh_mixture[c];
            }
            cout << endl;
        }
    }
//    cout << 100*(double(sum_nmix)/nmixture)/npattern << "% computation necessary" << endl;
    delete [] id_mixture;
    delete [] sorted_lh_mixture;
    delete [] lh_mixture;
    return static_cast<int>(num_best_mixture);
}

double PhyloTree::computeLogLVariance(double *ptn_lh, double tree_lh) {
    intptr_t nptn = getAlnNPattern();
    size_t nsite = getAlnNSite();
    double *pattern_lh = ptn_lh;
    if (!ptn_lh) {
        pattern_lh = new double[nptn];
        computePatternLikelihood(pattern_lh);
    }
    IntVector pattern_freq;
    aln->getPatternFreq(pattern_freq);
    if (tree_lh == 0.0) {
        for (intptr_t i = 0; i < nptn; ++i)
            tree_lh += pattern_lh[i] * pattern_freq[i];
    }
    double avg_site_lh = tree_lh / nsite;
    double variance = 0.0;
    for (intptr_t i = 0; i < nptn; ++i) {
        double diff = (pattern_lh[i] - avg_site_lh);
        variance += diff * diff * pattern_freq[i];
    }
    if (!ptn_lh) {
        delete[] pattern_lh;
    }
    if (nsite <= 1)
        return 0.0;
    return variance * ((double) nsite / (nsite - 1.0));
}

double PhyloTree::computeLogLDiffVariance(double *pattern_lh_other, double *ptn_lh) {
    intptr_t nptn = getAlnNPattern();
    size_t nsite = getAlnNSite();
    double *pattern_lh = ptn_lh;
    if (!ptn_lh) {
        pattern_lh = new double[nptn];
        computePatternLikelihood(pattern_lh);
    }
    IntVector pattern_freq;
    aln->getPatternFreq(pattern_freq);

    double avg_site_lh_diff = 0.0;
    for (intptr_t i = 0; i < nptn; ++i) {
        avg_site_lh_diff += (pattern_lh[i] - pattern_lh_other[i]) * pattern_freq[i];
    }
    avg_site_lh_diff /= nsite;
    double variance = 0.0;
    for (intptr_t i = 0; i < nptn; ++i) {
        double diff = (pattern_lh[i] - pattern_lh_other[i] - avg_site_lh_diff);
        variance += diff * diff * pattern_freq[i];
    }
    if (!ptn_lh)
        delete[] pattern_lh;
    if (nsite <= 1)
        return 0.0;
    return variance * ((double) nsite / (nsite - 1.0));
}

double PhyloTree::computeLogLDiffVariance(PhyloTree *other_tree, double *pattern_lh) {
    double *pattern_lh_other = new double[getAlnNPattern()];
    other_tree->computePatternLikelihood(pattern_lh_other);
    // BUG FIX found by Xcode analyze (use of memory after it is freed)
//    delete[] pattern_lh_other;
    double res = computeLogLDiffVariance(pattern_lh_other, pattern_lh);
    delete[] pattern_lh_other;
    return res;
}

void PhyloTree::getUnmarkedNodes(PhyloNodeVector& unmarkedNodes, PhyloNode* node, PhyloNode* dad) {
    if (!node) {
        node = getRoot();
    }
    if (markedNodeList.find(node->id) == markedNodeList.end()) {
        int numUnmarkedNei = 0;
        for (NeighborVec::iterator it = (node)->neighbors.begin(); it != (node)->neighbors.end(); it++) {
            if (markedNodeList.find((*it)->node->id) == markedNodeList.end()) {
                ++numUnmarkedNei;
            }
        }
        if (numUnmarkedNei == 1) {
            unmarkedNodes.push_back(node);
        }
    }
    FOR_EACH_ADJACENT_PHYLO_NODE(node, dad, it, child) {
        getUnmarkedNodes(unmarkedNodes, child, node);
    }
}

double PhyloTree::optimizeOneBranchLS(PhyloNode *node1, PhyloNode *node2) {
    if (!subTreeDistComputed) {
        if (params->ls_var_type == WLS_PAUPLIN) {
            computeNodeBranchDists();
        }
        computeSubtreeDists();
    }
    double A, B, C, D;
    A = B = C = D = 0;
    PhyloNode *nodeA = NULL, *nodeB = NULL, *nodeC = NULL, *nodeD = NULL;
    double lsBranch;

    // One of the node is a leaf
    if (node1->isLeaf() || node2->isLeaf()) {
        if (node1->isLeaf()) {
            // nodeA and nodeB are children of node2
            FOR_EACH_ADJACENT_PHYLO_NODE(node2, node1, it, node_A_or_B){
                if (A == 0) {
                    A = getNumTaxa(node_A_or_B, node2);
                    nodeA = node_A_or_B;
                } else {
                    B = getNumTaxa(node_A_or_B, node2);
                    nodeB = node_A_or_B;
                }
            }
            // nodeC is now node1
            nodeC = node1;
        } else {
            // nodeA and nodeB are children of node1
            FOR_EACH_ADJACENT_PHYLO_NODE(node1, node2, it, node_A_or_B) {
                if (A == 0) {
                    A = getNumTaxa(node_A_or_B, node1);
                    nodeA = node_A_or_B;
                } else {
                    B = getNumTaxa(node_A_or_B, node1);
                    nodeB = node_A_or_B;
                }
            }
            // nodeC is now node1
            nodeC = node2;
        }
        ASSERT(A != 0);
        ASSERT(B != 0);
        string keyAC = getBranchID(nodeA, nodeC);
        ASSERT(subTreeDists.count(keyAC));
        double distAC = subTreeDists[keyAC];
        double weightAC = subTreeWeights[keyAC];
        string keyBC = getBranchID(nodeB, nodeC);
        ASSERT(subTreeDists.count(keyBC));
        double distBC = subTreeDists[keyBC];
        double weightBC = subTreeWeights[keyBC];
        string keyAB = getBranchID(nodeA, nodeB);
        ASSERT(subTreeDists.count(keyAB));
        double distAB = subTreeDists[keyAB];
        double weightAB = subTreeWeights[keyAB];
        if (params->ls_var_type == OLS/* || params->ls_var_type == FIRST_TAYLOR || params->ls_var_type == FITCH_MARGOLIASH
                || params->ls_var_type == SECOND_TAYLOR*/) {
            lsBranch = 0.5 * (distAC / A + distBC / B - distAB / (A * B));
        } /*else if (params->ls_var_type == PAUPLIN) {
            // TODO: Chua test bao gio
            outError("Paulin formula not supported yet");
            lsBranch = 0.5 * (distAC + distBC) - 0.5 * distAB;
        }*/ else {
            // weighted least square
            lsBranch = 0.5*(distAC/weightAC + distBC/weightBC - distAB/weightAB);
        }
    } else { // Both node are internal node
        FOR_EACH_ADJACENT_PHYLO_NODE(node1, node2, it, node_A_or_B) {
            if (A == 0) {
                A = getNumTaxa(node_A_or_B, node1);
                nodeA = node_A_or_B;
            } else {
                B = getNumTaxa(node_A_or_B, node1);
                nodeB = node_A_or_B;
            }
        }

        FOR_EACH_ADJACENT_PHYLO_NODE(node2, node1, it, node_A_or_B) {
            if (C == 0) {
                C = getNumTaxa(node_A_or_B, node2);
                nodeC = node_A_or_B;
            } else {
                D = getNumTaxa(node_A_or_B, node2);
                nodeD = node_A_or_B;
            }
        }

        string keyAC = getBranchID(nodeA, nodeC);
        ASSERT(subTreeDists.count(keyAC));
        double distAC = subTreeDists[keyAC];
        double weightAC = subTreeWeights[keyAC];

        string keyBD = getBranchID(nodeB, nodeD);
        ASSERT(subTreeDists.count(keyBD));
        double distBD = subTreeDists[keyBD];
        double weightBD = subTreeWeights[keyBD];

        string keyBC = getBranchID(nodeB, nodeC);
        ASSERT(subTreeDists.count(keyBC));
        double distBC = subTreeDists[keyBC];
        double weightBC = subTreeWeights[keyBC];

        string keyAD = getBranchID(nodeA, nodeD);
        ASSERT(subTreeDists.count(keyAD));
        double distAD = subTreeDists[keyAD];
        double weightAD = subTreeWeights[keyAD];

        string keyAB = getBranchID(nodeA, nodeB);
        ASSERT(subTreeDists.count(keyAB));
        double distAB = subTreeDists[keyAB];
        double weightAB = subTreeWeights[keyAB];

        string keyCD = getBranchID(nodeC, nodeD);
        ASSERT(subTreeDists.count(keyCD));
        double distCD = subTreeDists[keyCD];
        double weightCD = subTreeWeights[keyCD];

        /*if (params->ls_var_type == PAUPLIN) {
            // this distance has a typo as also seen in Mihaescu & Pachter 2008
            //lsBranch = 0.25 * (distAC + distBD + distAD + distBC) - 0.5 * (distAB - distCD);
            outError("Paulin formula not supported yet");
            lsBranch = 0.25 * (distAC + distBD + distAD + distBC) - 0.5 * (distAB + distCD);
        } else*/ if (params->ls_var_type == OLS) {
            double gamma = (B * C + A * D) / ((A + B)*(C + D));
            lsBranch = 0.5 * (gamma * (distAC / (A * C) + distBD / (B * D))
                    + (1 - gamma) * (distBC / (B * C) + distAD / (A * D))
                    - distAB / (A * B) - distCD / (C * D));
        } else {
            // weighted least square
            double K = 1.0/weightAC + 1.0/weightBD + 1.0/weightAD + 1.0/weightBC;
            lsBranch =
                    ((distAC/weightAC+distBD/weightBD)*(weightAD+weightBC)/(weightAD*weightBC)+
                    (distAD/weightAD+distBC/weightBC)*(weightAC+weightBD)/(weightAC*weightBD))/K
                    - distAB/weightAB - distCD/weightCD;
            lsBranch = 0.5*lsBranch;
        }
    }
    return lsBranch;
}

void PhyloTree::updateSubtreeDists(const NNIMove &nnimove) {
    ASSERT(subTreeDistComputed);
    PhyloNode *nodeA = NULL, *nodeB = NULL, *nodeC = NULL, *nodeD = NULL;
    PhyloNode *node1 = nnimove.node1;
    PhyloNode *node2 = nnimove.node2;
    NeighborVec::iterator node1Nei_it = nnimove.node1Nei_it;
    NeighborVec::iterator node2Nei_it = nnimove.node2Nei_it;
    Neighbor *node1Nei = *(node1Nei_it);
    Neighbor *node2Nei = *(node2Nei_it);

    // ((A,C),(B,D))
    // C and D are the 2 subtree that get swapped
    FOR_EACH_ADJACENT_PHYLO_NODE(node1, node2, it, node_A_or_C) {
        if ((*it)->id != node1Nei->id) {
            nodeA = node_A_or_C;
        } else {
            nodeC = node_A_or_C;
        }
    }

    ASSERT(nodeA);
    ASSERT(nodeC);

    FOR_EACH_ADJACENT_PHYLO_NODE(node2, node1, it, node_B_or_D) {
        if ((*it)->id != node2Nei->id) {
            nodeB = node_B_or_D;
        } else {
            nodeD = node_B_or_D;
        }
    }

    ASSERT(nodeB);
    ASSERT(nodeD);

    NodeVector nodeListA, nodeListB, nodeListC, nodeListD;
    getAllNodesInSubtree(nodeA, node1, nodeListA);
    getAllNodesInSubtree(nodeC, node1, nodeListC);
    getAllNodesInSubtree(nodeB, node2, nodeListB);
    getAllNodesInSubtree(nodeD, node2, nodeListD);

    for (NodeVector::iterator it = nodeListA.begin(); it != nodeListA.end(); ++it) {
        string key = getBranchID((*it), node2);
        double distB = subTreeDists.find(getBranchID((*it), nodeB))->second;
        double distD = subTreeDists.find(getBranchID((*it), nodeD))->second;
        double newDist = distB + distD;
        StringDoubleMap::iterator dist_it = subTreeDists.find(key);
        ASSERT(dist_it != subTreeDists.end());
        dist_it->second = newDist;
    }

    for (NodeVector::iterator it = nodeListB.begin(); it != nodeListB.end(); ++it) {
        string key = getBranchID((*it), node1);
        double distC = subTreeDists.find(getBranchID((*it), nodeC))->second;
        double distA = subTreeDists.find(getBranchID((*it), nodeA))->second;
        double newDist = distC + distA;
        StringDoubleMap::iterator dist_it = subTreeDists.find(key);
        ASSERT(dist_it != subTreeDists.end());
        dist_it->second = newDist;
    }

    for (NodeVector::iterator it = nodeListC.begin(); it != nodeListC.end(); ++it) {
        string key = getBranchID((*it), node2);
        double distD = subTreeDists.find(getBranchID((*it), nodeD))->second;
        double distB = subTreeDists.find(getBranchID((*it), nodeB))->second;
        double newDist = distD + distB;
        StringDoubleMap::iterator dist_it = subTreeDists.find(key);
        ASSERT(dist_it != subTreeDists.end());
        dist_it->second = newDist;
    }

    for (NodeVector::iterator it = nodeListD.begin(); it != nodeListD.end(); ++it) {
        string key = getBranchID((*it), node1);
        double distA = subTreeDists.find(getBranchID((*it), nodeA))->second;
        double distC = subTreeDists.find(getBranchID((*it), nodeC))->second;
        double newDist = distA + distC;
        StringDoubleMap::iterator dist_it = subTreeDists.find(key);
        ASSERT(dist_it != subTreeDists.end());
        dist_it->second = newDist;
    }

    double distAB = subTreeDists.find(getBranchID(nodeA, nodeB))->second;
    double distAD = subTreeDists.find(getBranchID(nodeA, nodeD))->second;
    double distCB = subTreeDists.find(getBranchID(nodeC, nodeB))->second;
    double distCD = subTreeDists.find(getBranchID(nodeC, nodeD))->second;

    subTreeDists.find(getBranchID(node1, node2))->second = distAB + distAD + distCB + distCD;
}

void PhyloTree::computeSubtreeDists() {
    PhyloNodeVector unmarkedNodes;
    subTreeDists.clear();
    subTreeWeights.clear();
    do {
        // Generate a list of unmarked node that is adjacent to exactly one unmarked nodes
        // Here we will work up the tree in a bottom up manner
        unmarkedNodes.clear();
        getUnmarkedNodes(unmarkedNodes);
        if (unmarkedNodes.size() == 0)
            break;

        for (auto it = unmarkedNodes.begin(); it != unmarkedNodes.end(); ++it) {
            // if the node is an internal node then all of its child nodes should be marked
            // source_nei1 and source_nei2 are the 2 marked child node
            // nextNode is the other node, used for traversal
            PhyloNode* source_nei1 = NULL;
            PhyloNode* source_nei2 = NULL;
            PhyloNode* nextNode;
            if (!(*it)->isLeaf()) {
                // select the 2 marked child nodes
                for (NeighborVec::iterator it2 = (*it)->neighbors.begin(); it2 != (*it)->neighbors.end(); ++it2) {
                    PhyloNode* child_node = (PhyloNode*)(*it2)->node;
                    if (markedNodeList.find((*it2)->node->id) != markedNodeList.end()) {
                        if (!source_nei1) {
                            source_nei1 = child_node;
                        } else {
                            source_nei2 = child_node;
                        }
                    } else {
                        nextNode = child_node;
                    }
                }
                ASSERT(source_nei1);
                ASSERT(source_nei2);
            } else {
                nextNode = (*it)->firstNeighbor()->getNode();
            }
            // warning: 'nextNode' may be used uninitialized in this function
            computeAllSubtreeDistForOneNode((*it), source_nei1, source_nei2, (*it), nextNode);
            markedNodeList.insert(IntPhyloNodeMap::value_type((*it)->id, (*it)));
        }
    } while (true);
    markedNodeList.clear();
    subTreeDistComputed = true;
}

void PhyloTree::computeAllSubtreeDistForOneNode(PhyloNode* source, PhyloNode* source_nei1,
                                                PhyloNode* source_nei2,
                                                PhyloNode* node, PhyloNode* dad) {
    string key = getBranchID(source, dad);
    double dist, weight;
    if (markedNodeList.find(dad->id) != markedNodeList.end()) {
        return;
    } else if (source->isLeaf() && dad->isLeaf()) {
        ASSERT(dist_matrix);
        size_t nseq = aln->getNSeq();
        dist        = dist_matrix[dad->id * nseq + source->id];
        weight      = 1.0;
        if (params->ls_var_type != OLS && 0.0 != dist) {
            // this will take into account variances, also work for OLS since var = 1
            weight = 1.0/dist/dist;
            dist  *= weight;
        }
        subTreeDists.insert(StringDoubleMap::value_type(key, dist));
        subTreeWeights.insert(StringDoubleMap::value_type(key, weight));
    } else if (!source->isLeaf() && dad->isLeaf()) {
        ASSERT(source_nei1);
        ASSERT(source_nei2);
        string key1 = getBranchID(source_nei1, dad);
        ASSERT(subTreeDists.find(key1) == subTreeDists.end());
        double dist1 = subTreeDists.find(key1)->second;
        double weight1 = subTreeWeights.find(key1)->second;
        string key2 = getBranchID(source_nei2, dad);
        ASSERT(subTreeDists.find(key2) == subTreeDists.end());
        double dist2 = subTreeDists.find(key2)->second;
        double weight2 = subTreeWeights.find(key2)->second;
        dist = dist1 + dist2;
        weight = weight1 + weight2;
        subTreeDists.insert(StringDoubleMap::value_type(key, dist));
        subTreeWeights.insert(StringDoubleMap::value_type(key, weight));
    } else {
        PhyloNode* dad_nei1 = nullptr;
        PhyloNode* dad_nei2 = nullptr;
        FOR_EACH_ADJACENT_PHYLO_NODE(dad, node, it, grandma) {
            if (dad_nei1==nullptr) {
                dad_nei1 = grandma;
            } else {
                dad_nei2 = grandma;
            }
        }
        ASSERT(dad_nei1);
        ASSERT(dad_nei2);
        computeAllSubtreeDistForOneNode(source, source_nei1, source_nei2, dad, dad_nei1);
        computeAllSubtreeDistForOneNode(source, source_nei1, source_nei2, dad, dad_nei2);
        string key1 = getBranchID(source, dad_nei1);
        string key2 = getBranchID(source, dad_nei2);
        ASSERT(subTreeDists.find(key1) != subTreeDists.end());
        ASSERT(subTreeDists.find(key2) != subTreeDists.end());
        double dist1 = subTreeDists.find(key1)->second;
        double weight1 = subTreeWeights.find(key1)->second;
        double dist2 = subTreeDists.find(key2)->second;
        double weight2 = subTreeWeights.find(key2)->second;
        dist = dist1 + dist2;
        weight = weight1 + weight2;
        subTreeDists.insert(StringDoubleMap::value_type(key, dist));
        subTreeWeights.insert(StringDoubleMap::value_type(key, weight));
    }
}

set<int> PhyloTree::computeNodeBranchDists(Node *node, Node *dad) {
    set<int>::iterator i, j;
    if (!nodeBranchDists) {
        cout << "nodeNum = " << nodeNum << endl;
        nodeBranchDists = new int[nodeNum*nodeNum];
    }
    if (!node) {
        memset(nodeBranchDists, 0, sizeof(int)*nodeNum*nodeNum);
        ASSERT(root->isLeaf());
        dad = root;
        node = dad->neighbors[0]->node;
        set<int> res = computeNodeBranchDists(node, dad);
        for (i = res.begin(); i != res.end(); i++)
            nodeBranchDists[(*i)*nodeNum + dad->id] = nodeBranchDists[(dad->id)*nodeNum + (*i)] =
                nodeBranchDists[(*i)*nodeNum + node->id] + 1;
        // sanity check that all distances are filled
        for (int x = 0; x < nodeNum; x++)
            for (int y = 0; y < nodeNum; y++)
                if (x != y)
                    ASSERT(nodeBranchDists[x*nodeNum+y] != 0);
                else
                    ASSERT(nodeBranchDists[x*nodeNum+y] == 0);
        return res;
    }
    if (node->isLeaf()) {
        set<int> res;
        res.insert(node->id);
        return res;
    }
    ASSERT(node->degree() == 3);
    Node *left = NULL, *right = NULL;
    FOR_NEIGHBOR_IT(node, dad, it) {
        if (!left) left = (*it)->node; else right = (*it)->node;
    }
    ASSERT(left!=nullptr);
    ASSERT(right!=nullptr);
    if (left==nullptr || right==nullptr) {
        set<int> res;
        res.insert(node->id);
        return res;
    }
    set<int> resl = computeNodeBranchDists(left, node);
    set<int> resr = computeNodeBranchDists(right, node);
    for (i = resl.begin(); i != resl.end(); i++)
        nodeBranchDists[(*i)*nodeNum + node->id] = nodeBranchDists[(node->id)*nodeNum + (*i)] =
            nodeBranchDists[(*i)*nodeNum + left->id] + 1;
    for (i = resr.begin(); i != resr.end(); i++)
        nodeBranchDists[(*i)*nodeNum + node->id] = nodeBranchDists[(node->id)*nodeNum + (*i)] =
            nodeBranchDists[(*i)*nodeNum + right->id] + 1;
    for (i = resl.begin(); i != resl.end(); i++)
        for (j = resr.begin(); j != resr.end(); j++)
            nodeBranchDists[(*i)*nodeNum + (*j)] = nodeBranchDists[(*j)*nodeNum+(*i)] =
                nodeBranchDists[(*i)*nodeNum+node->id]+nodeBranchDists[(*j)*nodeNum+node->id];
    resl.insert(resr.begin(), resr.end());
    resl.insert(node->id);
    return resl;
}

/*
    b0: initial guess for the maximum
*/
double PhyloTree::approxOneBranch(PhyloNode *node, PhyloNode *dad, double b0) {
    double b_max, ddl, b1, b2, std, seqlen;
    double t1, t3, t5, t11, t18, t21, t26, t29, t30, t32, t44, t46, t48;
    double beps = 1/DBL_MAX;

    /* TODO: insert call to get sequence length */
    seqlen = static_cast<double>(getAlnNSite());

    /* use a robust first order approximation to the variance */
    std = sqrt(b0/seqlen);

    /* determine neighbour points */
    b1 = b0 - std;
    if (b1<=0) b1 = beps; /* only happens for b<=1 with small seq. len. */
    b2 = b0 + std;

    /* TODO: insert calls to log-likelihood function */
    PhyloNeighbor *dad_nei  = dad->findNeighbor(node);
    PhyloNeighbor *node_nei = node->findNeighbor(dad);
    double old_len = dad_nei->length;
    dad_nei->length = node_nei->length = b0;
    double l0 = computeLikelihoodBranch(dad_nei, dad, tree_buffers);
    dad_nei->length = node_nei->length = b1;
    double l1 = computeLikelihoodBranch(dad_nei, dad, tree_buffers);
    dad_nei->length = node_nei->length = b2;
    double l2 = computeLikelihoodBranch(dad_nei, dad, tree_buffers);
    dad_nei->length = node_nei->length = old_len;

    t1 = sqrt(b0);
    t3 = sqrt(b2);
    t5 = sqrt(b1);
    t11 = pow(-t1*l2+t3*l0+t5*l2+t1*l1-t5*l0-t3*l1,2.0);
    t18 = -b0*l2+b2*l0+b1*l2+b0*l1-b1*l0-b2*l1;
    t21 = t1-t5;
    t26 = -t1*t3+t1*t5+b2-t5*t3;
    t29 = t18*t18;
    t30 = 1/t11;
    t32 = sqrt(t29*t30);
    ddl = -2.0*t11/t18/t21/t26/t32;

    if (ddl > 0) {
        /* the analytic extremum is a minimum,
           so the maximum is at the lower bound */
        b_max = 0;
    } else {
        t44 = pow(-t1*b2+t5*b2-t5*b0+t3*b0-t3*b1+t1*b1,2.0);
        t46 = t21*t21;
        t48 = t26*t26;
        b_max = t29*t44/t46/t48*t30/4.0;
    }

    return(b_max);
}

void PhyloTree::approxAllBranches(PhyloNode *node, PhyloNode *dad) {
    if (!node) {
        node = getRoot();
    }

    if (dad) {
        PhyloNeighbor *node_dad_nei = node->findNeighbor(dad);
        PhyloNeighbor *dad_node_nei = dad->findNeighbor(node);
        double len = approxOneBranch(node, dad, dad_node_nei->length);
        node_dad_nei->length = len;
        dad_node_nei->length = len;
    }
    FOR_EACH_ADJACENT_PHYLO_NODE(node, dad, it, child) {
        approxAllBranches(child, node);
    }
}

/*
 void PhyloTree::computeAllSubtreeDists(PhyloNode* node, PhyloNode* dad) {
 if (!node) {
    node = getRoot();
 }

 if (dad) {
 // This function compute all pairwise subtree distance between subtree rooted at dad and others

 computeSubtreeDists(node, dad);
 }

 FOR_EACH_ADJACENT_PHYLO_NODE(node, dad, it, child) {

 computeAllSubtreeDists(child, node);
 }
 }

 void PhyloTree::computeSubtreeDists(PhyloNode* node, PhyloNode* dad) {
 // if both nodes are leaf then it is trivial, just retrieve the values from the distance matrix
 if (dad->isLeaf() && node->isLeaf()) {
 string key = nodePair2String(dad, node);
 assert(dist_matrix);
 size_t nseq = aln->getNSeq();
 double dist = dist_matrix[dad->id * nseq + node->id];
 interSubtreeDistances.insert(StringDoubleMap::value_type(key, dist));
 } else if (!dad->isLeaf() && node->isLeaf()) {

 FOR_EACH_ADJACENT_PHYLO_NODE(node, dad, it, child) {

 computeSubtreeDists(dad, child); //Todo: This looks wrong to me. -James B.
 }

 }
 }
 */

double PhyloTree::computeBayesianBranchLength(PhyloNeighbor*          dad_branch,
                                              PhyloNode*              dad) {
    double         obsLen      = 0.0;
    PhyloNode*     node        = dad_branch->getNode();
    PhyloNeighbor* node_branch = node->findNeighbor(dad);
    
    //PlacementTraversalInfo pti(*this, tree_buffers, nullptr, nullptr);
    if (dad->isLeaf()) {
        std::swap(dad,        node);
        std::swap(dad_branch, node_branch);
    }
    computeLikelihoodBranch(dad_branch, dad, tree_buffers);
    //pti.computePartialLikelihood(dad_branch,  dad);
    
    // now combine likelihood at the branch
    const int        nstates = aln->num_states;
    const int        numCat  = site_rate->getNRate();
    const size_t     line    = nstates * vector_size;
    const size_t     block   = numCat  * line;
    const intptr_t   nptn    = aln->size();
    //const intptr_t max_orig_nptn = roundUpToMultiple(nptn, vector_size);
    //const size_t   tip_mem_size  = max_orig_nptn * nstates;
    //const bool     site_model    = model_factory && model_factory->model &&
    //                               model_factory->model->isSiteSpecificModel();

    double* tmp_state_freq            = new double[nstates];
    double* tmp_ancestral_state_prob1 = new double[nstates];
    double* tmp_ancestral_state_prob2 = new double[nstates];

    model->getStateFrequency(tmp_state_freq);

    for (intptr_t ptn = 0; ptn < nptn; ++ptn) {
        // Compute the probability of each state for the current site
        double  sum_prob1        = 0.0;
        double  sum_prob2        = 0.0;
        size_t  offset           = (ptn / vector_size) * block +
                                   (ptn % vector_size);
        double* partial_lh_site  = node_branch->partial_lh + (offset);
        double* partial_lh_child = dad_branch->partial_lh  + (offset);
        size_t  leaf_state = node->isLeaf() ? aln->at(ptn)[node->id] : 0;
        double  no_match   = (leaf_state == aln->STATE_UNKNOWN) ? ( 1.0 / nstates ) : 0.0;
        
        for (size_t state = 0; state < nstates; state++) {
            tmp_ancestral_state_prob1[state] = 0.0;
            tmp_ancestral_state_prob2[state] = 0.0;
            if (node->isInterior()) {
                for (size_t cat_offset = 0; cat_offset < block; cat_offset += line) {
                    tmp_ancestral_state_prob1[state]
                        += partial_lh_site [cat_offset + state * vector_size];
                    tmp_ancestral_state_prob2[state]
                        += partial_lh_child[cat_offset + state * vector_size];
                }
            }
            else {
                //double* tip_partial_lh_child = (site_model)
                //   ? &tip_partial_lh[node->id*tip_mem_size + ptn*nstates]
                //   : pti.partial_lh_leaves;
                for (size_t cat_offset = 0; cat_offset < block; cat_offset += line) {
                    tmp_ancestral_state_prob1[state]
                        += partial_lh_site [cat_offset + state * vector_size];
                }
                tmp_ancestral_state_prob2[state] = (leaf_state == state) ? 1.0 : no_match;
            }
            tmp_ancestral_state_prob1[state] *= tmp_state_freq[state];
            tmp_ancestral_state_prob2[state] *= tmp_state_freq[state];
            sum_prob1 += tmp_ancestral_state_prob1[state];
            sum_prob2 += tmp_ancestral_state_prob2[state];
        }
        bool sameState = false;
        int  state1 = 0;
        int  state2 = 0;
        double cutoff = 1.0/nstates;
        for (int state = 0; state < nstates; state++) {
            tmp_ancestral_state_prob1[state] /= sum_prob1;
            tmp_ancestral_state_prob2[state] /= sum_prob2;
            if (tmp_ancestral_state_prob1[state] > tmp_ancestral_state_prob1[state1])
                state1 = state;
            if (tmp_ancestral_state_prob2[state] > tmp_ancestral_state_prob2[state2])
                state2 = state;
            if (tmp_ancestral_state_prob1[state] > cutoff &&
                tmp_ancestral_state_prob2[state] > cutoff) {
                sameState = true;
            }
        }
        sameState = sameState || (state1 == state2);
        if (!sameState) {
            obsLen += aln->at(ptn).frequency;
        }
    }
    obsLen /= getAlnNSite();
    if (obsLen < params->min_branch_length) {
        obsLen = params->min_branch_length;
    }
    delete[] tmp_ancestral_state_prob2;
    delete[] tmp_ancestral_state_prob1;
    delete[] tmp_state_freq;

#if (0)
    LOG_LINE(VerboseMode::VB_MIN, "Was " << node_branch->length << ","
             << " now " << obsLen
             << " node was leaf " << node->isLeaf() );
#endif
    return obsLen;
}

double PhyloTree::correctBranchLengthF81(double observedBran, double alpha) const {
    if (!model) {
        return jukesCantorCorrection(observedBran, alpha);
    }
    double H = 0.0;
    double correctedBranLen;
    for (int i = 0; i < model->num_states; i++) {
        H += model->state_freq[i] * (1 - model->state_freq[i]);
    }
    observedBran = 1.0 - observedBran / H;
    // no gamma
    if (observedBran <= 0.0) {
        return params->max_branch_length;
    }
    if (alpha <= 0.0) {
        correctedBranLen = -H * log(observedBran);
    } else {
        //if (verbose_mode >= VerboseMode::VB_MAX) cout << "alpha: " << alpha << endl;
        correctedBranLen = H * alpha * (pow(observedBran, -1 / alpha) - 1);
    }
    // Branch lengths under PoMo are #events, which is ~N^2 * #substitutions
    if (aln->seq_type == SeqType::SEQ_POMO) {
        correctedBranLen *= aln->virtual_pop_size * aln->virtual_pop_size;
    }
    if (correctedBranLen < params->min_branch_length) {
        correctedBranLen = params->min_branch_length;
    }
    if (correctedBranLen > params->max_branch_length) {
        correctedBranLen = params->max_branch_length;
    }
    return correctedBranLen;
}

double PhyloTree::computeCorrectedBayesianBranchLength(PhyloNeighbor* dad_branch,
                                                       PhyloNode*     dad) {
    double observedBran = computeBayesianBranchLength(dad_branch, dad);
    return correctBranchLengthF81(observedBran, site_rate->getGammaShape());
}

void PhyloTree::computeAllBayesianBranchLengths() {
    if (!progress_display::getProgressDisplay()) {
        cout << "Computing Bayesian branch lengths..." << endl;
    }
    PhyloNodeVector nodes, nodes2;
    computeBestTraversal(nodes, nodes2);
    initProgress((double)nodes.size(),
                 "Computing bayesian branch lengths", "", "", true);
    DoubleVector lengths;
    for (int j = 0; j < nodes.size(); j++) {
        PhyloNeighbor* nei = nodes[j]->findNeighbor(nodes2[j]);
        double bayesian_length = computeBayesianBranchLength(nei, nodes[j]);
        lengths.push_back ( bayesian_length );
        if ( (j % 100) == 99) {
            trackProgress(100.0);
        }
    }
    for (int j = 0; j < nodes.size(); j++) {
        PhyloNeighbor* nei = nodes[j]->findNeighbor(nodes2[j]);
        PhyloNeighbor* back_nei = nodes2[j]->findNeighbor(nodes[j]);
        nei->length = back_nei->length = lengths[j];
    }
    trackProgress(static_cast<double>(nodes.size() % 100));
    doneProgress();
}

double PhyloTree::computeLikelihoodZeroBranch(PhyloNeighbor *dad_branch, PhyloNode *dad) {
    double lh_zero_branch;
    double         saved_len   = dad_branch->length;
    PhyloNeighbor* node_branch = dad_branch->getNode()->findNeighbor(dad);
    dad_branch->length         = 0.0;
    node_branch->length        = 0.0;
    lh_zero_branch      = computeLikelihoodBranch(dad_branch, dad, tree_buffers);
    // restore branch length
    dad_branch->length  = saved_len;
    node_branch->length = saved_len;

    return lh_zero_branch;
}


/****************************************************************************
 Branch length optimization by maximum likelihood
 ****************************************************************************/

const double TOL_TREE_LENGTH_SCALE = 0.001;

double PhyloTree::optimizeTreeLengthScaling(double min_scaling, double &scaling, double max_scaling, double gradient_epsilon) {
    is_opt_scaling = true;
    current_scaling = scaling;
    double negative_lh, ferror;
    // 2018-08-20: make sure that max and min branch lengths do not go over bounds
    vector<DoubleVector> brlens;
    brlens.resize(branchNum);
    getBranchLengths(brlens);
    double min_brlen = params->max_branch_length;
    double max_brlen = params->min_branch_length;
    for (auto brlenvec = brlens.begin(); brlenvec != brlens.end(); brlenvec++) {
        for (auto brlen = brlenvec->begin(); brlen != brlenvec->end(); brlen++) {
            max_brlen = max(max_brlen, *brlen);
            min_brlen = min(min_brlen, *brlen);
        }
    }
    
    if (min_brlen <= 0.0) min_brlen = params->min_branch_length;
    if (max_scaling > 10.0*params->max_branch_length / max_brlen)
        max_scaling = 10.0*params->max_branch_length / max_brlen;
    if (min_scaling < 0.1*params->min_branch_length / min_brlen)
        min_scaling = 0.1*params->min_branch_length / min_brlen;
    
    scaling = minimizeOneDimen(min(scaling, min_scaling), scaling, max(max_scaling, scaling), max(TOL_TREE_LENGTH_SCALE, gradient_epsilon), &negative_lh, &ferror);
    if (scaling != current_scaling) {
        scaleLength(scaling / current_scaling);
        current_scaling = scaling;
        clearAllPartialLH();
    }
    is_opt_scaling = false;
    return computeLikelihood();
}

void PhyloTree::printTreeLengthScaling(const char *filename) {
//    double treescale = 1.0;
//    
//    cout << "Optimizing tree length scaling ..." << endl;
//    
//    double lh = optimizeTreeLengthScaling(MIN_BRLEN_SCALE, treescale, MAX_BRLEN_SCALE, 0.001);
//    
//    cout << "treescale: " << treescale << " / LogL: " << lh << endl;
    
    Checkpoint *saved_checkpoint = getModelFactory()->getCheckpoint();
    Checkpoint *new_checkpoint = new Checkpoint;
    new_checkpoint->setFileName(filename);
    new_checkpoint->setCompression(false);
    new_checkpoint->setHeader("IQ-TREE scaled tree length and model parameters");
    new_checkpoint->put("treelength", treeLength());
    saved_checkpoint->put("treelength", treeLength()); // also put treelength into current checkpoint
    
    getModelFactory()->setCheckpoint(new_checkpoint);    
    getModelFactory()->saveCheckpoint();
    new_checkpoint->dump();
    
    getModelFactory()->setCheckpoint(saved_checkpoint);
}

double PhyloTree::computeFunction(double value) {
    if (!is_opt_scaling) {
        current_it->length = value;
        current_it_back->length = value;
        return -computeLikelihoodBranch(current_it, current_it_back->getNode(), tree_buffers);
    } else {
        if (value != current_scaling) {
            scaleLength(value / current_scaling);
            current_scaling = value;
            clearAllPartialLH();
        }
        return -computeLikelihood();
    }
}

void PhyloTree::computeFuncDerv(double value, double &df, double &ddf) {
    current_it->length      = value;
    current_it_back->length = value;
    computeLikelihoodDerv(current_it, current_it_back->getNode(), &df, &ddf, tree_buffers);
    df  = -df;
    ddf = -ddf;
}

void PhyloTree::optimizePatternRates(DoubleVector &pattern_rates) {
    intptr_t nptn = aln->getNPattern();
    pattern_rates.resize(nptn, 1.0);
#pragma omp parallel for
    for (intptr_t ptn = 0; ptn < nptn; ptn++) {
        Alignment *paln = new Alignment;
        IntVector ptn_id;
        ptn_id.push_back(static_cast<int>(ptn));
        paln->extractPatterns(aln, ptn_id);
        PhyloTree *tree = new PhyloTree;
        tree->copyPhyloTree(this, false); //Local alignment, so tree can't "borrow" the summary of this
        tree->setParams(params);
        tree->setAlignment(paln);
        tree->prepareToComputeDistances();
        tree->sse = sse;
        tree->setNumThreads(1);
        // initialize model
        tree->setModelFactory(getModelFactory());

        // main optimization
        tree->optimizeTreeLengthScaling(MIN_SITE_RATE, pattern_rates[ptn], MAX_SITE_RATE, 0.0001);
        
        tree->setModelFactory(NULL);
        tree->doneComputingDistances();
        delete tree;
        delete paln;
    }
}

int PhyloTree::getNBranchParameters(int brlen_type) {
    if (params->fixed_branch_length || brlen_type == BRLEN_FIX)
        return 0;

    int df = 0;

    if (brlen_type == BRLEN_OPTIMIZE) {
        df = branchNum - (int)rooted;
    // If model is Lie-Markov, and is in fact time reversible, one of the
    // degrees of freedom is illusary. (Of the two edges coming from the
    // root, only sum of their lenghts affects likelihood.)
    // So correct for this. Without this correction, K2P and RY2.2b
    // would not be synonymous, for example.

//    string className(typeid(*model).name());
//    if (className.find("ModelLieMarkov")!=string::npos && model->isReversible())
//        df--;

        // BQM 2017-04-28, alternatively, check if there is a virtual_root and model is reversible
        if (rooted && model && model->isReversible())
            df--;

    } else if (brlen_type == BRLEN_SCALE)
        df = 1;
    return df;
}

void PhyloTree::optimizeOneBranch(PhyloNode *node1, PhyloNode *node2,
                                    bool clearLH, int maxNRStep) {
    if (rooted && (node1 == root || node2 == root)) {
        return; // does not optimize virtual branch from root
    }
    current_it      = node1->findNeighbor(node2);
    current_it_back = node2->findNeighbor(node1);

    ASSERT(current_it);
    ASSERT(current_it_back);

    double original_len = current_it->length;
    ASSERT(original_len >= 0.0);
    tree_buffers.theta_computed = false;
    double original_lh = computeLikelihoodBranch(current_it, node1, tree_buffers);
    double new_len;
    if (optimize_by_newton) {
        // Newton-Raphson method
        double derivative_of_likelihood_wrt_length = 0;
        new_len = minimizeNewton(params->min_branch_length, original_len,
            params->max_branch_length, params->min_branch_length,
            derivative_of_likelihood_wrt_length, maxNRStep);
    }
    else {
        // Brent method
        double negative_lh = 0;
        double ferror;
        new_len = minimizeOneDimen(params->min_branch_length, original_len,
            params->max_branch_length, params->min_branch_length,
            &negative_lh, &ferror);
    }
    current_it->length      = new_len;
    current_it_back->length = new_len;
    curScore = computeLikelihoodFromBuffer();
    if (curScore != original_lh) {
        LOG_LINE(VerboseMode::VB_MAX, "  branch=" << current_it->id
            << ", old_len=" << original_len << ", new_len=" << new_len
            << ", old_lh=" << original_lh << ", new_lh=" << curScore
            << ", delta=" << (curScore - original_lh));
    }
    if (optimize_by_newton && new_len > params->max_branch_length*0.95 && !isSuperTree()) {
        if (original_lh > curScore) {
            current_it->length      = original_len;
            current_it_back->length = original_len;
            new_len                 = original_len;
            curScore                = original_lh;
        }
    }
    if (clearLH && original_len != new_len) {
        node1->clearReversePartialLh(node2);
        node2->clearReversePartialLh(node1);
    }
}

double PhyloTree::optimizeChildBranches(PhyloNode *node, PhyloNode *dad) {
    FOR_EACH_ADJACENT_PHYLO_NODE(node, dad, it, child) {
        optimizeOneBranch(node, child);
    }
    return computeLikelihoodFromBuffer();
}

void PhyloTree::optimizeAllBranchesLS(PhyloNode *node, PhyloNode *dad) {
    if (!node) {
        node = getRoot();
    }
    if (dad) {
        double lsBran = optimizeOneBranchLS(node, dad);
        PhyloNeighbor *node_dad_nei = node->findNeighbor(dad);
        PhyloNeighbor *dad_node_nei = dad->findNeighbor(node);
        node_dad_nei->length = lsBran;
        dad_node_nei->length = lsBran;
    }
    FOR_EACH_ADJACENT_PHYLO_NODE(node, dad, it, child) {
        optimizeAllBranchesLS(child, node);
    }
}

void PhyloTree::optimizeAllBranches(PhyloNode *node, PhyloNode *dad, int maxNRStep) {
    if (!node) {
        node = getRoot();
    }
    FOR_EACH_ADJACENT_PHYLO_NODE(node, dad, it, child) {
        optimizeAllBranches(child, node, maxNRStep);
    }
    if (dad) {
        optimizeOneBranch(node, dad, true, maxNRStep); // BQM 2014-02-24: true was missing
    }
}

void PhyloTree::computeBestTraversal(NodeVector &nodes, NodeVector &nodes2) {
    PhyloNode *farleaf = findFarthestLeaf();
    // double call to farthest leaf to find the longest path on the tree
    findFarthestLeaf(farleaf);
    if (verbose_mode >= VerboseMode::VB_MAX) {
        cout << "Tree diameter: " << farleaf->height << endl;
    }
    getPreOrderBranches(nodes, nodes2, farleaf);
}

double PhyloTree::optimizeAllBranches(int my_iterations, double tolerance,
                                      int maxNRStep, bool were_lengths_consistent,
                                      PhyloTree* report_to_tree) {
    if (report_to_tree==nullptr) {
        report_to_tree = this;
    }
    TREE_LOG_LINE(*report_to_tree, VerboseMode::VB_MAX,
                  "Optimizing branch lengths (max " << my_iterations << " loops)...");
    PhyloNodeVector nodes, nodes2;
    computeBestTraversal(nodes, nodes2);
    PhyloNode*     firstNode      = nodes[0];
    PhyloNeighbor* firstNeighbor  = firstNode->findNeighbor(nodes2[0]);
    double         previous_score = computeLikelihoodBranch(firstNeighbor, firstNode,
                                                            tree_buffers);
    TREE_LOG_LINE(*report_to_tree, VerboseMode::VB_MAX,
                  "Initial tree log-likelihood: " << previous_score);
    DoubleVector lenvec;
    double work_estimate = (double)my_iterations * (double)nodes.size();
    report_to_tree->initProgress(work_estimate, "Optimizing branch lengths", "", "", true);
    for (int i = 0; i < my_iterations; i++) {
        LOG_LINE(VerboseMode::VB_MAX, 
                 "Likelihood before iteration " << i + 1 << " : " << previous_score);
        saveBranchLengths(lenvec);
        for (int j = 0; j < nodes.size(); j++) {
            optimizeOneBranch(nodes[j], nodes2[j]);
            //LOG_LINE(VerboseMode::VB_MAX, 
            //         "Branch " << nodes[j]->id << " " << nodes2[j]->id << 
            //         ": " << computeLikelihoodFromBuffer());
            if ( (j % 100) == 99) {
                report_to_tree->trackProgress(100.0);
            }
        }        
        report_to_tree->trackProgress(static_cast<double>(nodes.size() % 100));

        curScore = computeLikelihoodFromBuffer();
        TREE_LOG_LINE(*report_to_tree, VerboseMode::VB_MAX,
                      "Likelihood after iteration " << i + 1 << " : " << curScore);

        if (were_lengths_consistent && curScore < previous_score - tolerance*0.1) {
            // IN RARE CASE: tree log-likelihood decreases, revert the branch length and stop
            TREE_LOG_LINE(*report_to_tree, VerboseMode::VB_MED,
                          "NOTE: Restoring branch lengths"
                          << " as tree log-likelihood decreases"
                          << " after branch length optimization: "
                          << previous_score << " -> " << curScore);

            clearAllPartialLH();
            restoreBranchLengths(lenvec);

            double max_delta_lh = 1.0;
            // Increase max delta with PoMo because log likelihood is very much lower.
            if (aln->seq_type == SeqType::SEQ_POMO) max_delta_lh = 3.0;
            // Different max delta if (aln->seq_type == SeqType::SEQ_CODON) ?!
            curScore = computeLikelihood();
            if (fabs(curScore - previous_score) > max_delta_lh) {
                report_to_tree->hideProgress();
                printTree(cout);
                cout << endl;
                report_to_tree->showProgress();
                TREE_LOG_LINE(*report_to_tree, VerboseMode::VB_QUIET,
                              "new_tree_lh: " << curScore <<
                              " previous_tree_lh: " << previous_score);
                if (!params->ignore_any_errors) {
                    ASSERT(fabs(curScore - previous_score) < max_delta_lh);
                }
            }
            break;
        }
        // only return if the new_tree_lh >= tree_lh!
        // (in rare case that likelihood decreases, continue the loop)
        if (previous_score <= curScore && curScore <= previous_score + tolerance) {
            break;
        }
        LOG_LINE(VerboseMode::VB_DEBUG, 
                 "Prev " << previous_score << ", Curr " << curScore);
        previous_score = curScore;
        were_lengths_consistent = true;
    }
    report_to_tree->doneProgress();
    return curScore;
}

void PhyloTree::moveRoot(Node *node1, Node *node2) {
    // unplug root from tree
    PhyloNode* root_dad  = getRoot()->firstNeighbor()->getNode();
    PhyloNode* root_nei1 = nullptr;
    PhyloNode* root_nei2 = nullptr;
    double len = 0.0;
    FOR_EACH_ADJACENT_PHYLO_NODE(root_dad, root, it, node) {
        if (root_nei1==nullptr) {
            root_nei1 = node;
        }
        else if (root_nei2==nullptr) {
            root_nei2 = node;
        }
        else {
            outError("Cannot move multifurcating root branch");
        }
        len += (*it)->length;
    }
    root_nei1->updateNeighbor(root_dad, root_nei2, len);
    root_nei2->updateNeighbor(root_dad, root_nei1, len);

    // plug root to new branch
    len = node1->findNeighbor(node2)->length / 2.0;
    root_dad->updateNeighbor(root_nei1, node1, len);
    node1->updateNeighbor(node2, root_dad, len);
    root_dad->updateNeighbor(root_nei2, node2, len);
    node2->updateNeighbor(node1, root_dad, len);
    
    if (isSuperTree()) {
        ((PhyloSuperTree*) this)->mapTrees();
    }
    if (Params::getInstance().pll) {
        pllReadNewick(getTreeString());
    }
    resetCurScore();
    if (Params::getInstance().fixStableSplits || 
        Params::getInstance().adaptPertubation) {
        buildNodeSplit();
    }
    current_it = current_it_back = NULL;
    clearBranchDirection();
    computeBranchDirection();
}

double PhyloTree::optimizeRootPosition(int root_dist, bool write_info, double logl_epsilon) {
    if (!rooted) {
        return curScore;
    }
    NodeVector nodes1, nodes2;
    getBranches(root_dist+1, nodes1, nodes2);
    int i;
    Node *root_dad = root->neighbors[0]->node;

    double best_score = curScore;
    string best_tree = getTreeString();

    StrVector trees;

    // ignore branches directly descended from root branch
    for (i = 0; i != nodes1.size(); ) {
        if (nodes1[i] == root_dad || nodes2[i] == root_dad) {
            nodes1[i] = nodes1[nodes1.size()-1];
            nodes2[i] = nodes2[nodes2.size()-1];
            nodes1.pop_back();
            nodes2.pop_back();
        } else {
            ++i;
        }
    }

    // get all trees
    for (i = 0; i != nodes1.size(); i++) {
        moveRoot(nodes1[i], nodes2[i]);
        trees.push_back(getTreeString());
    }

    // optimize branch lengths for all trees
    for (auto t = trees.begin(); t != trees.end(); t++) {
        readTreeString(*t);
        optimizeAllBranches(100, logl_epsilon);
        if (verbose_mode >= VerboseMode::VB_MED) {
            cout << "Root pos " << (t - trees.begin())+1 << ": " << curScore << endl;
            if (verbose_mode >= VerboseMode::VB_DEBUG) {
                drawTree(cout);
            }
        }
        if (curScore > best_score + logl_epsilon) {
            if (verbose_mode >= VerboseMode::VB_MED || write_info) {
                cout << "Better root: " << curScore << endl;
            }
            best_score = curScore;
            best_tree = getTreeString();
        }
    }
    readTreeString(best_tree);
    curScore = computeLikelihood();

    ASSERT(fabs(curScore-best_score) < logl_epsilon);

    return curScore;
}

double PhyloTree::testRootPosition(bool write_info, double logl_epsilon) {
    if (!rooted) {
        return curScore;
    }    
    BranchVector branches;
    getBranches(branches);
    int i;
    Node *root_nei = root->neighbors[0]->node;
    ASSERT(root_nei->degree() == 3);
    Branch root_br = {NULL, NULL};
    FOR_NEIGHBOR_IT(root_nei, root, it) {
        if (root_br.first == NULL)
            root_br.first = (*it)->node;
        else
            root_br.second = (*it)->node;
    }
    
    double best_score = curScore, orig_score = curScore;
    
    multimap<double, string> logl_trees;
    
    // ignore branches directly descended from root branch
    for (i = 0; i != branches.size(); )
        if (branches[i].first == root_nei || branches[i].second == root_nei) {
            branches[i] = branches[branches.size()-1];
            branches.pop_back();
        } else {
            ++i;
        }
    branches.push_back(root_br);
        
    // get all trees
//    for (i = 0; i != nodes1.size(); i++) {
//        moveRoot(nodes1[i], nodes2[i]);
//        trees.push_back(getTreeString());
//    }
//
    // optimize branch lengths for all trees
    for (i = 0; i != branches.size(); i++) {
//    for (auto t = trees.begin()+1; t != trees.end(); t++) {
        moveRoot(branches[i].first, branches[i].second);
//        readTreeString(*t);
        optimizeAllBranches(100, logl_epsilon);
        stringstream ss;
        printTree(ss);
        logl_trees.insert({curScore, ss.str()});
        if (verbose_mode >= VerboseMode::VB_MED) {
            cout << "Root pos " << i+1 << ": " << curScore << endl;
            if (verbose_mode >= VerboseMode::VB_DEBUG)
                drawTree(cout);
        }
        if (curScore > best_score + logl_epsilon) {
            if (verbose_mode >= VerboseMode::VB_MED || write_info)
                cout << "Better root: " << curScore << endl;
            best_score = curScore;
        }
    }
    
//    readTreeString(best_tree);
//    curScore = computeLikelihood();
    
    ASSERT(curScore > orig_score - 0.1);
    if (curScore > orig_score)
        cout << "UPDATE BEST SCORE: " << curScore << endl;
    
    ofstream out;
    string out_file = (string)params->out_prefix + ".rooted_trees";
    out.open(out_file);
    out.precision(10);
    for (auto lt = logl_trees.rbegin(); lt != logl_trees.rend(); lt++) {
        out << "[ lh=" << lt->first << " ]" << lt->second << endl;
    }
    out.close();
    cout << "Rooted trees with log-likelihoods printed to " << out_file << endl;
    if (params->treeset_file.empty())
        params->treeset_file = out_file;

    // convert logL to weight based on the best score
//    ASSERT(logLs.size() == nodes1.size());
//    for (i = 0; i < logLs.size(); i++) {
//        double weight = exp(logLs[i] - best_score);
//        nodes1[i]->name = convertDoubleToString(weight);
//    }
    
    return curScore;
}

void PhyloTree::growTreeML(Alignment *alignment) {
    cout << "Stepwise addition using ML..." << endl;
    aln = alignment;
    size_t size = aln->getNSeq();
    if (size < 3) {
        outError(ERR_FEW_TAXA);
    }
    root = newNode();

    // create initial tree with 3 taxa
    for (leafNum = 0; leafNum < 3; leafNum++) {
        cout << "Add " << aln->getSeqName(leafNum) << " to the tree" << endl;
        Node* new_taxon = newNode(leafNum, aln->getSeqName(leafNum).c_str());
        root->addNeighbor(new_taxon, 1.0);
        new_taxon->addNeighbor(root, 1.0);
    }
    root = findNodeID(0);
    optimizeAllBranches();

    // stepwise adding the next taxon
    for (leafNum = 3; leafNum < size; leafNum++) {
        cout << "Add " << aln->getSeqName(leafNum) << " to the tree" << endl;
        // allocate a new taxon and a new ajedcent internal node
        PhyloNode* new_taxon  = newNode(leafNum, aln->getSeqName(leafNum).c_str());
        PhyloNode* added_node = newNode();
        added_node->addNeighbor(new_taxon, 1.0);
        new_taxon->addNeighbor(added_node, 1.0);

        // preserve two neighbors
        added_node->addNeighbor(DUMMY_NODE_1, 1.0);
        added_node->addNeighbor(DUMMY_NODE_2, 1.0);

        PhyloNode *target_node = nullptr;
        PhyloNode *target_dad  = nullptr;
        double dummy1, dummy2, dummy3;
        addTaxonML(new_taxon, added_node,
                   getRoot()->firstNeighbor()->getNode(), getRoot(),
                   false, target_node, target_dad,dummy1, dummy2, dummy3);
        // now insert the new node in the middle of the branch node-dad
        double len = target_dad->findNeighbor(target_node)->length;
        target_node->updateNeighbor(target_dad, added_node, len / 2.0);
        target_dad->updateNeighbor(target_node, added_node, len / 2.0);
        added_node->updateNeighbor(DUMMY_NODE_1, target_node, len / 2.0);
        added_node->updateNeighbor(DUMMY_NODE_2, target_dad, len / 2.0);
        // compute the likelihood
        added_node->clearReversePartialLh(new_taxon);
        optimizeAllBranches(10, TOL_LIKELIHOOD, 10, false); 
    }

    nodeNum = 2 * leafNum - 2;
}

/****************************************************************************
 Precalculation of "flattened" structure to speed determination of distance functions
 ****************************************************************************/

void PhyloTree::prepareToComputeDistances() {
#ifdef _OPENMP
    int threads = omp_get_max_threads();
#else
    int threads = 1;
#endif
    distanceProcessors.reserve(threads);
    for (threads -= static_cast<int>(distanceProcessors.size()); 0 < threads; --threads) {
        distanceProcessors.push_back(new AlignmentPairwise(this));
    }
    if (summary!=nullptr && !isSummaryBorrowed) {
        delete summary;
        summary = nullptr;
    }
    if (params->use_alignment_summary_for_distance_calculation
        && !isSummaryBorrowed) {
        summary = new AlignmentSummary(aln, true, true);
        summary->constructSequenceMatrix(false);
    }
}

bool PhyloTree::hasMatrixOfConvertedSequences() const {
    return summary!=nullptr && summary->hasSequenceMatrix();
}

size_t PhyloTree::getConvertedSequenceLength() const {
    if (summary==nullptr) {
        return 0;
    }
    return summary->getSequenceLength();
}

const char* PhyloTree::getConvertedSequenceByNumber(int seq1) const {
    if (!hasMatrixOfConvertedSequences()) {
        return nullptr;
    }
    return summary->getSequence(seq1);
}

const int* PhyloTree::getConvertedSequenceFrequencies() const {
    if (summary==nullptr) {
        return nullptr;
    }
    return summary->getSiteFrequencies().data();
}

const int* PhyloTree::getConvertedSequenceNonConstFrequencies() const {
    if (summary==nullptr) {
        return nullptr;
    }
    return summary->getNonConstSiteFrequencies().data();
}

size_t  PhyloTree::getSumOfFrequenciesForSitesWithConstantState(int state) const {
    if (summary==nullptr) {
        return 0;
    }
    return summary->getSumOfConstantSiteFrequenciesForState(state);
}

void PhyloTree::doneComputingDistances() {
    int p=0;
    for ( auto it = distanceProcessors.begin()
         ; it!=distanceProcessors.end(); ++it, ++p) {
        if (verbose_mode >= VerboseMode::VB_MAX) {
            if ( 0 < (*it)->costCalculationCount) {
                double ratio = (double) ((*it)->derivativeCalculationCount) /
                (double) ((*it)->costCalculationCount);
                std::cout << "Processor " << p << " processed "
                    << (*it)->pairCount << " pairs, evaluating cost "
                    << (*it)->costCalculationCount << " times, and finding "
                    << (*it)->derivativeCalculationCount << " derivatives "
                    << "( ratio " << ratio << " )"
                    << endl;
            }
        }
        delete (*it);
    }
    distanceProcessors.clear();
    if (!isSummaryBorrowed) {
        delete summary;
    }
    summary = nullptr;
}

/****************************************************************************
Distance function
****************************************************************************/
double PhyloTree::computeDist(int seq1, int seq2, double initial_dist, double &d2l) {
    // if no model or site rate is specified, return JC distance
    if (initial_dist == 0.0) {
        if (params->compute_obs_dist) {
            initial_dist = aln->computeObsDist(seq1, seq2);
            return initial_dist;
        } else {
            initial_dist = aln->computeDist(seq1, seq2);
        }
    }
    if (!model_factory || !site_rate) {
        return initial_dist; // MANUEL: here no d2l is return
    }
    // now optimize the distance based on the model and site rate
    AlignmentPairwise aln_pair ( this, seq1, seq2 );
    double dist = aln_pair.optimizeDist(initial_dist, d2l);
    return dist;
}

double PhyloTree::computeDist(int seq1, int seq2, double initial_dist) {
    double var;
    return computeDist(seq1, seq2, initial_dist, var);
}

double PhyloTree::correctDist(double *dist_mat) {
    size_t n = aln->getNSeq();
    size_t nsqr = n * n;
    // use Floyd algorithm to find shortest path between all pairs of taxa
    for (size_t k = 0; k < n; ++k) {
        for (size_t i = 0, pos = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j, ++pos) {
                double tmp = dist_mat[i * n + k] + dist_mat[k * n + j];
                if (dist_mat[pos] > tmp) {
                    dist_mat[pos] = tmp;
                }
            }
        }
    }
    double longest_dist = 0.0;
    for (size_t i = 0; i < nsqr; ++i) {
        if (dist_mat[i] > longest_dist) {
            longest_dist = dist_mat[i];
        }
    }
    return longest_dist;
}

template <class L, class F> double computeDistanceMatrix
    ( LEAST_SQUARE_VAR vartype
    , L unknown, const L* sequenceMatrix, intptr_t nseqs, size_t seqLen
    , double denominator, const F* frequencyVector
    , bool uncorrected, double num_states
    , double *dist_mat)
{
    //
    //L is the character type
    //sequenceMatrix is nseqs rows of seqLen characters
    //dist_mat is as in computeDist
    //F is the frequency count type
    //
    
    DoubleVector rowMaxDistance;
    rowMaxDistance.resize(nseqs, 0.0);
    double z = num_states / (num_states - 1.0);
    //Compute the upper-triangle of the distance matrix
    //(and write the row maximum onto the firt cell in the row)
    
    bool count_unknown_as_different = Params::getInstance().count_unknown_as_different;
        //Unknown counts as 50% different.
    
    #if USE_PROGRESS_DISPLAY
    progress_display progress(nseqs*(nseqs-1)/2, "Calculating observed distances");
    #else
    double progress = 0;
    #endif
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (intptr_t seq1 = 0; seq1<nseqs; ++seq1 ) {
        //Scanning is from bottom to top so that if "uneven execution"
        //results in the last few rows being allocated to some worker thread
        //just before the others finish... it won't be running
        //"all by itsef" for as long.
        intptr_t rowOffset     = nseqs * seq1;
        double*  distRow       = dist_mat       + rowOffset;
        const L* thisSequence  = sequenceMatrix + seq1 * seqLen;
        const L* otherSequence = thisSequence   + seqLen;
        double maxDistanceInRow = 0.0;
        for (intptr_t seq2 = seq1 + 1; seq2 < nseqs; ++seq2) {
            double distance = distRow[seq2];
            if ( 0.0 == distance ) {
                double unknownFreq = 0;
                double hamming =
                    hammingDistance ( unknown, thisSequence, otherSequence
                                    , seqLen, frequencyVector, unknownFreq );
                if (count_unknown_as_different) {
                    distance = ( hamming + unknownFreq * 0.75 ) / denominator;
                    if (0<distance && !uncorrected) {
                        double x      = (1.0 - z * distance);
                        distance      = (x<=0) ? MAX_GENETIC_DIST : ( -log(x) / z );
                    }
                } else if (0<hamming && unknownFreq < denominator) {
                    distance = hamming / (denominator - unknownFreq);
                    if (!uncorrected) {
                        double x      = (1.0 - z * distance);
                        distance      = (x<=0) ? MAX_GENETIC_DIST : ( -log(x) / z );
                    }
                }
                distRow[seq2] = distance;
            }
            if ( maxDistanceInRow < distance )
            {
                maxDistanceInRow = distance;
            }
            otherSequence += seqLen;
        }
        rowMaxDistance[seq1] = maxDistanceInRow;
        progress += (nseqs - 1 - seq1);
    }
    #if USE_PROGRESS_DISPLAY
    progress.done();
    #endif
    
    //
    //Determine the longest distance
    //Todo: Decide if it is worth figuring this out by
    //      compare-exchange first half, second half,
    //      (and shrink by half) repeatedly.
    //      Because each compare-exchange cycle can have
    //      its own #pragma omp parallel for.
    //Note: I don't think it is worth it.  Lots of tricky
    //      code for an O(n) step in an O(L.N^2) problem.
    //
    double longest_dist = 0.0;
    for ( intptr_t seq1 = 0; seq1 < nseqs; ++seq1  ) {
        if ( longest_dist < rowMaxDistance[seq1] ) {
            longest_dist = rowMaxDistance[seq1];
        }
    }

    //
    //Copy upper-triangle into lower-triangle and write
    //zeroes to the diagonal.
    //
    
    #if USE_PROGRESS_DISPLAY
    const char* taskReflect = "Determining distance matrix lower triangle from upper";
    if (verbose_mode < VerboseMode::VB_MED) {
        taskReflect = "";
    }
    double work_to_do = (double)nseqs * ((double)nseqs - 1.0) * 0.5;
    progress_display progress2(work_to_do, taskReflect);
    #else
    double progress2 = 0.0;
    #endif
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for ( intptr_t back_seq = 0; back_seq < nseqs; ++back_seq ) {
        intptr_t seq1      = nseqs - 1 - back_seq;
        intptr_t rowOffset = nseqs * seq1;
        double*  distRow   = dist_mat + rowOffset;
        double*  distCol   = dist_mat + seq1; //current entries in the columns
        for ( intptr_t seq2 = 0; seq2 < seq1; ++seq2, distCol+=nseqs ) {
            distRow [ seq2 ] = *distCol;
        }
        distRow [ seq1 ] = 0.0;
        progress2 += seq1;
    }
    #if USE_PROGRESS_DISPLAY
    progress2.done();
    #endif
    return longest_dist;
}

#define EX_START    double baseTime = getRealTime()
#define EX_TRACE(x) if (verbose_mode < VerboseMode::VB_MED) {} \
                    else cout << (getRealTime()-baseTime) << "s " << x << endl

double PhyloTree::computeDistanceMatrix_Experimental() {
    EX_START;
    PhyloTreeThreadingContext context(*this, params->distance_uses_max_threads);
    if (model_factory!=nullptr && site_rate!=nullptr) {
        return computeDistanceMatrix();
    }
    bool uncorrected = params->compute_obs_dist;
        //Use uncorrected (observed) distances
    intptr_t seqCount = aln->getNSeq();
    bool workToDo = false;
    cout.precision(6);
    EX_TRACE("Checking if distances already calculated...");
    //Check if there's any work to do.
    //If there are no zeroes off the diagonal, the distance
    //matrix (ick) has already been initialized and there's no
    //point doing all the following.
    {
        DoubleVector rowMaxDistance;
        rowMaxDistance.resize(seqCount, 0.0);
        #pragma omp parallel for
        for (intptr_t seq1=0; seq1<seqCount; ++seq1) {
            if (!workToDo) {
                const double* distRow   = dist_matrix + seq1 * seqCount;
                double maxDistanceInRow = 0;
                for (intptr_t seq2=0; seq2<seqCount; ++seq2) {
                    double distance = distRow[seq2];
                    if (0.0 == distance && ( seq1 != seq2 ) ) {
                        workToDo = true;
                        break;
                    }
                    if ( maxDistanceInRow < distance )
                    {
                        maxDistanceInRow = distance;
                    }
                }
                rowMaxDistance[seq1] = maxDistanceInRow;
            }
        }
        if (!workToDo) {
            EX_TRACE("No work to do");
            double longest_dist = 0.0;
            for ( intptr_t seq1 = 0; seq1 < seqCount; ++seq1  ) {
                if ( longest_dist < rowMaxDistance[seq1] ) {
                    longest_dist = rowMaxDistance[seq1];
                }
            }
            return longest_dist;
        }
    }
    EX_TRACE("Summarizing...");
    AlignmentSummary s(aln, false, false);
    int maxDistance = 0;
    
    EX_TRACE("Summarizing found " << s.getSequenceLength()
        << " sites with variation (and non-zero frequency),"
        << " and a state range of " << ( s.getStateCount()));
    intptr_t seq_len = s.getSequenceLength();
    for (intptr_t i=0; i<seq_len; ++i) {
        maxDistance += s.getSiteFrequencies()[i];
    }
    //Todo: Shouldn't this be totalFrequency, rather than
    //      totalFrequencyOfNonConst sites?
    double denominator = (double)s.getTotalFrequencyOfNonConstSites()
        + (double)s.getTotalFrequency() - (double)aln->num_variant_sites;
    EX_TRACE("Maximum possible uncorrected length "
        << ((double)maxDistance / (double)s.getTotalFrequencyOfNonConstSites())
        << " Denominator " << denominator);
    if ( 256 < s.getStateCount()  ) {
        EX_TRACE("Falling back to stock distance calculation");
        double longest = computeDistanceMatrix();
        EX_TRACE("Done stock distance calculation");
        return longest; //computeDist(dist_mat, var_mat);
    }
    EX_TRACE("Constructing sequence-major matrix of states"
        << " at " << s.getSequenceLength() << " varying sites"
        << " for " << s.getSequenceCount() << " sequences");
    s.constructSequenceMatrix(true);
    EX_TRACE("Determining distance matrix");
    const int* frequencies = s.getSiteFrequencies().data();
    double longest = ::computeDistanceMatrix
        ( params->ls_var_type, static_cast<char>(aln->STATE_UNKNOWN)
         , s.getSequenceMatrix(), s.getSequenceCount(), s.getSequenceLength()
         , denominator, frequencies, aln->num_states
         , uncorrected, dist_matrix);
    EX_TRACE("Longest distance was " << longest);
    return longest;
}

double PhyloTree::computeDistanceMatrix() {
    prepareToComputeDistances();
    int  nseqs = aln->getNSeq32();
    double longest_dist = 0.0;
    cout.precision(6);
    double baseTime = getRealTime();
    
    #if USE_PROGRESS_DISPLAY
        double work_estimate = static_cast<double>(nseqs) * static_cast<double>(nseqs) * 0.5;
        progress_display progress(work_estimate, "Calculating distance matrix");
    #else
        double progress = 0.0;
    #endif
    
    //compute the upper-triangle of distance matrix
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (int seq1 = 0; seq1 < nseqs; ++seq1) {
        #ifdef _OPENMP
            int threadNum = omp_get_thread_num();
        #else
            int threadNum = 0;
        #endif
        AlignmentPairwise* processor = distanceProcessors[threadNum];
        size_t rowStartPos = seq1 * nseqs;
        for (int seq2=seq1+1; seq2 < nseqs; ++seq2) {
            size_t sym_pos = rowStartPos + seq2;
            double d2l = 0;
            dist_matrix[sym_pos] = processor->recomputeDist(seq1, seq2, dist_matrix[sym_pos], d2l);
        }
        progress += static_cast<double>(nseqs - seq1 - 1);
    }
    #if USE_PROGRESS_DISPLAY
    progress.done();
    const char* taskReflect = "Determining distance matrix lower triangle from upper";
    if (verbose_mode < VerboseMode::VB_MED) {
        taskReflect = "";
    }
    progress_display progress2(work_estimate, taskReflect);
    #else
    double progress2 = 0.0;
    #endif

    //cout << (getRealTime()-baseTime) << "s Copying to lower triangle" << endl;
    //copy upper-triangle into lower-triangle and set diagonal = 0
    for (size_t seq1 = 1; seq1 < nseqs; ++seq1) {
        size_t rowStartPos = seq1 * nseqs;
        size_t rowStopPos  = rowStartPos + seq1;
        size_t colPos = seq1;
        for (size_t rowPos = rowStartPos; rowPos<rowStopPos; ++rowPos, colPos+=nseqs) {
            auto d = dist_matrix[colPos];
            dist_matrix [ rowPos ] = d;
            if (d > longest_dist) {
                longest_dist = d;
            }
        }
        dist_matrix [ rowStopPos] = 0.0;
        progress2 += seq1;
    }
    #if USE_PROGRESS_DISPLAY
    progress2.done();
    #endif
    doneComputingDistances();

    /*
     if (longest_dist > MAX_GENETIC_DIST * 0.99)
     outWarning("Some distances are saturated. Please check your alignment again");*/
    // NOTE: Bionj does handle long distances already (thanks Manuel)
    //return correctDist(dist_mat);
    EX_TRACE("Longest distance was " << longest_dist);
    return longest_dist;
}

void PhyloTree::decideDistanceFilePath(Params& params) {
    dist_file = params.out_prefix;
    if (!model_factory) {
        if (params.compute_obs_dist)
            dist_file += ".obsdist";
        else
            dist_file += ".mldist";
    } else
        dist_file += ".mldist";
}

void PhyloTree::ensureDistanceMatrixAllocated(size_t minimum_rank) {
    if (this->dist_matrix_rank < minimum_rank) {
        delete[] dist_matrix;
        dist_matrix = nullptr;
    }
    intptr_t nSquared = static_cast<size_t>(minimum_rank) 
                      * static_cast<size_t>(minimum_rank);
    if (this->dist_matrix == nullptr) {
        dist_matrix_rank = minimum_rank;
        dist_matrix = new double[nSquared];
        memset(dist_matrix, 0, sizeof(double) * nSquared);
        is_dist_file_read = false;
    }
}

double PhyloTree::computeDistanceMatrix(Params &params, Alignment *alignment) {
    this->params = &params;
    double longest_dist = 0.0;
    aln = alignment;

    size_t n = alignment->getNSeq();
    ensureDistanceMatrixAllocated(n);
    if (params.dist_file && !is_dist_file_read) {
        longest_dist = alignment->readDist(params.dist_file, params.incremental, dist_matrix);
        dist_file = params.dist_file;
        is_dist_file_read = true;
    }
    else {
        double begin_time = getRealTime();
        longest_dist = (params.use_alignment_summary_for_distance_calculation)
            ? computeDistanceMatrix_Experimental()
            : computeDistanceMatrix();
        if (verbose_mode >= VerboseMode::VB_MED) {
            cout << "Distance calculation time: "
            << getRealTime() - begin_time << " seconds" << endl;
        }
    }
    return longest_dist;
}

void PhyloTree::printDistanceFile() const {
    aln->printDist(params->dist_format, params->dist_compression_level,
                   dist_file.c_str(), dist_matrix);
    distanceFileWritten = dist_file.c_str();
}

double PhyloTree::computeObservedDistanceMatrix() {
    int nseqs = aln->getNSeq32();
    ensureDistanceMatrixAllocated(nseqs);
    double longest_dist = 0.0;
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (int seq1 = 0; seq1 < nseqs; ++seq1) {
        size_t pos       = static_cast<size_t>(seq1) 
                         * static_cast<size_t>(nseqs);
        size_t upper_pos = static_cast<size_t>(seq1);
        for (int seq2 = 0; seq2 < nseqs; 
            ++seq2, ++pos, upper_pos+=static_cast<size_t>(nseqs)) {
            if (seq1 == seq2)
                dist_matrix[pos] = 0.0;
            else if (seq1 < seq2) {
                dist_matrix[pos] = aln->computeObsDist(seq1, seq2);
            }
            else {
                //copy from upper triangle into lower triangle
                dist_matrix[pos] = dist_matrix[upper_pos];
            }
            if (dist_matrix[pos] > longest_dist) {
                longest_dist = dist_matrix[pos];
            }
        }
    }
    return longest_dist;
}

double PhyloTree::computeObservedDistanceMatrix(Params &params, Alignment *alignment) {
    double longest_dist = 0.0;
    aln = alignment;
    dist_file = params.out_prefix;
    dist_file += ".obsdist";

    longest_dist = computeObservedDistanceMatrix();
    return longest_dist;
}

/****************************************************************************
 compute BioNJ tree, a more accurate extension of Neighbor-Joining
 ****************************************************************************/

void PhyloTree::computeBioNJ(Params &params) {
    deleteAllPartialLhAndParsimony();
    string bionj_file = params.out_prefix;
    bionj_file += ".bionj";
    this->decideDistanceFilePath(params);
    auto treeBuilder
        = StartTree::Factory::getTreeBuilderByName
            ( params.start_tree_subtype_name);
    bool wasDoneInMemory = false;
    double timeToWriteDistanceFile = 0.0;
    double timeToCalculateMatrix = 0.0;
#ifdef _OPENMP
    #ifdef CLANG_UNDER_VS
        omp_set_max_active_levels(1);
    #else
        omp_set_nested(true); //allow nested OpenMP parallelism
    #endif
    #pragma omp parallel num_threads(2)
    {
        int thread = omp_get_thread_num();
#else
    for (int thread=0; thread<2; ++thread) {
#endif
        if (thread==0) {
            if (!params.dist_file) {
                //This will take longer
                double write_begin_time = getRealTime();
                printDistanceFile();
                if (verbose_mode >= VerboseMode::VB_MED) {
                    //Don't log it yet! It might mess up progress reporting!
                    timeToWriteDistanceFile = getRealTime() - write_begin_time;
                }
            }
        } else if (this->dist_matrix!=nullptr) {
            double start_time = getRealTime();
            wasDoneInMemory = treeBuilder->constructTreeInMemory
            ( this->aln->getSeqNames(), dist_matrix, bionj_file);
            if (wasDoneInMemory) {
                //Don't log it yet!
                timeToCalculateMatrix = getRealTime() - start_time;
            }
        }
    }
    #ifdef _OPENMP
        #pragma omp barrier
        #ifdef CLANG_UNDER_VS
            omp_set_max_active_levels(0);
        #else
            omp_set_nested(false); // don't allow nested OpenMP parallelism
        #endif
    #endif
    if (timeToWriteDistanceFile!=0.0) {
        //This information is logged *afterwards* because, if
        //timeToWriteDistanceFile was logged *during* the
        //calculation of the tree it could mess up the formatting
        //of the progress reporting (if any) done by the
        //distance matrix algorithm.
        LOG_LINE(VerboseMode::VB_MED, "Time taken to write distance file: "
                 << timeToWriteDistanceFile << " seconds " );
    }
    if (wasDoneInMemory) {
        cout << "Computing " << treeBuilder->getName() << " tree"
            << " (from in-memory) distance matrix took "
            << timeToCalculateMatrix << " sec." << endl;
    }
    if (!wasDoneInMemory) {
        double start_time = getRealTime();
        treeBuilder->constructTree(dist_file, bionj_file);
        if (verbose_mode >= VerboseMode::VB_MED) {
            cout << "Constructing " << treeBuilder->getName() << " tree"
                << " (from distance file " << dist_file << ") took "
                << (getRealTime()-start_time) << " sec." << endl;
        }
    }
    double tree_load_start_time = getRealTime();
    readTreeFile(bionj_file.c_str());
    if (verbose_mode >= VerboseMode::VB_MED) {
        cout << "Loading tree (from file " << bionj_file << ") took "
             << (getRealTime()-tree_load_start_time) << " sec." << endl;
    }
}

int PhyloTree::setNegativeBranch(bool force, double newlen, Node *node, Node *dad) {
    if (!node) node = root;
    int fixed = 0;
    FOR_NEIGHBOR_IT(node, dad, it) {
        if ((*it)->length < 0.0 || force) { // negative branch length detected
            (*it)->length = newlen;
            // set the backward branch length
            (*it)->node->findNeighbor(node)->length = (*it)->length;
            ++fixed;
        }
        fixed += setNegativeBranch(force, newlen, (*it)->node, node);
    }
    return fixed;
}

void PhyloTree::fixOneNegativeBranch(double branch_length, Neighbor *dad_branch, Node *dad) {
    dad_branch->length = branch_length;
    // set the backward branch length
    dad_branch->node->findNeighbor(dad)->length = branch_length;
}

int PhyloTree::fixNegativeBranch(bool force, PhyloNode *node, PhyloNode *dad) {

    // 2019-02-05: fix crash when no variant sites found
    if (aln->num_variant_sites == 0) {
        return setNegativeBranch(force, params->min_branch_length, root, NULL);
    }
    if (!node) {
        node = getRoot();
        // 2015-11-30: if not bifurcating, initialize unknown branch lengths with 0.1
        if (!isBifurcating()) {
            return setNegativeBranch(force, 0.1, root, NULL);
        }
    }
    int fixed = 0;
    if (force && !isUsingSankoffParsimony()) {
        return setParsimonyBranchLengths();
    }
    double alpha = (site_rate) ? site_rate->getGammaShape() : 1.0;
    
    FOR_EACH_PHYLO_NEIGHBOR(node, dad, it, nei){
        if (nei->length < 0.0 || force) { // negative branch length detected
            int branch_subst;
            computeParsimonyBranch(nei, node, &branch_subst);
            if (branch_subst < 1) {
                branch_subst = 1;
            }
            double observed_parsimony_distance = (double) branch_subst / getAlnNSite();
            double branch_length = correctBranchLengthF81(observed_parsimony_distance, alpha);
            if (branch_length <= 0) {
                branch_length = params->min_branch_length;
            }
            fixOneNegativeBranch(branch_length, nei, node);
            ++fixed;
        }
        if (nei->length <= 0.0 && (!rooted || node != root)) {
            nei->length = params->min_branch_length;
            nei->getNode()->findNeighbor(node)->length = nei->length;
        }
        fixed += fixNegativeBranch(force, nei->getNode(), node);
    }
    return fixed;
}

//int PhyloTree::assignRandomBranchLengths(bool force, Node *node, Node *dad) {
//
//    if (!node)
//        node = root;
//    int fixed = 0;
//
//    FOR_NEIGHBOR_IT(node, dad, it){
//        if ((*it)->length < 0.0 || force) { // negative branch length detected
//            if (verbose_mode >= VerboseMode::VB_DEBUG)
//            cout << "Negative branch length " << (*it)->length << " was set to ";
//            (*it)->length = random_double() + 0.1;
//            if (verbose_mode >= VerboseMode::VB_DEBUG)
//            cout << (*it)->length << endl;
//            // set the backward branch length
//            (*it)->node->findNeighbor(node)->length = (*it)->length;
//            ++fixed;
//        }
//        if ((*it)->length <= 0.0) {
//            (*it)->length = 1e-6;
//            (*it)->node->findNeighbor(node)->length = (*it)->length;
//        }
//        fixed += assignRandomBranchLengths(force, (*it)->node, node);
//    }
//    return fixed;
//}

/****************************************************************************
 Nearest Neighbor Interchange by maximum likelihood
 ****************************************************************************/

/*
void PhyloTree::doOneRandomNNI(Branch branch) {
    assert(isInnerBranch(branch.first, branch.second));

    if (((PhyloNeighbor*)branch.first->findNeighbor(branch.second))->direction == TOWARD_ROOT) {
        // swap node1 and node2 if the direction is not right, only for nonreversible models
        std::swap(branch.first, branch.second);
    }

    NNIMove nni;
    nni.node1 = (PhyloNode*) branch.first;
    nni.node2 = (PhyloNode*) branch.second;
    FOR_NEIGHBOR_IT(branch.first, branch.second, node1NeiIt)
    if (((PhyloNeighbor*)*node1NeiIt)->direction != TOWARD_ROOT)
    {
        nni.node1Nei_it = node1NeiIt;
        break;
    }
    int randInt = random_int(branch.second->neighbors.size()-1);
    int cnt = 0;
    FOR_NEIGHBOR_IT(branch.second, branch.first, node2NeiIt) {
        // if this loop, is it sure that direction is away from root because node1->node2 is away from root
        if (cnt == randInt) {
            nni.node2Nei_it = node2NeiIt;
            break;
        } else {
            ++cnt;
        }
    }
    assert(*nni.node1Nei_it != NULL && *nni.node2Nei_it != NULL);
    assert(((PhyloNeighbor*)*nni.node1Nei_it)->direction != TOWARD_ROOT && ((PhyloNeighbor*)*nni.node2Nei_it)->direction != TOWARD_ROOT);

    if (constraintTree.isCompatible(nni))
        doNNI(nni, true);
}
*/
    
NNIMove PhyloTree::getRandomNNI(Branch &branch) {
    ASSERT(isInnerBranch(branch.first, branch.second));
    // for rooted tree
    if (((PhyloNeighbor*)branch.first->findNeighbor(branch.second))->direction == TOWARD_ROOT) {
        // swap node1 and node2 if the direction is not right, only for nonreversible models
        swap(branch.first, branch.second);
    }
    NNIMove nni;
    nni.node1 = (PhyloNode*) branch.first;
    nni.node2 = (PhyloNode*) branch.second;

    FOR_NEIGHBOR_IT(branch.first, branch.second, node1NeiIt)
        if (((PhyloNeighbor*)*node1NeiIt)->direction != TOWARD_ROOT) {
            nni.node1Nei_it = node1NeiIt;
            break;
        }
    int randInt = random_int(static_cast<int>(branch.second->neighbors.size())-1);
    int cnt = 0;
    FOR_NEIGHBOR_IT(branch.second, branch.first, node2NeiIt) {
        // if this loop, is it sure that direction is away from root because node1->node2 is away from root
        if (cnt == randInt) {
            nni.node2Nei_it = node2NeiIt;
            break;
        } else {
            ++cnt;
        }
    }
    ASSERT(*nni.node1Nei_it != NULL && *nni.node2Nei_it != NULL);
    ASSERT(((PhyloNeighbor*)*nni.node1Nei_it)->direction != TOWARD_ROOT && ((PhyloNeighbor*)*nni.node2Nei_it)->direction != TOWARD_ROOT);
    nni.newloglh = 0.0;
    return nni;
}

void PhyloTree::doNNI(const NNIMove &move, bool clearLH) {
    PhyloNode*            node1       = move.node1;
    PhyloNode*            node2       = move.node2;
    NeighborVec::iterator node1Nei_it = move.node1Nei_it;
    NeighborVec::iterator node2Nei_it = move.node2Nei_it;
    PhyloNeighbor*        node1Nei    = (PhyloNeighbor*)*(node1Nei_it);
    PhyloNeighbor*        node2Nei    = (PhyloNeighbor*)*(node2Nei_it);

    ASSERT(node1->degree() == 3 && node2->degree() == 3);

    PhyloNeighbor* node12_it = node1->findNeighbor(node2); 
    PhyloNeighbor* node21_it = node2->findNeighbor(node1); 

    // reorient partial_lh before swap
    if (params->compute_likelihood && !isSuperTree()) {
        reorientPartialLh(node12_it, node1);
        reorientPartialLh(node21_it, node2);
    }
    
    // do the NNI swap
    LOG_LINE(VerboseMode::VB_DEBUG, "  Swapping branches " << node2Nei->id << " and " << node1Nei->id);
    node1->updateNeighbor(node1Nei_it, node2Nei);
    node2Nei->node->updateNeighbor(node2, node1);
    PhyloNeighbor* backNei2 = node2Nei->getNode()->findNeighbor(node1);//JB
    backNei2->length        = node2Nei->length;                        //JB

    node2->updateNeighbor(node2Nei_it, node1Nei);
    node1Nei->node->updateNeighbor(node1, node2);
    PhyloNeighbor* backNei1 = node1Nei->getNode()->findNeighbor(node2); //JB
    backNei1->length        = node1Nei->length;                         //JB

    ASSERT(backNei2->id == node2Nei->id);
    ASSERT(backNei1->id == node1Nei->id);
    
    PhyloNeighbor *nei12 = node1->findNeighbor(node2); // return neighbor of node1 which points to node 2
    PhyloNeighbor *nei21 = node2->findNeighbor(node1); // return neighbor of node2 which points to node 1

    if (clearLH) {
        if (params->compute_likelihood) {
            // clear partial likelihood vector
            nei12->clearPartialLh();
            nei21->clearPartialLh();
            node2->clearReversePartialLh(node1);
            node1->clearReversePartialLh(node2);
        }
        nei12->clearPartialParsimony();
        nei21->clearPartialParsimony();
        node2->clearReversePartialParsimony(node1);
        node1->clearReversePartialParsimony(node2);
    }

    if (params->leastSquareNNI) {
        updateSubtreeDists(move);
    }

    // update split store in node
    if (nei12->split != NULL || nei21->split != NULL) {
        delete nei12->split;
        nei12->split = new Split(leafNum);
        delete nei21->split;
        nei21->split = new Split(leafNum);

        FOR_NEIGHBOR_IT(nei12->node, node1, it)
                *(nei12->split) += *((*it)->split);

        FOR_NEIGHBOR_IT(nei21->node, node2, it)
                *(nei21->split) += *((*it)->split);
    }
}

void PhyloTree::changeNNIBrans(const NNIMove &nnimove) {
    PhyloNode*     node1           = nnimove.node1;
    PhyloNode*     node2           = nnimove.node2;
    PhyloNeighbor* node1_node2_nei = node1->findNeighbor(node2);
    PhyloNeighbor* node2_node1_nei = node2->findNeighbor(node1);
    node1_node2_nei->setLength(nnimove.newLen[0]);
    node2_node1_nei->setLength(nnimove.newLen[0]);
    if (params->nni5) {
        int i = 1;
        FOR_EACH_PHYLO_NEIGHBOR(node1, node2, it, nei) {
            PhyloNeighbor* nei_back = (node1)->findNeighbor(nei->node);
            nei->setLength(nnimove.newLen[i]);
            nei_back->setLength(nnimove.newLen[i]);
            ++i;
        }
        FOR_EACH_PHYLO_NEIGHBOR(node2, node1, it, nei) {
            PhyloNeighbor* nei_back = (node2)->findNeighbor(nei->node);
            nei->setLength(nnimove.newLen[i]);
            nei_back->setLength(nnimove.newLen[i]);
            ++i;
        }
    }
}

void PhyloTree::clearInwardViewsFromNeighbors(PhyloNode* node1, PhyloNode* node2) {
    node1->findNeighbor(node2)->clearPartialLh();
    node2->findNeighbor(node1)->clearPartialLh();
    if (!params->nni5) {
        return;
    }
    FOR_EACH_ADJACENT_PHYLO_NODE(node1, nullptr, it, node_X)
    {
        PhyloNeighbor* nei = node_X->findNeighbor(node1);
        nei->clearPartialLh();
    }
    FOR_EACH_ADJACENT_PHYLO_NODE(node2, nullptr, it, node_Y) {
        PhyloNeighbor* nei = node_Y->findNeighbor(node2);
        nei->clearPartialLh();
    }
}

NNIMove PhyloTree::getBestNNIForBran(PhyloNode *node1, PhyloNode *node2, NNIMove* nniMoves) {
    ASSERT(!node1->isLeaf() && !node2->isLeaf());
    ASSERT(node1->degree() == 3 && node2->degree() == 3);
    PhyloNeighbor* nei1 = (node1->findNeighbor(node2));
    if (nei1->direction == TOWARD_ROOT) {
        // swap node1 and node2 if the direction is not right, only for nonreversible models
        std::swap(node1, node2);
    }
    int branch_id = nei1->id;

    NNIContext context(this, node1, node2);
    NNIMove    localNNIMoves[2];
    if (!nniMoves) {
        nniMoves = localNNIMoves;
        //NNIMove constructor now sets node1, node2, and ptnlh to nullptr (James B. 06-Aug-2020)
    }
    
    if (nniMoves[0].node1 != nullptr) {
        // assuming that node1Nei_it and node2Nei_it is defined in nniMoves structure
        for (int cnt = 0; cnt < 2; cnt++) {
            // sanity check
            if (!node1->findNeighbor((*nniMoves[cnt].node1Nei_it)->node)) outError(__func__);
            if (!node2->findNeighbor((*nniMoves[cnt].node2Nei_it)->node)) outError(__func__);
        }
    } else {
        int cnt = 0;
        FOR_EACH_PHYLO_NEIGHBOR(node1, node2, node1_it, nei) {
            if (nei->direction != TOWARD_ROOT)
            {
                cnt = 0;
                FOR_NEIGHBOR_IT(node2, node1, node2_it) {
                    //   Initialize the 2 NNI moves
                    nniMoves[cnt].node1Nei_it = node1_it;
                    nniMoves[cnt].node2Nei_it = node2_it;
                    ++cnt;
                }
                break;
            }
        }
        ASSERT(cnt == 2);
    }

    // Initialize node1 and node2 in nniMoves
    nniMoves[0].central_branch_id = nniMoves[1].central_branch_id = branch_id;
    nniMoves[0].node1    = nniMoves[1].node1    = node1;
    nniMoves[0].node2    = nniMoves[1].node2    = node2;
    nniMoves[0].newloglh = nniMoves[1].newloglh = -DBL_MAX;

    for (int cnt = 0; cnt < 2; cnt++) {
        if (constraintTree.isCompatible(nniMoves[cnt])) {
            NNIMove& move = nniMoves[cnt];
            optimizeOneBranch(move.node1, move.node2, false, NNI_MAX_NR_STEP);
            double old_log_lh = computeLikelihoodFromBuffer();
            move.doSwap(this);
            int nni5_num_eval = max(params->nni5_num_eval, getMixlen());
            double new_log_lh = move.optimizeNNIBranches(this, params->nni5, nni5_num_eval);
            if (old_log_lh<new_log_lh) {
                LOG_LINE(VerboseMode::VB_DEBUG, "cbi=" << move.central_branch_id
                         << ", cnt=" << cnt
                         << ", old_log_lh=" << old_log_lh
                         << ", new_log_lh=" << new_log_lh
                         << ", delta=" << (new_log_lh - old_log_lh));
            }
            move.getLengths(params->nni5);
            if (save_all_trees == 2) {
                saveCurrentTree(move.newloglh); // BQM: for new bootstrap
            }
            move.doSwap(this);
        }
    }
    context.restore();

    LOG_LINE(VerboseMode::VB_MAX, "NNI scores were " << nniMoves[0].newloglh 
             << " and " << nniMoves[1].newloglh);
    
    bool     takeFirst = nniMoves[0].newloglh > nniMoves[1].newloglh;
    NNIMove& res       = nniMoves[ takeFirst ? 0 : 1 ];
    LOG_LINE(VerboseMode::VB_MAX, 
             "  Node1<->Node2 branch has id " << res.node1->findNeighbor(res.node2)->id 
             << ", new len " << res.newLen);
    LOG_LINE(VerboseMode::VB_MAX, 
             "  Swapped branch 1 has id " << (*res.node1Nei_it)->id 
             << " and length " << (*res.node1Nei_it)->length);
    LOG_LINE(VerboseMode::VB_MAX, 
             "  Swapped branch 2 has id " << (*res.node2Nei_it)->id 
             << " and length " << (*res.node2Nei_it)->length);

    #if (0)
        double saved_score = curScore;
        double rescore     = computeLikelihood();
        LOG_LINE(VerboseMode::VB_MAX, "  Rescore after restore was " << rescore
            << " versus previous score of " << saved_score );
    #endif

    return res;
}



void PhyloTree::pruneSubtree(PhyloNode *node, PhyloNode *dad, PruningInfo &info) {

    bool first = true;
    info.node = node;
    info.dad = dad;

    FOR_EACH_PHYLO_NEIGHBOR(dad, node, it, nei){
        if (first) {
            info.dad_it_left  = it;
            info.dad_nei_left = nei;
            info.dad_lh_left  = nei->partial_lh;
            info.left_node    = nei->node;
            info.left_len     = nei->length;
            first = false;
        } else {
            info.dad_it_right  = it;
            info.dad_nei_right = nei;
            info.dad_lh_right  = nei->partial_lh;
            info.right_node    = nei->node;
            info.right_len     = nei->length;
        }
    }
    info.left_it   = info.left_node->findNeighborIt(dad);
    info.right_it  = info.right_node->findNeighborIt(dad);
    info.left_nei  = (*info.left_it);
    info.right_nei = (*info.right_it);

    info.left_node->updateNeighbor(info.left_it, info.dad_nei_right);
    info.right_node->updateNeighbor(info.right_it, info.dad_nei_left);
    ((PhyloNeighbor*) info.dad_nei_right)->partial_lh = newPartialLh();
    ((PhyloNeighbor*) info.dad_nei_left)->partial_lh = newPartialLh();
    //James B. 07-Oct-2020.  Look to see if these are leaked?
}

void PhyloTree::regraftSubtree(PruningInfo &info, PhyloNode *in_node, PhyloNode *in_dad) {
    NeighborVec::iterator in_node_it = in_node->findNeighborIt(in_dad);
    NeighborVec::iterator in_dad_it = in_dad->findNeighborIt(in_node);
    Neighbor *in_dad_nei = (*in_dad_it);
    Neighbor *in_node_nei = (*in_node_it);
    info.dad->updateNeighbor(info.dad_it_right, in_dad_nei);
    info.dad->updateNeighbor(info.dad_it_left, in_node_nei);
}

/****************************************************************************
 Approximate Likelihood Ratio Test with SH-like interpretation
 ****************************************************************************/

/*void PhyloTree::computeNNIPatternLh(double cur_lh, double &lh2, double *pattern_lh2, double &lh3, double *pattern_lh3,
        PhyloNode *node1, PhyloNode *node2) {

    assert(node1->degree() == 3 && node2->degree() == 3);

    // recompute pattern scaling factors if necessary
    PhyloNeighbor *node12_it = node1->findNeighbor(node2);
    PhyloNeighbor *node21_it = node2->findNeighbor(node1);
    NeighborVec::iterator it;
    const int IT_NUM = 6;

    NeighborVec::iterator saved_it[IT_NUM];
    int id = 0;

    FOR_NEIGHBOR(node1, node2, it)
    {
        saved_it[id++] = (*it)->node->findNeighborIt(node1);
    } else {

        saved_it[id++] = it;
    }

    FOR_NEIGHBOR(node2, node1, it)
    {
        saved_it[id++] = (*it)->node->findNeighborIt(node2);
    } else {
        saved_it[id++] = it;
    }
    assert(id == IT_NUM);

    Neighbor * saved_nei[IT_NUM];
    // save Neighbor and allocate new Neighbor pointer
    for (id = 0; id < IT_NUM; id++) {
        saved_nei[id] = (*saved_it[id]);
        // NOTE BUG DOWN HERE!
        *saved_it[id] = new PhyloNeighbor(saved_nei[id]->node, saved_nei[id]->length); // BUG for PhyloSuperTree!
        ((PhyloNeighbor*) (*saved_it[id]))->partial_lh = newPartialLh();
        ((PhyloNeighbor*) (*saved_it[id]))->scale_num = newScaleNum();
    }

    // get the Neighbor again since it is replaced for saving purpose
    node12_it = node1->findNeighbor(node2);
    node21_it = node2->findNeighbor(node1);

    // PhyloNeighbor *node2_lastnei = NULL;

    // save the first found neighbor of node 1 (excluding node2) in node1_it
    FOR_NEIGHBOR_DECLARE(node1, node2, node1_it)

        break;
    Neighbor *node1_nei = *node1_it;

    bool first = true;

    FOR_NEIGHBOR_IT(node2, node1, node2_it) {
        // do the NNI swap
        Neighbor *node2_nei = *node2_it;
        node1->updateNeighbor(node1_it, node2_nei);
        node2_nei->node->updateNeighbor(node2, node1);

        node2->updateNeighbor(node2_it, node1_nei);
        node1_nei->node->updateNeighbor(node1, node2);

        // re-optimize five adjacent branches
        double old_score = -INFINITY, new_score = old_score;

        // clear partial likelihood vector
        node12_it->clearPartialLh();
        node21_it->clearPartialLh();
        int i;
        for (i = 0; i < 2; i++) {

            new_score = optimizeOneBranch(node1, node2, false);

            FOR_EACH_PHYLO_NEIGHBOR(node1, node2, it, nei) {
                //for (id = 0; id < IT_NUM; id++)
                //((PhyloNeighbor*)(*saved_it[id]))->clearPartialLh();
                nei->getNode()->findNeighbor(node1)->clearPartialLh();
                new_score = optimizeOneBranch(node1, nei->getNode(), false);
            }

            node21_it->clearPartialLh();

            FOR_EACH_PHYLO_NEIGHBOR(node2, node1, it, nei) {
                //for (id = 0; id < IT_NUM; id++)
                //((PhyloNeighbor*)(*saved_it[id]))->clearPartialLh();
                nei->getNode()->findNeighbor(node2)->clearPartialLh();
                new_score = optimizeOneBranch(node2, nei->getNode(), false);
                //node2_lastnei = nei;
            }
            node12_it->clearPartialLh();
            if (new_score < old_score + TOL_LIKELIHOOD) break;
            old_score = new_score;
        }
        saveCurrentTree(new_score); // BQM: for new bootstrap

        //new_score = optimizeOneBranch(node1, node2, false);
        if (new_score > cur_lh + TOL_LIKELIHOOD)
        cout << "Alternative NNI shows better likelihood " << new_score << " > " << cur_lh << endl;
        double *result_lh;
        if (first) {
            result_lh = pattern_lh2;
            lh2 = new_score;
        } else {
            result_lh = pattern_lh3;
            lh3 = new_score;
        }
        old_score = new_score;
        computePatternLikelihood(result_lh);
        // swap back and recover the branch lengths
        node1->updateNeighbor(node1_it, node1_nei);
        node1_nei->node->updateNeighbor(node2, node1);
        node2->updateNeighbor(node2_it, node2_nei);
        node2_nei->node->updateNeighbor(node1, node2);
        first = false;
    }

// restore the Neighbor*
    for (id = 0; id < IT_NUM; id++) {

        aligned_free( ((PhyloNeighbor*) *saved_it[id])->scale_num);
        aligned_free( ((PhyloNeighbor*) *saved_it[id])->partial_lh);
        delete (*saved_it[id]);
        (*saved_it[id]) = saved_nei[id];
    }

    // restore the length of 4 branches around node1, node2
    FOR_NEIGHBOR(node1, node2, it)
        (*it)->length = (*it)->node->findNeighbor(node1)->length;
    FOR_NEIGHBOR(node2, node1, it)
        (*it)->length = (*it)->node->findNeighbor(node2)->length;
}*/

void PhyloTree::computeNNIPatternLh(double cur_lh, double &lh2, double *pattern_lh2, double &lh3, double *pattern_lh3,
        PhyloNode *node1, PhyloNode *node2) {
    NNIMove nniMoves[2];
    nniMoves[0].ptnlh = pattern_lh2;
    nniMoves[1].ptnlh = pattern_lh3;
    bool nni5 = params->nni5;
    params->nni5 = true; // always optimize 5 branches for accurate SH-aLRT
    //NNIMove constructor now sets node1 and node2 members to nullptr (James B. 06-Aug-2020)
    getBestNNIForBran(node1, node2, nniMoves);
    params->nni5 = nni5;
    lh2 = nniMoves[0].newloglh;
    lh3 = nniMoves[1].newloglh;
    if (max(lh2,lh3) > cur_lh + TOL_LIKELIHOOD)
        cout << "Alternative NNI shows better log-likelihood " << max(lh2,lh3) << " > " << cur_lh << endl;
}

void PhyloTree::resampleLh(double **pat_lh, double *lh_new, int *rstream) {
    intptr_t nptn = getAlnNPattern();
    memset(lh_new, 0, sizeof(double) * 3);
    int *boot_freq = aligned_alloc<int>(getAlnNPattern());
    aln->createBootstrapAlignment(boot_freq, params->bootstrap_spec, rstream);
    for (intptr_t i = 0; i < nptn; i++) {
        lh_new[0] += boot_freq[i] * pat_lh[0][i];
        lh_new[1] += boot_freq[i] * pat_lh[1][i];
        lh_new[2] += boot_freq[i] * pat_lh[2][i];
    }
    aligned_free(boot_freq);
}

/*********************************************************/
/** THIS FUNCTION IS TAKEN FROM PHYML source code alrt.c
* Convert an aLRT statistic to a none parametric support
* param in: the statistic
*/

double Statistics_To_Probabilities(double in)
{
  double rough_value=0.0;
  double a=0.0;
  double b=0.0;
  double fa=0.0;
  double fb=0.0;

  if(in>=0.000000393 && in<0.00000157)
    {
      a=0.000000393;
      b=0.00000157;
      fa=0.0005;
      fb=0.001;
    }
  else if(in>=0.00000157 && in<0.0000393)
    {
      a=0.00000157;
      b=0.0000393;
      fa=0.001;
      fb=0.005;
    }
  else if(in>=0.0000393 && in<0.000157)
    {
      a=0.0000393;
      b=0.000157;
      fa=0.005;
      fb=0.01;
    }
  else if(in>=0.000157 && in<0.000982)
    {
      a=0.000157;
      b=0.000982;
      fa=0.01;
      fb=0.025;
    }
  else if(in>0.000982 && in<0.00393)
    {
      a=0.000982;
      b=0.00393;
      fa=0.025;
      fb=0.05;
    }
  else if(in>=0.00393 && in<0.0158)
    {
      a=0.00393;
      b=0.0158;
      fa=0.05;
      fb=0.1;
    }
  else if(in>=0.0158 && in<0.0642)
    {
      a=0.0158;
      b=0.0642;
      fa=0.1;
      fb=0.2;
    }
  else if(in>=0.0642 && in<0.148)
    {
      a=0.0642;
      b=0.148;
      fa=0.2;
      fb=0.3;
    }
  else if(in>=0.148 && in<0.275)
    {
      a=0.148;
      b=0.275;
      fa=0.3;
      fb=0.4;
    }
  else if(in>=0.275 && in<0.455)
    {
      a=0.275;
      b=0.455;
      fa=0.4;
      fb=0.5;
    }
  else if(in>=0.455 && in<0.708)
    {
      a=0.455;
      b=0.708;
      fa=0.5;
      fb=0.6;
    }
  else if(in>=0.708 && in<1.074)
    {
      a=0.708;
      b=1.074;
      fa=0.6;
      fb=0.7;
    }
  else if(in>=1.074 && in<1.642)
    {
      a=1.074;
      b=1.642;
      fa=0.7;
      fb=0.8;
    }
  else if(in>=1.642 && in<2.706)
    {
      a=1.642;
      b=2.706;
      fa=0.8;
      fb=0.9;
    }
  else if(in>=2.706 && in<3.841)
    {
      a=2.706;
      b=3.841;
      fa=0.9;
      fb=0.95;
    }
  else if(in>=3.841 && in<5.024)
    {
      a=3.841;
      b=5.024;
      fa=0.95;
      fb=0.975;
    }
  else if(in>=5.024 && in<6.635)
    {
      a=5.024;
      b=6.635;
      fa=0.975;
      fb=0.99;
    }
  else if(in>=6.635 && in<7.879)
    {
      a=6.635;
      b=7.879;
      fa=0.99;
      fb=0.995;
    }
  else if(in>=7.879 && in<10.828)
    {
      a=7.879;
      b=10.828;
      fa=0.995;
      fb=0.999;
    }
  else if(in>=10.828 && in<12.116)
    {
      a=10.828;
      b=12.116;
      fa=0.999;
      fb=0.9995;
    }
  if (in>=12.116)
    {
      rough_value=0.9999;
    }
  else if(in<0.000000393)
    {
      rough_value=0.0001;
    }
  else
    {
      rough_value=(b-in)/(b-a)*fa + (in - a)/(b-a)*fb;
    }
  rough_value=rough_value+(1.0-rough_value)/2.0;
  rough_value=rough_value*rough_value*rough_value;
  return rough_value;
}

// Implementation of testBranch follows Guindon et al. (2010)

double PhyloTree::testOneBranch(double best_score, double *pattern_lh, int reps, int lbp_reps,
        PhyloNode *node1, PhyloNode *node2, double &lbp_support, double &aLRT_support, double &aBayes_support) {
    const int NUM_NNI = 3;
    double lh[NUM_NNI];
    double *pat_lh[NUM_NNI];
    lh[0] = best_score;
    pat_lh[0] = pattern_lh;
    intptr_t nptn = getAlnNPattern();
    pat_lh[1] = new double[nptn];
    pat_lh[2] = new double[nptn];
    int tmp = save_all_trees;
    save_all_trees = 0;
    computeNNIPatternLh(best_score, lh[1], pat_lh[1], lh[2], pat_lh[2], node1, node2);
    save_all_trees = tmp;
    double aLRT;
    if (lh[1] > lh[2])
        aLRT = (lh[0] - lh[1]);
    else
        aLRT = (lh[0] - lh[2]);

    // compute parametric aLRT test support
    double aLRT_stat = 2*aLRT;
    aLRT_support = 0.0;
    if (aLRT_stat >= 0) {
        aLRT_support = Statistics_To_Probabilities(aLRT_stat);
    }

    aBayes_support = 1.0 / (1.0 + exp(lh[1]-lh[0]) + exp(lh[2]-lh[0]));

    int SH_aLRT_support = 0;
    int lbp_support_int = 0;

    int times = max(reps, lbp_reps);

    if (max(lh[1],lh[2]) == -DBL_MAX) {
        SH_aLRT_support = times;
        outWarning("Branch where both NNIs violate constraint tree will show 100% SH-aLRT support");
    } else
#ifdef _OPENMP
#pragma omp parallel
    {
        int *rstream;
        init_random(params->ran_seed + omp_get_thread_num(), false, &rstream);
#pragma omp for reduction(+: lbp_support_int, SH_aLRT_support)
#endif
    for (int i = 0; i < times; i++) {
        double lh_new[NUM_NNI];
        // resampling estimated log-likelihood (RELL)
#ifdef _OPENMP
        resampleLh(pat_lh, lh_new, rstream);
#else
        resampleLh(pat_lh, lh_new, randstream);
#endif
        if (lh_new[0] > lh_new[1] && lh_new[0] > lh_new[2]) {
            ++lbp_support_int;
        }
        double cs[NUM_NNI], cs_best, cs_2nd_best;
        cs[0] = lh_new[0] - lh[0];
        cs[1] = lh_new[1] - lh[1];
        cs[2] = lh_new[2] - lh[2];
        if (cs[0] >= cs[1] && cs[0] >= cs[2]) {
            cs_best = cs[0];
            if (cs[1] > cs[2])
                cs_2nd_best = cs[1];
            else
                cs_2nd_best = cs[2];
        } else if (cs[1] >= cs[2]) {
            cs_best = cs[1];
            if (cs[0] > cs[2])
                cs_2nd_best = cs[0];
            else
                cs_2nd_best = cs[2];
        } else {
            cs_best = cs[2];
            if (cs[0] > cs[1])
                cs_2nd_best = cs[0];
            else
                cs_2nd_best = cs[1];
        }
        if (aLRT > (cs_best - cs_2nd_best) + 0.05) {
            ++SH_aLRT_support;
        }
    }
#ifdef _OPENMP
    finish_random(rstream);
    }
#endif
    delete[] pat_lh[2];
    delete[] pat_lh[1];
    
    lbp_support = 0.0;
    if (times > 0)
        lbp_support = ((double)lbp_support_int) / times;

    if (times > 0)
        return ((double) SH_aLRT_support) / times;
    else
        return 0.0;
}

int PhyloTree::testAllBranches(int threshold, double best_score, double *pattern_lh,
                               int reps, int lbp_reps, bool aLRT_test, bool aBayes_test,
                               PhyloNode *node, PhyloNode *dad) {
    int num_low_support = 0;
    if (!node) {
        node = getRoot();
        node->firstNeighbor()->getNode()->name = "";
        if (isSuperTree()) {
            int tmp = save_all_trees;
            save_all_trees = 2;
            bool nni5 = params->nni5;
            params->nni5 = true; // always optimize 5 branches for accurate SH-aLRT
            initPartitionInfo();
            params->nni5 = nni5;
            save_all_trees = tmp;
        }
    }
    if (dad && !node->isLeaf() && !dad->isLeaf()) {
        double lbp_support, aLRT_support, aBayes_support;
        double SH_aLRT_support = (testOneBranch(best_score, pattern_lh, reps, lbp_reps,
            node, dad, lbp_support, aLRT_support, aBayes_support) * 100);
        ostringstream ss;
        ss.precision(3);
        ss << node->name;
        if (!node->name.empty())
            ss << "/";
        if (reps)
            ss << SH_aLRT_support;
        if (lbp_reps)
            ss << "/" << lbp_support * 100;
        if (aLRT_test)
            ss << "/" << aLRT_support;
        if (aBayes_test)
            ss << "/" << aBayes_support;
        node->name = ss.str();
        if (SH_aLRT_support < threshold)
            num_low_support = 1;
        if (( node->findNeighbor(dad))->partial_pars) {
            ( node->findNeighbor(dad))->partial_pars[0] = static_cast<int>(round(SH_aLRT_support));
            ( dad->findNeighbor(node))->partial_pars[0] = static_cast<int>(round(SH_aLRT_support));
        }
    }
    FOR_EACH_ADJACENT_PHYLO_NODE(node, dad, it, child)
        num_low_support += testAllBranches(threshold, best_score, pattern_lh, reps,
                                           lbp_reps, aLRT_test, aBayes_test, child, node);

    return num_low_support;
}

/****************************************************************************
 Collapse stable (highly supported) clades by one representative
 ****************************************************************************/

void PhyloTree::deleteLeaf(Node *leaf) {

    Node *near_node = leaf->neighbors[0]->node;
    ASSERT(leaf->isLeaf() && near_node->degree() == 3);
    Node *node1 = NULL;
    Node *node2 = NULL;
    double sum_len = 0.0;

    FOR_NEIGHBOR_IT(near_node, leaf, it){
    sum_len += (*it)->length;
    if (!node1)
    node1 = (*it)->node;

    else
    node2 = (*it)->node;
}
// make sure that the returned node1 and node2 are correct
    ASSERT(node1 && node2);
    // update the neighbor
    node1->updateNeighbor(near_node, node2, sum_len);
    node2->updateNeighbor(near_node, node1, sum_len);
}

void PhyloTree::reinsertLeaf(Node *leaf, Node *node, Node *dad) {

    bool first = true;
    Node *adjacent_node = leaf->neighbors[0]->node;
    Neighbor *nei = node->findNeighbor(dad);
    //double len = nei->length;
    double len = max(nei->length, params->min_branch_length * 2);
    // to avoid too small branch length when reinserting leaf

    FOR_NEIGHBOR_IT(adjacent_node, leaf, it){
        if (first) {
            (*it)->node = node;
            (*it)->length = len / 2;
            node->updateNeighbor(dad, adjacent_node, len / 2);
        } else {
            (*it)->node = dad;
            (*it)->length = len / 2;
            dad->updateNeighbor(node, adjacent_node, len / 2);
        }
        first = false;
    }
}

bool PhyloTree::isSupportedNode(PhyloNode* node, int min_support) {
    FOR_EACH_PHYLO_NEIGHBOR(node, nullptr, it, nei) {
        if (!nei->node->isLeaf()) {
            if ( nei->partial_pars[0] < (UINT)min_support) {
                return false;
            }
        }
    }
    return true;
}

int PhyloTree::collapseStableClade(int min_support, NodeVector &pruned_taxa,
                                   StrVector &linked_name, double* &dist_mat) {
    PhyloNodeVector taxa;
    getTaxa(taxa);
    IntVector linked_taxid;
    linked_taxid.resize(leafNum, -1);
    int num_pruned_taxa; // global num of pruned taxa
    int ntaxa = leafNum;
    do {
        num_pruned_taxa = 0;
        for (auto tax_it = taxa.begin(); tax_it != taxa.end(); tax_it++)
            if (linked_taxid[(*tax_it)->id] < 0) {
                PhyloNode *taxon = (*tax_it);
                PhyloNode *near_node = taxon->firstNeighbor()->getNode();
                Node *adj_taxon = NULL;
                FOR_NEIGHBOR_DECLARE(near_node, taxon, it)
                    if ((*it)->node->isLeaf()) {
                        adj_taxon = (*it)->node;
                        break;
                    }
                // if it is not a cherry
                if (!adj_taxon) {
                    continue;
                }
                ASSERT(linked_taxid[adj_taxon->id] < 0);
                PhyloNeighbor * near_nei = NULL;
                FOR_EACH_PHYLO_NEIGHBOR(near_node, taxon, it2, nei)
                    if (nei->node != adj_taxon) {
                        near_nei = nei;
                        break;
                    }
                ASSERT(near_nei);
                // continue if the cherry is not stable, or distance between two taxa is near ZERO
                if (!isSupportedNode(near_nei->getNode(), min_support)
                        && dist_mat[taxon->id * ntaxa + adj_taxon->id] > 2e-6)
                    continue;
                // now do the taxon pruning
                Node * pruned_taxon = taxon, *stayed_taxon = adj_taxon;
                // prune the taxon that is far away
                if (adj_taxon->neighbors[0]->length > taxon->neighbors[0]->length) {
                    pruned_taxon = adj_taxon;
                    stayed_taxon = taxon;
                }
                deleteLeaf(pruned_taxon);
                linked_taxid[pruned_taxon->id] = stayed_taxon->id;
                pruned_taxa.push_back(pruned_taxon);
                linked_name.push_back(stayed_taxon->name);
                ++num_pruned_taxa;
                // do not prune more than n-4 taxa
                if (pruned_taxa.size() >= ntaxa - 4) {
                    break;
                }
            }
    } while (num_pruned_taxa && pruned_taxa.size() < ntaxa - 4);

    if (pruned_taxa.empty())
        return 0;

    if (verbose_mode >= VerboseMode::VB_MED) {
        auto tax_it    = pruned_taxa.begin();
        auto linked_it = linked_name.begin();
        for (; tax_it != pruned_taxa.end(); tax_it++, linked_it++) {
            cout << "Delete " << (*tax_it)->name << " from " << (*linked_it) << endl;
        }
    }

    // set root to the first taxon which was not deleted
    for (auto tax_it = taxa.begin(); tax_it != taxa.end(); tax_it++)
        if (linked_taxid[(*tax_it)->id] < 0) {
            root = (*tax_it);
            break;
        }
    // extract the sub alignment
    IntVector stayed_id;
    int i, j;
    for (i = 0; i < taxa.size(); i++) {
        if (linked_taxid[i] < 0) {
            stayed_id.push_back(i);
        }
    }
    ASSERT(stayed_id.size() + pruned_taxa.size() == leafNum);
    Alignment * pruned_aln = new Alignment();
    pruned_aln->extractSubAlignment(aln, stayed_id, 2); // at least 2 informative characters
    nodeNum = leafNum = static_cast<int>(stayed_id.size());
    initializeTree();
    setAlignment(pruned_aln);

    double *pruned_dist = new double[leafNum * leafNum];
    for (i = 0; i < static_cast<int>(leafNum); i++) {
        for (j = 0; j < static_cast<int>(leafNum); j++) {
            pruned_dist[i * leafNum + j] = dist_mat[stayed_id[i] * ntaxa + stayed_id[j]];
        }
    }
    dist_mat = pruned_dist;
    return static_cast<int>(pruned_taxa.size());
}

int PhyloTree::restoreStableClade(Alignment *original_aln, NodeVector &pruned_taxa, StrVector &linked_name) {
    //int num_inserted_taxa;
    NodeVector::reverse_iterator tax_it;
    StrVector::reverse_iterator linked_it;
    tax_it = pruned_taxa.rbegin();
    linked_it = linked_name.rbegin();
    for (; tax_it != pruned_taxa.rend(); tax_it++, linked_it++) {
        //cout << "Reinsert " << (*tax_it)->name << " to " << (*linked_it) << endl;
        Node *linked_taxon = findNodeName((*linked_it));
        ASSERT(linked_taxon);
        ASSERT(linked_taxon->isLeaf());
        ++leafNum;
        reinsertLeaf((*tax_it), linked_taxon, linked_taxon->neighbors[0]->node);
    }
    ASSERT(leafNum == original_aln->getNSeq());
    nodeNum = leafNum;
    initializeTree();
    setAlignment(original_aln);
    setRootNode(params->root);
    //if (verbose_mode >= VerboseMode::VB_MED) drawTree(cout);

    return 0;
}

bool PhyloTree::checkEqualScalingFactor(double &sum_scaling, PhyloNode *node, PhyloNode *dad) {
    if (!node) {
        node = getRoot();
    }
    if (dad) {
        double scaling = node->findNeighbor(dad)->lh_scale_factor
                + dad->findNeighbor(node)->lh_scale_factor;
        if (sum_scaling > 0)
            sum_scaling = scaling;
        if (fabs(sum_scaling - scaling) > 1e-6) {
            cout << sum_scaling << " " << scaling << endl;
            return false;
        }
    }
    FOR_EACH_ADJACENT_PHYLO_NODE(node, dad, it, child) {
        if (!checkEqualScalingFactor(sum_scaling, child, node)) {
            return false;
        }
    }

    return true;
}

void PhyloTree::randomizeNeighbors(Node *node, Node *dad) {

    if (!node)
        node = root;
    FOR_NEIGHBOR_IT(node, dad, it)randomizeNeighbors((*it)->node, node);

    my_random_shuffle(node->neighbors.begin(), node->neighbors.end());
}

void PhyloTree::printTransMatrices(Node *node, Node *dad) {
    if (!node)
        node = root;
    int nstates = aln->num_states;

    if (dad) {
        double *trans_cat = new double[nstates * nstates];
        model_factory->computeTransMatrix(dad->findNeighbor(node)->length * site_rate->getRate(0), trans_cat);
        cout << "Transition matrix " << dad->name << " to " << node->name << endl;
        for (int i = 0; i < nstates; i++) {
            for (int j = 0; j < nstates; j++) {
                cout << "\t" << trans_cat[i * nstates + j];
            }
            cout << endl;
        }
        delete[] trans_cat;
    }
    FOR_NEIGHBOR_IT(node, dad, it)printTransMatrices((*it)->node, node);
}

void PhyloTree::removeIdenticalSeqs(Params &params) {
    // commented out because it also works for SuperAlignment now!
    Alignment *new_aln;
    // 2017-03-31: always keep two identical sequences no matter if -bb or not, to avoid conflict between 2 subsequent runs
    if (params.root) {
        new_aln = aln->removeIdenticalSeq((string)params.root, true, removed_seqs, twin_seqs);
    } else {
        new_aln = aln->removeIdenticalSeq("", true, removed_seqs, twin_seqs);
    }
    if (removed_seqs.size() > 0) {
        cout << "NOTE: " << removed_seqs.size() << " identical sequences (see below) will be ignored for subsequent analysis" << endl;
        for (int i = 0; i < removed_seqs.size(); i++) {
            if (!params.suppress_duplicate_sequence_warnings) {
                cout << "NOTE: " << removed_seqs[i] << " (identical to " << twin_seqs[i] << ") is ignored but added at the end" << endl;
            }
        }
        delete aln;
        aln = new_aln;
    }
    if (!constraintTree.empty()) {
        int count = constraintTree.removeTaxa(removed_seqs, true, "");
        if (count) {
            cout << count << " taxa removed from constraint tree" << endl;
        }
    }
}

void PhyloTree::reinsertIdenticalSeqs(Alignment *orig_aln) {
    if (removed_seqs.empty()) return;

    //delete all partial likelihood and partial parsimony 
    //vectors, which will be automatically recreated later
    //James B. 03-Dec-2020.  Better to do this *before* 
    //inserting taxa, rather than *after*.
    deleteAllPartialLhAndParsimony();
    
    insertTaxa(removed_seqs, twin_seqs);
    setAlignment(orig_aln);
}

void PhyloTree::computeSeqIdentityAlongTree(Split &sp, Node *node, Node *dad) {
    ASSERT(node && !node->isLeaf());
    // recursive
    FOR_NEIGHBOR_IT(node, dad, it) {
        if ((*it)->node->isLeaf()) {
            sp.addTaxon((*it)->node->id);
        } else {
            Split newsp(leafNum);
            computeSeqIdentityAlongTree(newsp, (*it)->node, node);
            sp += newsp;
        }
    }
    if (!dad) return;
    // now going along alignment to compute seq identity
    int ident = 0;
    int nseqs = aln->getNSeq32();
    for (Alignment::iterator it = aln->begin(); it != aln->end(); it++) {
        char state = aln->STATE_UNKNOWN;
        bool is_const = true;
        for (int i = 0; i < nseqs; ++i) {
            if (static_cast<int>((*it)[i]) < aln->num_states && sp.containTaxon(i)) {
                if (state < aln->num_states && state != (*it)[i]) {
                    is_const = false;
                    break;
                }
                state = (*it)[i];
            }
        }
        if (is_const) {
            ident += it->frequency;
        }
    }
    ident = static_cast<int>((ident*100)/aln->getNSite());
    if (node->name == "")
        node->name = convertIntToString(ident);
    else
        node->name += "/" + convertIntToString(ident);
}

void PhyloTree::computeSeqIdentityAlongTree() {
    Split sp(leafNum);
    if (root->isLeaf())
        computeSeqIdentityAlongTree(sp, root->neighbors[0]->node);
    else
        computeSeqIdentityAlongTree(sp, root);
}

void PhyloTree::generateRandomTree(TreeGenType tree_type) {
    if (!constraintTree.empty() && tree_type != YULE_HARDING)
        outError("Only Yule-Harding ramdom tree supported with constraint tree");
    ASSERT(aln);
    int orig_size = params->sub_size;
    params->sub_size = aln->getNSeq32();
    MExtTree ext_tree;
    if (constraintTree.empty()) {
        switch (tree_type) {
        case YULE_HARDING: 
            ext_tree.generateYuleHarding(*params);
            break;
        case UNIFORM:
            ext_tree.generateUniform(params->sub_size);
            break;
        case CATERPILLAR:
            ext_tree.generateCaterpillar(params->sub_size);
            break;
        case BALANCED:
            ext_tree.generateBalanced(params->sub_size);
            break;
        case STAR_TREE:
            ext_tree.generateStarTree(*params);
            break;
        default:
            break;
        }
        NodeVector taxa;
        ext_tree.getTaxa(taxa);
        ASSERT(taxa.size() == aln->getNSeq());
        for (NodeVector::iterator it = taxa.begin(); it != taxa.end(); it++)
            (*it)->name = aln->getSeqName((*it)->id);
    } else {
        ext_tree.generateConstrainedYuleHarding(*params, &constraintTree, aln->getSeqNames());
    }
    params->sub_size = orig_size;
    stringstream str;
    ext_tree.printTree(str);
    PhyloTree::readTreeStringSeqName(str.str());
}

void PhyloTree::computeBranchDirection(PhyloNode *node, PhyloNode *dad) {
    if (!node) {
        node = getRoot();
    }
    if (dad) {
        node->findNeighbor(dad)->direction = TOWARD_ROOT;
    }
    FOR_EACH_PHYLO_NEIGHBOR(node, dad, it, nei) {
        // do not update if direction was already computed
        ASSERT(nei->direction != TOWARD_ROOT);
        if (nei->direction != UNDEFINED_DIRECTION) {
            continue;
        }
        // otherwise undefined.
        nei->direction = AWAYFROM_ROOT;
        computeBranchDirection(nei->getNode(), node);
    }
}

void PhyloTree::clearBranchDirection(PhyloNode *node, PhyloNode *dad) {
    if (!node) {
        node = getRoot();
    }
    if (dad) {
        node->findNeighbor(dad)->direction = UNDEFINED_DIRECTION;
    }
    FOR_EACH_PHYLO_NEIGHBOR(node, dad, it, nei) {
        nei->direction = UNDEFINED_DIRECTION;
        clearBranchDirection(nei->getNode(), node);
    }

}

/*
void PhyloTree::sortNeighborBySubtreeSize(PhyloNode *node, PhyloNode *dad) {

    // already sorted, return
    PhyloNeighbor *nei = dad->findNeighbor(node);
    if (nei->size >= 1)
        return;

    if (dad && node->isLeaf()) {
        nei->size = 1;
        return;
    }

    nei->size = 0;
    FOR_NEIGHBOR_DECLARE(node, dad, it) {
        sortNeighborBySubtreeSize((PhyloNode*)(*it)->node, node);
        nei->size += ((PhyloNeighbor*)*it)->size;
    }
    
    // sort neighbors in descending order of sub-tree size
    FOR_NEIGHBOR(node, dad, it)
        for (NeighborVec::iterator it2 = it+1; it2 != node->neighbors.end(); it2++)
            if ((*it2)->node != dad && ((PhyloNeighbor*)*it)->size < ((PhyloNeighbor*)*it2)->size) {
                Neighbor *nei;
                nei = *it;
                *it = *it2;
                *it2 = nei;
            }
}
*/

void PhyloTree::convertToRooted() {
    ASSERT(leafNum == aln->getNSeq());
    Node *node, *dad;
    double node_len, dad_len;
    if (params->root) {
        string name = params->root;
        node = findNodeName(name);
        if (!node)
            outError("Cannot find leaf with name " + name);
        ASSERT(node->isLeaf());
        dad = node->neighbors[0]->node;
        node_len = dad_len = node->neighbors[0]->length*0.5;
    } else {
        //midpoint rooting
        Node *node1, *node2;
        double longest = root->longestPath2(node1, node2);
        longest *= 0.5;
        double curlen = 0.0;
        for (; node1 != node2 && curlen + node1->highestNei->length < longest; node1 = node1->highestNei->node) {
            curlen += node1->highestNei->length;
        }
        // place root on node1->heigestNei
        node = node1;
        dad = node1->highestNei->node;
        node_len = longest - curlen;
        dad_len = node1->highestNei->length - node_len;
        ASSERT(dad_len >= 0.0);
        /*
        // place root to the longest branch
        NodeVector nodes1, nodes2;
        getBranches(nodes1, nodes2);
        double max_length = -1.0;
        for (int i = 0; i < nodes1.size(); i++)
            if (max_length < nodes1[i]->findNeighbor(nodes2[i])->length) {
                max_length = nodes1[i]->findNeighbor(nodes2[i])->length;
                node = nodes1[i];
                dad = nodes2[i];
            }
         */
    }
    rooted = true;
    root = newNode(leafNum, ROOT_NAME);
    Node *root_int = newNode();
    root->addNeighbor(root_int, 0.0);
    root_int->addNeighbor(root, 0.0);
    ++leafNum;
    //double newlen = node->findNeighbor(dad)->length/2.0;
    node->updateNeighbor(dad, root_int, node_len);
    root_int->addNeighbor(node, node_len);
    dad->updateNeighbor(node, root_int, dad_len);
    root_int->addNeighbor(dad, dad_len);
    initializeTree();
    computeBranchDirection();
    current_it = current_it_back = NULL;
}

void PhyloTree::convertToUnrooted() {
    ASSERT(rooted);
    if (aln)
        ASSERT(leafNum == aln->getNSeq()+1);
    ASSERT(root);
    ASSERT(root->isLeaf() && root->id == leafNum-1);
    Node *node = root->neighbors[0]->node;
    Node *taxon = findFirstTaxon();

    rooted = false;
    leafNum--;

    // delete root node
    if (node->degree() == 3) {
        // delete and join adjacent branches
        Node *node1 = NULL, *node2 = NULL;
        double len = 0.0;
        FOR_NEIGHBOR_IT(node, root, it) {
            if (!node1) node1 = (*it)->node; else node2 = (*it)->node;
            len += (*it)->length;
        }
        node1->updateNeighbor(node, node2, len);
        node2->updateNeighbor(node, node1, len);
        delete node;
    } else {
        // only delete root node
        auto it = node->findNeighborIt(root);
        delete *it;
        node->neighbors.erase(it);

    }

    delete root;
    // set a temporary taxon so that tree traversal works
    root = taxon;
    
    if (params) {
        setRootNode(params->root);
    }

    initializeTree();
    clearBranchDirection();
//    computeBranchDirection();
}

void PhyloTree::reorientPartialLh(PhyloNeighbor* dad_branch,
                                  PhyloNode *dad) {
    ASSERT(!isSuperTree());
    if ( dad_branch->partial_lh != nullptr ) {
        return;
    }
    PhyloNode* node = dad_branch->getNode();
    FOR_EACH_PHYLO_NEIGHBOR(node, dad, it, nei) {
        PhyloNeighbor *backnei = nei->getNode()->findNeighbor(node);
        if (backnei->partial_lh) {
            mem_slots.takeover(dad_branch, backnei);
            break;
        }
    }
    if ( dad_branch->partial_lh != nullptr ) {
        return;
    }
    if (params->lh_mem_save == LM_PER_NODE) {
        hideProgress();
        ASSERT(dad_branch->partial_lh && "partial_lh is not re-oriented");
        showProgress();
    }
}

/****************************************************************************
        helper functions for computing tree traversal
 ****************************************************************************/

bool PhyloTree::computeTraversalInfo(PhyloNeighbor* dad_branch,
                                     PhyloNode* dad,
                                     double* &buffer) {
    size_t     nstates = aln->num_states;
    PhyloNode* node    = dad_branch->getNode();

    if (dad_branch->isLikelihoodComputed() || node->isLeaf()) {
        return mem_slots.lock(dad_branch);
    }
#ifndef _MSC_VER
    bool locked[node->degree()];
#else
    boost::scoped_array<bool> locked(new bool[node->degree()]);
#endif
    size_t num_leaves = 0;
    for (int i=node->degree()-1; 0 <= i; --i) {
        locked[i] = false;
    }

    // sort neighbor in descending size order (with a selection sort!)
    PhyloNeighborVec neivec(node->neighbors);
    for (size_t i = 0; i < neivec.size(); ++i) {
        for (size_t j = i+1; j < neivec.size(); ++j) {
            if (neivec[i]->size < neivec[j]->size) {
                PhyloNeighbor* t = neivec[i];
                neivec[i] = neivec[j];
                neivec[j] = t;
            }
        }
    }

    // recursive
    for (int i = 0; i < static_cast<int>(neivec.size()); ++i) {
        PhyloNode* child = neivec[i]->getNode();
        if (child != dad) {
            locked[i]   = computeTraversalInfo(neivec[i], node, buffer);
            num_leaves += child->isLeaf() ? 1 : 0;
        }
    }
    dad_branch->setLikelihoodComputed(true);

    // prepare information for this branch
    TraversalInfo info(dad_branch, dad);

    // re-orient partial_lh
    reorientPartialLh(dad_branch, dad);

    if (!dad_branch->partial_lh || mem_slots.locked(dad_branch)) {
        // still no free entry found, memory saving technique
        int slot_id = mem_slots.allocate(dad_branch);
        if (slot_id < 0) {
            cout << "traversal order:";
            for (auto it = traversal_info.begin(); it != traversal_info.end(); it++) {
                it->dad_branch->node->name = convertIntToString(it->dad_branch->size);
                cout << "  ";
                if (it->dad->isLeaf())
                    cout << it->dad->name;
                else
                    cout << it->dad->id;
                cout << "->";
                if (it->dad_branch->node->isLeaf())
                    cout << it->dad_branch->node->name;
                else
                    cout << it->dad_branch->node->id;
                if (params->lh_mem_save == LM_MEM_SAVE) {
                    if (it->dad_branch->isLikelihoodComputed())
                        cout << " [";
                    else
                        cout << " (";
                    cout << mem_slots.findNei(it->dad_branch) - mem_slots.begin();
                    if (it->dad_branch->isLikelihoodComputed())
                        cout << "]";
                    else
                        cout << ")";
                }
            }
            cout << endl;
            drawTree(cout);
            ASSERT(0 && "No free/unlocked mem slot found!");
        }
    } else {
        mem_slots.update(dad_branch);
    }

    if (params->lh_mem_save == LM_MEM_SAVE) {
        if (verbose_mode >= VerboseMode::VB_MED) {
            int slot_id = static_cast<int>(mem_slots.findNei(dad_branch) - mem_slots.begin());
            node->name = convertIntToString(slot_id);
            //cout << "Branch " << dad->id << "-" << node->id << " assigned slot " << slot_id << endl;
        }
        for (int i = 0; i < static_cast<int>(neivec.size()); ++i) {
            PhyloNode* child = neivec[i]->getNode();
            if (child != dad) {
                if (child->isInterior() && locked[i]) {
                    mem_slots.unlock(neivec[i]);
                }
            }
        }
    }

    if (!model->isSiteSpecificModel() && !Params::getInstance().buffer_mem_save) {
        //------- normal model -----
        info.echildren = buffer;
        size_t block   = nstates * ((model_factory->fused_mix_rate)
                                    ? site_rate->getNRate()
                                    : site_rate->getNRate()*model->getNMixtures());
        buffer += get_safe_upper_limit(block*nstates*(node->degree()-1));
        if (0 < num_leaves) {
            info.partial_lh_leaves = buffer;
            buffer += get_safe_upper_limit((aln->STATE_UNKNOWN+1)*block*num_leaves);
        }
    }
    traversal_info.push_back(info);
    return mem_slots.lock(dad_branch);
}

void PhyloTree::writeSiteLh(ostream &out, SiteLoglType wsl, int partid) {
    // error checking
    if (!getModel()->isMixture()) {
        if (wsl != WSL_RATECAT) {
            outWarning("Switch now to '-wslr' as it is the only option for non-mixture model");
            wsl = WSL_RATECAT;
        }
    } else {
        // mixture model
        if (wsl == WSL_MIXTURE_RATECAT && getModelFactory()->fused_mix_rate) {
            outWarning("-wslmr is not suitable for fused mixture model, switch now to -wslm");
            wsl = WSL_MIXTURE;
        }
    }
    size_t nsites = getAlnNSite();
    size_t ncat = getNumLhCat(wsl);
    double *pattern_lh, *pattern_lh_cat;
    pattern_lh = aligned_alloc<double>(getAlnNPattern());
    pattern_lh_cat = aligned_alloc<double>(getAlnNPattern()*ncat);
    computePatternLikelihood(pattern_lh, NULL, pattern_lh_cat, wsl);
    for (size_t i = 0; i < nsites; ++i) {
        if (partid >= 0) {
            out << partid << "\t";
        }
        intptr_t ptn = aln->getPatternID(i);
        out << i+1 << "\t" << pattern_lh[ptn];
        for (int j = 0; j < ncat; j++) {
            out << "\t" << pattern_lh_cat[ptn*ncat+j];
        }
        out << endl;
    }
    aligned_free(pattern_lh_cat);
    aligned_free(pattern_lh);
}

void PhyloTree::writeSiteRates(ostream &out, bool bayes, int partid) {
    DoubleVector pattern_rates;
    IntVector pattern_cat;
    int ncategory = 1;
    
    if (bayes)
        ncategory = site_rate->computePatternRates(pattern_rates, pattern_cat);
    else
        optimizePatternRates(pattern_rates);
    
    if (pattern_rates.empty()) return;
    size_t nsite = aln->getNSite();
    
    out.setf(ios::fixed,ios::floatfield);
    out.precision(5);
    //cout << __func__ << endl;
    IntVector count;
    count.resize(ncategory, 0);
    for (size_t i = 0; i < nsite; ++i) {
        int ptn = aln->getPatternID(i);
        if (partid >= 0)
            out << partid << "\t";
        out << i+1 << "\t";
        if (pattern_rates[ptn] >= MAX_SITE_RATE) out << "100.0"; else out << pattern_rates[ptn];
        //cout << i << " "<< ptn << " " << pattern_cat[ptn] << endl;
        if (!pattern_cat.empty()) {
            int site_cat;
            double cat_rate;
            if (site_rate->getPInvar() == 0.0) {
                site_cat = pattern_cat[ptn]+1;
                cat_rate = site_rate->getRate(pattern_cat[ptn]);
            } else {
                site_cat = pattern_cat[ptn];
                if (site_cat == 0)
                    cat_rate = 0.0;
                else
                    cat_rate = site_rate->getRate(pattern_cat[ptn]-1);
            }
            out << "\t" << site_cat << "\t" << cat_rate;
            ++(count[pattern_cat[ptn]]);
        }
        out << endl;
    }
    if (bayes) {
        cout << "Empirical proportions for each category:";
        for (size_t i = 0; i < count.size(); ++i) {
            cout << " " << ((double)count[i])/nsite;
        }
        cout << endl;
    }
}

void PhyloTree::writeBranches(ostream &out) {
    outError("Please only use this feature with partition model");
}

const string& PhyloTree::getDistanceFileWritten() const {
    return distanceFileWritten;
}
    
void PhyloTree::initProgress(double size, std::string name, 
                             const char* verb, const char* noun, 
                             bool isUpperBound) {
    ++progressStackDepth;
    if (progressStackDepth==1 && !isShowingProgressDisabled) {
        #if USE_PROGRESS_DISPLAY
        progress = new progress_display(size, name.c_str(), verb, noun);
        if (isUpperBound) {
            progress->setIsEstimateABound(true);
        }
        #endif
    }
}
    
void PhyloTree::trackProgress(double amount) {
    #if USE_PROGRESS_DISPLAY
    if (!isShowingProgressDisabled && progressStackDepth==1 && progress!=nullptr) {
        (*progress) += amount;
    }
    #endif
}

void PhyloTree::hideProgress() const {
    #if USE_PROGRESS_DISPLAY
    if (!isShowingProgressDisabled && progressStackDepth>0 && progress!=nullptr) {
        progress->hide();
    }
    #endif
}

void PhyloTree::showProgress() const {
    #if USE_PROGRESS_DISPLAY
    if (!isShowingProgressDisabled && progressStackDepth>0 && progress!=nullptr) {
        progress->show();
    }
    #endif
}

void PhyloTree::doneProgress() {
    --progressStackDepth;
    if (progressStackDepth==0) {
        #if USE_PROGRESS_DISPLAY
        if (progress!=nullptr) {
            progress->done();
        }
        delete progress;
        #endif
        progress = nullptr;
    }
}
    
void PhyloTree::showNoProgress() {
    #if USE_PROGRESS_DISPLAY
    if (progress!=nullptr) {
        delete progress;
        progress = nullptr;
    }
    #endif
    isShowingProgressDisabled = true;
}

    
PhyloNode* PhyloTree::findFarthestLeaf(PhyloNode *node,
                                       PhyloNode *dad) {
    return (PhyloNode*) super::findFarthestLeaf(node, dad);
}

void PhyloTree::computePatternPacketBounds(int v_size, int threads, int packets,
                                           intptr_t elements, vector<intptr_t> &limits) {
    //It is assumed that threads divides packets evenly
    limits.reserve(packets+1);
    elements = roundUpToMultiple(elements, v_size);
    size_t block_start = 0;
    
    for (int wave = packets/threads; wave>=1; --wave) {
        intptr_t elementsThisWave = (elements-block_start);
        if (1<wave) {
            elementsThisWave = (elementsThisWave * 3) / 4;
        }
        elementsThisWave = roundUpToMultiple(elementsThisWave, v_size);
        size_t stopElementThisWave = block_start + elementsThisWave;
        for (int threads_to_go=threads; 1<=threads_to_go; --threads_to_go) {
            limits.push_back(block_start);
            size_t block_size = (stopElementThisWave - block_start)/threads_to_go;
            block_size = roundUpToMultiple(block_size, v_size);
            block_start += block_size;
        }
    }
    limits.push_back(elements);
    
    if (limits.size() != packets+1) {
        if (Params::getInstance().num_threads == 0)
            outError("Too many threads may slow down analysis [-nt option]. Reduce threads");
        else
            outError("Too many threads may slow down analysis [-nt option]. Reduce threads or use -nt AUTO to automatically determine it");
    }
}

double PhyloTree::getCurScore() {
    return curScore;
}

void PhyloTree::setCurScore(double score) {
    this->curScore = score;
}

void PhyloTree::resetCurScore(double score) {
    if (score != 0.0) {
        curScore = score;
    }
    else {
        curScore = -DBL_MAX;
    }
    if (model && this->tree_buffers.buffer_partial_lh != nullptr) {
        initializeAllPartialLh();
    }
}

double PhyloTree::getGammaShape() const {
    return site_rate ? site_rate->getGammaShape() : 1.0;
}

double* PhyloTree::getPatternLhCatPointer() {
    return tree_buffers._pattern_lh_cat;
}

void PhyloTree::initPartitionInfo() {
}
