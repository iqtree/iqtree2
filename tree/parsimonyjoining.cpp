//
//  parsimonyjoining.cpp
//  This is a Neighbour Joining algorithm, that uses parsimony
//  cost (for joining parsimony subtrees).
//  Created by James Barbetti on 31-Oct-2020.
//

#include "phylotree.h"
#include <placement/parallelparsimonycalculator.h>
#include <placement/blockallocator.h>
#include <placement/taxontoplace.h>
#include "phylotreethreadingcontext.h"
#include <utils/rapidnj.h>
#include <utils/auctionmatrix.h>
#include <utils/timekeeper.h>
#include <sprng/sprng.h>

class ParsimonyMatrix: public StartTree::NJMatrix<NJFloat> {
protected:
    typedef StartTree::NJMatrix<NJFloat> super;
    typedef NJFloat T;
    using   super::rows;
    using   super::row_count;
    using   super::clusters;
    using   super::rowToCluster;
    using   super::rowTotals;
    
    PhyloTree*       tree;
    PhyloNeighborVec topOfCluster;
    PhyloNode*       last_interior_node; //The last interior node
                                         //connected to the tree.
    PhyloNode*       true_root;          //The node for sequence zero
    UINT*            next_partial_pars;
    
    T                rowTotalMultiplier;
    
public:
    ParsimonyMatrix(): tree(nullptr), last_interior_node(nullptr)
                     , true_root(nullptr), next_partial_pars(nullptr)
                     , rowTotalMultiplier(2.0) {
        //rowTotalMultiplier of 0.0: ignore row totals (so
        //basically a parsimony equivalent of Kruskal's Algorithm)
        //rowTotalMultiplier of 1.0 is like NJ
        //Higher multipliers settle the placement of outgroup sequences
        //sooner.
    }
    virtual std::string getAlgorithmName() const {
        return "PJ";
    }
    void setTree(PhyloTree* treeToUse) {
        tree = treeToUse;
        tree->ensureCentralPartialParsimonyIsAllocated(0);
        last_interior_node = nullptr;
        true_root          = nullptr;
        next_partial_pars  = tree->central_partial_pars;
        ASSERT(next_partial_pars && next_partial_pars < tree->tip_partial_pars);
    }
    virtual void calculateLeafParsimonies() {
        ParallelParsimonyCalculator calculator(*tree, false);
        int nseq = static_cast<int>(tree->aln->getNSeq());
        for (int i=0; i<nseq; ++i) {
            auto           leaf_name = tree->aln->getSeqName(i);
            PhyloNode*     leafNode  = tree->newNode(i, leaf_name.c_str());
            PhyloNeighbor* topNei    = new PhyloNeighbor(leafNode, -1);
                //Todo: should come via tree.newNeighbor(), but that doesn't
                //      exist.
            allocateParsimonyFor(topNei);
            calculator.schedulePartialParsimony(topNei, DUMMY_NODE_1);
            topOfCluster.emplace_back(topNei);
        }
        true_root = topOfCluster[0]->getNode();
        calculator.calculate(0, "Calculating leaf parsimony vectors");
    }
    virtual void calculateLeafParsimonyDistances() {
        #if USE_PROGRESS_DISPLAY
        double work_estimate  = (double)row_count*(double)(row_count-1)/2.0;
        const char* task_name = "Constructing leaf-leaf"
                                " parsimony distance matrix";
        progress_display progress(work_estimate, task_name);
        #else
        double progress = 0.0;
        #endif
        #if _OPENMP
        #pragma omp parallel for schedule(dynamic)
        #endif
        for (intptr_t r=0; r<row_count; ++r) {
            auto currRow   = rows[r];
            auto rowVector = topOfCluster[r]->partial_pars;
            for (intptr_t c=r+1; c<row_count; ++c) {
                auto colVector = topOfCluster[c]->partial_pars;
                int score;
                tree->computeParsimonyOutOfTree( rowVector, colVector, &score );
                currRow[c] = static_cast<T>(score);
            }
            progress += (row_count - r);
        }
        #if USE_PROGRESS_DISPLAY
        progress.done();
        #endif
        
        #if _OPENMP
        #pragma omp parallel for
        #endif
        for (intptr_t r=1; r<row_count; ++r) {
            auto currRow    = rows[r];
            for (intptr_t c=0; c<r; ++c) {
                currRow[c] = rows[c][r];
            }
        }
    }
    void allocateParsimonyFor(PhyloNeighbor* nei) {
        if (nei->partial_pars==nullptr) {
            //Note: Not thread-safe; only call from main thread
            ASSERT( next_partial_pars < tree->tip_partial_pars );
            nei->partial_pars    = next_partial_pars;
            nei->setParsimonyComputed(false);
            next_partial_pars += tree->pars_block_size;
        }
    }
    void allocateParsimonyForAll(PhyloNode* near_node, PhyloNode* back) {
        FOR_EACH_PHYLO_NEIGHBOR(near_node, back, it, nei) {
            allocateParsimonyFor(nei);
            allocateParsimonyForAll(nei->getNode(), near_node);
            PhyloNeighbor* backNei = nei->getNode()->findNeighbor(near_node);
            allocateParsimonyFor(backNei);
        }
    }
    virtual void cluster(intptr_t a, intptr_t b) {
        auto aRow       = rows[a];
        auto bRow       = rows[b];
        T cTotal        = 0;

        PhyloNode*     topNodeInCluster = tree->newNode();
        last_interior_node->addNeighbor(topNodeInCluster, -1);
        PhyloNeighbor* neighToCluster   = last_interior_node->firstNeighbor();
        topOfCluster.emplace_back(neighToCluster);
        last_interior_node->neighbors.resize(0);
        {
            allocateParsimonyFor(neighToCluster);
            
            PhyloNeighbor* neighborA = topOfCluster[rowToCluster[a]];
            topNodeInCluster->neighbors.emplace_back(neighborA);
            neighborA->getNode()->addNeighbor(topNodeInCluster, -1);

            PhyloNeighbor* neighborB = topOfCluster[rowToCluster[b]];
            topNodeInCluster->neighbors.emplace_back(neighborB);
            neighborB->getNode()->addNeighbor(topNodeInCluster, -1);
            
            tree->computePartialParsimonyOutOfTree
                ( neighborA->partial_pars, neighborB->partial_pars,
                  neighToCluster->partial_pars );
        }
        auto abVector = neighToCluster->partial_pars;
        auto abScore  = tree->getSubTreeParsimony(neighToCluster);
        
        #ifdef _OPENMP
        #pragma omp parallel for reduction(+:cTotal)
        #endif
        for (intptr_t i=0; i<row_count; ++i) {
            if (i!=a && i!=b) {
                T Dai         = aRow[i];
                T Dbi         = bRow[i];

                //Parsimony distance
                auto iCluster = topOfCluster[rowToCluster[i]];
                auto iVector  = iCluster->partial_pars;
                int  iScore   = tree->getSubTreeParsimony(iCluster);
                int score     = 0;
                tree->computeParsimonyOutOfTree( abVector, iVector, &score );
                score         = score - abScore - iScore;                 

                T Dci         = static_cast<T>(score);

                aRow[i]       = Dci;
                rows[i][a]    = Dci;
                rowTotals[i] += Dci - Dai - Dbi;
                                //JB2020-06-18 Adjust row totals on fly
                cTotal       += Dci;
            }
        }
        //Negative lengths will be corrected later
        clusters.addCluster ( rowToCluster[a], -1,
                              rowToCluster[b], -1);
        rowTotals[a]    = cTotal;
        rowToCluster[a] = clusters.size()-1;
        rowToCluster[b] = rowToCluster[row_count-1];
        removeRowAndColumn(b);
    }
    virtual void finishClustering() {
        ASSERT( row_count == 3);
        for (size_t i=0;i<3;++i) {
            PhyloNeighbor* clusterTopNei = topOfCluster[rowToCluster[i]];
            last_interior_node->neighbors.emplace_back( clusterTopNei );
            clusterTopNei->getNode()->addNeighbor(last_interior_node, -1);
        }
        clusters.addCluster
            ( rowToCluster[0], -1
            , rowToCluster[1], -1
            , rowToCluster[2], -1);
        row_count      = 0;
        tree->leafNum  = tree->aln->getNSeq32();
        tree->root     = true_root;
        tree->rooted   = false;
        tree->initializeTree();
        tree->setAlignment(tree->aln);
        allocateParsimonyForAll(last_interior_node, nullptr);
    }
    virtual void calculateScaledRowTotals() const {
        super::calculateScaledRowTotals();
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (intptr_t r=0; r<row_count; ++r) {
            scaledRowTotals[r] *= rowTotalMultiplier;
        }
    }
    virtual void join() {
        last_interior_node = tree->newNode();
        size_t n           = tree->aln->getNSeq();
        auto   names       = tree->aln->getSeqNames();
        setSize(n);
        for (auto it = names.begin(); it != names.end(); ++it) {
            clusters.addCluster(*it);
        }
        calculateLeafParsimonies();
        calculateLeafParsimonyDistances();
        calculateRowTotals();
        constructTree();
        last_interior_node = nullptr;
    }
};

class RapidParsimonyMatrix: public StartTree::BoundingMatrix<NJFloat, ParsimonyMatrix>  {
public:
    typedef StartTree::BoundingMatrix<NJFloat, ParsimonyMatrix> super;
    RapidParsimonyMatrix(PhyloTree& tree) {
        setTree(&tree);
    }
    virtual std::string getAlgorithmName() const {
        return "RapidPJ";
    }
};

class AuctionParsimonyMatrix:public StartTree::AuctionMatrix<NJFloat, ParsimonyMatrix> {
public:
    AuctionParsimonyMatrix(PhyloTree& tree) {
        setTree(&tree);
    }
    virtual std::string getAlgorithmName() const {
        return "AuctionPJ";
    }
};

typedef RapidParsimonyMatrix ParsimonyJoiningMatrixType;

class RoutedTaxon: public TaxonToPlace {
public:
    int best_branch;
    int best_score;
    bool operator < (const RoutedTaxon& r) const {
        return best_branch < r.best_branch;
    }
    bool operator <= (const RoutedTaxon& r) const {
        return best_branch <= r.best_branch;
    }
};

class ParsimonyRouter {
    //Todo: allow caller to pass in rand_stream
public:
    PhyloTree&                 tree;
    PhyloTreeThreadingContext& context;
    int*                       rand_stream;
    bool                       owner_of_stream;
    PhyloBranchVector          branches;     //branches added to tree
    IntVector                  sample;       //indicates indices of branches in sample

    ParsimonyRouter(PhyloTree& tree_to_route,
                    PhyloTreeThreadingContext &context_to_use)
        : tree(tree_to_route), context(context_to_use) {
        rand_stream = init_sprng(0, 1, tree.params->ran_seed,
                                 SPRNG_DEFAULT);
        owner_of_stream = true;
    }
    ParsimonyRouter(PhyloTree& tree_to_route,
                    PhyloTreeThreadingContext &context_to_use,
                    int* stream)
        : tree(tree_to_route), context(context_to_use)
        , rand_stream(stream), owner_of_stream(false) {
    }
    virtual ~ParsimonyRouter() {
        if (owner_of_stream) {
            if (rand_stream!=nullptr) {
                free_sprng(rand_stream);
                rand_stream = nullptr;
            }
            owner_of_stream = false;
        }
    }
    void selectSample(int sample_size, int branch_count) {
        if (sample_size==branch_count) {
            sample.resize(branch_count);
            for (int i=0; i<branch_count; ++i) {
                sample[i] = i;
            }
            return;
        }
        sample.resize(sample_size);
        int g = 0;
        do {
            //This is a junk way to choose a sample
            //but it's adequate if sample_size is
            //significantly smaller than branch_count.
            for (size_t i = g; i<sample_size; ++i) {
                sample[i] = random_int(branch_count);
            }
            std::sort(sample.begin(), sample.end());
            //keep unique values, and count them (with g)
            g = 1;
            for (int h=1; h<sample_size; ++h) {
                if (sample[h-1]<sample[h]) {
                    sample[g]=sample[h];
                    ++g;
                }
            }
            //g now, index after last non-duplicate value
        } while (g<sample_size);
    }
    void bestSampleBranch(RoutedTaxon& candidate,
                          UINT* scratch_vector) const {
        auto     pars        = candidate.getParsimonyBlock();
        intptr_t sample_size = sample.size();
        for (intptr_t s = 0; s<sample_size; ++s) {
            PhyloBranch    b(branches[sample[s]]);
            PhyloNeighbor* nei1 = b.getLeftNeighbor();
            PhyloNeighbor* nei2 = b.getRightNeighbor();
            tree.computePartialParsimonyOutOfTree(nei1->get_partial_pars(),
                                                  nei2->get_partial_pars(),
                                                  scratch_vector);
            int branch_score = 0;
            tree.computeParsimonyOutOfTree(pars, scratch_vector,
                                           &branch_score);
            if ( s == 0 || branch_score<candidate.best_score) {
                candidate.best_score  = branch_score;
                candidate.best_branch = sample[s];
            }
        }
        /*TREE_LOG_LINE(tree, VerboseMode::VB_MIN, "C Taxon " << candidate.new_leaf->id
                      << " @ Branch " << best_branch
                      << " with score " << best_branch_score );*/
    }
    inline static bool doesFirstTouch(const PhyloBranch& b1,
                                      const PhyloBranch& b2) {
        return (b1.first == b2.first || b1.first == b2.second);
    }
    int routeCandidate(RoutedTaxon& taxon, UINT* scratch_vector) const {
        int prev_branch = -1;
        int dest_branch_id = taxon.best_branch;
        
        auto pars = taxon.getParsimonyBlock();
        for (;;) {
            int best_branch_id = dest_branch_id;
            PhyloBranch target = branches[dest_branch_id];
            IntVector neighboring_branch_ids;
            FOR_EACH_PHYLO_NEIGHBOR(target.first,  target.second, it, nei) {
                if (nei->id != prev_branch) {
                    neighboring_branch_ids.push_back(nei->id);
                }
            }
            FOR_EACH_PHYLO_NEIGHBOR(target.second, target.first,  it, nei) {
                if (nei->id != prev_branch) {
                    neighboring_branch_ids.push_back(nei->id);
                }
            }
            int equal_scores       = 0;
            int best_subtree_score = -1;
            for (auto id : neighboring_branch_ids) {
                PhyloBranch    b    = branches[id];
                PhyloNeighbor* nei1 = b.getLeftNeighbor();
                PhyloNeighbor* nei2 = b.getRightNeighbor();
                tree.computePartialParsimonyOutOfTree(nei1->get_partial_pars(),
                                                      nei2->get_partial_pars(),
                                                      scratch_vector);
                int branch_score = 0;
                tree.computeParsimonyOutOfTree(pars, scratch_vector, &branch_score);
                auto touching_nei = doesFirstTouch(b, target) ? nei1 : nei2;
                int subtree_score = tree.getSubTreeParsimony(touching_nei);
                if (branch_score<taxon.best_score) {
                    taxon.best_score   = branch_score;
                    best_subtree_score = subtree_score;
                    best_branch_id     = id;
                    equal_scores       = 0;
                    /*TREE_LOG_LINE(tree, VerboseMode::VB_MIN, "I Taxon " << taxon.new_leaf->id
                                  << " @ Branch " << best_branch_id
                                  << " with score " << best_score );*/
                } else if (branch_score==taxon.best_score) {
                    if (subtree_score < best_subtree_score ||
                        best_subtree_score == -1) {
                        best_branch_id     = id;
                        best_subtree_score = subtree_score;
                    }
                    ++equal_scores;
                }
            }
            if (best_branch_id == dest_branch_id || 1 < equal_scores) {
                return best_branch_id;
            }
            prev_branch    = dest_branch_id;
            dest_branch_id = best_branch_id;
        }
    }
    void updateBranch(intptr_t branch_id) {
        PhyloBranch    branch = branches[branch_id];
        PhyloNeighbor* front  = branch.first->findNeighbor(branch.second);
        front->id             = static_cast<int>(branch_id);
        PhyloNeighbor* back   = branch.second->findNeighbor(branch.first);
        back->id              = static_cast<int>(branch_id);
    }
    void insertCandidate(RoutedTaxon& taxon,
                         BlockAllocator& block_allocator) {
        PhyloBranch target = branches[taxon.best_branch];
        taxon.new_interior->addNeighbor(target.first, -1);
        taxon.new_interior->addNeighbor(target.second, -1);
        target.first->updateNeighbor(target.second, taxon.new_interior, -1);
        target.second->updateNeighbor(target.first, taxon.new_interior, -1);
        //
        // A
        //  \
        //   C--D  (where C is taxon.new_interior, D is taxon.new_leaf)
        //  /
        // B
        //
        PhyloNeighbor* AC = target.first->findNeighbor(taxon.new_interior);   //used to be AB
        PhyloNeighbor* CB = taxon.new_interior->findNeighbor(target.second);
        block_allocator.allocateMemoryFor(CB);
        std::swap(AC->partial_pars, CB->partial_pars);
        CB->setParsimonyComputed(AC->isParsimonyComputed());

        PhyloNeighbor* BC = target.second->findNeighbor(taxon.new_interior); //used to be BA
        PhyloNeighbor* CA = taxon.new_interior->findNeighbor(target.first);
        block_allocator.allocateMemoryFor(CA);
        std::swap(BC->partial_pars, CA->partial_pars);
        CA->setParsimonyComputed(BC->isParsimonyComputed());
        
        PhyloNeighbor* DC = taxon.new_leaf->findNeighbor(taxon.new_interior);
        block_allocator.allocateMemoryFor(DC);
        DC->setParsimonyComputed(false);

        taxon.new_interior->clearReversePartialParsimony(taxon.new_leaf);
        taxon.new_leaf->clearReversePartialParsimony(taxon.new_interior);
                
        intptr_t next_branch_id = static_cast<intptr_t>(branches.size());
        branches.emplace_back(target.first,  taxon.new_interior);
        updateBranch(next_branch_id);
        
        ++next_branch_id;
        branches.emplace_back(target.second, taxon.new_interior);
        updateBranch(next_branch_id);
        
        target.first  = taxon.new_leaf;
        target.second = taxon.new_interior;
        branches[taxon.best_branch] = target;
        updateBranch(taxon.best_branch);
        
        tree.nodeNum   += 2;
        tree.branchNum += 2;
        ++tree.leafNum;
    }
    void constructTree() {
        int         nseq          = tree.aln->getNSeq32();
        double      work_estimate = (double)nseq * sqrt((nseq) + 12);
        const char* task          = "Constructing tree with Parsimony Routing";
        tree.initProgress(work_estimate, task, "", "");

        TimeKeeper initializing ("Initializing");
        initializing.start();
        
        tree.freeNode();
        tree.deleteAllPartialLhAndParsimony();
        IntVector taxon_order;
        tree.create3TaxonTree(taxon_order, rand_stream);
        tree.initializeTree();
        tree.setParsimonyKernel(tree.params->SSE);
        tree.ensureCentralPartialParsimonyIsAllocated(tree.num_threads);
        int index_parsimony = 0;
        tree.initializeAllPartialPars(index_parsimony);
        BlockAllocator block_allocator(tree, index_parsimony);
        std::vector<UINT*> buffer;
        block_allocator.allocateVectorOfParsimonyBlocks(tree.num_threads, buffer);

        double parsimony_score = tree.computeParsimony("Computing pre-PR parsimony",
                                                       true, false);
        
        PhyloNodeVector v1, v2;
        tree.getBranchesInIDOrder(v1, v2);
        intptr_t branch_count = v1.size();
        for (intptr_t i=0; i<branch_count; ++i) {
            branches.emplace_back(v1[i], v2[i]);
        }
        
        TypedTaxaToPlace<RoutedTaxon> candidates;
        candidates.resize(nseq);
        int first_new_interior_id = tree.renumberInternalNodes();
        for (int i=3; i<nseq; ++i) {
            int         taxonId   = taxon_order[i];
            std::string taxonName = tree.aln->getSeqName(taxonId);
            candidates[i].initialize(&block_allocator, first_new_interior_id+i-3, taxonId, taxonName, false);
            candidates[i].new_interior->id = static_cast<int>(i + nseq - 2);
        }

        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (intptr_t i=3; i<nseq; ++i) {
            candidates[i].computeParsimony(&tree);
        }
        
        initializing.stop();

        TimeKeeper sampling ("Sampling branches");
        TimeKeeper rescoring("Rescoring parsimony");
        TimeKeeper searching("Searching for taxa insertion positions");
        TimeKeeper inserting("Inserting taxa");

        int sample_count = 2;
        double total_work = 0;
        for (int stop_batch=3; stop_batch<nseq; ) {
            
            sampling.start();
            ++sample_count;
            selectSample(sample_count, (int)branches.size());
            sampling.stop();
            
            int start_batch = stop_batch;
            stop_batch += sample_count/4 + 1;
            if (nseq < stop_batch) {
                stop_batch = nseq;
            }

            searching.start();
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for (int i=start_batch;i<stop_batch;++i) {
                candidates[i].computeParsimony(&tree);
                int t = context.getThreadNumber();
                bestSampleBranch(candidates[i], buffer[t]);
                routeCandidate(candidates[i], buffer[t]);
            }
            searching.stop();

            inserting.start();
            std::sort(candidates.begin()+start_batch,
                      candidates.begin()+stop_batch);
            
            for (int j=start_batch; j<stop_batch; ++j) {
                insertCandidate(candidates[j], block_allocator);
            }
            inserting.stop();
            double work_done = (double)(stop_batch - start_batch)
                             * (double) (sample_count+2);
            tree.trackProgress(work_done);
            total_work += work_done;
            
            rescoring.start();
            parsimony_score = tree.computeParsimony("Computing post-batch parsimony",
                                                    true, false);
            rescoring.stop();
        }
        
        TREE_LOG_LINE(tree, VerboseMode::VB_MIN, 
                      "Score " << parsimony_score
                      << " total_work " << total_work);
        if (VerboseMode::VB_MED <= verbose_mode) {
            tree.hideProgress();
            std::cout.precision(4);
            initializing.report();
            sampling.report();
            rescoring.report();
            searching.report();
            inserting.report();
            tree.showProgress();
        }
        tree.doneProgress();
    }
};

int PhyloTree::joinParsimonyTree(const char *out_prefix,
                                 Alignment *alignment) {
    aln = alignment;
    size_t nseq = aln->getNSeq();
    if (nseq < 3) {
        outError(ERR_FEW_TAXA);
    }
    PhyloTreeThreadingContext context(*this, params->parsimony_uses_max_threads);

    if (1) {
        ParsimonyJoiningMatrixType pjm(*this);
        pjm.join();
    } else {
        ParsimonyRouter pr(*this, context);
        pr.constructTree();
    }
    deleteAllPartialParsimony();
    initializeAllPartialPars();

    /* how long does this take?! */
    double parsimony_score = computeParsimony("Computing post PJ parsimony",
                                              true, false);
    setAllBranchLengthsFromParsimony(false, parsimony_score);
    
    // convert to rooted tree if originally so
    if (out_prefix) {
        string file_name = out_prefix;
        file_name += ".parstree";
        printTree(file_name.c_str(), WT_NEWLINE + WT_BR_LEN);
    }
    return static_cast<int>(floor(parsimony_score));
}
