//
//  parsimonyjoining.cpp
//  alignment
//
//  Created by James Barbetti on 31/10/20.
//

#include "phylotree.h"
#include <placement/parallelparsimonycalculator.h>
#include <utils/rapidnj.h>

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
    PhyloNode*       dummyRoot;
    UINT*            next_partial_pars;
    
    double rowTotalMultiplier;
    
public:
    ParsimonyMatrix(): tree(nullptr), dummyRoot(nullptr), next_partial_pars(nullptr)
                    , rowTotalMultiplier(2.0){
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
        dummyRoot = nullptr;
        next_partial_pars = tree->central_partial_pars;
        ASSERT(next_partial_pars && next_partial_pars < tree->tip_partial_pars);
    }
    virtual void calculateLeafParsimonies() {
        ParallelParsimonyCalculator calculator(*tree);
        for (size_t i=0; i<tree->aln->getNSeq(); ++i) {
            std::string&   leaf_name = tree->aln->getSeqName(i);
            PhyloNode*     leafNode  = tree->newNode(i, leaf_name.c_str());
            PhyloNeighbor* topNei    = new PhyloNeighbor(leafNode, -1);
                //Todo: should come via tree.newNeighbor(), but that doesn't
                //      exist.
            allocateParsimonyFor(topNei);
            calculator.computePartialParsimony(topNei, DUMMY_NODE_1);
            topOfCluster.emplace_back(topNei);
        }
        calculator.calculate(0, "Calculating leaf parsimony vectors");
    }
    virtual void calculateLeafParsimonyDistances() {
        progress_display progress(row_count*(row_count-1)/2,
                                  "Constructing leaf-leaf parsimony distance matrix");
        #if _OPENMP
        #pragma omp parallel for schedule(dynamic)
        #endif
        for (size_t r=0; r<row_count; ++r) {
            auto currRow   = rows[r];
            auto rowVector = topOfCluster[r]->partial_pars;
            for (size_t c=r+1; c<row_count; ++c) {
                auto colVector = topOfCluster[c]->partial_pars;
                int score;
                tree->computeParsimonyOutOfTree( rowVector, colVector, &score );
                currRow[c] = score;
            }
            progress += (row_count - r);
        }
        progress.done();
        #if _OPENMP
        #pragma omp parallel for
        #endif
        for (size_t r=1; r<row_count; ++r) {
            auto currRow    = rows[r];
            for (size_t c=0; c<r; ++c) {
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
    virtual void cluster(size_t a, size_t b) {
        auto aRow       = rows[a];
        auto bRow       = rows[b];
        T cTotal        = 0;

        PhyloNode*     topNodeInCluster = tree->newNode();
        dummyRoot->addNeighbor(topNodeInCluster, -1);
        PhyloNeighbor* neighToCluster   = dummyRoot->firstNeighbor();
        topOfCluster.emplace_back(neighToCluster);
        dummyRoot->neighbors.resize(0);
        {
            allocateParsimonyFor(neighToCluster);
            
            PhyloNeighbor* neighborA = topOfCluster[rowToCluster[a]];
            topNodeInCluster->neighbors.emplace_back(neighborA);
            neighborA->getNode()->addNeighbor(topNodeInCluster, -1);

            PhyloNeighbor* neighborB = topOfCluster[rowToCluster[b]];
            topNodeInCluster->neighbors.emplace_back(neighborB);
            neighborB->getNode()->addNeighbor(topNodeInCluster, -1);
            
            tree->computePartialParsimonyOutOfTree
                ( neighborA->partial_pars, neighborB->partial_pars, neighToCluster->partial_pars );
        }
        auto abVector = neighToCluster->partial_pars;
        
        #ifdef _OPENMP
        #pragma omp parallel for reduction(+:cTotal)
        #endif
        for (size_t i=0; i<row_count; ++i) {
            if (i!=a && i!=b) {
                T Dai         = aRow[i];
                T Dbi         = bRow[i];

                //Parsimony distance
                auto iVector  = topOfCluster[rowToCluster[i]]->partial_pars;
                int score     = 0;
                tree->computeParsimonyOutOfTree( abVector, iVector, &score );
                T Dci         = score;

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
            dummyRoot->neighbors.emplace_back( clusterTopNei );
            clusterTopNei->getNode()->addNeighbor(dummyRoot, -1);
        }
        clusters.addCluster
            ( rowToCluster[0], -1
            , rowToCluster[1], -1
            , rowToCluster[2], -1);
        row_count      = 0;
        tree->leafNum  = tree->aln->getNSeq();
        tree->root     =  topOfCluster[0]->getNode();
        tree->rooted   = false;
        tree->initializeTree();
        tree->setAlignment(tree->aln);
        allocateParsimonyForAll(dummyRoot, nullptr);
    }
    virtual void calculateScaledRowTotals() const {
        super::calculateScaledRowTotals();
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (size_t r=0; r<row_count; ++r) {
            scaledRowTotals[r] *= rowTotalMultiplier;
        }
    }
    virtual void join() {
        dummyRoot  = tree->newNode();
        size_t n   = tree->aln->getNSeq();
        auto names = tree->aln->getSeqNames();
        setSize(n);
        for (auto it = names.begin(); it != names.end(); ++it) {
            clusters.addCluster(*it);
        }
        calculateLeafParsimonies();
        calculateLeafParsimonyDistances();
        calculateRowTotals();
        constructTree();
        dummyRoot = nullptr;
    }
};

class RapidParsimonyMatrix: public StartTree::BoundingMatrix<NJFloat, ParsimonyMatrix>  {
public:
    typedef StartTree::BoundingMatrix<NJFloat, ParsimonyMatrix> super;
    RapidParsimonyMatrix(PhyloTree& tree) {
        setTree(&tree);
    }
};

class AuctionParsimonyMatrix: public StartTree::BoundingMatrix<NJFloat, ParsimonyMatrix>  {
    //
    //A reverse auction algorithm: each row "bids" to be clustered
    //(lowest bid wins) based on the difference between the best join
    //involving it (and another lower-numbered cluster) and the
    //second-best.  if there aren't two lower-numbered clusters,
    //no bid is placed (infiniteDistance is returned).
    //
    //This is a little faster than actually looking at all of the
    //parsimony distances!  And it practice, yields answers that aren't
    //that much worse.
    //
public:
    typedef StartTree::BoundingMatrix<NJFloat, ParsimonyMatrix> super;
    AuctionParsimonyMatrix(PhyloTree& tree) {
        setTree(&tree);
    }
    StartTree::Position<T> getRowMinimum(size_t row, T maxTot, T qBest) const {
        StartTree::Position<T> pos(row, 0, infiniteDistance, 0);
        const T*   rowData   = entriesSorted.rows[row];
        const int* toCluster = entryToCluster.rows[row];
        size_t i = 0;
        int   bestOtherRow;
        for (i=0; rowData[i]<infiniteDistance; ++i) {
            bestOtherRow = clusterToRow[toCluster[i]];
            if (bestOtherRow != notMappedToRow) {
                break;
            }
        }
        if (rowData[i] < infiniteDistance) {
            T     bestV = rowData[i];
            int   secondBestOtherRow;
            for (; rowData[i]<infiniteDistance; ++i) {
                secondBestOtherRow = clusterToRow[toCluster[i]];
                if (secondBestOtherRow != notMappedToRow) {
                    break;
                }
            }
            if (rowData[i] < infiniteDistance) {
                pos.row   = bestOtherRow;
                pos.value = bestV - rowData[secondBestOtherRow];
            }
            if (1<i) {
                purgeRow(row);
            }
        }
        return pos;
    }
};

typedef RapidParsimonyMatrix ParsimonyJoiningMatrixType;

int PhyloTree::joinParsimonyTree(const char *out_prefix,
                                 Alignment *alignment) {
    aln = alignment;
    size_t nseq = aln->getNSeq();
    if (nseq < 3) {
        outError(ERR_FEW_TAXA);
    }

    ParsimonyJoiningMatrixType pjm(*this);
    pjm.join();

    /* how long does this take?! */
    fixNegativeBranch(true);
    // convert to rooted tree if originally so
    if (out_prefix) {
        string file_name = out_prefix;
        file_name += ".parstree";
        printTree(file_name.c_str(), WT_NEWLINE + WT_BR_LEN);
    }
    
    return computeParsimony();
}


