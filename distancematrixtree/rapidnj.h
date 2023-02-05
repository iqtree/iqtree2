//
//  rapidnj.h - RapidNJ and RapidBIONJ distance matrix tree construction.
//
//  BoundingNJ implementation loosely based on ideas from
//        https://birc.au.dk/software/rapidnj/.
//        Paper: "Inference of Large Phylogenies using Neighbour-Joining."
//               Martin Simonsen, Thomas Mailund, Christian N. S. Pedersen.
//               Communications in Computer and Information Science
//               (Biomedical Engineering Systems and Technologies:
//               3rd International Joint Conference, BIOSTEC 2010,
//               Revised Selected Papers), volume 127, pages 334-344,
//               Springer Verlag, 2011.
//        Tag:  [SMP2011].
//        (but, optionally, using a variance matrix, as in BIONJ, and
//        keeping the distance and variance matrices square -
//        they're not triangular because
//                  (i) *read* memory access patterns are more favourable
//                 (ii) *writes* don't require conditional transposition
//                      of the row and column coordinates (but their
//                      access patterns aren't as favourable, but
//                (iii) reads vastly outnumber writes)
//        (there's no code yet for removing duplicated rows either;
//        those that have distance matrix rows identical to earlier rows;
//        Rapid NJ "hates" them) (this is also covered in section 2.5)
//
//        The BoundingMatrix class adds branch-and-bound optimization
//        to *other* distance matrix implementations.  In this file, only
//        to NJ and BIONJ, via the RapidNJ and RapidBIONJ classes.
//        
//        It sets up auxiliary S and I matrices.  In each row:
//        S = unadjusted distances to clusters that were in play
//            when this row was set up (ascending order)
//        I = the cluster indices that corresponded to each of
//            the cells in S.
//
//        Rows of S and I are sorted via a mergesort.
//        (During matrix set up, a sequential mergesort, because
//         row construction is parallelized; later, during
//         clustering, a parallel one).
//        (S is sorted, I is permuted to match).
//        (S and I are the names of these matrices in [SMP2011]).
//
//Notes:  1.An SI matrix, of pair<T,int> would probably be better,
//          as that could be sorted faster.
//        2.An adaptive row-sorting routine could be used,
//          particularly if new SI rows (after cluster joins) were 
//          constructed (not merely in cluster index order) in an 
//          order "suggested" by the content of one (or both?) of the
//          *existing* SI rows.
//          (but... row sorting is only ~10% of running time)
//
//  This file, created by James Barbetti on 31-Oct-2020.
//  (But the bulk of the code in it was from bionj2.cpp,
//  which dates back to 18-Jun-2020).
//
//  LICENSE:
//* This program is free software; you can redistribute it and/or modify
//* it under the terms of the GNU General Public License as published by
//* the Free Software Foundation; either version 2 of the License, or
//* (at your option) any later version.
//*
//* This program is distributed in the hope that it will be useful,
//* but WITHOUT ANY WARRANTY; without even the implied warranty of
//* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//* GNU General Public License for more details.
//*
//* You should have received a copy of the GNU General Public License
//* along with this program; if not, write to the
//* Free Software Foundation, Inc.,
//* 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//

#ifndef rapidnj_h
#define rapidnj_h

#include "nj.h"
#include <utils/parallel_mergesort.h>
#include <utils/timeutil.h>  //for getRealTime

namespace StartTree
{

/**
 * @brief  Applies branch-and-bound trickery (as per the 
 *         RapidNJ paper - see details at the top of this file)
 * @tparam T the distance type
 * @tparam SUPER the superclass, that is to be "munged" into a branch-and-bound
 *         version that will (on average), run a lot faster. It is assumed that 
 *         SUPER is NJMatrix, inherits from NJMatrix, or is 
 *         "duck-type-equivalent" to NJMatrix.
 * @note   The mutable members are calculated repeatedly, from others, in 
 *         member functions marked as const. They're declared at the class 
 *         level so that they don't need to be reallocated over and over again.
 *         It's all about amortizing allocation and deallocation cost (the
 *         price paid for that is a slight reduction in readability. -James B)
 * @note   Mapping members to the RapidNJ papers:
 *         rows           is the D matrix
 *         entriesSorted  is the S matrix
 *         entryToCluster is the I matrix
 * @note   scaledMaxEarlierClusterTotal[c] is the largest row total for a 
 *         cluster with a lower number, than cluster c (where c is the index
 *         and the cluster number) (if c indicates a cluster for which there 
 *         are still rows in the distance matrix: call this an in-play cluster).
 *         This is a tighter bound, when searching for the minimum Qij... 
 *         and processing distances from cluster c to earlier clusters, 
 *         than the largest row total for ALL the live clusters.
 *         See section 2.5 of Simonsen, Mailund, Pedersen [2011].
 * @note   rowOrderChosen is a vector of int rather than bool because 
 *         simultaneous non-overlapping random-access writes to elements of 
 *         std::vector<bool> *can* interfere with each other (because 
 *         std::vector<bool> maps multiple nearby elements onto bitfields, 
 *         so the writes... *do* overlap) (ouch!).
 * @note   If _OPENMP is defined, it is assumed that the threadcount,
 *         as returned by omp_get_max_threads() won't change during
 *         execution! Or at least, it won't change between the construction
 *         of a BoundingMatrix<T,S> instance, and its destruction.
 * @note   several "hide the OMP" functions, getThreadCount() and
 *         getThreadNumber() should perhaps move to NJMatrix<T>.
 */
template <class T=NJFloat, class SUPER=BIONJMatrix<T>>
class BoundingMatrix: public SUPER
{
public:
    typedef SUPER super;
    typedef MirrorMergeSorter<T, int> Sorter;
protected:
    using super::rows;
    using super::row_count;
    using super::rowMinima;
    using super::rowTotals;
    using super::rowToCluster;
    using super::clusters;
    using super::silent;
    using super::isRooted;
    using super::finishClustering;
    using super::clusterDuplicates;
    using super::prepareToConstructTree;

    std::vector<intptr_t> clusterToRow;   //Maps clusters to their rows (-1 means not mapped)
                                          //for each cluster
    std::vector<T>        clusterTotals;  //"Row" totals indexed by cluster

    mutable std::vector<intptr_t> clusterPartners;//Counts the number of possible partners,
    mutable std::vector<T>        scaledClusterTotals;   //The same, multiplied by
                                                        //(1.0 / (n-2)).
    mutable std::vector<T>        scaledMaxEarlierClusterTotal;
    mutable std::vector<int>      rowOrderChosen; //Indicates if a row's scanning
                                                 //order chosen has been chosen.
                                                 //Only used in... getRowScanningOrder().
    mutable std::vector<size_t>   rowScanOrder;   //Order in which rows are to be scanned
                                                 //Only used in... getRowMinima().
    
    SquareMatrix<T>   entriesSorted; //The S matrix: Entries in distance matrix
                                     //(each row sorted by ascending value)
    SquareMatrix<int> entryToCluster;//The I matrix: for each entry in S, which
                                     //cluster the row (that the entry came from)
                                     //was mapped to (at the time).
    int                 threadCount;
    std::vector<Sorter> sorters;
    
public:
    typedef T distance_type;
    BoundingMatrix() : super() {
        #ifdef _OPENMP
            threadCount = omp_get_max_threads();
        #else
            threadCount = 1;
        #endif
        sorters.resize(threadCount);
    }
    /**
     * @brief  Get the (current) thread count (returns 1 if _OPENMP is
     *         not defined)
     * @return intptr_t - the current rhead count
     */
    intptr_t getThreadCount() const {
        intptr_t t=1;
        #ifdef _OPENMP
            #pragma omp parallel
            {
                if (omp_get_thread_num()==0) {
                    t = omp_get_num_threads();
                }
            }
        #endif
        return t;
    }
    /**
     * @brief  Get the (current) threads thread number, between 0 
     *         inclusive and getThreadCount() exclusive (well, so we
     *         hope! If the number of threads hasn't changed since 
     *         we last called getThreadCount()! -James B)
     *         (returns 0 if _OPENMP is not defined)
     * @return intptr_t - the current rhead count
     */
    intptr_t getThreadNumber() const {
        #ifdef _OPENMP
            return omp_get_thread_num();
        #else
            return 0;
        #endif
    }
    virtual std::string getAlgorithmName() const {
        return "Rapid" + super::getAlgorithmName();
    }

    /**
     * @brief  Sets up auxiliary arrays, and infers a phylogenetic tree.
     * @return true 
     * @return false 
     */
    virtual bool constructTree() {
        //0. Ensure variance matrix is initialized
        prepareToConstructTree(); //Need this if SUPER is BIONJMatrix

        //1. Set up vectors indexed by cluster number,
        clusterToRow.resize(row_count);
        clusterTotals.resize(row_count);
        clusterPartners.resize(row_count*2-2);
        for (intptr_t r=0; r<row_count; ++r) {
            clusterToRow[r]  = static_cast<int>(r);
            clusterTotals[r] = rowTotals[r];
        }
        
        //2. Set up "scratch" vectors used in getRowMinima
        //   so that it won't be necessary to reallocate them
        //   for each call.
        scaledClusterTotals.resize(row_count);
        scaledMaxEarlierClusterTotal.resize(row_count);
        rowOrderChosen.resize(row_count);
        rowScanOrder.resize(row_count);

        {
            #if USE_PROGRESS_DISPLAY
            const char* taskName = silent
                ? "" :  "Setting up auxiliary I and S matrices";
            progress_display setupProgress(row_count, taskName, "sorting", "row");
            #endif
            //2. Set up the matrix with row sorted by distance
            //   And the matrix that tracks which distance is
            //   to which cluster (the S and I matrices, in the
            //   RapidNJ papers).
            entriesSorted.setSize(row_count);
            entryToCluster.setSize(row_count);
            #ifdef _OPENMP
            #pragma omp parallel num_threads(threadCount)
            #endif
            {
                int threadNum = getThreadNumber();
                for (intptr_t r=threadNum; r<row_count; r+=threadCount) {
                    sortRow(r,r,false,sorters[threadNum]);
                    //copies the "left of the diagonal" portion of
                    //row r from the D matrix and sorts it
                    //into ascending order.
                    #if USE_PROGRESS_DISPLAY
                        ++setupProgress;
                    #endif
                }
            }
        }
        clusterDuplicates();
        {
            intptr_t nextPurge = (row_count+row_count)/3;
            #if USE_PROGRESS_DISPLAY
            std::string taskName;
            if (!silent) {
                taskName = "Constructing " + getAlgorithmName() + " tree";
            }
            double triangle = row_count * (row_count + 1.0) * 0.5;
            progress_display show_progress(triangle, taskName.c_str(), "", "");
            #endif
            intptr_t degree_of_root = isRooted ? 2 : 3;
            while (degree_of_root<row_count) {
                Position<T> best;
                super::getMinimumEntry(best);
                cluster(best.column, best.row);
                if ( row_count == nextPurge ) {
                    #ifdef _OPENMP
                    #pragma omp parallel num_threads(threadCount)
                    #endif
                    {
                        for (intptr_t r=getThreadNumber(); 
                            r<row_count; r+=threadCount) {
                            purgeRow(r);
                        }
                    }
                    nextPurge = (row_count + row_count)/3;
                }
                #if USE_PROGRESS_DISPLAY
                show_progress+=row_count;
                #endif
            }
            #if USE_PROGRESS_DISPLAY
            show_progress.done();
            #endif
            finishClustering();
        }
        return true;
    }    
    /**
     * @brief Calculate row r of the S and I matrices, from row r
     *        of the D matrix.
     * @param r        - the row number
     * @param c        - upper bound on cluster index (but NOT a cluster index!)
     * @param parallel - true if a the sort should be parallel (when running
     *                   in code that was multi-threaded higher up, this should
     *                   be false.  When running in code that wasn't this could
     *                   be true; it is up to the caller).
     * @param sorter   - a sorter (a sorting algorithm implementation, perhaps
     *                   with memory allocated for auxiliary arrays that it uses
     *                   to sort), that is "owned" by the current thread.
     * @note 
     */
    void sortRow(intptr_t r /*row index*/, int c /*upper bound on cluster index*/
        ,  bool parallel, Sorter& sorter) {
        //1. copy data from a row of the D matrix into the S matrix
        //   (and write the cluster identifiers that correspond to
        //    the values in the D row into the same-numbered
        //    row in the I matrix), for distances between the cluster
        //    in that row, and other live clusters (up to, but not including c).
        T*       sourceRow      = rows[r];
        T*       values         = entriesSorted.rows[r];
        int*     clusterIndices = entryToCluster.rows[r];
        intptr_t w = 0;
        for (intptr_t i=0; i<row_count; ++i) {
            values[w]         = sourceRow[i];
            clusterIndices[w] = static_cast<int>(rowToCluster[i]);
            if ( i != r && clusterIndices[w] < c ) {
                ++w;
            }
        }
        values[w]         = infiniteDistance; //sentinel value, to stop row search
        clusterIndices[w] = static_cast<int>(rowToCluster[r]);
            //Always room for this, because distance to self
            //was excluded via the i!=r check above.
        
        //2. Sort the row in the S matrix and mirror the sort
        //   on the same row of the I matrix.
        if (parallel) {
            sorter.parallel_mirror_sort(values, w, clusterIndices);
        } else {
            sorter.single_thread_mirror_sort(values, w, clusterIndices);
        }
        clusterPartners[c] = w;
    }

    /**
     * @brief Scan a row of the I matrix, so as to remove entries that refer 
     *        to clusters that are no longer in-play. Remove the corresponding
     *        values in the same row of the S matrix, at the same time.
     * @param r the row number
     */
    void purgeRow(intptr_t r /*row index*/) const {
        T*    values         = entriesSorted.rows[r];
        int*  clusterIndices = entryToCluster.rows[r];
        int   c              = rowToCluster[r];
        intptr_t w = 0;
        intptr_t i = 0;
        for (; i<row_count ; ++i ) {
            values[w]         = values[i];
            clusterIndices[w] = clusterIndices[i];
            if ( infiniteDistance <= values[i] ) {
                break;
            }
            bool liveRow = ( clusterToRow[clusterIndices[i]] != 
                             notMappedToRow);
            w += (liveRow ? 1 : 0);
        }
        if (w<row_count) {
            values[w] = infiniteDistance;
        }
        clusterPartners[c] = w;
    }
    /**
     * @brief Join the clusters (clusterA and clusterB, below) that correspond
     *        to rows a and b of the working distance matrix (0<=a<b<n, where
     *        n is the current rank of the working distance matrix).
     * @param a the lower-numbered row
     * @param b the higher-numbered row
     * @note  It is assumed a<b, but this is NOT checked.
     * @note  It is assumed that distances, for the merged cluster, 
     *        will be written over the top of row a of the D 
     *        matrix and the D matrix's rank (row count and column count),
     *        reduced by one, by super::cluster().
     *        Yes, this IS tight long-range coupling. -James B.
     */
    virtual void cluster(intptr_t a, intptr_t b) {
        size_t clusterA         = rowToCluster[a];
        size_t clusterB         = rowToCluster[b];
        size_t clusterMoved     = rowToCluster[row_count-1];
        clusterToRow[clusterA]  = notMappedToRow;
        clusterTotals[clusterA] = -infiniteDistance;
        clusterToRow[clusterB]  = notMappedToRow;
        clusterTotals[clusterB] = -infiniteDistance;
        size_t clusterC         = clusters.size(); //cluster # of new cluster

        super::cluster(a,b);

        if (b<row_count) {
            clusterToRow[clusterMoved] = static_cast<int>(b);
        }
        clusterToRow.emplace_back(a);
        clusterTotals.emplace_back(rowTotals[a]);
        scaledClusterTotals.emplace_back(rowTotals[a] / (T)( row_count - 1.0 ) );
        scaledMaxEarlierClusterTotal.emplace_back((T)0.0);
        //Mirror row rearrangement done on the D (distance) matrix
        //(and possibly also on the V (variance estimate) matrix),
        //onto the S and I matrices.
        entriesSorted.removeRowOnly(b);
        entryToCluster.removeRowOnly(b);
        
        //Recalculate cluster totals.
        for (size_t wipe = 0; wipe<clusterC; ++wipe) {
            clusterTotals[wipe] = -infiniteDistance;
            //A trick.  This way we don't need to check if clusters
            //are still "live" in the inner loop of getRowMinimum().
            //When we are "subtracting" cluster totals to calculate
            //entries in Q, they will come out so big they won't be
            //considered as candidates for neighbour join.
            //If we didn't do this we'd have to check, all the time,
            //when calculating entries in Q, if clusters are still
            //"live" (have corresponding rows in the D matrix).
        }
        for (intptr_t r = 0; r<row_count; ++r) {
            size_t cluster_num = rowToCluster[r];
            clusterTotals[cluster_num] = rowTotals[r];
        }
        sortRow(a, clusterC, true, sorters[0]);
    }
    /**
     * @brief Derive a bound (on the best adjusted distance), by looking
     *        at the first entry (for a live cluster), in each row of 
     *        the S matrix.
     * @param qBest 
     * @note  Since we always have to check these entries when we process
     *        the row, why not process them up front, hoping to get a
     *        better bound on min(V) (and perhaps even "rule" entire rows
     *        "out of consideration", using that bound)? -James B).
     */
    void deriveBoundFromFirstColumn(T& qBest) const {
        intptr_t rSize = rowMinima.size();
        std::vector<T> qLocalBestVector;
        qLocalBestVector.resize( threadCount, qBest);
        T* qLocalBest =  qLocalBestVector.data();

        #ifdef _OPEN_MP
        #pragma omp parallel for num_threads(threadCount)
        #endif
        for (int b=0; b<threadCount; ++b) {
            T      qBestForThread = qBest;
            int    b_plus_1       = (b + 1);
            size_t rStart         = (b*rSize)        / threadCount;
            size_t rStop          = (b_plus_1*rSize) / threadCount;
            for (size_t r=rStart; r < rStop
                    && rowMinima[r].value < infiniteDistance; ++r) {
                intptr_t rowA     = rowMinima[r].row;
                intptr_t rowB     = rowMinima[r].column;
                if (rowA < row_count && rowB < row_count ) {
                    size_t clusterA = rowToCluster[rowA];
                    size_t clusterB = rowToCluster[rowB];
                    T qHere = this->rows[rowA][rowB]
                            - scaledClusterTotals[clusterA]
                            - scaledClusterTotals[clusterB];
                    if (qHere < qBestForThread) {
                        qBestForThread = qHere;
                    }
                }
            }
            qLocalBest[b] = qBestForThread;
        }
        for (int b=0; b<threadCount; ++b) {
            if ( qLocalBest[b] < qBest ) {
                qBest = qLocalBest[b];
            }
        }
    }

    /**
     * @brief Partially sort the row minima (in parallel)
     * @note  The algorithm used here is roughly to compare and perhaps 
     *        reorder rows in the first half of the input (rounded up)
     *        with rows in the second (rounded down) (in parallel!), 
     *        then repeat...) to find the minimum.
     * @note  the point is that each compare-perhaps-reorder operation
     *        in a given "reduce the problem by half" pass is entirely
     *        independent of all the others, and they can be done in
     *        parallel.
     * @note  A "conditional-move" version might be faster than 
     *        swapping conditionally (if (rowMinima[j] < rowMinima[i])).
     *        (I don't remember whether I tried that. -James B)
     */
    void partiallySortRowMinima() const {
        intptr_t rSize     = rowMinima.size();
        #ifdef _OPENMP
        int      threshold = threadCount << 7; /* multiplied by 128*/
        #endif
        //Note, rowMinima might have size 0 (the first time this member
        //function is called during processing of a distance matrix)
        //Or it might have a size of n+1 (later times), but it won't be n.
        for ( intptr_t len = rSize; 1<len; len=(len+1)/2 ) {
            intptr_t halfLen = len/2; //rounded down
            intptr_t gap     = len-halfLen;
            #ifdef _OPENMP
            #pragma omp parallel for if(threshold<halfLen) num_threads(threadCount)
            #endif
            for ( intptr_t i = 0; i<halfLen; ++i) {
                intptr_t j = i + gap;
                if ( rowMinima[j] < rowMinima[i] ) {
                    std::swap(rowMinima[i], rowMinima[j]);
                }
            }
        }
    }

    /**
     * @brief Rig the order in which rows are scanned based on which rows 
     *        (might) have the lowest row minima based on what we saw last 
     *        time.
     * @param qBest 
     * @note  The original RapidNJ puts the second-best row from last time 
     *        first. And, apart from that, goes in (cyclical) row order.
     *        But rows in the D, S, and I matrices are (all) shuffled
     *        in memory, so why not do all the rows in ascending order
     *        of their best Q-values from the last iteration?
     *        Or, better yet... From this iteration?!
     *        The idea is to (hopefully) find better bounds sooner.
     */
    void decideOnRowScanningOrder(T& qBest) const {

        deriveBoundFromFirstColumn(qBest);
             
        partiallySortRowMinima();

        #ifdef _OPENMP
        #pragma omp parallel for if((threadCount<<7)<row_count) num_threads(threadCount)
        #endif
        for (intptr_t i = 0; i < row_count; ++i) {
            rowOrderChosen[i]=0; //Not chosen yet
        }
        
        intptr_t rSize = rowMinima.size();
        intptr_t w = 0;
        for (intptr_t r=0; r < rSize
             && rowMinima[r].row < row_count
             && rowMinima[r].value < infiniteDistance; ++r) {
            intptr_t rowA   = rowMinima[r].row;
            intptr_t rowB   = rowMinima[r].column;
            size_t clusterA = (rowA < row_count)    ? rowToCluster[rowA] : 0;
            size_t clusterB = (rowB < row_count)    ? rowToCluster[rowB] : 0;
            size_t row      = (clusterA < clusterB) ? rowA : rowB;
            if (row < (size_t)row_count) {
                rowScanOrder[w] = row;
                w += rowOrderChosen[row] ? 0 : 1;
                rowOrderChosen[row] = 1; //Chosen
            }
        }
        
        //The weird-looking middle term in the for loop is as
        //intended: when w reaches n all of the rows (0..n-1)
        //must be in rowScanOrder, so there's no need to continue
        //until row==n.
        for (intptr_t row=0; w < row_count ; ++row) {
            rowScanOrder[w] = row;
            w += ( rowOrderChosen[row] ? 0 : 1 );
        }
    }

    /**
     * @brief Get the Row Minima (that is, the minimum adjusted 
     *        inter-cluster distance, for each row, r), by looking up
     *        (for each row, r), the associated cluster, y, and 
     *        in the S row for y, finding the adjusted distance 
     *        for each other in-play cluster, x, that might be less
     *        than our current upper bound on the global minimum
     *        adjusted distance.
     * @note  Rather than multiplying distances by (n-2) repeatedly, 
     *        it is cheaper to work with scaled cluster totals,
     *        multiplied by (1.0/(T)(n-2)).
     *        Better n multiplications than 0.5*n*(n-1).
     */
    virtual void getRowMinima() const {
        size_t c           = clusters.size();
        T      nless2      = (T)( row_count - 2 );
        T      tMultiplier = ( row_count <= 2 ) ? (T)0.0 : ((T)1.0 / nless2);
        T      maxTot      = -infiniteDistance; //maximum row total divided by (n-2)
        for (size_t i=0; i<c; ++i) {
            scaledClusterTotals[i] = clusterTotals[i] * tMultiplier;
            scaledMaxEarlierClusterTotal[i] = maxTot;
            if ( clusterToRow[i] != notMappedToRow ) {
                if (maxTot < scaledClusterTotals[i] ) {
                    maxTot=scaledClusterTotals[i];
                }
            }
        }
        
        T qBest = infiniteDistance;
            //upper bound on minimum Q[row,col]
            //  = D[row,col] - R[row]*tMultipler - R[col]*tMultiplier
            //

        decideOnRowScanningOrder(qBest);
        rowMinima.resize(row_count);
        
        #ifdef _OPENMP
        #pragma omp parallel num_threads(threadCount)
        #endif
        {
            for (intptr_t r=getThreadNumber(); r<row_count; r+=threadCount) {
                size_t row             = rowScanOrder[r];
                size_t cluster_num     = rowToCluster[row];
                T      maxEarlierTotal = scaledMaxEarlierClusterTotal[cluster_num];
                rowMinima[r]           = getRowMinimum(row, maxEarlierTotal, qBest);
                if (rowMinima[r].value < qBest) {
                    qBest = rowMinima[r].value;
                }
            }
        }
    }

    /**
     * @brief Get the cluster merge description (Position<T>), with
     *        the minimum adjusted distance between the two merged clusters,
     *        that will merge a cluster, x, with a cluster number less 
     *        than that of the cluster, y, the cluster correpsonding to 
     *        row (row) in the current distance matrix, with y.
     * @param row    - row number in the current distance matrix
     * @param maxTot - the maximum (scaled) cluster distance total, 
     *                 for any cluster, with a cluster number less 
     *                 than y (used to get a slightly tighter bound,
     *                 on how bad a raw distance, read from S can be,
     *                 if it is still to be possible for the adjusted
     *                 distance, to the same cluster the raw distance
     *                 is to, to be equal to or better than qBest).
     * @param qBest  - a copy of an upper bound on the global minimum
     *                 (possibly found earlier).
     * @return Position<T> 
     */
    virtual Position<T> getRowMinimum(intptr_t row, T maxTot, T qBest) const {
        T nless2      = (T)( row_count - 2 );
        T tMultiplier = ( row_count <= 2 ) ? (T)0.0 : ( (T)1.0 / nless2 );
        auto    tot   = scaledClusterTotals.data();
        T rowTotal    = rowTotals[row] * tMultiplier; //scaled by (1/(n-2)).
        T rowBound    = qBest + maxTot + rowTotal;
                //Upper bound for distance, in this row, that
                //could (after row totals subtracted) provide a
                //better min(Q).

        Position<T> pos(row, 0, infiniteDistance, 0);
        const T*   rowData   = entriesSorted.rows[row];
        const int* toCluster = entryToCluster.rows[row];
        for (size_t i=0; ; ++i) {
            T Drc = rowData[i];
            if (rowBound<Drc && 0<i) {
                break;
            }
            size_t  cluster_num = toCluster[i];
                //The cluster associated with this distance
                //The c in Qrc and Drc.
            T Qrc = Drc - tot[cluster_num] - rowTotal;
            if (Qrc < pos.value) {
                intptr_t otherRow = clusterToRow[cluster_num];
                if (otherRow != notMappedToRow) {
                    pos.column = (otherRow < row ) ? otherRow : row;
                    pos.row    = (otherRow < row ) ? row : otherRow;
                    pos.value  = Qrc;
                    if (Qrc < qBest ) {
                        qBest    = Qrc;
                        rowBound = qBest + maxTot + rowTotal;
                    }
                }
            }
        }
        return pos;
    }
};

typedef BoundingMatrix<NJFloat,   NJMatrix<NJFloat>>    RapidNJ;
typedef BoundingMatrix<NJFloat,   BIONJMatrix<NJFloat>> RapidBIONJ;


#if USE_VECTORCLASS_LIBRARY
/**
 * @brief  Munging class, that adds vectorization to a BoundingMatrix<T> class
 *         (or to another class that is duck-type equivalent to
 *          BoundingMatrix<T>)
 * @tparam T - the distance type
 * @tparam V - the vector type (for a continguous blocks of T, that can 
 *             be Single Instruction Multiple Data (SMD) instructions)
 * @tparam VB - the corresponding boolean vector type
 * @tparam M     - the matrix type
 * @tparam SUPER - the immediate super class
 * @note   binary search - see indexOfFirstGreater() - is used, to figure
 *         out how many distances will need to be read from the row of S
 *         (the hope is that log(w) unpredictable distance comparisons will 
 *          be cheaper than the w predictable comparisons, that would
 *          otherwise be needed) 
 *         (in practice this hope seems to be realized. -James B)
 */
template <class T=NJFloat, class V=FloatVector, 
          class VB=FloatBoolVector,
          class M=BIONJMatrix<T>, class SUPER=BoundingMatrix<T,M> >
class VectorizedBoundingMatrix: public SUPER
{
protected:
    const intptr_t  block_size;
    std::vector<T>  vector_storage;
    std::vector<T*> per_thread_vector_storage;

public:
    typedef SUPER super;
    using super::threadCount;
    using super::getThreadNumber;
    using super::getThreadCount;

    using super::row_count;
    using super::rowTotals;
    using super::rowToCluster;

    using super::scaledClusterTotals;
    using super::clusterToRow;
    using super::clusterPartners;

    using super::entriesSorted;
    using super::entryToCluster;

    VectorizedBoundingMatrix(): super(), block_size(VB().size()),
        vector_storage  (block_size*4) {
        #ifdef _OPENMP
        vector_storage.resize(threadCount*block_size*4);
        #endif
        for (int i=0; i<threadCount; ++i) {
            T* ptr = vector_storage.data() + block_size * i * 4;
            per_thread_vector_storage.push_back( ptr ) ;
        }
    }
    /**
     * @brief  Binary search
     * @param  data  pointer to (sorted) block of (distance) T
     * @param  hi    number of distances in the block
     * @param  bound the bound to search for
     * @return size_t - the index (in data) of the first element in the
     *         block that is greater than bound. Or, if there isn't one, hi.
     * @note   this is a binary search after the style of Bottenbruch
     *         (see his 1962 paper).
     */
    size_t indexOfFirstGreater(const T* data, size_t hi, T bound) const {
        size_t lo = 0;
        while (lo<hi) {
            size_t guess = lo + (hi-lo) / 2;
            if (data[guess]<=bound) {
                lo = guess + 1;
            } else {
                hi = guess;
            }
        }
        return lo;
    }

    /**
     * @brief  Get the Postion<T> description of the best cluster merge,
     *         for which the higher-numbered cluster, y, is the one that
     *         corresponds to row (row) of the current distance matrix,
     *         by reading row (row) of both the S and I matrices.
     * @param  row    - the row number 
     * @param  maxTot - the (deep breath) maximum, of the scaled cluster
     *                  distance totals, for all of the clusters, that have
     *                  a cluster number, less than y.
     * @param  qBest  - a *copy* of the current upper bound on the global 
     *                  minimum adjusted distance, between clusters.
     * @return Position<T> - describes the best cluster merge 
     *         (if there is one!)
     */
    virtual Position<T> getRowMinimum(intptr_t row, T maxTot, T qBest) const override {
        T nless2        = (T)( row_count - 2 );
        T tMultiplier   = ( row_count <= 2 ) ? (T)0.0 : ( (T)1.0 / nless2 );
        auto tot        = scaledClusterTotals.data();
        T rowTotal      = rowTotals[row] * tMultiplier; //scaled by (1/(n-2)).
        T rowBound      = qBest + maxTot + rowTotal;
                //Upper bound for distance, in this row, that
                //could (after row totals subtracted) provide a
                //better min(Q).

        Position<T> pos(row, 0, infiniteDistance, 0);
        const T*   rowData   = entriesSorted.rows[row];
        const int* toCluster = entryToCluster.rows[row];
        const int  cluster   = rowToCluster[row];

        size_t partners   = indexOfFirstGreater
                            (rowData, clusterPartners[cluster], rowBound);
        size_t v_partners = (partners/block_size)*block_size;
        auto thread_num   = getThreadNumber();
        T*   blockRawDist = per_thread_vector_storage[thread_num];
        T*   blockCluster = blockRawDist + block_size;
        T*   blockIndex   = blockCluster + block_size;
        T*   blockHCDist  = blockIndex   + block_size;

        V    best_hc_vector ((T)(infiniteDistance));
        V    best_ix_vector ((T)(-1));
        V    raw(0);
        V    c_tot(0);
        V    ix(0);

        for (size_t i=0; i<v_partners; i+=block_size) {
            for (int j=0; j<block_size; ++j) {
                int k = toCluster[i+j];
                blockCluster[j] = tot[k];
                blockIndex[j]   = (T)k;
            }
            raw.load(rowData+i);
            c_tot.load(blockCluster);
            ix.load(blockIndex);
            V  hc   = raw - c_tot; //subtract cluster totals to get half-cooked distances?
            VB less = hc < best_hc_vector; //which are improvements?
            best_hc_vector = select(less, hc, best_hc_vector);
            best_ix_vector = select(less, ix, best_ix_vector);
        }

        if (0<v_partners) {
            best_hc_vector.store(blockHCDist);
            best_ix_vector.store(blockIndex);
            for (int j=0; j<block_size; ++j) {
                T Qrc = blockHCDist[j] - rowTotal;
                if (Qrc < pos.value) {
                    size_t cluster_num = (size_t)blockIndex[j];
                    intptr_t otherRow = clusterToRow[cluster_num];
                    if (otherRow != notMappedToRow) {
                        pos.column = (otherRow < row ) ? otherRow : row;
                        pos.row    = (otherRow < row ) ? row : otherRow;
                        pos.value  = Qrc;
                        if (Qrc < qBest ) {
                            qBest    = Qrc;
                            rowBound = qBest + maxTot + rowTotal;
                        }
                    }
                }
            }
        }
        adjustRowMinimumFromLeftovers(row, rowData, rowBound, rowTotal, maxTot, qBest, 
                                      toCluster, tot, v_partners, partners, pos);
        return pos;
    }

    /**
     * @brief Consider the last partners-v_partners items in the S and I
     *        rows, "in" getRowMinimum() (where v_partners is the number of 
     *        items that were already handled by the vectorized code), 
     *        using NON-vectorized code.
     * @note  this function only exists, seperately from getRowMinimum(), 
     *        to trick Lizard into *not* complaining about the Cyclomatic 
     *        Complexity of getRowMinimum(). This function only exists, 
     *        seperately from getRowMinimum(), to trick Lizard into *not* 
     *        complaining about the Cyclomatic Complexity of getRowMinimum().
     * @note  All the parameters are merely copies, of parameters to, 
     *        or local variables declared in, getRowMinimum().
     */
    //note: 
    inline void adjustRowMinimumFromLeftovers
        (   intptr_t row, const T* rowData, T rowBound, T rowTotal, T maxTot, T qBest,
            const int* toCluster, T* tot, size_t v_partners, size_t partners, 
            Position<T>& pos ) const {
        for (size_t i=v_partners; i<partners; ++i) {
            T Drc = rowData[i];
            if (rowBound<Drc && 0<i) {
                break;
            }
            size_t  cluster_num = toCluster[i];
                //The cluster associated with this distance
                //The c in Qrc and Drc.
            T Qrc = Drc - tot[cluster_num] - rowTotal;
            if (Qrc < pos.value) {
                intptr_t otherRow = clusterToRow[cluster_num];
                if (otherRow != notMappedToRow) {
                    pos.column = (otherRow < row ) ? otherRow : row;
                    pos.row    = (otherRow < row ) ? row : otherRow;
                    pos.value  = Qrc;
                    if (Qrc < qBest ) {
                        qBest    = Qrc;
                        rowBound = qBest + maxTot + rowTotal;
                    }
                }
            }
        }
    }
};

typedef VectorizedBoundingMatrix
        <NJFloat, FloatVector, FloatBoolVector, NJMatrix<NJFloat> >
        Vectorized_RapidNJ;
        
typedef VectorizedBoundingMatrix
        <NJFloat, FloatVector, FloatBoolVector, BIONJMatrix<NJFloat> >
        Vectorized_RapidBIONJ;

#endif //USE_VECTORCLASS_LIBRARY

}

#endif /* rapidnj_h */
