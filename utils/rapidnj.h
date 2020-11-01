//
//  rapidnj.h
//  iqtree
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
//        those that has distance matrix rows identical to earlier rows;
//        Rapid NJ "hates" them) (this is also covered in section 2.5)
//        See the BoundingBIONJMatrix class.
//
//  Created by James Barbetti on 31/10/20.
//

#ifndef rapidnj_h
#define rapidnj_h

#include "nj.h"
#include "heapsort.h"                //for mirroredHeapsort, used to sort

namespace StartTree
{

template <class T=NJFloat, class super=BIONJMatrix<T>>
class BoundingMatrix: public super
{
protected:
    using super::rows;
    using super::row_count;
    using super::rowMinima;
    using super::rowTotals;
    using super::rowToCluster;
    using super::clusters;
    using super::silent;
    using super::finishClustering;
    //
    //Note 1: mutable members are calculated repeatedly, from
    //        others, in member functions marked as const.
    //        They're declared at the class level so that they
    //        don't need to be reallocated over and over again.
    //Note 2: Mapping members to the RapidNJ papers:
    //        rows           is the D matrix
    //        entriesSorted  is the S matrix
    //        entryToCluster is the I matrix
    //Note 3: scaledMaxEarlierClusterTotal[c] is the largest row total
    //        for a cluster with a lower number, than cluster c
    //        (if c indicates a cluster for which there are still rows
    //        in the distance matrix: call this a live cluster).
    //        This is a tighter bound, when searching for
    //        the minimum Qij... and processing distances from
    //        cluster c to earlier clusters, than the largest
    //        row total for ALL the live clusters.
    //        See section 2.5 of Simonsen, Mailund, Pedersen [2011].
    //Note 4: rowOrderChosen is a vector of int rather than bool
    //        because simultaneous non-overlapping random-access writes
    //        to elements of std::vector<bool> *can* interfere with each other
    //        (because std::vector<bool> maps multiple nearby elements onto
    //         bitfields, so the writes... *do* overlap) (ouch!).
    //
    std::vector<int> clusterToRow;   //Maps clusters to their rows (-1 means not mapped)
    std::vector<T>   clusterTotals;  //"Row" totals indexed by cluster

    mutable std::vector<T>       scaledClusterTotals;   //The same, multiplied by
                                                        //(1.0 / (n-2)).
    mutable std::vector<T>       scaledMaxEarlierClusterTotal;
    mutable std::vector<int>     rowOrderChosen; //Indicates if a row's scanning
                                                 //order chosen has been chosen.
                                                 //Only used in... getRowScanningOrder().
    mutable std::vector<size_t>  rowScanOrder;   //Order in which rows are to be scanned
                                                 //Only used in... getRowMinima().
    
    SquareMatrix<T>   entriesSorted; //The S matrix: Entries in distance matrix
                                     //(each row sorted by ascending value)
    SquareMatrix<int> entryToCluster;//The I matrix: for each entry in S, which
                                     //cluster the row (that the entry came from)
                                     //was mapped to (at the time).
    double rowSortingTime;
    int    threadCount;
    
public:
    typedef T distance_type;
    BoundingMatrix() : super(), rowSortingTime(0) {
        #ifdef _OPENMP
            threadCount = omp_get_max_threads();
        #else
            threadCount = 1;
        #endif
    }
    virtual std::string getAlgorithmName() const {
        return "Rapid" + super::getAlgorithmName();
    }
    virtual bool constructTree() {
        //1. Set up vectors indexed by cluster number,
        clusterToRow.resize(row_count);
        clusterTotals.resize(row_count);
        for (size_t r=0; r<row_count; ++r) {
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
            const char* taskName = silent ? "" : "Setting up auxiliary I and S matrices";
            progress_display setupProgress(row_count, taskName, "sorting", "row");
            //2. Set up the matrix with row sorted by distance
            //   And the matrix that tracks which distance is
            //   to which cluster (the S and I matrices, in the
            //   RapidNJ papers).
            entriesSorted.setSize(row_count);
            entryToCluster.setSize(row_count);
            #ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic)
            #endif
            for (size_t r=0; r<row_count; ++r) {
                sortRow(r,r,false);
                ++setupProgress;
                //copies the "left of the diagonal" portion of
                //row r from the D matrix and sorts it
                //into ascending order.
            }
        }
        {
            size_t nextPurge = (row_count+row_count)/3;
            std::string taskName = "Constructing " + getAlgorithmName() + " tree";
            if (silent) {
                taskName = "";
            }
            double triangle = row_count * (row_count + 1.0) * 0.5;
            progress_display show_progress(triangle, taskName.c_str(), "", "");
            while (3<row_count) {
                Position<T> best;
                super::getMinimumEntry(best);
                cluster(best.column, best.row);
                if ( row_count == nextPurge ) {
                    #ifdef _OPENMP
                    #pragma omp parallel for
                    #endif
                    for (size_t r=0; r<row_count; ++r) {
                        purgeRow(r);
                    }
                    nextPurge = (row_count + row_count)/3;
                }
                show_progress+=row_count;
            }
            show_progress.done();
            finishClustering();
        }
        return true;
    }
    void sortRow(size_t r /*row index*/, size_t c /*upper bound on cluster index*/
        , bool inParallel) {
        //1. copy data from a row of the D matrix into the S matrix
        //   (and write the cluster identifiers that correspond to
        //    the values in the D row into the same-numbered
        //    row in the I matrix), for distances between the cluster
        //    in that row, and other live clusters (up to, but not including c).
        T*     sourceRow      = rows[r];
        T*     values         = entriesSorted.rows[r];
        int*   clusterIndices = entryToCluster.rows[r];
        size_t w = 0;
        for (size_t i=0; i<row_count; ++i) {
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
        if ( inParallel ) {
            double now = getRealTime();
            mirroredHeapsort(values, 0, w, clusterIndices);
            rowSortingTime += (getRealTime() - now);
        } else {
            mirroredHeapsort(values, 0, w, clusterIndices);
        }
    }
    void purgeRow(size_t r /*row index*/) const {
        //Scan a row of the I matrix, so as to remove
        //entries that refer to clusters that are no longer
        //being processed. Remove the corresponding values
        //in the same row of the S matrix.
        T*    values         = entriesSorted.rows[r];
        int*  clusterIndices = entryToCluster.rows[r];
        size_t w = 0;
        size_t i = 0;
        for (; i<row_count ; ++i ) {
            values[w]         = values[i];
            clusterIndices[w] = clusterIndices[i];
            if ( infiniteDistance <= values[i] ) {
                break;
            }
            if ( clusterToRow[clusterIndices[i]] != notMappedToRow ) {
                ++w;
            }
        }
        if (w<row_count) {
            values[w] = infiniteDistance;
        }
    }
    virtual void cluster(size_t a, size_t b) {
        size_t clusterA         = rowToCluster[a];
        size_t clusterB         = rowToCluster[b];
        size_t clusterMoved     = rowToCluster[row_count-1];
        clusterToRow[clusterA]  = notMappedToRow;
        clusterTotals[clusterA] = -infiniteDistance;
        clusterToRow[clusterB]  = notMappedToRow;
        clusterTotals[clusterB] = -infiniteDistance;
        size_t clusterC = clusters.size(); //cluster # of new cluster
        super::cluster(a,b);
        if (b<row_count) {
            clusterToRow[clusterMoved] = static_cast<int>(b);
        }
        clusterToRow.emplace_back(a);
        clusterTotals.emplace_back(rowTotals[a]);
        scaledClusterTotals.emplace_back(rowTotals[a] / (T)( row_count - 1.0 ) );
        scaledMaxEarlierClusterTotal.emplace_back(0.0);
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
        for (size_t r = 0; r<row_count; ++r) {
            size_t cluster = rowToCluster[r];
            clusterTotals[cluster] = rowTotals[r];
        }
        sortRow(a, clusterC, true);
    }
    void decideOnRowScanningOrder(T& qBest) const {
        size_t rSize = rowMinima.size();
        //
        //Rig the order in which rows are scanned based on
        //which rows (might) have the lowest row minima
        //based on what we saw last time.
        //
        //The original RapidNJ puts the second-best row from last time first.
        //And, apart from that, goes in row order.
        //But rows in the D, S, and I matrices are (all) shuffled
        //in memory, so why not do all the rows in ascending order
        //of their best Q-values from the last iteration?
        //Or, better yet... From this iteration?!
        //
        
        #define DERIVE_BOUND_FROM_FIRST_COLUMN 1
        #if (DERIVE_BOUND_FROM_FIRST_COLUMN)
        {
            //
            //Since we always have to check these entries when we process
            //the row, why not process them up front, hoping to get a
            //better bound on min(V) (and perhaps even "rule" entire rows
            //"out of consideration", using that bound)? -James B).
            //
            std::vector<T> qLocalBestVector;
            qLocalBestVector.resize( threadCount, qBest);
            T* qLocalBest =  qLocalBestVector.data();

            #ifdef _OPEN_MP
            #pragma omp parallel for
            #endif
            for (size_t b=0; b<threadCount; ++b) {
                T      qBestForThread = qBest;
                size_t rStart         = b*rSize / threadCount;
                size_t rStop          = (b+1)*rSize / threadCount;
                for (size_t r=rStart; r < rStop
                     && rowMinima[r].value < infiniteDistance; ++r) {
                    size_t rowA       = rowMinima[r].row;
                    size_t rowB       = rowMinima[r].column;
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
            for (size_t b=0; b<threadCount; ++b) {
                if ( qLocalBest[b] < qBest ) {
                    qBest = qLocalBest[b];
                }
            }
        }
        #endif
        
        int threshold = threadCount << 7; /* multiplied by 128*/
        //Note, rowMinima might have size 0 (the first time this member
        //function is called during processing of a distance matrix)
        //Or it might have a size of n+1 (later times), but it won't be n.
        for ( size_t len = rSize; 1<len; len=(len+1)/2 ) {
            size_t halfLen = len/2; //rounded down
            size_t gap     = len-halfLen;
            #ifdef _OPENMP
            #pragma omp parallel for if(threshold<halfLen)
            #endif
            for ( size_t i=0; i<halfLen; ++i) {
                size_t j = i + gap;
                if ( rowMinima[j] < rowMinima[i] ) {
                    std::swap(rowMinima[i], rowMinima[j]);
                }
            }
        }
        #ifdef _OPENMP
        #pragma omp parallel for if(threshold<row_count)
        #endif
        for (size_t i = 0; i < row_count; ++i) {
            rowOrderChosen[i]=0; //Not chosen yet
        }
        
        size_t w = 0;
        for (size_t r=0; r < rSize
             && rowMinima[r].row < row_count
             && rowMinima[r].value < infiniteDistance; ++r) {
            size_t rowA     = rowMinima[r].row;
            size_t rowB     = rowMinima[r].column;
            size_t clusterA = (rowA<row_count) ? rowToCluster[rowA] : 0;
            size_t clusterB = (rowB<row_count) ? rowToCluster[rowB] : 0;
            size_t row      = (clusterA<clusterB) ? rowA : rowB;
            if (row < row_count) {
                rowScanOrder[w] = row;
                w += rowOrderChosen[row] ? 0 : 1;
                rowOrderChosen[row] = 1; //Chosen
            }
        }
        
        //The weird-looking middle term in the for loop is as
        //intended: when w reaches n all of the rows (0..n-1)
        //must be in rowScanOrder, so there's no need to continue
        //until row==n.
        for (size_t row=0; w < row_count ; ++row) {
            rowScanOrder[w] = row;
            w += ( rowOrderChosen[row] ? 0 : 1 );
        }
    }
    virtual void getRowMinima() const {
        //
        //Note: Rather than multiplying distances by (n-2)
        //      repeatedly, it is cheaper to work with cluster
        //      totals multiplied by (1.0/(T)(n-2)).
        //      Better n multiplications than 0.5*n*(n-1).
        //Note 2: Note that these are indexed by cluster number,
        //      and *not* by row number.
        //
        size_t c           = clusters.size();
        T      nless2      = ( row_count - 2 );
        T      tMultiplier = ( row_count <= 2 ) ? 0 : (1 / nless2);
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
        #pragma omp parallel for
        #endif
        for (size_t r=0; r<row_count ; ++r) {
            size_t row             = rowScanOrder[r];
            size_t cluster         = rowToCluster[row];
            T      maxEarlierTotal = scaledMaxEarlierClusterTotal[cluster];
            //Note: Older versions of RapidNJ used maxTot rather than
            //      maxEarlierTotal here...
            rowMinima[r]           = getRowMinimum(row, maxEarlierTotal, qBest);
        }
    }
    Position<T> getRowMinimum(size_t row, T maxTot, T qBest) const {
        T nless2      = ( row_count - 2 );
        T tMultiplier = ( row_count <= 2 ) ? 0 : ( 1.0 / nless2 );
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
            size_t  cluster = toCluster[i];
                //The cluster associated with this distance
                //The c in Qrc and Drc.
            T Qrc = Drc - tot[cluster] - rowTotal;
            if (Qrc < pos.value) {
                int otherRow = clusterToRow[cluster];
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

}

#endif /* rapidnj_h */
