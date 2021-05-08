//
//  upgma.h
//  UPGMA_Matrix template class.
//  Implementation of the UPGMA algorithm of Rober R. Sokal, and
//  Charles D. Michener (1958), "Evaluating Systematic Relationships"
//  (in the University of Kansas Science Bulletin).
//  UPGMA is (slightly) simpler than NJ,BIONJ, and UNJ.
//  Created by James Barbetti on 31/10/20.
//
//  UPGMA_Matrix extends SquareMatrix like so:
//  1. It maintains a mapping between row numbers (the rows
//     for clusters still being considered) and cluster numbers,
//     in its rowToCluster member.  That's initialized in setSize().
//  2. It keeps track of the clusters that have been created
//     thus far, in its cluster member. Each single taxon is
//     considered a cluster (and to begin with, row i corresponds
//     to cluster i, for each of the rank rows in the V matrix).
//     The first rank clusters are added to the vector in setSize().
//  3. It keeps track of the best candidate "join" found looking
//     at each row in the V matrix.  In a rowMinima vector.
//  4. It defines a number of public member functions that are
//     overridden in its subclasses:
//     (a) loadMatrixFromFile
//     (b) loadMatrix
//     (c) constructTree
//     (d) writeTreeFile
//  5. It defines a number of protected member functions that are
//     overridden in its subclasses:
//     (a) getMinimumEntry() - identify the the row an column that
//                             correspond to the next two clusters
//                             to be joined.
//     (b) getRowMinima()    - find, for each row in the matrix,
//                             which column (corresponding to another
//                             clusters) corresponds to the cluster
//                             that is most "cheaply" joined with the
//                             cluster corresponding to the row.
//                             Write the answers in rowMinima.
//     (c) getImbalance      - determine, for two clusters that might
//                             be joined "how out of balance" the sizes
//                             of the clusters are.  This is used for
//                             tie-breaking, and to try to avoid
//                             degenerate trees when many taxa are
//                             identical.
//     (d) cluster           - given two row/column numbers a and b
//                             (where a is less), for rows that
//                             correspond to clusters to be joined,
//                             record that they have been joined,
//                             calculate a new row for the joined
//                             cluster, write that over the top of row a,
//                             and remove row b via removeRowAndColumn
//                             (which writes the content of the last row
//                              in the matrix over the top of b, and then
//                              removes the last row from the matrix).
//     (e) finishClustering  - join up the last three clusters
//
//  Notes:
//  A. rowMinima could be defined in constructTree() and passed
//     down to getMinimumEntry() and rowMinima(), but declaring it as a
//     member function of the class makes it easier to look at it in a
//     debugger (and saves on some passing around of pointers between
//     member functions).
//  B. The convention is that column numbers are less than row numbers.
//     (it is assumed that the matrix is symmetric around its diagonal).
//  C. Rows are *swapped* (and the last row/column removed from the
//     matrix, because this approach avoids keeping track of which rows
//     or columns are "out of use") (all are in use, all the time!),
//     and reduces the number of memory accesses by a factor of about 3
//     because, asymptotically, the sum of the squares of the numbers
//     up to N is N*(N+1)*(2*N+1)/6.  But the real benefit is avoiding the
//     pipeline stalls that would result from mispredicted branches for
//     *if* statements that would otherwise be required, for the checks
//     whether a given row is in use.  Row processing is also more easily
//     vectorized but, in terms of performance, that matters less
//     (Vectorization is... x2 or so, avoiding the ifs is... x5 or more).
//

#ifndef upgma_h
#define upgma_h

#include "distancematrix.h"          //for Matrix template class
#include "vectortypes.h"             //for StrVector and IntVector
#include "clustertree.h"             //for ClusterTree template class
#include <vector>                    //for std::vector
#include <string>                    //sequence names stored as std::string
#include <functional>                //for std::hash
#include <algorithm>                 //for std::sort
#include "progress.h"                //for progress_display
#include "my_assert.h"               //for ASSERT macro

#if (!USE_PROGRESS_DISPLAY)
typedef double progress_display;
#endif

typedef float    NJFloat;
const   NJFloat  infiniteDistance = (NJFloat)(1e+36);
const   intptr_t notMappedToRow = -1;

#ifdef   USE_VECTORCLASS_LIBRARY
#include <vectorclass/vectorclass.h> //for Vec4d and Vec4db vector classes
typedef  Vec8f    FloatVector;
typedef  Vec8fb   FloatBoolVector;
#endif

namespace StartTree
{
template <class T=NJFloat> struct Position
{
    //A position (row, column) in an UPGMA or NJ matrix
    //Note that column should be strictly less than row.
    //(Because that is the convention in RapidNJ).
public:
    intptr_t row;
    intptr_t column;
    T        value;
    size_t   imbalance;
    Position() : row(0), column(0), value(0), imbalance(0) {}
    Position(size_t r, size_t c, T v, size_t imbalance)
        : row(r), column(c), value(v) {}
    Position& operator = (const Position &rhs) {
        row       = rhs.row;
        column    = rhs.column;
        value     = rhs.value;
        imbalance = rhs.imbalance;
        return *this;
    }
    bool operator< ( const Position& rhs ) const {
        return value < rhs.value
        || (value == rhs.value && imbalance < rhs.imbalance);
    }
    bool operator<= ( const Position& rhs ) const {
        return value < rhs.value
        || (value == rhs.value && imbalance <= rhs.imbalance);
    }
};

template <class T> class Positions : public std::vector<Position<T>>
{
};


template <class T=NJFloat> class UPGMA_Matrix: public SquareMatrix<T> {
    //UPGMA_Matrix is a D matrix (a matrix of distances).
public:
    typedef SquareMatrix<T> super;
    using super::rows;
    using super::setSize;
    using super::loadDistancesFromFlatArray;
    using super::calculateRowTotals;
    using super::removeRowAndColumn;
    using super::row_count;
    using super::column_count;
protected:
    std::vector<size_t>  rowToCluster; //*not* initialized by setSize
    ClusterTree<T>       clusters;     //*not* touched by setSize
    bool                 isOutputToBeZipped;
    mutable Positions<T> rowMinima;    //*not* touched by setSize
    bool                 silent;
public:
    UPGMA_Matrix():super(), isOutputToBeZipped(false), silent(false) {
    }
    virtual std::string getAlgorithmName() const {
        return "UPGMA";
    }
    virtual void setSize(intptr_t rank) {
        super::setSize(rank);
        rowToCluster.clear();
        for (intptr_t r=0; r<row_count; ++r) {
            rowToCluster.emplace_back(r);
        }
    }
    virtual void addCluster(const std::string &name) {
        clusters.addCluster(name);
    }
    virtual bool loadMatrixFromFile(const std::string &distanceMatrixFilePath) {
        bool rc = loadDistanceMatrixInto(distanceMatrixFilePath, true, *this);
        calculateRowTotals();
        return rc;
    }
    virtual bool loadMatrix(const StrVector& names,
                            const double* matrix) {
        //Assumptions: 2 < names.size(), all names distinct
        //  matrix is symmetric, with matrix[row*names.size()+col]
        //  containing the distance between taxon row and taxon col.
        setSize(static_cast<intptr_t>(names.size()));
        clusters.clear();
        for (auto it = names.begin(); it != names.end(); ++it) {
            clusters.addCluster(*it);
        }
        loadDistancesFromFlatArray(matrix);
        calculateRowTotals();
        return true;
    }
    virtual bool constructTree() {
        clusterDuplicates();

        Position<T> best;
        #if USE_PROGRESS_DISPLAY
        std::string taskName = "Constructing " + getAlgorithmName() + " tree";
        if (silent) {
            taskName="";
        }
        double triangle = row_count * (row_count + 1.0) * 0.5;
        progress_display show_progress(triangle, taskName.c_str(), "", "");
        #else
        double show_progress = 0;
        #endif
        while ( 3 < row_count ) {
            getMinimumEntry(best);
            cluster(best.column, best.row);
            show_progress += row_count;
        }
        finishClustering();
        #if USE_PROGRESS_DISPLAY
        show_progress.done();
        #endif
        return true;
    }
    virtual void setZippedOutput(bool zipIt) {
        isOutputToBeZipped = zipIt;
    }
    virtual void beSilent() {
        silent = true;
    }
    bool writeTreeFile(int precision, const std::string &treeFilePath) const {
        return clusters.writeTreeFile(isOutputToBeZipped, precision, treeFilePath);
    }
protected:
    void getMinimumEntry(Position<T> &best) {
        getRowMinima();
        best.value = infiniteDistance;
        for (intptr_t r=0; r<row_count; ++r) {
            Position<T> & here = rowMinima[r];
            if (here.value < best.value && here.row != here.column) {
                best = here;
            }
        }
    }
    virtual void getRowMinima() const
    {
        rowMinima.resize(row_count);
        rowMinima[0].value = infiniteDistance;
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic)
        #endif
        for (intptr_t row=1; row<row_count; ++row) {
            T      bestVrc    = (T)infiniteDistance;
            size_t bestColumn = 0;
            const  T* rowData = rows[row];
            for (intptr_t col=0; col<row; ++col) {
                T    v      = rowData[col];
                bool better = ( v < bestVrc );
                if (better) {
                    bestColumn = col;
                    bestVrc = v;
                }
            }
            rowMinima[row] = Position<T>(row, bestColumn, bestVrc, getImbalance(row, bestColumn));
        }
    }
    virtual void finishClustering() {
        //But:  The formula is probably wrong. Felsenstein [2004] chapter 11 only
        //      covers UPGMA for rooted trees, and I don't know what
        //      the right formula is.
        //
        ASSERT( row_count == 3);
        T weights[3];
        T denominator = (T)0.0;
        for (size_t i=0; i<3; ++i) {
            weights[i]   = (T)clusters[rowToCluster[i]].countOfExteriorNodes;
            denominator += weights[i];
        }
        for (size_t i=0; i<3; ++i) {
            weights[i] /= ((T)2.0 * denominator);
        }
        clusters.addCluster
            ( rowToCluster[0], weights[1]*rows[0][1] + weights[2]*rows[0][2]
            , rowToCluster[1], weights[0]*rows[0][1] + weights[2]*rows[1][2]
            , rowToCluster[2], weights[0]*rows[0][2] + weights[1]*rows[1][2]);
        row_count = 0;
    }
    virtual void cluster(intptr_t a, intptr_t b) {
        T      aLength = rows[b][a] * (T)0.5;
        T      bLength = aLength;
        size_t aCount  = clusters[rowToCluster[a]].countOfExteriorNodes;
        size_t bCount  = clusters[rowToCluster[b]].countOfExteriorNodes;
        size_t tCount  = aCount + bCount;
        T      lambda  = (T)aCount / (T)tCount;
        T      mu      = (T)1.0 - lambda;
        auto rowA = rows[a];
        auto rowB = rows[b];
        for (intptr_t i=0; i<row_count; ++i) {
            if (i!=a && i!=b) {
                T Dai      = rowA[i];
                T Dbi      = rowB[i];
                T Dci      = lambda * Dai + mu * Dbi;
                rowA[i]    = Dci;
                rows[i][a] = Dci;
            }
        }
        clusters.addCluster ( rowToCluster[a], aLength,
                              rowToCluster[b], bLength);
        rowToCluster[a] = clusters.size()-1;
        rowToCluster[b] = rowToCluster[row_count-1];
        removeRowAndColumn(b);
    }
    struct HashRow {
    public:
        size_t hash;
        size_t row_num;
        int compare(const HashRow& rhs, UPGMA_Matrix<T>* me) const {
            if (hash!=rhs.hash) {
                return (hash < rhs.hash) ? -1 : 1;
            }
            T* rowA = me->rows[row_num];
            T* rowB = me->rows[rhs.row_num];
            for (intptr_t i=0; i<me->row_count; ++i) {
                if (rowA[i]!=rowB[i]) {
                    return (rowA[i] < rowB[i]) ? -1 : 1;
                }
            }
            return 0;
        }
    };
    void clusterDuplicates() {
        #if (USE_PROGRESS_DISPLAY)
        std::string taskName = "Identifying identical (and nearly identical) taxa";
        if (silent) {
            taskName="";
        }
        progress_display show_progress(row_count*2, taskName.c_str(), "", "");
        #else
        progress_display show_progress = 0;
        #endif

        std::vector<HashRow> hashed_rows;
        calculateRowHashes(hashed_rows, show_progress);
        std::vector< std::vector< intptr_t > > vvc;
        identifyDuplicateClusters(hashed_rows, vvc, show_progress);
        
        #if (USE_PROGRESS_DISPLAY)
        size_t dupes_clustered = joinUpDuplicateClusters(vvc, show_progress);
        show_progress.done();
        if (0<dupes_clustered && !silent) {
            std::cout << "Clustered " << dupes_clustered
                      << " identical (or near-identical) taxa." << std::endl;
        }
        #else
            joinUpDuplicateClusters(vvc, show_progress);
        #endif
    }
    void calculateRowHashes(std::vector<HashRow>& hashed_rows,
                            progress_display& show_progress) {
        //1. Calculate row hashes
        hashed_rows.resize(row_count);
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (intptr_t i=0; i<row_count; ++i) {
            T* row_i = this->rows[i];
            size_t row_hash = 0;
            for (intptr_t col=0; col<row_count;++col ) {
                row_hash ^= std::hash<double>()(row_i[col])
                    + 0x9e3779b9 + (row_hash<<6) + (row_hash>>2);
            }
            hashed_rows[i].hash     = row_hash;
            hashed_rows[i].row_num  = i;
            if ((i%1000)==999) {
                show_progress += (double)(1000.0);
            }
        }
        //2. Sort rows by hash (and tiebreak on row content)
        auto me = this;
        std::sort(hashed_rows.begin(), hashed_rows.end(),
              [me](const HashRow& a, const HashRow& b) {
                return a.compare(b, me) < 0;
            });
        show_progress += (double)(row_count%1000);
    }
    void identifyDuplicateClusters(const std::vector<HashRow>& hashed_rows,
                                   std::vector< std::vector< intptr_t > >& vvc,
                                   progress_display& show_progress) {
        //3. Now, any duplicate rows are next to each other in
        //   hashed_rows.  Construct a vector of vectors of
        //   cluster numbers (vvc).
        std::vector< intptr_t> vc; //vector of cluster #s
        for (intptr_t i=1; i<row_count; ++i) {
            if (hashed_rows[i].compare(hashed_rows[i-1], this)!=0) {
                //Not a duplicate of the previous row.
                if (!vc.empty()) {
                    //Add vector of the clusters in previous group of
                    //duplicates, to vvc.  And start a new one (empty).
                    vvc.emplace_back(vc);
                    vc.clear();
                }
            } else { //duplicate of previous row
                if (vc.empty()) {
                    //This row and the previous one are to be the first
                    //two rows whose clusters belong in the current
                    //vector of duplicate clusters, which is empty.
                    //Add the cluster that corresponds to the previous row.
                    vc.push_back(rowToCluster[hashed_rows[i-1].row_num]);
                }
                vc.push_back(rowToCluster[hashed_rows[i].row_num]);
            }
        }
        if (!vc.empty()) {
            vvc.emplace_back(vc);
        }
    }
    size_t joinUpDuplicateClusters(std::vector< std::vector< intptr_t > >& vvc,
                                 progress_display& show_progress) {
        if (vvc.empty()) {
            show_progress += (double)row_count;
            return 0; //Nothing to do!
        }
        //4. Set up an array to map cluster numbers to row numbers
        //   (step 5 will maintain this, as it goes)
        std::vector< intptr_t > cluster_to_row;
        cluster_to_row.resize(clusters.size(), -1 );
        for (int i=0; i<row_count; ++i) {
            cluster_to_row[rowToCluster[i]] = i;
        }

        double dupes = 0.0;
        for (const std::vector<intptr_t>& vc : vvc) {
            dupes += vc.size();
        }
        double work_per_dupe = (double)row_count / dupes;
        
        //5. Join up any duplicate clusters
        double work_done = 0.0;
        size_t dupes_removed = 0;
        for (std::vector<intptr_t> vc : vvc) {
            double work_here = static_cast<double>(vc.size()) * work_per_dupe;
            dupes_removed += vc.size();
            --dupes_removed;
            while (vc.size()>1) {
                intptr_t first_half = vc.size() / 2;
                intptr_t second_half = vc.size()-first_half;
                for (intptr_t i=0; i<first_half; ++i) {
                    intptr_t cluster_a = vc[i];
                    intptr_t row_a     = cluster_to_row[cluster_a];
                    intptr_t cluster_b = vc[i+second_half];
                    intptr_t row_b     = cluster_to_row[cluster_b];
                    intptr_t cluster_c = clusters.size();
                    if (row_b<row_a) {
                        std::swap(row_a, row_b);
                    }
                    cluster(row_a, row_b);
                    vc[i] = cluster_c;
                    cluster_to_row.push_back(row_a);
                    if (row_b < row_count) {
                        cluster_to_row[rowToCluster[row_b]] = row_b;
                    }
                }
                vc.resize(second_half);
                //Not first_half (rounded down), second_half (rounded up), 
                //because, if there was an odd cluster, this keeps it
            }
            work_done += work_here;
            if (work_done > 1000.0) {
                show_progress += 1000.0;
                work_done -= 1000.0;
            }
        }
        show_progress += work_done;
        return dupes_removed;
    }
    size_t getImbalance(size_t rowA, size_t rowB) const {
        size_t clusterA = rowToCluster[rowA];
        size_t clusterB = rowToCluster[rowB];
        size_t sizeA    = clusters[clusterA].countOfExteriorNodes;
        size_t sizeB    = clusters[clusterB].countOfExteriorNodes;
        return (sizeA<sizeB) ? (sizeB-sizeA) : (sizeA-sizeB);
    }
};

#ifdef USE_VECTORCLASS_LIBRARY
template <class T=NJFloat, class V=FloatVector, class VB=FloatBoolVector>
class VectorizedUPGMA_Matrix: public UPGMA_Matrix<T>
{
protected:
    typedef UPGMA_Matrix<T> super;
    using super::rowMinima;
    using super::rows;
    using super::row_count;
    using super::calculateRowTotals;
    using super::getImbalance;
    const intptr_t blockSize;
    mutable std::vector<T> scratchColumnNumbers;
public:
    VectorizedUPGMA_Matrix() : super(), blockSize(VB().size()) {
    }
    virtual std::string getAlgorithmName() const {
        return "Vectorized-" + super::getAlgorithmName();
    }
    virtual void calculateRowTotals() const {
        size_t fluff = MATRIX_ALIGNMENT / sizeof(T);
        scratchColumnNumbers.resize(row_count + fluff, 0.0);
    }
    virtual void getRowMinima() const
    {
        T* nums = matrixAlign ( scratchColumnNumbers.data() );
        rowMinima.resize(row_count);
        rowMinima[0].value = infiniteDistance;
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic)
        #endif
        for (intptr_t row=1; row<row_count; ++row) {
            Position<T> pos(row, 0, infiniteDistance, 0);
            const T* rowData = rows[row];
            intptr_t col;
            V minVector  = infiniteDistance;
            V ixVector   = -1 ;

            for (col=0; col+blockSize<row; col+=blockSize) {
                V  rowVector; rowVector.load_a(rowData+col);
                VB less      = rowVector < minVector;
                V  numVector; numVector.load_a(nums+col);
                ixVector     = select(less, numVector, ixVector);
                minVector    = select(less, rowVector, minVector);
            }
            //Extract minimum and column number
            for (int c=0; c<blockSize; ++c) {
                if (minVector[c] < pos.value) {
                    pos.value  = minVector[c];
                    pos.column = (size_t)ixVector[c];
                }
            }
            for (; col<row; ++col) {
                T dist = rowData[col];
                if (dist < pos.value) {
                    pos.column = col;
                    pos.value  = dist;
                }
            }
            pos.imbalance  = getImbalance(pos.row, pos.column);
            rowMinima[row] = pos;
        }
    }
};
#endif //USE_VECTORCLASS_LIBRARY

}

#endif /* upgma_h */
