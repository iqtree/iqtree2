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
//                             (almost) identical.
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

#include <utils/vectortypes.h>       //for StrVector and IntVector
#include <utils/progress.h>          //for progress_display
#include <utils/my_assert.h>         //for ASSERT macro
#include <vector>                    //for std::vector
#include <string>                    //sequence names stored as std::string
#include <functional>                //for std::hash
#include <algorithm>                 //for std::sort

#include "distancematrix.h"          //for Matrix template class
#include "hashrow.h"                 //for HashRow template class
#include "clustertree.h"             //for ClusterTree template class
#include "utils/parallel_mergesort.h"

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

/**
 * @brief  A position (row, column) in an UPGMA or NJ matrix
 *         Note that column should be strictly less than row,
 *         as that is the convention in RapidN. 
 *         A "nullatory" position has row and column 
 *         both set to 0.
 * @tparam T the floating-point type used to represent 
 *           distances
 */
template <class T=NJFloat> struct Position
{
    //
public:
    intptr_t row;       //should be < column
    intptr_t column;    //should be > row (unless both zero)
    T        value;     //
    size_t   imbalance; //indicates how out-of-balance a join
                        //between the clusters (indicated by row
                        //and column) would be. Typically the 
                        //absolute value of the differences
                        //between the counts of taxa.
                        
    Position() : row(0), column(0), value(0), imbalance(0) {}
    Position(size_t r, size_t c, T v, size_t imbalance_val)
        : row(r), column(c), value(v), imbalance(imbalance_val) {}
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

/**
 * @brief  A matrix, of distances, between clusters used by the UPGMA
 *         algorithm.
 * @tparam T the floating point type that represents a distance.
 */
template <class T=NJFloat> class UPGMA_Matrix: public SquareMatrix<T> {
    //UPGMA_Matrix is a D matrix (a matrix of distances).
public:
    typedef SquareMatrix<T> super;
    using super::rows;
    using super::setSize;
    using super::loadDistancesFromFlatArray;
    using super::calculateRowTotals;
    using super::removeRowAndColumn;  //reduces matrix rank by 1
    using super::row_count;
    using super::column_count;        
protected:
    std::vector<size_t>  rowToCluster;
    ClusterTree<T>       clusters;     //*not* touched by setSize
    mutable Positions<T> rowMinima;    //*not* touched by setSize
    bool                 silent;
    bool                 isOutputToBeZipped;
    bool                 isOutputToBeAppended;
    bool                 isRooted;
    bool                 subtreeOnly;
public:
    UPGMA_Matrix(): super(), silent(false)
                  , isOutputToBeZipped(false), isRooted(false)
                  , subtreeOnly(false) {
    }
    /**
     * @brief  return the name of the algorithm this isntance
     *         implements (subclasses override it).
     * @return std::string - the name of the algorithm
     */
    virtual std::string getAlgorithmName() const {
        return "UPGMA";
    }
    /**
     * @brief Set the number of rows and columns
     * @param rank 
     * @note  initializes rowToCluster 
     * @note  does not write anything into clusters
     * @note  does not initialize rowMinima
     */
    virtual void setSize(intptr_t rank) override {
        super::setSize(rank);
        rowToCluster.clear();
        for (intptr_t r=0; r<row_count; ++r) {
            rowToCluster.emplace_back(r);
        }
    }
    /**
     * @brief Add a cluster (with the specified name) to
     *        the clusters.  The name is what will be 
     *        displayed in a newick format tree.       
     * @param name 
     * @note  This is normally called only for leaf taxa.
     */
    virtual void addCluster(const std::string &name) override {
        clusters.addCluster(name);
    }

    /**
     * @brief  Load sequence names, and distances
     *         from the phylip-format distance matrix
     *         file indicated by the supplied file path.
     * @param  distanceMatrixFilePath 
     * @return true  - if the distance matrix file is loaded
     * @return false - if an error prevents the loading of
     *         the sequence names and distances, from the file. 
     */
    virtual bool loadMatrixFromFile(const std::string &distanceMatrixFilePath) {
        bool rc = loadDistanceMatrixInto(distanceMatrixFilePath, true, *this);
        calculateRowTotals();
        return rc;
    }

    /**
     * @brief  Load sequence names, and a distance matrix
     * @param  names  an n-item vector of leaf taxon names
     * @param  matrix a flat n*n square matrix, stored in row major order.
     * @return true   always.
     * @note   Assumes that 2<names.size(), all names are distinct
     * @note   Assumes matrix is symmetric, with 
     *         both matrix[row*names.size()+col]
     *         *and* matrix[col*names.size()+row]
     *         containing the distance between taxon row and taxon col.
     * @note   Assumes that the entries down the diagonal are all zero
     *         (If they're not, they'll be blindly copied anyway).
     */
    virtual bool loadMatrix(const StrVector& names,
                            const double* matrix) {
        setSize(static_cast<intptr_t>(names.size()));
        clusters.clear();
        for (auto it = names.begin(); it != names.end(); ++it) {
            clusters.addCluster(*it);
        }
        loadDistancesFromFlatArray(matrix);
        calculateRowTotals();
        return true;
    }

    /**
     * @brief  Sets a flag that indicates whether we are building a
     *         subtree (where the "last" node of the tree has degree 2)
     *         Usually these routines are used to construct unrooted
     *         trees (where the "last" node of the tree has degree 3).
     * @param  rootIt true, if the "last" node is to be of degree 2.
     *                false, if the "last" node is to be of degree 3.
     * @return true - indicating that this implementation allows rooted true
     */
    virtual bool setIsRooted(bool rootIt) {
        isRooted = rootIt;
        return true;
    }

    /**
     * @brief  Controls whether the children of the "last" node will
     *         be output as ( node1 node2 [node3] );
     *         (which is the default, false), or whether they will 
     *         be as node1 node [node3] ... without the enclosing (, )
     *         and the terminating semi-colon.
     * @param  wantSubtree true if the (, ), and ; are to be omitted.
     * @return true - indicating this implementation can be asked to
     *         generate a subtree (stuff hanging off a node, inside a
     *         larger tree). wantSubtree=true is used for doing neighbour 
     *         joining of subtrees, to generate a tree that is
     *         consistent with a constraint tree.
     */
    virtual bool setSubtreeOnly(bool wantSubtree) {
        subtreeOnly = wantSubtree;
        return true;
    }

    /**
     * @brief Called before tree construction begins (can be overridden
     *        in subclasses that must carry out additional preparation
     *        before it can construct a tree).
     * @note  BIONJ implementations use this to ensure that their
     *        variance matrix is properly initialized.
     */
    virtual void prepareToConstructTree() {
    }

    /**
     * @brief  Construct a tree 
     * @return true - always
     * @note   clusterDuplicates() handles taxa that are identical
     *         (effectively, taxa which correspond to rows, in the
     *         distance matrix, that are identical, are considered
     *         to be identical).  It could... do nothing. That's okay.
     * @note   the heavy lifting is done by the getMinimumEntry() and
     *         cluster() member functions which subclasses may override.
     */
    virtual bool constructTree() {
        prepareToConstructTree();
        clusterDuplicates();

        Position<T> best;
        #if USE_PROGRESS_DISPLAY
        std::string taskName = "Constructing " + getAlgorithmName() + " tree";
        if (silent) {
            taskName="";
        }
        double triangle = row_count * (row_count + 1.0) * 0.5;
        progress_display show_progress(triangle, taskName.c_str(), "", "");
        #endif

        int degree_of_root = isRooted ? 2 : 3;
        while ( degree_of_root < row_count ) {
            getMinimumEntry(best);
            cluster(best.column, best.row);
            #if USE_PROGRESS_DISPLAY
            show_progress += row_count;
            #endif
        }
        finishClustering();
        #if USE_PROGRESS_DISPLAY
        show_progress.done();
        #endif
        return true;
    }

    /**
     * @brief  Called to indicate that output is to be compressed.
     *         This won't have any real effect unless compression
     *         libraries are used and linked (unless USE_GZSTREAM
     *         is defined and is non-zero)
     * @param  zipIt - true if compression is rquested
     * @return true - always returns true
     */
    virtual bool setZippedOutput(bool zipIt) {
        isOutputToBeZipped = zipIt;
        return true;
    }

    /**
     * @brief  Indicate whether the output file is, if it already exists,
     *         to be appended (true) or overwritten (false)
     * @param  appendIt true for append, false for overwrite
     * @return true - always returns true (this implementation supports 
     *         appending output files).
     */
    virtual bool setAppendFile(bool appendIt) {
        isOutputToBeAppended = appendIt;
        return true;
    }

    /**
     * @brief Called to suppress logging and displaying of progress.
     */
    virtual void beSilent() {
        silent = true;
    }

    /**
     * @brief  Writes a tree (or subtree) to an open file
     * @param  stream the file to write to
     * @return true if it could be written, false if an error occurs
     * @note   if something goes wrong, an error message ought have
     *         been logged to std::cerr before this returns control
     *         to the caller.
     * @note   it is assumed that the precision (the number of digits
     *         to display after a decimal point, for a number), will
     *         already have been set before this member function is
     *         called.
     */
    bool writeTreeToOpenFile(std::iostream& stream) const {
        return clusters.writeTreeToOpenFile(subtreeOnly, stream);
    }

    /**
     * @brief  Overwrites, or appends, the file with specified path.
     * @param  precision    how many digits of precision
     * @param  treeFilePath the name of the file path
     * @return true  if the file is overwritten or appended successfully.
     * @return false if an error occurs
     * @note   member variables that control behaviour:
     *         isOutputToBeZipped   - if true, file compressed in gzip format
     *         isOutputToBeAppended - if true, file is to be appended
     *         subtreeOnly - if true, the leading (, trailing ) and ; for
     *         the "last" or "top" node are to be omitted.
     */
    bool writeTreeFile(int precision, const std::string &treeFilePath) const {
        return clusters.writeTreeFile(isOutputToBeZipped, precision,
                                      treeFilePath, isOutputToBeAppended,
                                      subtreeOnly);
    }
    
    /**
     * @brief  Calculate the root-mean-square of the matrix of 
     *         differences between the tree-distances (between taxa) 
     *         implied by the edge lengths in (clusters), and the
     *         input distances (in a distance matrix, which should be
     *         a copy of the input that was provided to the phylogenetic
     *         inference algorithm).
     * @param  matrix a flat matrix (rank*rank) of doubles in row-major 
     *                order.
     * @param  rank   the rank of the flat matrix.
     * @param  rms    the root mean square will be written here.
     * @return true if the calculation succeeded
     * @note   it is assumed that rank is equal to the number of leaf
     *         taxa in clusters, and that the rows and columns in
     *         the distance matrix (pointed to by matrix) correspond
     *         one-to-one and in order to the first (rank) clusters.
     * @note   it is assumed that the distance matrix is symmetric.
     *         it is assumed that the diagnal entries of the distance 
     *         matrix are all zeroes.
     */
    bool calculateRMSOfTMinusD(const double* matrix, intptr_t rank, double& rms) {
        return clusters.calculateRMSOfTMinusD(matrix, rank, rms);
    }

protected:
    /**
     * @brief find the pair of clusters, indicated by a Position<T>
     *        that should be joined next.
     * @param best - reference to a Position<T> into which the answer
     *        is to be returned.
     * @note  assumes that 0 <= best.col < best.row < row_count
     */
    void getMinimumEntry(Position<T> &best) {
        getRowMinima();
        best.value = infiniteDistance;
        for (intptr_t r=0; r<row_count; ++r) {
            const Position<T> & here = rowMinima[r];
            if (here.value < best.value && here.row != here.column) {
                best = here;
            }
        }
    }

    /**
     * @brief  Determine, for each row, r, the cluster (in that row)
     *         that would be best used, and write information about
     *         that cluster into rowMinima[r].
     * @note   rowMinima[r].value  is the adjusted distance for the cluster
     *                             (lower is better)
     *         rowMinima[r].column is the lower-numbered row 
     *                             (should always be >=0 and less than r)
     *         rowMinima[r].row    is the higher-numbered row 
     *                             (should always be equal to r)
     * @note   The search is parallelized over the rows of the matrix
     *         (if _OPENMP is defined).
     */
    virtual void getRowMinima() const {
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

    /**
     * @brief Link together the last two or three nodes to a "last" or "top"
     *        node of the tree. 
     * @note  The formula might be wrong in the two-node case.
     *        Felsenstein [2004] chapter 11 only covers UPGMA 
     *        for rooted trees, and I don't know what
     *        the right formula is for unrooted trees.
     */
    virtual void finishClustering() {
        ASSERT( row_count == 2 || row_count == 3);
        std::vector<T> weights(row_count);
        T denominator = (T)0.0;
        for (intptr_t i=0; i<row_count; ++i) {
            weights[i]   = (T)clusters[rowToCluster[i]].countOfExteriorNodes;
            denominator += weights[i];
        }
        for (intptr_t i=0; i<row_count; ++i) {
            weights[i] /= ((T)2.0 * denominator);
        }
        if (row_count==3) {
            //Unrooted tree. Last cluster has degree 3.
            clusters.addCluster
                ( rowToCluster[0], weights[1]*rows[0][1] + weights[2]*rows[0][2]
                , rowToCluster[1], weights[0]*rows[0][1] + weights[2]*rows[1][2]
                , rowToCluster[2], weights[0]*rows[0][2] + weights[1]*rows[1][2]);
        } else {
            //Rooted tree. Last cluster has degree 2.
            clusters.addCluster
                ( rowToCluster[0], weights[1]*rows[0][1]
                , rowToCluster[1], weights[0]*rows[0][1] );
        }
        row_count = 0;
    }

    /**
     * @brief Join two clusters (corresponding with column a and
     *        row b in the matrix).
     * @param a the column 
     * @param b the row (a<b)
     * @note  It is assumed throughout that 0<a<b<row_count.
     * @note  UPGMA weights the clusters based on the number of leaf taxa
     *        they contain.  In the language of BIONJ, lambda (the weight
     *        to give, in the calculations, to the cluster that currently
     *        corresponds to row a in the matrix, relative to the overall
     *        weight of both the a and b clusters) is:
     *
     *        taxa(A) / (taxa(A) + taxa(B)).
     */
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

    /**
     * @brief Identify taxa that are effectively identical 
     *        (for the purposes of phylogenetic inference using a 
     *        neighbour-joining algorithm), because their rows in
     *        the distance matrix are identical.  For each "row 
     *        content equivalence class", construct balanced 
     *        subtrees, by joining them into clusters.
     * @note  hashing is used to speed up comparisons.
     * @note  calculateRowHashes() does most of the real work,
     *        and HashRow<T>::identifyDuplicateClusters doesn't
     *        really do very much.
     */
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

        std::vector<HashRow<T>> hashed_rows;
        calculateRowHashes(hashed_rows, show_progress);
        DuplicateTaxa vvc;
        HashRow<T>::identifyDuplicateClusters(hashed_rows, vvc);
 
        size_t dupes_clustered = joinUpDuplicateClusters(vvc, show_progress);
        #if (USE_PROGRESS_DISPLAY)
        show_progress.done();
        #endif
        if (0<dupes_clustered && !silent) {
            std::cout << "Clustered " << dupes_clustered
                      << " identical (or near-identical) taxa." << std::endl;
        }
    }

    /**
     * @brief Calculate hashes for all the rows in the distance matrix.
     *        Construct a vector of HashRow<T> (which is hash + row pointer)
     *        Sort it (first by hash, tie-breaking by comparing entry by
     *        entry in the rows).
     * @param hashed_rows   vector of HashRow<T> (output)
     * @param show_progress used to display progress report
     * @note  The hashing is parallelized over rows (if _OPENMP is defined
     *        and is non-zero).
     * @note  The sorting is NOT parallelize (the assumption is that it
     *        won't take very long). But perhaps it should be parallelized.
     */
    void calculateRowHashes(std::vector<HashRow<T>>& hashed_rows,
                            progress_display& show_progress) {
        hashed_rows.resize(row_count);
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (intptr_t i=0; i<row_count; ++i) {
            hashed_rows[i] = HashRow<T>(rowToCluster[i], rows[i], row_count);
            #if (USE_PROGRESS_DISPLAY)
                if ((i%1000)==999) {
                    show_progress += (double)(1000.0);
                }
            #endif
        }
        std::sort(hashed_rows.begin(), hashed_rows.end());
        #if (USE_PROGRESS_DISPLAY)
            show_progress += (double)(row_count%1000);
        #endif
    }
    
    /**
     * @brief Join up clusters, that have been identified as
     *        duplicates.
     * @param vvc           a vector of vectors of int
     *                      each item in vvc is a vector 
     *                      of row numbers.
     * @param show_progress used for displaying progress
     * @return size_t - how many rows of the matrix were
     *         gotten rid of by joining duplicate taxa,
     *         or clusters of same.
     * @note   it is assumed that, for any pair of duplicate clusters, a, b,
     *         the distance between other cluster i and the cluster, c,
     *         (the results of joining a and b) is the same as was the
     *         distance between i and a). If (only if!) a and b are duplicates, 
     *         this is true of NJ, UNJ, BIONJ, and UPGMA.
     */
    size_t joinUpDuplicateClusters(DuplicateTaxa& vvc,
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
        
        double work_done = 0.0;
        size_t dupes_removed = 0;
        for (std::vector<intptr_t>& vc : vvc) {
            double work_here = static_cast<double>(vc.size()) * work_per_dupe;
            dupes_removed += (vc.size()-1);
            while (vc.size()>1 && 3<row_count) {
                intptr_t first_half  = vc.size() / 2;          //half, rounded down
                intptr_t second_half = vc.size() - first_half; //half, rounded up
                for (intptr_t i=0; i<first_half && 3<row_count; ++i) {
                    intptr_t cluster_a = vc[i];
                    intptr_t row_a     = cluster_to_row[cluster_a];
                    intptr_t cluster_b = vc[i+second_half];
                    intptr_t row_b     = cluster_to_row[cluster_b];
                    intptr_t cluster_c = clusters.size();
                    if (row_b<row_a) {
                        std::swap(row_a, row_b);
                    }
                    intptr_t cluster_x = rowToCluster[row_count-1];
                    cluster(row_a, row_b);
                    vc[i] = cluster_c;
                    cluster_to_row.push_back(row_a);
                    cluster_to_row[cluster_x] = row_b;
                }
                vc.resize(second_half);
                //Not first_half (rounded down), but second_half (rounded up), 
                //because, if there was an odd cluster, it must be kept
                //in play.
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

    /**
     * @brief  Determine the imbalance (the absolute value of the 
     *         difference between the number of taxa in the cluster
     *         associated with row A, and the number of taxa in the
     *         cluster associated with row B) that will result if/when
     *         those clsuters are joined.
     * @param  rowA one row (0<=rowA<row_count)
     * @param  rowB another (0<=rowB<row_count, and rowB!=rowA)
     * @return size_t the imbalance.
     */
    size_t getImbalance(size_t rowA, size_t rowB) const {
        size_t clusterA = rowToCluster[rowA];
        size_t clusterB = rowToCluster[rowB];
        size_t sizeA    = clusters[clusterA].countOfExteriorNodes;
        size_t sizeB    = clusters[clusterB].countOfExteriorNodes;
        return (sizeA<sizeB) ? (sizeB-sizeA) : (sizeA-sizeB);
    }
};

/**
 * @brief Vectorized version of UPGMA_Matrix
 * @note  Overrides getRowMinima() with a vectorized implementation.
 */
#if USE_VECTORCLASS_LIBRARY
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

    virtual std::string getAlgorithmName() const override {
        return "Vectorized-" + super::getAlgorithmName();
    }

    virtual void calculateRowTotals() const override {
        size_t fluff = MATRIX_ALIGNMENT / sizeof(T);
        scratchColumnNumbers.resize(row_count + fluff, 0.0);
    }
    
    virtual void getRowMinima() const override {
        T* nums = matrixAlign ( scratchColumnNumbers.data() );
        rowMinima.resize(row_count);
        rowMinima[0].value = infiniteDistance;
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic)
        #endif
        for (intptr_t row=1; row<row_count; ++row) {
            Position<T> pos(row, 0, infiniteDistance, 0);
            const T*    rowData    = rows[row];
            intptr_t    col;
            V           minVector  = infiniteDistance;
            V           ixVector   = -1;

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
