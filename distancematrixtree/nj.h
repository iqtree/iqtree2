//
//  nj.h - Implementations of NJ, BIONJ, and UNJ algorithms.
//
//  NJ    implementation based on the same (but original NJ, without
//        a matrix of variance estimates (see NJMatrix).
//        Paper:  "The neighbor-joining method: a new method
//                for reconstructing phylogenetic trees",
//                Naurya Saitou and Masatoshi Nei (1987).
//        Tag:    [NS1987]
//        Class:  NJMatrix<T> (where T is the floating point type)
//        Notes:  Subclass of UPGMA_Matrix<T>.
//
//  BIONJ implementation based on http://www.lirmm.fr/~w3ifa/MAAS/BIONJ/BIONJ.html
//        (see BIONJMatrix).  Original authors: Olivier Gascuel
//        and Hoa Sien Cuong (the code for the Unix version).
//        Paper: "BIONJ: An Improved Version of the NJ Algorithm
//                Based on a Simple Model of Sequence Data" (2009).
//        Precis: Based on NJ. Adds a V (variance estimate) matrix
//                which is used to determine 
//                (doubles the memor requirements)
//        Tag:    [GAS2009].
//        Class:  BIONJMatrix<T>
//        Notes:  Subclass of NJMatrix.
//                Uses *square* matrices (rather than triangular matrices)
//                (square matrices make the program simpler, and gives its
//                 memory access patterns better spatial locality of reference,
//                 at the cost of doubling the memory requirements).
//                Calculations are "moved left of the summation" wherever possible
//                (the formulae used are equivialent to those in [GAS2009], but
//                with calculations moved outside of for-loops wherever possible,
//                and function calls avoided as much as was practicable).
//
//  UNJ   implementation based on a 1997 paper.
//        Paper: "Concerning the NJ algorithm and its unweighted version, UNJ",
//               Olivier Gascuel
//        TAG:   [GAS1997]
//        Class: UNJMatrix<T>
//        Note:  Subclass of NJMatrix.
//               UNJ's name ("unweighted...") is a bit misleading.
//               To recover the "unweighted" contributions of the taxa in
//               two clusters of sizes s1 and s2 it is necessary to reweight
//               contributions from both (e.g. d = d1*s1/(s1+s2) + d2*s2/(s1+s2)),
//               where d1 and d2 are distances from another cluster (cluster 3)
//               to clusters 1 and 2.
//               This makes UNJ slower than NJ, in partice.
//
// The vectorized implementations (of BIONJ and NJ) use Agner Fog's
// vectorclass library.
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

#ifndef nj_h
#define nj_h

#include "upgma.h"     //for UPGMA_Matrix template class
#include <utils/my_assert.h> //for ASSERT macro

namespace StartTree
{

/**
 * @brief  Neighbour Joining implementation
 * @tparam T the distance type
 * @note   The main structure is a D matrix 
 *         (a matrix of distances) and a U vector
 *         (a vector of row totals for that matrix)
 * @note   The textbook implementation uses formulae that feature 
 *         row totals, and multiplies distances by (n-2).  This 
 *         implementation calculates "scaled" row totals, which are
 *         row totals, U', multiplied by (1/(n-2)), so simpler
 *         formulae of the form Dij-U'i-U'j rather than (n-2)Dij-Ui-Uj,
 *         can be used (since we're just comparing magnitudes to determine
 *         which pair of clusters to join, it doesn't matter if those 
 *         magnitudes are all rescaled.
 *         This trick saves a multiplication on every evaluation, at the
 *         cost of row_count multiplications (to calculate U' from U)
 *         per iteration.
 */
template <class T=NJFloat> class NJMatrix: public UPGMA_Matrix<T>
{
public:
    typedef UPGMA_Matrix<T> super;
    using super::clusters;
    using super::rows;
    using super::row_count;
    using super::rowMinima;
    using super::rowTotals;
    using super::rowToCluster;
    using super::removeRowAndColumn;
    using super::calculateRowTotals;
    using super::getImbalance;
    using super::silent;
    using super::isRooted;
protected:
    mutable std::vector<T> scaledRowTotals; //used in getRowMinima
public:
    NJMatrix(): super() { }
    virtual std::string getAlgorithmName() const override {
        return "NJ";
    }
protected:
    /**
     * @brief Calculate scaled row totals, U' (in scaledRowTotals)
     *        by multiplying row totals, U (in rowTotals) by
     *        (1/(n-2)). If U' is used, it isn't necessarily to 
     *        multiply distances (from the D matrix) by (n-2)
     *        when comparing adjusted distances to determine which
     *        cluster to join.
     */
    virtual void calculateScaledRowTotals() const {
        scaledRowTotals.resize(row_count);
        T nless2      = (T)( row_count - 2 );
        T tMultiplier = ( row_count <= 2 ) ? (T)0.0 : ((T)1.0 / nless2);
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (intptr_t r=0; r<row_count; ++r) {
            scaledRowTotals[r] = rowTotals[r] * tMultiplier;
        }
    }
    /**
     * @brief Override of calculateRowTotals() - that ensures that the
     *        scaled row totals are also recalculated.
     */
    virtual void calculateRowTotals() const override {
        super::calculateRowTotals();
        calculateScaledRowTotals();
    }
    /**
     * @brief determine, for each row, r, which column,row pair of 
     *        clusters (with column c less than row r), would be the
     *        "best" choice for the next pair of clusters to join.
     * @note  scaled row totals are looked up "around" the  
     *        scaledRowTotals vector via a "naked" array pointer, 
     *        tot to avoid range-checking overhead.
     * @note  bestVrc is, initially, the "best" value of Drc - U'c 
     *        for the part of the current row in the lower-left 
     *        triangle (c between 0 and r-1 inclusive). U'r 
     *        is only subtracted *after* the minimum for (Drc - U'c)
     *        has been found.  We don't need to calculate 
     *        Drc - U'c - U'r except for the c for which Drc - U'c
     *        was minimized.
     */
    virtual void getRowMinima() const override {
        calculateScaledRowTotals();
        rowMinima.resize(row_count);
        rowMinima[0].value = infiniteDistance;
        const T* tot = this->scaledRowTotals.data();
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic)
        #endif
        for (intptr_t row=1; row<row_count; ++row) {
            T        bestVrc    = infiniteDistance;
            size_t   bestColumn = 0;
            const T* rowData    = rows[row];
            for (intptr_t col=0; col<row; ++col) {
                T    v      = rowData[col] - tot[col];
                bool better = ( v < bestVrc );
                if (better) {
                    bestColumn = col;
                    bestVrc = v;
                }
            }
            bestVrc       -= tot [row];
            rowMinima[row] = Position<T>(row, bestColumn, bestVrc
                                         , getImbalance(row, bestColumn));
        }
    }
    /**
     * @brief Join two clusters (corresponding with column a and
     *        row b in the matrix), and update row totals on the go.
     * @param a the column 
     * @param b the row (a<b)
     * @note  It is assumed throughout that 0<a<b<row_count.
     * @note  Most implementations of Neighbour Joining periodically
     *        recalculate row totals "from scratch" (the fear is that
     *        rounding error will otherwise accumulate).  Originally
     *        I had code to do that, periodically.  But it didn't 
     *        seem to make much real difference.  -James B.
     * @note  In the following c is the cluster that results by joining
     *        the clusters corresponding to rows a and b.  In practice,
     *        cluster c will (by the time this function returns) 
     *        correspond to the (a)th row/column in the matrix.
     *        And what *was* cluster (row_count-1) will correspond 
     *        to the (b)th row/column.
     * @note  If we think of NJ as a "variant" of BIONJ, it is as 
     *        though we always have lambda=0.5.  And mu (1-lambda) 
     *        is also 0.5.  That's why those variable names are used.      
     * @note  Perhaps scaledRowTotals[i] ought to be updated here, too.
     *        But then every subclass would have to do the same.      
     */
    virtual void cluster(intptr_t a, intptr_t b) override {
        T nless2        = (T)(row_count-2);
        T tMultiplier   = (row_count<3) ? (T)0.0 : ((T)0.5 / nless2);
        T lambda        = (T)0.5;
        T medianLength  = lambda * rows[a][b];
        T fudge         = (rowTotals[a] - rowTotals[b]) * tMultiplier;
        T aLength       = medianLength + fudge;
        T bLength       = medianLength - fudge;
        T mu            = (T)1.0 - lambda;
        T dCorrection   = - lambda * aLength - mu * bLength;
        auto aRow       = rows[a];
        auto bRow       = rows[b];
        T cTotal        = (T)0.0;
        #ifdef _OPENMP
        #pragma omp parallel for reduction(+:cTotal)
        #endif
        for (intptr_t i=0; i<row_count; ++i) {
            if (i!=a && i!=b) {
                T Dai         = aRow[i];
                T Dbi         = bRow[i];
                T Dci         = lambda * Dai + mu * Dbi + dCorrection;
                aRow[i]       = Dci;
                rows[i][a]    = Dci;
                rowTotals[i] += Dci - Dai - Dbi;
                cTotal       += Dci;
            }
        }
        clusters.addCluster ( rowToCluster[a], aLength,
                              rowToCluster[b], bLength);
        rowTotals[a]    = cTotal;
        rowToCluster[a] = clusters.size()-1;
        rowToCluster[b] = rowToCluster[row_count-1];
        removeRowAndColumn(b);
    }
    /**
     * @brief Finish the tree, by joining the last two or three clusters.
     *        (Conventionally there are three clusters, but if neighbour
     *        joining is being used to cluster a subtree, when honouring
     *        a constraint tree, depending on what the caller wants, the
     *        "top" node of the tree might have degree 2 rather than 3).
     */
    virtual void finishClustering() override {
        ASSERT( row_count == 2 || row_count == 3);
        T halfD01 = (T)0.5 * rows[0][1];
        if (row_count==3) {
            T halfD02 = (T)0.5 * rows[0][2];
            T halfD12 = (T)0.5 * rows[1][2];
            clusters.addCluster
                ( rowToCluster[0], halfD01 + halfD02 - halfD12
                , rowToCluster[1], halfD01 + halfD12 - halfD02
                , rowToCluster[2], halfD02 + halfD12 - halfD01);
            row_count = 0;
        }
        else {
            clusters.addCluster
                ( rowToCluster[0], halfD01
                , rowToCluster[1], halfD01);
            row_count = 0;
        }
    }
};

/**
 * @brief  Unweighted Neighbour Joining implementation
 * @tparam T the distance type.
 * @note   This is based on NJMatrix<T> (and why not, when the only
 *         difference between the two algorithms is how distances
 *         between newly joined clusters, and other existing clusters
 *         are calculated.
 */
template <class T=NJFloat> class UNJMatrix: public NJMatrix<T> {
protected:
    intptr_t original_n; //Needed in formulae used in distance formulae
                         //during clutering.
public:
    typedef NJMatrix<T> super;
    UNJMatrix(): super(), original_n(0) { }
    virtual std::string getAlgorithmName() const override {
        return "UNJ";
    }
    virtual void setSize(intptr_t rank) override {
        super::setSize(rank);
        original_n = rank;
    }

protected:
    using super::row_count;
    using super::rows;
    using super::rowTotals;
    using super::clusters;
    using super::rowToCluster;
    using super::removeRowAndColumn;
    
    /**
     * @brief Join two clusters (those corresponding to rows a and b
     *        in the working distance matrix).
     * @param a the lower-numbered row (or, if you prefer the column,
     *          in the lower-left triangle of the matrix)
     * @param b the higher-numbered row (a<b)
     * @note  It is assumed a<b.
     * @note  The clusters are weighted in terms of the number of 
     *        taxa they contain (aCount and bCount), as per [GAS1997].
     *        (Conceptually, the leaf nodes for the taxa are NOT weighted,
     *        and standard NJ is "downweighting" each taxon, with a weight 
     *        of 1/aCount in the cluster for row a, and by a weight 
     *        of 1/bCount in the cluster for row b; but in practice, 
     *        "undoing the effect of the weighting" requires more 
     *        multiplication, and more time).  
     *        So UNJMatrix runs slower than NJMatrix! The End. -James B.
     * @note  The greek letter lambda is given a different role
     *        in [GAS1997], but I wanted to use it in a fashion
     *        a bit more consistent with later BIONJ implementations
     *        (since BIONJ is much more widely used). -James B.
     * @note  mu is just (1-lambda).
     */
    virtual void cluster(intptr_t a, intptr_t b) override {
        size_t aCount     = clusters[rowToCluster[a]].countOfExteriorNodes;
        size_t bCount     = clusters[rowToCluster[b]].countOfExteriorNodes;
        size_t cCount     = aCount + bCount;
        T tMultiplier     = (row_count<3) ? (T)0.0 : ((T)0.5 / (T)( original_n - cCount));
        T lambda          = (T)aCount / (T)(aCount + bCount); //relative weight of a in u
        T mu              = (T)1.0 - lambda;                  //relative weight of b in u
        auto aRow         = rows[a];
        auto bRow         = rows[b];
        
        T aFudge          = 0.0;
        #ifdef _OPENMP
        #pragma omp parallel for reduction(+:aFudge)
        #endif
        for (intptr_t i=0; i<row_count; ++i) {
            if (i!=a && i!=b) {
                T iCount  = (T)(clusters[rowToCluster[i]].countOfExteriorNodes);
                aFudge   += iCount * (aRow[i] - bRow[i]);
                //reweight the contribution for the cluster in the ith row
                //according to the number (iCount) of leaf nodes it contains,
                //as per estimation formula (4) in [GAS1997].
            }
        }
        T abLength        = rows[a][b];
        T auLength        = (T)0.5 * abLength + aFudge * tMultiplier;
        T buLength        = abLength - auLength;
        T dCorrection     = - lambda * auLength - mu * buLength;
        T cTotal          = (T)0.0;
        #ifdef _OPENMP
        #pragma omp parallel for reduction(+:cTotal)
        #endif
        for (intptr_t i=0; i<row_count; ++i) {
            if (i!=a && i!=b) {
                T Dai         = aRow[i];
                T Dbi         = bRow[i];
                T Dci         = lambda * Dai + mu * Dbi + dCorrection;
                aRow[i]       = Dci;
                rows[i][a]    = Dci;
                rowTotals[i] += Dci - Dai - Dbi;
                                //JB2020-06-18 Adjust row totals on fly
                cTotal       += Dci;
            }
        }
        clusters.addCluster ( rowToCluster[a], auLength,
                              rowToCluster[b], buLength);
        rowTotals[a]    = cTotal;
        rowToCluster[a] = clusters.size()-1; //cluster u
        rowToCluster[b] = rowToCluster[row_count-1];
        removeRowAndColumn(b);
    }
};

/**
 * @brief  Implementation of the BIONJ neighbour joining variant
 * @tparam T the distance type
 * @note   The big difference between BIONJ and NJ is that BIONJ
 *         maintains an auxiliary V matrix, of "estimated variances"
 *         of "positions" of clusters in a notional genetic "space".
 *         The variance member is V.  Row and column rearrangements
 *         in V duplicate those in the primary distance matrix, D.
 *         Values in the V matrix are used for determining lambda
 *         (the "weight" to be assigned to the first cluster, when 
 *         merging two clusters, corresponding to rows a and b, a<b,
 *         in the D matrix).
 */
template <class T=NJFloat> class BIONJMatrix : public NJMatrix<T> {
public:
    typedef NJMatrix<T> super;
    using super::clusters;
    using super::row_count;
    using super::column_count;
    using super::rows;
    using super::rowToCluster;
    using super::rowTotals;
    using super::removeRowAndColumn;
    using super::silent;
protected:
    SquareMatrix<T>  variance;       //The V matrix
public:
    virtual std::string getAlgorithmName() const override {
        return "BIONJ";
    }
    virtual bool loadMatrixFromFile(const std::string &distanceMatrixFilePath) override
    {
        bool rc = super::loadMatrixFromFile(distanceMatrixFilePath);
        return rc;
    }
    virtual bool loadMatrix(const StrVector& names, 
                            const double* matrix) override {
        bool rc = super::loadMatrix(names, matrix);
        return rc;
    }
    /**
     * @brief Calculate lambda, the weight to give distances 
     *        to the cluster that currently corresponds to
     *        row a of the distance matrix, if it is joined
     *        with the custer that currently corresponds to
     *        row b of the distance matrix (given row numbers
     *        a and b, 0<=a<b<row_count).
     * @param a   the lower row number  (0<=a)
     * @param b   the higher row number (a<b<row_count)
     * @param Vab the estimated variance of the genetic 
     *            distance between the clusters that currently
     *            correspond to rows a and b of the distance matrix.
     * @return T  lambda, the weight to give distances to the a cluster,
     *            between 0 and 1.  mu, the weight to give distances
     *            to the b cluster, is always (1-lambda).
     * @note      The code that adds for all i, and then subtracts
     *            for a and b, runs faster than would code that,
     *            for all i, checked that i!=a && i!=b before deciding
     *            whether to add.  Yes, there's a slight increase in 
     *            rounding error, but the performance benefit appears
     *            to be worth it. -James B.
     */
    inline T chooseLambda(intptr_t a, intptr_t b, T Vab) {
        T lambda = 0;
        if (Vab==0.0) {
            return 0.5; //special case (to avoid division by zero)
        }
        auto vRowA = variance.rows[a];
        auto vRowB = variance.rows[b];
        for (intptr_t i=0; i<column_count; ++i) {
            lambda += vRowB[i] - vRowA[i];
        }
        if (a<column_count) {
            lambda += vRowA[a] - vRowB[a];
        }
        if (a!=b && b<column_count) {
            lambda += vRowA[b] - vRowB[b];
        }
        lambda = (T)0.5 + lambda / ((T)2.0*((T)row_count-2)*Vab);
        if (1.0<lambda) lambda=(T)1.0;
        if (lambda<0.0) lambda=(T)0.0;
        return lambda;
    }
    /**
     * @brief Join the clusters that currently correspond to rows
     *        a and b of the distance matrix.  Update row totals
     *        (but not scaled row totals), and variances, too.      
     * @param a the lower row number  (0<=a<b)
     * @param b the higher row number (a<b<row_count)
     * @note  The bits of this that differ, in a meaningful way, from 
     *        the corresponding code in NJMatrix are tagged with BIO.
     *        In NJ, lambda is always 0.5, and elements in V for the new
     *        cluster (which overwrite what is in row a in V), need to
     *        be calcualted.  And lastly, what is done to the D matrix
     *        to shrink it, also has to be "echoed" on the V matrix.
     */
    virtual void cluster(intptr_t a, intptr_t b) override {
        //Assumed 0<=a<b<n
        //Bits that differ from super::cluster tagged BIO
        T nless2          = (T)(row_count - 2);
        T tMultiplier     = ( row_count < 3 ) ? (T)0.0 : ( (T)0.5 / nless2 );
        T medianLength    = (T)0.5 * rows[b][a];
        T fudge           = (rowTotals[a] - rowTotals[b]) * tMultiplier;
        T aLength         = medianLength + fudge;
        T bLength         = medianLength - fudge;
        T Vab             = variance.rows[b][a];     //BIO
        T lambda          = chooseLambda(a, b, Vab); //BIO
        T mu              = (T)1.0 - lambda;
        T dCorrection     = - lambda * aLength - mu * bLength;
        T vCorrection     = - lambda * mu * Vab;
        auto rowA         = rows[a];
        auto rowB         = rows[b];
        auto varianceRowA = variance.rows[a];
        auto varianceRowB = variance.rows[b];
        T cTotal          = (T)0.0;
        #ifdef _OPENMP
        #pragma omp parallel for reduction(+:cTotal)
        #endif
        for (intptr_t i=0; i<row_count; ++i) {
            if (i!=a && i!=b) {
                //Dci as per reduction 4 in [GAS2009]
                T Dai         = rowA[i];
                T Dbi         = rowB[i];
                T Dci         = lambda * Dai + mu * Dbi + dCorrection;
                rowA[i]       = Dci;
                rows[i][a]    = Dci;
                rowTotals[i] += Dci - Dai - Dbi; //JB2020-06-18 Adjust row totals
                cTotal       += Dci;
                
                //BIO begin (Reduction 10 on variance estimates) in [GAS2009]
                T Vci   = lambda * varianceRowA[i]
                        + mu * varianceRowB[i]
                        + vCorrection;
                varianceRowA[i] = Vci;
                variance.rows[i][a] = Vci;
                //BIO finish
            }
        }
        clusters.addCluster ( rowToCluster[a], aLength, rowToCluster[b], bLength);
        rowTotals[a]    = cTotal;
        rowToCluster[a] = clusters.size()-1;
        rowToCluster[b] = rowToCluster[row_count-1];
        removeRowAndColumn(b);
        variance.removeRowAndColumn(b); //BIO
    }
    virtual void prepareToConstructTree() override {
        variance = *this;
    }
};

#if USE_VECTORCLASS_LIBRARY
/**
 * @brief 
 * 
 * @tparam T the distance type
 * @tparam SUPER the superclass (in practice, NJMatrix<T>, 
 *         UNJMatrix<T>, or BIONJMatrix<T>) this class is to
 *         vectorize.  
 * @tparam V  the vector type
 * @tparam VB the boolean vector type
 * @note   This works by overriding getRowMinima() with 
 *         a vectorized version.
 * @note   This can't subclass UPGMAMatrix<T> (so don't try
 *         using it with SUPER=UPGMAMatrix<T>!), because 
 *         that has a rather different getRowMinima().
 *         There's a separate class for vectorizing that.
 */
template <class T=NJFloat, class SUPER=BIONJMatrix<T>, 
          class V=FloatVector, class VB=FloatBoolVector>
class VectorizedMatrix: public SUPER
{
    typedef SUPER super;
    using super::rows;
    using super::row_count;
    using super::rowMinima;
    using super::rowTotals;
    using super::calculateRowTotals;
    using super::getImbalance;
protected:
    mutable std::vector<T> scratchTotals;
    mutable std::vector<T> scratchColumnNumbers;
    const intptr_t  blockSize;
public:
    VectorizedMatrix() : super(), blockSize(VB().size()) {
    }
    virtual std::string getAlgorithmName() const {
        return "Vectorized-" + super::getAlgorithmName();
    }
    /**
     * @brief Override for calculateRowTotals().
     *        Make sure the scaled total and row number
     *        scratch vectors are large enough.
     */
    virtual void calculateRowTotals() const {
        super::calculateRowTotals();
        size_t fluff = MATRIX_ALIGNMENT / sizeof(T);
        scratchTotals.resize(row_count + fluff, 0.0);
        scratchColumnNumbers.resize(row_count + fluff, 0.0);
    }
    /**
     * @brief Get the Row Minima object
     * @note  tot is an *aligned* pointer to 
     *        elements in a "scratch" vector of rescaled
     *        row totals, calculated each time this member
     *        is called.
     * @note  nums is an *aligned* pointer to elements in
     *        a "scratch" vector, which is used to keep 
     *        track of *where* (which is to say, which)
     *        column, each adjusted distance in the current
     *        row comes from.
     */
    virtual void getRowMinima() const {
        T nless2      = (T)( row_count - 2 );
        T tMultiplier = ( row_count <= 2 ) ? (T)0.0 : ((T)1.0 / nless2);
        T* tot  = matrixAlign ( scratchTotals.data() );
        T* nums = matrixAlign ( scratchColumnNumbers.data() );
        for (intptr_t r=0; r<row_count; ++r) {
            tot[r]  = rowTotals[r] * tMultiplier;
            nums[r] = (T)r;
        }
        rowMinima.resize(row_count);
        rowMinima[0].value = infiniteDistance;
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic)
        #endif
        for (intptr_t row=1; row<row_count; ++row) {
            Position<T> pos(row, 0, infiniteDistance, 0);
            const T* rowData   = rows[row];
            intptr_t col;
            V        minVector = infiniteDistance;
                     //The minima of columns with indices
                     //"congruent modulo blockSize"
                     //For example, if blockSize is 4,
                     //minVector[1] holds the minimum of
                     //columns 1,5,9,13,17,...
            V        ixVector  = -1;
                     //For each entry in minVector, the column 
                     //from which that value came.
            
            intptr_t colStop = row - blockSize;
            V  rowVector(0);
            V  totVector(0);
            V  numVector(0);
            for (col=0; col<colStop; col+=blockSize) {
                rowVector.load_a(rowData+col);
                totVector.load_a(tot+col);
                V  adjVector = rowVector - totVector;
                VB less      = adjVector < minVector;
                numVector.load_a(nums+col);
                ixVector     = select(less, numVector, ixVector);
                minVector    = select(less, adjVector, minVector);
            }
            //Extract minimum and column number
            for (int c=0; c<blockSize; ++c) {
                if (minVector[c] < pos.value) {
                    pos.value  = minVector[c];
                    pos.column = (size_t)(ixVector[c]);
                }
            }
            for (; col<row; ++col) {
                T v = rowData[col] - tot[col];
                if (v < pos.value) {
                    pos.column = col;
                    pos.value  = v;
                }
            }
            pos.value     -= tot [row];
            pos.imbalance  = getImbalance(pos.row, pos.column);
            rowMinima[row] = pos;
        }
    }
};//end of class

    typedef VectorizedMatrix<NJFloat, NJMatrix<NJFloat>>    VectorNJ;
    typedef VectorizedMatrix<NJFloat, BIONJMatrix<NJFloat>> VectorBIONJ;
#endif
}
#endif /* nj_h */
