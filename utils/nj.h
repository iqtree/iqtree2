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
//               To recover the "unweighted" contributions of the taxons in
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
#include "my_assert.h" //for ASSERT macro

namespace StartTree
{
template <class T=NJFloat> class NJMatrix: public UPGMA_Matrix<T>
    //NJMatrix is a D matrix (a matrix of distances).
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
protected:
    mutable std::vector<T> scaledRowTotals; //used in getRowMinima
public:
    NJMatrix(): super() { }
    virtual std::string getAlgorithmName() const {
        return "NJ";
    }
protected:
    virtual void calculateScaledRowTotals() const {
        //
        //Note: Rather than multiplying distances by (n-2)
        //      repeatedly, it is cheaper to work with row
        //      totals multiplied by (1/(T)(n-2)).
        //      Better n multiplications than n*(n-1)/2.
        //
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
    virtual void calculateRowTotals() const {
        super::calculateRowTotals();
        calculateScaledRowTotals();
    }
    virtual void getRowMinima() const {
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
    virtual void cluster(intptr_t a, intptr_t b) {
        //Cluster two active rows, identified by row indices a and b).
        //Assumed 0<=a<b<n
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
                                //JB2020-06-18 Adjust row totals on fly
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
    virtual void finishClustering() {
        ASSERT( row_count == 3);
        T halfD01 = (T)0.5 * rows[0][1];
        T halfD02 = (T)0.5 * rows[0][2];
        T halfD12 = (T)0.5 * rows[1][2];
        clusters.addCluster
            ( rowToCluster[0], halfD01 + halfD02 - halfD12
            , rowToCluster[1], halfD01 + halfD12 - halfD02
            , rowToCluster[2], halfD02 + halfD12 - halfD01);
        row_count = 0;
    }
};

template <class T=NJFloat> class UNJMatrix: public NJMatrix<T> {
protected:
    intptr_t original_n;
public:
    typedef NJMatrix<T> super;
    UNJMatrix(): super(), original_n(0) { }
    virtual std::string getAlgorithmName() const {
        return "UNJ";
    }
    virtual void setSize(intptr_t rank) {
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
    
    virtual void cluster(intptr_t a, intptr_t b) {
        //Cluster two active rows, identified by row indices a and b).
        //Assumes: 0<=a<b<n
        //
        //The clusters are weighted in terms of the number of taxa they
        //contain (aCount and bCount), as per [GAS1997].
        //(Conceptually, the leaf nodes for the taxa are NOT weighted,
        // and standard NJ is "downweighting" each taxon, with a weight of
        // 1/aCount in the cluster for row a, and by a weight of
        // 1/bCount in the cluster for row b; but in practice, "undoing
        // the effect of the weighting" requires more multiplication, and
        // more time).
        //
        //Note: The greek letter lambda is given a different role
        //      in [GAS1997], but I wanted to use it in a fashion
        //      a bit more consistent with later BIONJ implementations.
        //      -James B. 19-Oct-2020.
        //
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
        #pragma omp parallel for
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
    virtual std::string getAlgorithmName() const {
        return "BIONJ";
    }
    virtual bool loadMatrixFromFile(const std::string &distanceMatrixFilePath)
    {
        bool rc = super::loadMatrixFromFile(distanceMatrixFilePath);
        variance = *this;
        return rc;
    }
    virtual bool loadMatrix(const StrVector& names, 
                            const double* matrix) {
        bool rc = super::loadMatrix(names, matrix);
        variance = *this;
        return rc;
    }
    inline T chooseLambda(intptr_t a, intptr_t b, T Vab) {
        //Assumed 0<=a<b<n
        T lambda = 0;
        if (Vab==0.0) {
            return 0.5;
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
    virtual void cluster(intptr_t a, intptr_t b) {
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
};

#ifdef USE_VECTORCLASS_LIBRARY
template <class T=NJFloat, class SUPER=BIONJMatrix<T>, 
          class V=FloatVector, class VB=FloatBoolVector>
class VectorizedMatrix: public SUPER
{
    typedef SUPER super;
    //
    //Note: this is a first attempt at hand-vectorizing
    //      BIONJMatrix::getRowMinima (via Agner Fog's
    //      vectorclass library).
    //      It can subclass either NJMatrix or BIONJMatrix.
    //      It cannot subclass UPGMA_Matrix.
    //
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
    virtual void calculateRowTotals() const {
        super::calculateRowTotals();
        size_t fluff = MATRIX_ALIGNMENT / sizeof(T);
        scratchTotals.resize(row_count + fluff, 0.0);
        scratchColumnNumbers.resize(row_count + fluff, 0.0);
    }
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
            for (col=0; col<colStop; col+=blockSize) {
                V  rowVector; rowVector.load_a(rowData+col);
                V  totVector; totVector.load_a(tot+col);
                V  adjVector = rowVector - totVector;
                VB less      = adjVector < minVector;
                V  numVector; numVector.load_a(nums+col);
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
