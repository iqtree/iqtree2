//
//  nj.h
//  iqtree
//
//  BIONJ implementation based on http://www.lirmm.fr/~w3ifa/MAAS/BIONJ/BIONJ.html
//        (see BIONJMatrix).  Original authors: Olivier Gascuel
//        and Hoa Sien Cuong (the code for the Unix version).
//        Paper: "BIONJ: An Improved Version of the NJ Algorithm
//                Based on a Simple Model of Sequence Data" (2009).
//        Tag:   [GAS2009].
//  NJ    implementation based on the same (but original NJ, without
//        a matrix of variance estimates (see NJMatrix).
//        Paper: "The neighbor-joining method: a new method
//               for reconstructing phylogenetic trees",
//               Naurya Saitou and Masatoshi Nei (1987).
//        Tag:   [NS1987]
//  UNJ   implementation based on
//        Paper: "Concerning the NJ algorithm and its unweighted version, UNJ",
//               Olivier Gascuel
//        TAG:   [GAS1997]
//
// The vectorized implementations (of BIONJ and NJ) use Agner Fog's
// vectorclass library.

//  Created by James Barbetti on 31/10/20.
//

#ifndef nj_h
#define nj_h

#include "upgma.h"
#include "tools.h" //for ASSERT macro

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
        T nless2      = ( row_count - 2 );
        T tMultiplier = ( row_count <= 2 ) ? 0 : (1 / nless2);
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (size_t r=0; r<row_count; ++r) {
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
        for (size_t row=1; row<row_count; ++row) {
            float  bestVrc    = infiniteDistance;
            size_t bestColumn = 0;
            const T* rowData = rows[row];
            for (size_t col=0; col<row; ++col) {
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
    virtual void cluster(size_t a, size_t b) {
        //Cluster two active rows, identified by row indices a and b).
        //Assumed 0<=a<b<n
        T nless2        = row_count-2;
        T tMultiplier   = (row_count<3) ? 0 : (0.5 / nless2);
        T lambda        = 0.5;
        T medianLength  = lambda * rows[a][b];
        T fudge         = (rowTotals[a] - rowTotals[b]) * tMultiplier;
        T aLength       = medianLength + fudge;
        T bLength       = medianLength - fudge;
        T mu            = 1.0 - lambda;
        T dCorrection   = - lambda * aLength - mu * bLength;
        auto aRow       = rows[a];
        auto bRow       = rows[b];
        T cTotal        = 0;
        #ifdef _OPENMP
        #pragma omp parallel for reduction(+:cTotal)
        #endif
        for (size_t i=0; i<row_count; ++i) {
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
        T halfD01 = 0.5 * rows[0][1];
        T halfD02 = 0.5 * rows[0][2];
        T halfD12 = 0.5 * rows[1][2];
        clusters.addCluster
            ( rowToCluster[0], halfD01 + halfD02 - halfD12
            , rowToCluster[1], halfD01 + halfD12 - halfD02
            , rowToCluster[2], halfD02 + halfD12 - halfD01);
        row_count = 0;
    }
};

template <class T=NJFloat> class UNJMatrix: public NJMatrix<T> {
protected:
    size_t original_n;
public:
    typedef NJMatrix<T> super;
    UNJMatrix(): super(), original_n(0) { }
    virtual std::string getAlgorithmName() const {
        return "UNJ";
    }
    virtual void setSize(size_t rank) {
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
    
    virtual void cluster(size_t a, size_t b) {
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
        T aCount          = clusters[rowToCluster[a]].countOfExteriorNodes;
        T bCount          = clusters[rowToCluster[b]].countOfExteriorNodes;
        T cCount          = aCount + bCount;
        T tMultiplier     = (row_count<3) ? 0 : (0.5 / ( original_n - cCount));
        T lambda          = aCount / (aCount + bCount); //relative weight of a in u
        T mu              = 1.0 - lambda;               //relative weight of b in u
        auto aRow         = rows[a];
        auto bRow         = rows[b];
        
        T aFudge          = 0.0;
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (size_t i=0; i<row_count; ++i) {
            if (i!=a && i!=b) {
                T iCount  = clusters[rowToCluster[i]].countOfExteriorNodes;
                aFudge   += iCount * (aRow[i] - bRow[i]);
                //reweight the contribution for the cluster in the ith row
                //according to the number (iCount) of leaf nodes it contains,
                //as per estimation formula (4) in [GAS1997].
            }
        }
        T abLength        = rows[a][b];
        T auLength        = 0.5 * abLength + aFudge * tMultiplier;
        T buLength        = abLength - auLength;
        T dCorrection     = - lambda * auLength - mu * buLength;
        T cTotal          = 0.0;
        #ifdef _OPENMP
        #pragma omp parallel for reduction(+:cTotal)
        #endif
        for (size_t i=0; i<row_count; ++i) {
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
    virtual bool loadMatrix(const std::vector<std::string>& names, const double* matrix) {
        bool rc = super::loadMatrix(names, matrix);
        variance = *this;
        return rc;
    }
    inline T chooseLambda(size_t a, size_t b, T Vab) {
        //Assumed 0<=a<b<n
        T lambda = 0;
        if (Vab==0.0) {
            return 0.5;
        }
        auto vRowA = variance.rows[a];
        auto vRowB = variance.rows[b];
        for (size_t i=0; i<column_count; ++i) {
            lambda += vRowB[i] - vRowA[i];
        }
        if (a<column_count) {
            lambda += vRowA[a] - vRowB[a];
        }
        if (a!=b && b<column_count) {
            lambda += vRowA[b] - vRowB[b];
        }
        lambda = 0.5 + lambda / (2.0*((T)row_count-2)*Vab);
        if (1.0<lambda) lambda=1.0;
        if (lambda<0.0) lambda=0.0;
        return lambda;
    }
    virtual void cluster(size_t a, size_t b) {
        //Assumed 0<=a<b<n
        //Bits that differ from super::cluster tagged BIO
        T nless2          = row_count - 2 ;
        T tMultiplier     = ( row_count < 3 ) ? 0 : ( 0.5 / nless2 );
        T medianLength    = 0.5 * rows[b][a];
        T fudge           = (rowTotals[a] - rowTotals[b]) * tMultiplier;
        T aLength         = medianLength + fudge;
        T bLength         = medianLength - fudge;
        T Vab             = variance.rows[b][a];     //BIO
        T lambda          = chooseLambda(a, b, Vab); //BIO
        T mu              = 1.0 - lambda;
        T dCorrection     = - lambda * aLength - mu * bLength;
        T vCorrection     = - lambda * mu * Vab;
        auto rowA         = rows[a];
        auto rowB         = rows[b];
        auto varianceRowA = variance.rows[a];
        auto varianceRowB = variance.rows[b];
        T cTotal          = 0;
        #ifdef _OPENMP
        #pragma omp parallel for reduction(+:cTotal)
        #endif
        for (size_t i=0; i<row_count; ++i) {
            if (i!=a && i!=b) {
                //Dci as per reduction 4 in [GAS2009]
                T Dai         = rowA[i];
                T Dbi         = rowB[i];
                T Dci         = lambda * Dai + mu * Dbi + dCorrection;
                rowA[i]       = Dci;
                rows[i][a]    = Dci;
                rowTotals[i] += Dci - Dai - Dbi; //JB2020-06-18 Adjust row totals
                cTotal       += Dci;
                
                //BIO begin (Reduction 10 on variance estimates)
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

template <class T=NJFloat, class super=BIONJMatrix<T>, class V=FloatVector, class VB=FloatBoolVector>
    class VectorizedMatrix: public super
{
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
    const size_t  blockSize;
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
        T nless2      = ( row_count - 2 );
        T tMultiplier = ( row_count <= 2 ) ? 0 : (1 / nless2);
        T* tot  = matrixAlign ( scratchTotals.data() );
        T* nums = matrixAlign ( scratchColumnNumbers.data() );
        for (size_t r=0; r<row_count; ++r) {
            tot[r]  = rowTotals[r] * tMultiplier;
            nums[r] = r;
        }
        rowMinima.resize(row_count);
        rowMinima[0].value = infiniteDistance;
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic)
        #endif
        for (size_t row=1; row<row_count; ++row) {
            Position<T> pos(row, 0, infiniteDistance, 0);
            const T* rowData   = rows[row];
            size_t   col;
            V        minVector = infiniteDistance;
                     //The minima of columns with indices
                     //"congruent modulo blockSize"
                     //For example, if blockSize is 4,
                     //minVector[1] holds the minimum of
                     //columns 1,5,9,13,17,...
            V        ixVector  = -1;
                     //For each entry in minVector, the column from which
                     //that value came.
            
            for (col=0; col+blockSize<row; col+=blockSize) {
                V  rowVector; rowVector.load_a(rowData+col);
                V  totVector; totVector.load_a(tot+col);
                V  adjVector = rowVector - totVector;
                VB less      = adjVector < minVector;
                V  numVector; numVector.load_a(nums+col);
                ixVector  = select(less, numVector, ixVector);
                minVector = select(less, adjVector, minVector);
            }
            //Extract minimum and column number
            for (int c=0; c<blockSize; ++c) {
                if (minVector[c] < pos.value) {
                    pos.value  = minVector[c];
                    pos.column = ixVector[c];
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
}

#endif /* nj_h */
