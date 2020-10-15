//
//  bionj2.cpp - Implementations of NJ and BIONJ algorithms
//               (that work in terms of .mldist inputs and
//                NEWICK outputs).
//
//  Copyright (C) 2020, James Barbetti.
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
// The vectorized implementations (of BIONJ and NJ) use Agner Fog's
// vectorclass library.
//
// Short names used for matrices and vectors (all indices start at 0).
//   D distance matrix (input, read from an .mldist file)
//   V estimated variance matrix (used in BIONJ, but not in NJ)
//   S bottom left triangle of the distance matrix,
//     with each row sorted in ascending order (see [SMP2011])
//   I index matrix, indicating (for each row of S, which cluster
//     each distance was calculated for).
//   Q connection-cost matrix (doesn't actually exist, cells of it
//     are calculated).
//   U vector of row totals (each row, the sum of the corresponding
//     row of the D matrix).  BIONJ implementations also use vectors
//     (indexed by cluster) of cluster totals.
//
// Created by James Barbetti on 18/6/2020.
// Note: Before 12-Aug, the Matrix and ClusterTree template
//       classes were declared in this file (rather than in
//       separate headers)
//

#include "starttree.h"
#include "distancematrix.h"          //for Matrix template class
#include "clustertree.h"             //for ClusterTree template class
#include "heapsort.h"                //for mirroredHeapsort, used to sort
                                     //rows of the S and I matrices
                                     //See [SMP2011], section 2.5.
#include <vector>                    //for std::vector
#include <string>                    //sequence names stored as std::string
#include <vectorclass/vectorclass.h> //for Vec4d and Vec4db vector classes
#include "progress.h"                //for progress_display

typedef float   NJFloat;
typedef Vec8f   FloatVector;
typedef Vec8fb  FloatBoolVector;
const   NJFloat infiniteDistance = 1e+36;
const   int     notMappedToRow = -1;

namespace StartTree
{

template <class T=NJFloat> struct Position
{
    //A position (row, column) in an NJ matrix
    //Note that column is always less than row.
    //(Because that's the convention in RapidNJ).
public:
    size_t  row;
    size_t  column;
    T       value;
    size_t  imbalance;
    Position() : row(0), column(0), value(0), imbalance(0) {}
    Position(size_t r, size_t c, T v, size_t imbalance)
        : row(r), column(c), value(v) {}
    Position& operator = (const Position &rhs) {
        row       = rhs.row;
        column    = rhs.column;
        value     = rhs.value;
        imbalance = rhs.value;
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
    virtual void setSize(size_t rank) {
        super::setSize(rank);
        rowToCluster.clear();
        for (int r=0; r<row_count; ++r) {
            rowToCluster.emplace_back(r);
        }
    }
    virtual void addCluster(const std::string &name) {
        clusters.addCluster(name);
    }
    bool loadMatrixFromFile(const std::string &distanceMatrixFilePath) {
        bool rc = loadDistanceMatrixInto(distanceMatrixFilePath, true, *this);
        calculateRowTotals();
        return rc;
    }
    virtual bool loadMatrix(const std::vector<std::string>& names, const double* matrix) {
        //Assumptions: 2 < names.size(), all names distinct
        //  matrix is symmetric, with matrix[row*names.size()+col]
        //  containing the distance between taxon row and taxon col.
        setSize(names.size());
        clusters.clear();
        for (auto it = names.begin(); it != names.end(); ++it) {
            clusters.addCluster(*it);
        }
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (size_t row=0; row<row_count; ++row) {
            const double* sourceStart = matrix + row * column_count;
            const double* sourceStop  = sourceStart + column_count;
            T*            dest        = rows[row];
            for (const double* source=sourceStart;
                 source<sourceStop; ++source, ++dest ) {
                *dest = (T) *source;
            }
        }
        calculateRowTotals();
        return true;
    }
    virtual bool constructTree() {
        Position<T> best;
        std::string taskName = "Constructing " + getAlgorithmName() + " tree";
        if (silent) {
            taskName="";
        }
        double triangle = row_count * (row_count + 1.0) * 0.5;
        progress_display show_progress(triangle, taskName.c_str(), "", "");
        while ( 3 < row_count ) {
            getMinimumEntry(best);
            cluster(best.column, best.row);
            show_progress += row_count;
        }
        finishClustering();
        show_progress.done();
        return true;
    }
    virtual void setZippedOutput(bool zipIt) {
        isOutputToBeZipped = zipIt;
    }
    virtual void beSilent() {
        silent = true;
    }
    bool writeTreeFile(const std::string &treeFilePath) const {
        return clusters.writeTreeFile(isOutputToBeZipped, treeFilePath);
    }
protected:
    void getMinimumEntry(Position<T> &best) {
        getRowMinima();
        best.value = infiniteDistance;
        for (size_t r=0; r<row_count; ++r) {
            Position<T> & here = rowMinima[r];
            if (here.value < best.value) {
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
        for (size_t row=1; row<row_count; ++row) {
            float  bestVrc    = infiniteDistance;
            size_t bestColumn = 0;
            const  T* rowData = rows[row];
            for (size_t col=0; col<row; ++col) {
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
        //Assumes: n is always exactly 3
        //But:  The formula is probably wrong. Felsenstein [2004] chapter 11 only
        //      covers UPGMA for rooted trees, and I don't know what
        //      the right formula is.
        //
        T weights[3];
        T denominator = 0;
        for (size_t i=0; i<3; ++i) {
            weights[i] = clusters[rowToCluster[i]].countOfExteriorNodes;
            denominator += weights[i];
        }
        for (size_t i=0; i<3; ++i) {
            weights[i] /= (2.0 * denominator);
        }
        clusters.addCluster
            ( rowToCluster[0], weights[1]*rows[0][1] + weights[2]*rows[0][2]
            , rowToCluster[1], weights[0]*rows[0][1] + weights[2]*rows[1][2]
            , rowToCluster[2], weights[0]*rows[0][2] + weights[1]*rows[1][2]);
        row_count = 0;
    }
    virtual void cluster(size_t a, size_t b) {
        double aLength = rows[b][a] * 0.5;
        double bLength = aLength;
        size_t aCount  = clusters[rowToCluster[a]].countOfExteriorNodes;
        size_t bCount  = clusters[rowToCluster[b]].countOfExteriorNodes;
        size_t tCount  = aCount + bCount;
        double lambda  = (double)aCount / (double)tCount;
        double mu      = 1.0 - lambda;
        for (size_t i=0; i<row_count; ++i) {
            if (i!=a && i!=b) {
                T Dai      = rows[a][i];
                T Dbi      = rows[b][i];
                T Dci      = lambda * Dai + mu * Dbi;
                rows[a][i] = Dci;
                rows[i][a] = Dci;
            }
        }
        clusters.addCluster ( rowToCluster[a], aLength,
                              rowToCluster[b], bLength);
        rowToCluster[a] = clusters.size()-1;
        rowToCluster[b] = rowToCluster[row_count-1];
        removeRowAndColumn(b);
    }
    size_t getImbalance(size_t rowA, size_t rowB) const {
        size_t clusterA = rowToCluster[rowA];
        size_t clusterB = rowToCluster[rowB];
        size_t sizeA = clusters[clusterA].countOfExteriorNodes;
        size_t sizeB = clusters[clusterB].countOfExteriorNodes;
        return (sizeA<sizeB) ? (sizeB-sizeA) : (sizeA-sizeB);
    }
};

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
    using super::recalculateTotalForOneRow;
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
        //
        //Note: Rather than multiplying distances by (n-2)
        //      repeatedly, it is cheaper to work with row
        //      totals multiplied by (1/(T)(n-2)).
        //      Better n multiplications than n*(n-1)/2.
        //
        T nless2      = ( row_count - 2 );
        T tMultiplier = ( row_count <= 2 ) ? 0 : (1 / nless2);
        calculateScaledRowTotals();
        T* tot = scaledRowTotals.data();
        for (size_t r=0; r<row_count; ++r) {
            tot[r] = rowTotals[r] * tMultiplier;
        }
        rowMinima.resize(row_count);
        rowMinima[0].value = infiniteDistance;
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
        T medianLength  = 0.5 * rows[a][b];
        T fudge         = (rowTotals[a] - rowTotals[b]) * tMultiplier;
        T aLength       = medianLength + fudge;
        T bLength       = medianLength - fudge;
        T lambda        = 0.5;
        T mu            = 1.0 - lambda;
        T dCorrection   = - lambda * aLength - mu * bLength;
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (size_t i=0; i<row_count; ++i) {
            if (i!=a && i!=b) {
                T Dai   = rows[a][i];
                T Dbi   = rows[b][i];
                T Dci   = lambda * Dai + mu * Dbi + dCorrection;
                rows[a][i]    = Dci;
                rows[i][a]    = Dci;
                rowTotals[i] += Dci - Dai - Dbi;
                                //JB2020-06-18 Adjust row totals on fly
            }
        }
        recalculateTotalForOneRow(a,b);
        rowTotals[a] -= rows[a][b];
        clusters.addCluster ( rowToCluster[a], aLength,
                              rowToCluster[b], bLength);
        rowToCluster[a] = clusters.size()-1;
        rowToCluster[b] = rowToCluster[row_count-1];
        removeRowAndColumn(b);
    }
    virtual void finishClustering() {
        //Assumes that n is 3
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

template <class T=NJFloat> class BIONJMatrix : public NJMatrix<T> {
public:
    typedef NJMatrix<T> super;
    using super::clusters;
    using super::row_count;
    using super::column_count;
    using super::rows;
    using super::rowToCluster;
    using super::rowTotals;
    using super::recalculateTotalForOneRow;
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
        for (size_t i=0; i<a; ++i) {
            lambda += variance.rows[b][i] - variance.rows[a][i];
        }
        for (size_t i=a+1; i<b; ++i) {
            lambda += variance.rows[b][i] - variance.rows[a][i];
        }
        for (size_t i=b+1; i<column_count; ++i) {
            lambda += variance.rows[b][i] - variance.rows[a][i];
        }
        lambda = 0.5 + lambda / (2.0*((T)row_count-2)*Vab);
        if (1.0<lambda) lambda=1.0;
        if (lambda<0.0) lambda=0.0;
        return lambda;
    }
    virtual void cluster(size_t a, size_t b) {
        //Assumed 0<=a<b<n
        //Bits that differ from super::cluster tagged BIO
        T nless2        = row_count - 2 ;
        T tMultiplier   = ( row_count < 3 ) ? 0 : ( 0.5 / nless2 );
        T medianLength  = 0.5 * rows[b][a];
        T fudge         = (rowTotals[a] - rowTotals[b]) * tMultiplier;
        T aLength       = medianLength + fudge;
        T bLength       = medianLength - fudge;
        T Vab           = variance.rows[b][a];     //BIO
        T lambda        = chooseLambda(a, b, Vab); //BIO
        T mu            = 1.0 - lambda;
        T dCorrection   = - lambda * aLength - mu * bLength;
        T vCorrection   = - lambda * mu * Vab;
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (size_t i=0; i<row_count; ++i) {
            if (i!=a && i!=b) {
                //Dci as per reduction 4 in [Gascuel]
                T Dai         = rows[a][i];
                T Dbi         = rows[b][i];
                T Dci         = lambda * Dai + mu * Dbi + dCorrection;
                rows[a][i]    = Dci;
                rows[i][a]    = Dci;
                rowTotals[i] += Dci - Dai - Dbi; //JB2020-06-18 Adjust row totals
                
                //BIO begin (Reduction 10 on variance estimates)
                T Vci   = lambda * variance.rows[a][i]
                        + mu * variance.rows[b][i]
                        + vCorrection;
                variance.rows[a][i] = Vci;
                variance.rows[i][a] = Vci;
                //BIO finish
            }
        }
        recalculateTotalForOneRow(a,b);
        clusters.addCluster ( rowToCluster[a], aLength, rowToCluster[b], bLength);
        rowToCluster[a] = clusters.size()-1;
        rowToCluster[b] = rowToCluster[row_count-1];
        removeRowAndColumn(b);
        variance.removeRowAndColumn(b); //BIO
    }
};

template <class T=NJFloat, class super=BIONJMatrix<T>>
class BoundingMatrix: public super
{
    using super::rows;
    using super::row_count;
    using super::rowMinima;
    using super::rowTotals;
    using super::rowToCluster;
    using super::clusters;
    using super::silent;
    using super::finishClustering;
protected:
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
    void purgeRow(size_t r /*row index*/) {
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
                    size_t rowA     = rowMinima[r].row;
                    size_t rowB     = rowMinima[r].column;
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
        
        //The weird-looking middle term in the for loop is
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
            T      qBestForThread  = qBest;
            size_t row             = rowScanOrder[r];
            size_t cluster         = rowToCluster[row];
            T      maxEarlierTotal = scaledMaxEarlierClusterTotal[cluster];
            //Note: Older versions of RapidNJ used maxTot rather than
            //      maxEarlierTotal here...
            rowMinima[r]           = getRowMinimum(row, maxEarlierTotal, qBestForThread);
            T      qBestInRow      = rowMinima[row].value;
            if ( qBestInRow < qBestForThread ) {
                qBestForThread = qBestInRow;
            }
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
            const T* rowData = rows[row];
            size_t col;
            V minVector = infiniteDistance;
                //The minima of columns with indices
                //"congruent modulo blockSize"
                //For example, if blockSize is 4,
                //minVector[1] holds the minimum of
                //columns 1,5,9,13,17,...
            V ixVector = -1;
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
            pos.value -= tot [row];
            pos.imbalance = getImbalance(pos.row, pos.column);
            rowMinima[row] = pos;
        }
    }
};//end of class

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
    const size_t blockSize;
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
        for (size_t row=1; row<row_count; ++row) {
            Position<T> pos(row, 0, infiniteDistance, 0);
            const T* rowData = rows[row];
            size_t col;
            V minVector  = infiniteDistance;
            V ixVector   = -1 ;

            for (col=0; col+blockSize<row; col+=blockSize) {
                V  rowVector; rowVector.load_a(rowData+col);
                VB less      = rowVector < minVector;
                V  numVector; numVector.load_a(nums+col);
                ixVector  = select(less, numVector, ixVector);
                minVector = select(less, rowVector, minVector);
            }
            //Extract minimum and column number
            for (int c=0; c<blockSize; ++c) {
                if (minVector[c] < pos.value) {
                    pos.value  = minVector[c];
                    pos.column = ixVector[c];
                }
            }
            for (; col<row; ++col) {
                T dist = rowData[col];
                if (dist < pos.value) {
                    pos.column = col;
                    pos.value  = dist;
                }
            }
            pos.imbalance = getImbalance(pos.row, pos.column);
            rowMinima[row] = pos;
        }
    }
};

typedef BoundingMatrix<NJFloat, NJMatrix<NJFloat>>      RapidNJ;
typedef BoundingMatrix<NJFloat, BIONJMatrix<NJFloat>>   RapidBIONJ;
typedef VectorizedMatrix<NJFloat, NJMatrix<NJFloat>>    VectorNJ;
typedef VectorizedMatrix<NJFloat, BIONJMatrix<NJFloat>> VectorBIONJ;

void addBioNJ2020TreeBuilders(Factory& f) {
    f.advertiseTreeBuilder( new Builder<NJMatrix<NJFloat>>    ("NJ",      "Neighbour Joining (Saitou, Nei [1987])"));
    f.advertiseTreeBuilder( new Builder<RapidNJ>              ("NJ-R",    "Rapid Neighbour Joining (Simonsen, Mailund, Pedersen [2011])"));
    f.advertiseTreeBuilder( new Builder<VectorNJ>             ("NJ-V",    "Vectorized Neighbour Joining (Saitou, Nei [1987])"));
    f.advertiseTreeBuilder( new Builder<BIONJMatrix<NJFloat>> ("BIONJ",   "BIONJ (Gascuel, Cong [2009])"));
    f.advertiseTreeBuilder( new Builder<RapidBIONJ>  ("BIONJ-R", "Rapid BIONJ (Saitou, Nei [1987], Gascuel [2009], Simonson Mailund Pedersen [2011])"));
    f.advertiseTreeBuilder( new Builder<VectorBIONJ> ("BIONJ-V", "Vectorized BIONJ (Gascuel, Cong [2009])"));
    f.advertiseTreeBuilder( new Builder<UPGMA_Matrix<NJFloat>>("UPGMA",    "UPGMA (Sokal, Michener [1958])"));
    f.advertiseTreeBuilder( new Builder<VectorizedUPGMA_Matrix<NJFloat>>("UPGMA-V", "Vectorized UPGMA (Sokal, Michener [1958])"));
    f.advertiseTreeBuilder( new Builder<BoundingMatrix<double>> ("NJ-R-D", "Double precision Rapid Neighbour Joining"));
    const char* defaultName = "RapidNJ";
    f.advertiseTreeBuilder( new Builder<RapidNJ>                (defaultName, "Rapid Neighbour Joining (Simonsen, Mailund, Pedersen [2011]) (default)"));  //Default.
    f.setNameOfDefaultTreeBuilder(defaultName);
}
}; //end of namespace
