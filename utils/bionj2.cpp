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
//

#include "starttree.h"
#include "heapsort.h"                //for mirroredHeapsort, used to sort
                                     //rows of the S and I matrices
                                     //See [SMP2011], section 2.5.
#include "gzstream.h"                //for igzstream
#include <vector>                    //for std::vector
#include <string>                    //sequence names stored as std::string
#include <fstream>
#include <iostream>                  //for std::istream
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

template <class T=NJFloat> struct Link {
    //
    //Describes a link between an interior node and
    //a cluster (clusters are identified by index).
    //
public:
    size_t  clusterIndex;
    T       linkDistance;
    Link(size_t index, T distance) {
        clusterIndex = index;
        linkDistance = distance;
    }
};

template <class T=NJFloat> struct Cluster
{
    //
    //Describes a cluster (either a single exterior
    //node, with no links out from it), or an inerior
    //node, with links to clusters that were formed
    //earlier.
    //
public:
    size_t countOfExteriorNodes;
    std::string name;
    std::vector<Link<T>> links;
    Cluster(): countOfExteriorNodes(0) {
    }
    explicit Cluster(const std::string &taxon_name) {
        countOfExteriorNodes = 1;
        name = taxon_name;
    }
};

template <class T> class ClusterTree: public std::vector<Cluster<T>>
{
public:
    typedef std::vector<Cluster<T>> super;
    using super::at;
    using super::back;
    using super::emplace_back;
    using super::push_back;
    using super::size;
    Cluster<T>& addCluster(const std::string& taxon_name) {
        emplace_back(taxon_name);
        return back();
    }
    Cluster<T>& addCluster(size_t a, T aLength, size_t b, T bLength) {
        push_back(Cluster<T>());
        Cluster<T>& cluster = back();
        cluster.links.emplace_back(a, aLength);
        cluster.links.emplace_back(b, bLength);
        cluster.countOfExteriorNodes = at(a).countOfExteriorNodes + at(b).countOfExteriorNodes;
        return cluster;
    }
    Cluster<T>& addCluster
    ( size_t a, T aLength, size_t b, T bLength
     , size_t c, T cLength)  {
        Cluster<T>& cluster = addCluster(a, aLength, b, bLength);
        cluster.links.emplace_back(c, cLength);
        cluster.countOfExteriorNodes += at(c).countOfExteriorNodes;
        return cluster;
    }
    template <class F> bool writeTreeToFile(const std::string &treeFilePath, F& out) const {
        struct Place
        {
            //
            //Used for keep of tracking where we're up to when
            //we are writing out the description of a Cluster.
            //
        public:
            size_t clusterIndex;
            size_t linkNumber;
            Place(size_t ix, size_t num) {
                clusterIndex = ix;
                linkNumber = num;
            }
        };
        
        out.exceptions(std::ios::failbit | std::ios::badbit);
        try {
            out.open(treeFilePath.c_str(), std::ios_base::out);
            out.precision(8);
            
            std::vector<Place> stack;
            bool failed = false; //Becomes true if clusters
            //defines cycles (should never happen)
            //Indicates a fatal logic error
            size_t maxLoop = 3 * size();
            //More than this, and there must be
            //a cycle.  Or something.
            
            stack.emplace_back(size()-1, 0); //assumes: size is at least 1!
            do {
                --maxLoop;
                if (maxLoop==0) {
                    failed = true;
                    break;
                }
                Place here = stack.back();
                const Cluster<T>& cluster = at(here.clusterIndex);
                stack.pop_back();
                if (cluster.links.empty()) {
                    out << cluster.name;
                    continue;
                }
                if (here.linkNumber==0) {
                    out << "(";
                    stack.emplace_back(here.clusterIndex, 1);
                    stack.emplace_back(cluster.links[0].clusterIndex, 0);
                    continue;
                }
                size_t nextChildNum = here.linkNumber;
                const Link<T> & linkPrev = cluster.links[nextChildNum-1];
                out << ":" << linkPrev.linkDistance;
                if (nextChildNum<cluster.links.size()) {
                    out << ",";
                    const Link<T> & linkNext = cluster.links[nextChildNum];
                    stack.emplace_back(here.clusterIndex, nextChildNum+1);
                    stack.emplace_back(linkNext.clusterIndex, 0);
                } else {
                    out << ")";
                }
            } while (0 < stack.size());
            out << ";" << std::endl;
            out.close();
            return true;
        } catch (std::ios::failure &) {
            std::cerr << "IO error"
            << " opening/writing file: " << treeFilePath << std::endl;
            return false;
        } catch (const char *str) {
            std::cerr << "Writing newick file failed: " << str << std::endl;
            return false;
        } catch (std::string &str) {
            std::cerr << "Writing newick file failed: " << str << std::endl;
            return false;
        }
    }
    bool writeTreeFile(bool zipIt, const std::string &treeFilePath) const {
        if (zipIt) {
            ogzstream out;
            return writeTreeToFile(treeFilePath, out);
        } else {
            std::fstream out;
            return writeTreeToFile(treeFilePath, out);
        }
    }
};

#define MATRIX_ALIGNMENT 64
    //MUST be a power of 2 (else x & MATRIX_ALIGNMENT_MASK
    //would be no good and x % MATRIX_ALIGNMENT would be needed).
    //Assumed: sizeof(NJFloat) divides MATRIX_ALIGNMENT
    //Vectorized versions run faster (particularly on older
    //hardware), if rows are aligned.
#define MATRIX_ALIGNMENT_MASK (MATRIX_ALIGNMENT - 1)

template <class P> inline P* matrixAlign(P* p) {
    //If we've got an array that mighnt't be MATRIX_ALIGNMENT-byte
    //aligned, but we've got MATRIX_ALIGNMENT/sizeof(P) extra items
    //in it, we can point to the first item in the array that *is*
    //MATRIX_ALIGNMENT-byte aligned.  This function returns the
    //address of that item.
    uintptr_t address = reinterpret_cast<uintptr_t>(p);
    auto offset = address & MATRIX_ALIGNMENT_MASK;
    if (0<offset)
    {
        return p + (MATRIX_ALIGNMENT - offset)/sizeof(P);
    } else {
        return p;
    }
}

template <class T=NJFloat> class Matrix
{
    //Note 1: This is a separate class so that it can be
    //        used for square variance (V) and rectangular
    //        sorted distance (S) and index (I) matrices,
    //        not just square distance (D) matrices.
    //        Lines that access the upper-right triangle
    //        of the matrix are tagged with U-R.
    //Note 2: I resorted to declaring the data, rows, and
    //        rowTotals members public, because of problems
    //        I had accessing them from BoundingMatrix.
    //Note 3: Perhaps there should be separate SquareMatrix
    //        and RectangularMatrix classes?
public:
    size_t n;
    size_t shrink_n; //if n reaches *this*, pack the array
    T*     data;
    T**    rows;
    T*     rowTotals; //The U vector
    virtual void setSize(size_t rank) {
        clear();
        if (0==rank) {
            return;
        }
        try {
            size_t w    = widthNeededFor(rank);
            n           = rank;
            shrink_n    = (rank+rank)/3;
            if (shrink_n<100) {
                shrink_n=0;
            }
            data        = new T[n*w + MATRIX_ALIGNMENT/sizeof(T)];
            rows        = new T*[n];
            rowTotals   = new T[n];
            T *rowStart = matrixAlign(data);
            for (size_t r=0; r<n; ++r) {
                rows[r]      = rowStart;
                rowStart    += w;
                rowTotals[r] = 0.0;
            }
            #pragma omp parallel for
            for (size_t r=0; r<n; ++r) {
                zeroRow(r);
            }
        }
        catch (...) {
            clear();
            throw;
        }
    }
    size_t widthNeededFor(size_t width) {
        //
        //returns width, rounded up so that each row will
        //have a starting address that is MATRIX_ALIGNMENT-byte aligned.
        //
        if (MATRIX_ALIGNMENT<=sizeof(T))
        {
            return width;
        }
        size_t leftOver  = (width * sizeof(T)) & MATRIX_ALIGNMENT_MASK;
        if (leftOver==0) {
            return width;
        }
        return width + (MATRIX_ALIGNMENT-leftOver) / sizeof(T);
    }
    void zeroRow(size_t r) {
        T* rowStart = rows[r];
        T* rowStop  = rowStart + n;
        for (T* rowZap=rowStart; rowZap<rowStop; ++rowZap) {
            *rowZap = 0;
        }
    }
    void assign(const Matrix& rhs) {
        if (this==&rhs) {
            return;
        }
        setSize(rhs.n);
        #pragma omp parallel for
        for (size_t r=0; r<n; ++r) {
            T *             destRow = rows[r];
            T const *       sourceRow = rhs.rows[r];
            T const * const endSourceRow = sourceRow + n;
            for (; sourceRow<endSourceRow; ++destRow, ++sourceRow) {
                *destRow = *sourceRow;
            }
            rowTotals[r] = rhs.rowTotals[r];
        }
    }
    Matrix(): n(0), shrink_n(0), data(nullptr), rows(nullptr), rowTotals(nullptr) {
    }
    Matrix(const Matrix& rhs): data(nullptr), rows(nullptr), rowTotals(nullptr) {
        assign(rhs);
    }
    virtual ~Matrix() {
        clear();
    }
    void clear() {
        n = 0;
        delete [] data;
        delete [] rows;
        delete [] rowTotals;
        data = nullptr;
        rows = nullptr;
        rowTotals = nullptr;
    }
    Matrix& operator=(const Matrix& rhs) {
        assign(rhs);
        return *this;
    }
    virtual void calculateRowTotals() const {
        //Note: Although this isn't currently in use,
        //it's been kept, in case it is needed
        //(after, say, every 200 iterations of
        //neighbour-joining) to deal with accumulated
        //rounding error.  It might be.
        #pragma omp parallel for
        for (size_t r=0; r<n; ++r) {
            T total = 0;
            const T* rowData = rows[r];
            for (size_t c=0; c<r; ++c) {
                total += rowData[c];
            }
            for (size_t c=r+1; c<n; ++c) {
                total += rowData[c]; //U-R
            }
            rowTotals[r] = total;
        }
    }
    void recalculateTotalForOneRow(size_t a, size_t b) {
        //recalculate total for row, a, excluding
        //column b (a<=b).
        T replacementRowTotal = 0;
        for (size_t i=0; i<a; ++i) {
            replacementRowTotal += rows[a][i];
        }
        for (size_t i=a+1; i<b; ++i) {
            replacementRowTotal += rows[a][i];
        }
        for (size_t i=b+1; i<n; ++i) {
            replacementRowTotal += rows[a][i];
        }
        rowTotals[a] = replacementRowTotal;
    }
    void removeRowAndColumn(size_t rowNum)  {
        //Remove row (and matching column) from a
        //square matrix, by swapping the last row
        //(and column) into its place.
        #pragma omp parallel for
        for (size_t r=0; r<n; ++r) {
            if (r!=rowNum) {
              T* rowData = rows[r];
              rowData[rowNum] = rowData[n-1]; //U-R
            }
        }
        --n;
        rowTotals[rowNum] = rowTotals[n];
        //was rows[rowNum] = rows[n];... but let's copy
        //instead.  On average it seems (very slightly) faster.
        T*       destRow   = rows[rowNum];
        const T* sourceRow = rows[n];
        rows[n] = nullptr;
        if (destRow!=sourceRow) {
            #pragma omp parallel for
            for (size_t c=0; c<n; ++c) {
                destRow[c] = sourceRow[c];
            }
        }
        if ( n == shrink_n && 0 < shrink_n) {
            //Move the data in the array closer to the front.
            //This also helps (but: only very slightly. 5%ish?).
            size_t   w = widthNeededFor(n);
            T* destRow = data;
            for (size_t r=1; r<n; ++r) {
                destRow += w;
                const T* sourceRow = rows[r];
                #pragma omp parallel for
                for (size_t c=0; c<n; ++c) {
                    destRow[c] = sourceRow[c];
                }
                rows[r] = destRow;
            }
            shrink_n    = (n+n)/3;
            if (shrink_n<100) shrink_n=0;
        }
    }
    void removeRowOnly(size_t rowNum) {
        //Remove row from a rectangular matrix.
        //Don't touch the columns in the row
        //(Used for the S and I matrices in BIONJ).
        rowTotals[rowNum] = rowTotals[n-1];
        rows[rowNum]      = rows[n-1];
        rows[n-1]         = nullptr;
        --n;
    }
};

template <class T=NJFloat> class UPGMA_Matrix: public Matrix<T> {
    //UPGMA_Matrix is a D matrix (a matrix of distances).
public:
    typedef Matrix<T> super;
    using super::rows;
    using super::setSize;
    using super::n;
    using super::calculateRowTotals;
    using super::removeRowAndColumn;
protected:
    std::vector<size_t>  rowToCluster; //*not* initialized by setSize
    ClusterTree<T>       clusters;     //*not* touched by setSize
    mutable Positions<T> rowMinima;    //*not* touched by setSize
    bool isOutputToBeZipped;
public:
    UPGMA_Matrix():super(), isOutputToBeZipped(false) {
    }
    virtual std::string getAlgorithmName() const {
        return "UPGMA";
    }
    bool loadMatrixFromFile(const std::string &distanceMatrixFilePath) {
        size_t rank;
        igzstream in;
        try {
            in.exceptions(std::ios::failbit | std::ios::badbit);
            in.open(distanceMatrixFilePath.c_str(), std::ios_base::in);
            in >> rank;
            setSize(rank);
            progress_display progress(rank, "Loading distance matrix", "loaded", "row");
            for (size_t r=0; r<n; ++r) {
                std::string name;
                in >> name;
                clusters.addCluster(name);
                for (size_t c=0; c<n; ++c) {
                    in >> rows[r][c];
                    //Ensure matrix is symmetric (as it is read!)
                    if (c<r && rows[r][c] != rows[c][r]) {
                        T v = ( rows[r][c] + rows[c][r] ) * 0.5;
                        rows[c][r] = v; //U-R
                        rows[r][c] = v;
                    }
                }
                rowToCluster.emplace_back(r);
                ++progress;
            }
            in.close();
        } catch (std::ios::failure &) {
            std::cerr << "Load matrix failed: IO error"
                << " reading file: " << distanceMatrixFilePath << std::endl;
            return false;
        } catch (const char *str) {
            std::cerr << "Load matrix failed: " << str << std::endl;
            return false;
        } catch (std::string &str) {
            std::cerr << "Load matrix failed: " << str << std::endl;
            return false;
        }
        calculateRowTotals();
        //Note: The old code wrote a message to standard output,
        //      if the matrix was not symmetric.  This code doesn't.
        return true;
    }
    virtual bool loadMatrix(const std::vector<std::string>& names, double* matrix) {
        //Assumptions: 2 < names.size(), all names distinct
        //  matrix is symmetric, with matrix[row*names.size()+col]
        //  containing the distance between taxon row and taxon col.
        setSize(names.size());
        clusters.clear();
        for (auto it = names.begin(); it != names.end(); ++it) {
            clusters.addCluster(*it);
        }
        rowToCluster.resize(n, 0);
        for (size_t r=0; r<n; ++r) {
            rowToCluster[r]=r;
        }
        #pragma omp parallel for
        for (size_t row=0; row<n; ++row) {
            double* sourceStart = matrix + row * n;
            double* sourceStop  = sourceStart + n;
            T*      dest        = rows[row];
            for (double* source=sourceStart; source<sourceStop; ++source, ++dest ) {
                *dest = (T) *source;
            }
        }
        calculateRowTotals();
        return true;
    }
    virtual bool constructTree() {
        Position<T> best;
        std::string taskName = "Constructing " + getAlgorithmName() + " tree";
        progress_display show_progress(n*(n+1)/2, taskName.c_str(), "", "");
        while (3<n) {
            getMinimumEntry(best);
            cluster(best.column, best.row);
            show_progress+=n;
        }
        finishClustering();
        show_progress.done();
        return true;
    }
    virtual void setZippedOutput(bool zipIt) {
        isOutputToBeZipped = zipIt;
    }
    bool writeTreeFile(const std::string &treeFilePath) const {
        return clusters.writeTreeFile(isOutputToBeZipped, treeFilePath);
    }
protected:
    virtual void setSize(size_t rank) {
        super::setSize(rank);
        rowToCluster.clear();
    }
    void getMinimumEntry(Position<T> &best) {
        getRowMinima();
        best.value = infiniteDistance;
        for (size_t r=0; r<n; ++r) {
            Position<T> & here = rowMinima[r];
            if (here.value < best.value && here.row != here.column) {
                best = here;
            }
        }
    }
    virtual void getRowMinima() const
    {
        rowMinima.resize(n);
        rowMinima[0].value = infiniteDistance;
        #pragma omp parallel for schedule(dynamic)
        for (size_t row=1; row<n; ++row) {
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
    void finishClustering() {
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
        n = 0;
    }
    virtual void cluster(size_t a, size_t b) {
        double aLength = rows[b][a] * 0.5;
        double bLength = aLength;
        size_t aCount  = clusters[rowToCluster[a]].countOfExteriorNodes;
        size_t bCount  = clusters[rowToCluster[b]].countOfExteriorNodes;
        size_t tCount  = aCount + bCount;
        double lambda  = (double)aCount / (double)tCount;
        double mu      = 1.0 - lambda;
        for (size_t i=0; i<n; ++i) {
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
        rowToCluster[b] = rowToCluster[n-1];
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
    using super::n;
    using super::rows;
    using super::rowMinima;
    using super::rowTotals;
    using super::rowToCluster;
    using super::removeRowAndColumn;
    using super::calculateRowTotals;
    using super::recalculateTotalForOneRow;
    using super::getImbalance;
protected:
    mutable std::vector<T> scaledRowTotals; //used in getRowMinima
public:
    NJMatrix(): super() { }
    virtual std::string getAlgorithmName() const {
        return "NJ";
    }
protected:
    virtual void calculateScaledRowTotals() const {
        scaledRowTotals.resize(n);
        T nless2      = ( n - 2 );
        T tMultiplier = ( n <= 2 ) ? 0 : (1 / nless2);
        #pragma omp parallel for
        for (size_t r=0; r<n; ++r) {
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
        T nless2      = ( n - 2 );
        T tMultiplier = ( n <= 2 ) ? 0 : (1 / nless2);
        calculateScaledRowTotals();
        T* tot = scaledRowTotals.data();
        for (size_t r=0; r<n; ++r) {
            tot[r] = rowTotals[r] * tMultiplier;
        }
        rowMinima.resize(n);
        rowMinima[0].value = infiniteDistance;
        #pragma omp parallel for schedule(dynamic)
        for (size_t row=1; row<n; ++row) {
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
        T nless2        = n-2;
        T tMultiplier   = (n<3) ? 0 : (0.5 / nless2);
        T medianLength  = 0.5 * rows[a][b];
        T fudge         = (rowTotals[a] - rowTotals[b]) * tMultiplier;
        T aLength       = medianLength + fudge;
        T bLength       = medianLength - fudge;
        T lambda        = 0.5;
        T mu            = 1.0 - lambda;
        T dCorrection   = - lambda * aLength - mu * bLength;
        #pragma omp parallel for
        for (size_t i=0; i<n; ++i) {
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
        rowToCluster[b] = rowToCluster[n-1];
        removeRowAndColumn(b);
    }
    void finishClustering() {
        //Assumes that n is 3
        T halfD01 = 0.5 * rows[0][1];
        T halfD02 = 0.5 * rows[0][2];
        T halfD12 = 0.5 * rows[1][2];
        clusters.addCluster
            ( rowToCluster[0], halfD01 + halfD02 - halfD12
            , rowToCluster[1], halfD01 + halfD12 - halfD02
            , rowToCluster[2], halfD02 + halfD12 - halfD01);
        n = 0;
    }
};

template <class T=NJFloat> class BIONJMatrix : public NJMatrix<T> {
public:
    typedef NJMatrix<T> super;
    using super::clusters;
    using super::n;
    using super::rows;
    using super::rowToCluster;
    using super::rowTotals;
    using super::recalculateTotalForOneRow;
    using super::removeRowAndColumn;
protected:
    Matrix<T>  variance;       //The V matrix
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
    virtual bool loadMatrix(const std::vector<std::string>& names, double* matrix) {
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
        for (size_t i=b+1; i<n; ++i) {
            lambda += variance.rows[b][i] - variance.rows[a][i];
        }
        lambda = 0.5 + lambda / (2.0*((T)n-2)*Vab);
        if (1.0<lambda) lambda=1.0;
        if (lambda<0.0) lambda=0.0;
        return lambda;
    }
    virtual void cluster(size_t a, size_t b) {
        //Assumed 0<=a<b<n
        //Bits that differ from super::cluster tagged BIO
        T nless2        = n - 2 ;
        T tMultiplier   = ( n < 3 ) ? 0 : ( 0.5 / nless2 );
        T medianLength  = 0.5 * rows[b][a];
        T fudge         = (rowTotals[a] - rowTotals[b]) * tMultiplier;
        T aLength       = medianLength + fudge;
        T bLength       = medianLength - fudge;
        T Vab           = variance.rows[b][a];     //BIO
        T lambda        = chooseLambda(a, b, Vab); //BIO
        T mu            = 1.0 - lambda;
        T dCorrection   = - lambda * aLength - mu * bLength;
        T vCorrection   = - lambda * mu * Vab;
        #pragma omp parallel for
        for (size_t i=0; i<n; ++i) {
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
        rowToCluster[b] = rowToCluster[n-1];
        removeRowAndColumn(b);
        variance.removeRowAndColumn(b); //BIO
    }
};

template <class T=NJFloat, class super=BIONJMatrix<T>>
class BoundingMatrix: public super
{
    using super::n;
    using super::rows;
    using super::rowMinima;
    using super::rowTotals;
    using super::rowToCluster;
    using super::clusters;
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
    //
    std::vector<int> clusterToRow;   //Maps clusters to their rows (-1 means not mapped)
    std::vector<T>   clusterTotals;  //"Row" totals indexed by cluster

    mutable std::vector<T>       scaledClusterTotals;   //The same, multiplied by
                                                        //(1.0 / (n-2)).
    mutable std::vector<T>       scaledMaxEarlierClusterTotal;
    mutable std::vector<bool>    rowOrderChosen; //Indicates if row order chosen
    mutable std::vector<size_t>  rowScanOrder;   //Order in which rows are to be scanned
                                                 //Only used in... getRowMinima().
    
    Matrix<T>   entriesSorted; //The S matrix: Entries in distance matrix
                               //(each row sorted by ascending value)
    Matrix<int> entryToCluster;//The I matrix: for each entry in S, which
                               //cluster the row (that the entry came from)
                               //was mapped to (at the time).
    double rowSortingTime;
    
public:
    BoundingMatrix() : super() {
        rowSortingTime = 0;
    }
    virtual std::string getAlgorithmName() const {
        return "Rapid" + super::getAlgorithmName();
    }
    virtual bool constructTree() {
        //1. Set up vectors indexed by cluster number,
        clusterToRow.resize(n);
        clusterTotals.resize(n);
        for (size_t r=0; r<n; ++r) {
            clusterToRow[r]  = static_cast<int>(r);
            clusterTotals[r] = rowTotals[r];
        }
        
        //2. Set up "scratch" vectors used in getRowMinima
        //   so that it won't be necessary to reallocate them
        //   for each call.
        scaledClusterTotals.resize(n);
        scaledMaxEarlierClusterTotal.resize(n);
        rowOrderChosen.resize(n);
        rowScanOrder.resize(n);

        {
            progress_display setupProgress(n, "Setting up auxiliary I and S matrices", "sorting", "row");
            //2. Set up the matrix with row sorted by distance
            //   And the matrix that tracks which distance is
            //   to which cluster (the S and I matrices, in the
            //   RapidNJ papers).
            entriesSorted.setSize(n);
            entryToCluster.setSize(n);
#pragma omp parallel for schedule(dynamic)
            for (size_t r=0; r<n; ++r) {
                sortRow(r,r);
                ++setupProgress;
                //copies the "left of the diagonal" portion of
                //row r from the D matrix and sorts it
                //into ascending order.
            }
        }
        {
            size_t nextPurge = (n+n)/2;
            std::string taskName = "Constructing " + getAlgorithmName() + " tree";
            progress_display show_progress(n*(n+1)/2, taskName.c_str(), "", "");
            while (3<n) {
                Position<T> best;
                super::getMinimumEntry(best);
                cluster(best.column, best.row);
                if ( n == nextPurge ) {
                    #pragma omp parallel for
                    for (size_t r=0; r<n; ++r) {
                        purgeRow(r);
                    }
                    nextPurge = n*2/3;
                }
                show_progress+=n;
            }
            show_progress.done();
            super::finishClustering();
        }
        return true;
    }
    void sortRow(size_t r /*row index*/, size_t c /*upper bound on cluster index*/) {
        //1. copy data from a row of the D matrix into the S matrix
        //   (and write the cluster identifiers that correspond to
        //    the values in the D row into the same-numbered
        //    row in the I matrix), for distances between the cluster
        //    in that row, and other live clusters (up to, but not including c).
        T*     sourceRow      = rows[r];
        T*     values         = entriesSorted.rows[r];
        int*   clusterIndices = entryToCluster.rows[r];
        size_t w = 0;
        for (size_t i=0; i<n; ++i) {
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
        if ( n<=c) {
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
        for (; i<n ; ++i ) {
            values[w]         = values[i];
            clusterIndices[w] = clusterIndices[i];
            if ( infiniteDistance <= values[i] ) {
                break;
            }
            if ( clusterToRow[clusterIndices[i]] != notMappedToRow ) {
                ++w;
            }
        }
        if (w<n) {
            values[w] = infiniteDistance;
        }
    }
    virtual void cluster(size_t a, size_t b) {
        size_t clusterA         = rowToCluster[a];
        size_t clusterB         = rowToCluster[b];
        size_t clusterMoved     = rowToCluster[n-1];
        clusterToRow[clusterA]  = notMappedToRow;
        clusterTotals[clusterA] = -infiniteDistance;
        clusterToRow[clusterB]  = notMappedToRow;
        clusterTotals[clusterB] = -infiniteDistance;
        size_t clusterC = clusters.size(); //cluster # of new cluster
        super::cluster(a,b);
        if (b<n) {
            clusterToRow[clusterMoved] = static_cast<int>(b);
        }
        clusterToRow.emplace_back(a);
        clusterTotals.emplace_back(rowTotals[a]);
        scaledClusterTotals.emplace_back(rowTotals[a] / (T)( n - 1.0 ) );
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
        for (size_t r = 0; r<n; ++r) {
            size_t cluster = rowToCluster[r];
            clusterTotals[cluster] = rowTotals[r];
        }
        sortRow(a, clusterC);
    }
    void decideOnRowScanningOrder() const {
        //
        //Rig the order in which rows are scanned based on
        //which rows (might) have the lowest row minima
        //based on what we saw last time.
        //The original RapidNJ puts the second-best row from last time first.
        //And, apart from that, goes in row order.
        //But rows in the D, S, and I matrices are (all) shuffled
        //in memory, so why not do all the rows in ascending order
        //of their best Q-values from the last iteration?
        //
        //(I have my doubts about this.  I'm now thinking, it should be
        //better to "read off" the first column of the S matrix and
        //calculate the corresponding entry in Q for each live cluster;
        //use up-to-date rather than out-of-date information.
        //
        //Since we always have to check these entries when we process
        //the row, why not process them up front, hoping to get a
        //better bound on min(V) (and "rule out" entire rows with that
        //bound)? -James B).
        //
        
        //Note, rowMinima might have size 0 (the first time this member
        //function is called during processing of a distance matrix)
        //Or it might have a size of n+1 (later times), but it won't be n.
        for ( size_t len = rowMinima.size(); 1<len; len=(len+1)/2 ) {
            size_t halfLen = len/2; //rounded down
            size_t gap     = len-halfLen;
            //Although the following loop could in theory be parallelized,
            //using #pragma omp for on it did not seem to help performance any.
            for ( size_t i=0; i<halfLen; ++i) {
                size_t j = i + gap;
                if ( rowMinima[j] < rowMinima[i] ) {
                    std::swap(rowMinima[i], rowMinima[j]);
                }
            }
        }
        for (size_t i=0; i<n; ++i) {
            rowOrderChosen[i]=false;
        }
        size_t w = 0;
        for (size_t r=0; r < rowMinima.size()
             && rowMinima[r].value < infiniteDistance; ++r) {
            size_t rowA     = rowMinima[r].row;
            size_t rowB     = rowMinima[r].column;
            size_t clusterA = (rowA<n) ? rowToCluster[rowA] : 0;
            size_t clusterB = (rowB<n) ? rowToCluster[rowB] : 0;
            size_t row      = (clusterA<clusterB) ? rowB : rowA;
            rowScanOrder[w] = row;
            w += ( row < n && !rowOrderChosen[row] ) ? 1 : 0;
            rowOrderChosen[row] = true;
        }
        for (size_t r=0; r<n; ++r) {
            rowScanOrder[w] = r;
            w += ( rowOrderChosen[r] ? 0 : 1 );
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
        T      nless2      = ( n - 2 );
        T      tMultiplier = ( n <= 2 ) ? 0 : (1 / nless2);
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

        decideOnRowScanningOrder();
        rowMinima.resize(n);
        #pragma omp parallel for
        for (size_t r=0; r<n; ++r) {
            size_t row             = rowScanOrder[r];
            size_t cluster         = rowToCluster[row];
            T      maxEarlierTotal = scaledMaxEarlierClusterTotal[cluster];
            //Note: Older versions of RapidNJ used maxTot rather than
            //      maxEarlierTotal here...
            rowMinima[r]          = getRowMinimum(row, maxEarlierTotal, qBest);
            T      v              = rowMinima[r].value;
            {
                if ( v < qBest ) {
                    #pragma omp critical(checkmin)
                    if (v < qBest) {
                        qBest = v;
                    }
                }
            }
        }
    }
    Position<T> getRowMinimum(size_t row, T maxTot, T qBest) const {
        T nless2      = ( n - 2 );
        T tMultiplier = ( n <= 2 ) ? 0 : ( 1.0 / nless2 );
        auto    tot   = scaledClusterTotals.data();
        T rowTotal    = rowTotals[row] * tMultiplier; //scaled by (1/(n-2)).
        T rowBound    = qBest + maxTot + rowTotal;
                //Upper bound for distance, in this row, that
                //could (after row totals subtracted) provide a
                //better min(Q).

        Position<T> pos(row, 0, infiniteDistance, 0);
        const T*   rowData   = entriesSorted.rows[row];
        const int* toCluster = entryToCluster.rows[row];
        T Drc;
        for (size_t i=0; (Drc=rowData[i])<rowBound; ++i) {
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
    using super::n;
    using super::rows;
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
        scratchTotals.resize(n + fluff, 0.0);
        scratchColumnNumbers.resize(n + fluff, 0.0);
    }
    virtual void getRowMinima() const {
        T nless2      = ( n - 2 );
        T tMultiplier = ( n <= 2 ) ? 0 : (1 / nless2);
        T* tot  = matrixAlign ( scratchTotals.data() );
        T* nums = matrixAlign ( scratchColumnNumbers.data() );
        for (size_t r=0; r<n; ++r) {
            tot[r]  = rowTotals[r] * tMultiplier;
            nums[r] = r;
        }
        rowMinima.resize(n);
        rowMinima[0].value = infiniteDistance;
        #pragma omp parallel for schedule(dynamic)
        for (size_t row=1; row<n; ++row) {
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
    using super::n;
    using super::rowMinima;
    using super::rows;
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
        scratchColumnNumbers.resize(n + fluff, 0.0);
    }
    virtual void getRowMinima() const
    {
        T* nums = matrixAlign ( scratchColumnNumbers.data() );
        rowMinima.resize(n);
        rowMinima[0].value = infiniteDistance;
        #pragma omp parallel for schedule(dynamic)
        for (size_t row=1; row<n; ++row) {
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
