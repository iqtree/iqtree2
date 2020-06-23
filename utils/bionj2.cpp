//
//  bionj2.cpp - Implementations of NJ and BIONJ algorithms
//               (that work in terms of .mldist inputs and
//                NEWICK outputs).
//
//  BIONJ implementation based on http://www.lirmm.fr/~w3ifa/MAAS/BIONJ/BIONJ.html
//        (see BIONJMatrix).  Original authors: Olivier Gascuel
//        and Hoa Sien Cuong (the code for the Unix version).
//  NJ    implementation based on the same (but original NJ, without
//        a matrix of variance estimates (see NJMatrix).
//  BoundingNJ implementation loosely based on ideas from
//        https://birc.au.dk/software/rapidnj/ and
//        from: Inference of Large Phylogenies using Neighbour-Joining.
//              Martin Simonsen, Thomas Mailund, Christian N. S. Pedersen.
//              Communications in Computer and Information Science
//              (Biomedical Engineering Systems and Technologies:
//              3rd International Joint Conference, BIOSTEC 2010,
//              Revised Selected Papers), volume 127, pages 334-344,
//              Springer Verlag, 2011.
//         Tag: [SMP2011].
//         (but using a variance matrix, as in BIONJ, and keeping the
//          distance and variance matrices square - they're not triangular
//          because (i) *read* memory access patterns are more favourable
//                 (ii) *writes* don't require conditional transposition
//                      of the row and column coordinates (but their
//                      access patterns aren't as favourable, but
//                (iii) reads vastly outnumber writes)
//         (there's no code yet for removing duplicated rows either;
//          those that has distance matrix rows identical to earlier rows;
//          Rapid NJ "hates" them) (this is also covered in section 2.5)
//         See the BoundingBIONJMatrix class.
//
//  Created by James Barbetti on 18/6/2020.
//

#include "bionj2.h"
#include "heapsort.h"             //for mirroredHeapsort
#include <vector>
#include <string>
#include <fstream>
#include <iostream>               //for std::istream
#include <boost/scoped_array.hpp> //for boost::scoped_array
#include "utils/timeutil.h"       //for getRealTime()
#include <vectorclass/vectorclass.h> //for Vec4d

typedef double NJFloat;

const NJFloat infiniteDistance = 1e+300;

namespace
{

struct Position
{
    //A position (row, column) in an NJ matrix
    //Note that column is always less than row.
    //(Because that's the convention in RapidNJ).
public:
    size_t  row;
    size_t  column;
    NJFloat value;
    Position() : row(0), column(0), value(0) {}
    Position(size_t r, size_t c, NJFloat v)
        : row(r), column(c), value(v) {}
    Position& operator = (const Position &rhs) {
        row    = rhs.row;
        column = rhs.column;
        value  = rhs.value;
        return *this;
    }
    bool operator< ( const Position& rhs ) const {
        return value < rhs.value;
    }
    bool operator<= ( const Position& rhs ) const {
        return value <= rhs.value;
    }
};

typedef std::vector<Position> Positions;

struct Link {
    //
    //Describes a link between an interior node and
    //a cluster (clusters are identified by index).
    //
public:
    size_t  clusterIndex;
    NJFloat linkDistance;
    Link(size_t index, NJFloat distance) {
        clusterIndex = index;
        linkDistance = distance;
    }
};

struct Cluster
{
    //
    //Describes a cluster (either a single exterior
    //node, with no links out from it), or an inerior
    //node, with links to clusters that were formed
    //earlier.
    //
public:
    std::string name;
    std::vector<Link> links;
    explicit Cluster(const std::string &taxon_name) {
        name = taxon_name;
    }
    Cluster(size_t a, NJFloat aLength, size_t b, NJFloat bLength) {
        links.emplace_back(a, aLength);
        links.emplace_back(b, bLength);
    }
    Cluster
        ( size_t a, NJFloat aLength, size_t b, NJFloat bLength
        , size_t c, NJFloat cLength) {
        links.emplace_back(a, aLength);
        links.emplace_back(b, bLength);
        links.emplace_back(c, cLength);
    }
};

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
    T*     data;
    T**    rows;
    T*     rowTotals;
    void setSize(size_t rank) {
        clear();
        if (0==rank) {
            return;
        }
        try {
            n         = rank;
            data      = new T[n*n];
            rows      = new T*[n];
            rowTotals = new T[n];
            T *rowStart = data;
            for (int r=0; r<n; ++r) {
                rows[r] = rowStart;
                rowStart += n;
                rowTotals[r] = 0.0;
            }
            #pragma omp parallel for
            for (int r=0; r<n; ++r) {
                zeroRow(r);
            }
        }
        catch (...) {
            clear();
            throw;
        }
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
    Matrix(): n(0), data(0), rows(0), rowTotals(0) {
    }
    Matrix(const Matrix& rhs) {
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
    void calculateRowTotals() const {
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
    void removeRowAndColumn(size_t rowNum)  {
        //Remove row (and matching column) from a
        //square matrix, by swapping the last row
        //(and column) into its place.
        #pragma omp parallel for
        for (size_t r=0; r<n; ++r) {
            T* rowData = rows[r];
            rowData[rowNum] = rowData[n-1]; //U-R
        }
        rowTotals[rowNum] = rowTotals[n-1];
        rows[rowNum] = rows[n-1];
        rows[n-1] = nullptr;
        --n;
    }
    void removeRowOnly(size_t rowNum) {
        //Remove row from a rectangular matrix.
        rowTotals[rowNum] = rowTotals[n-1];
        rows[rowNum]      = rows[n-1];
        rows[n-1]         = nullptr;
        --n;
    }
};

class NJMatrix: public Matrix<NJFloat>
{
protected:
    std::vector<size_t>          rowToCluster;
    std::vector<Cluster>         clusters;
    mutable Positions            rowMinima;
    mutable std::vector<NJFloat> scaledRowTotals; //used in getRowMinima
public:
    explicit NJMatrix(const std::string &distanceMatrixFilePath) {
        size_t rank;
        std::fstream in;
        in.open(distanceMatrixFilePath, std::ios_base::in);
        in >> rank;
        setSize(rank);
        for (int r=0; r<n; ++r) {
            std::string name;
            in >> name;
            clusters.emplace_back(name);
            for (int c=0; c<n; ++c) {
                in >> rows[r][c];
                //Ensure matrix is symmetric (as it is read!)
                if (c<r && rows[r][c]<rows[c][r]) {
                    NJFloat v = ( rows[r][c] + rows[c][r] ) * 0.5;
                    rows[c][r] = v; //U-R
                    rows[r][c] = v;
                }
            }
            rowToCluster.emplace_back(r);
        }
        in.close();
        calculateRowTotals();
        scaledRowTotals.resize(n, 0.0);
        calculateScaledRowTotals();
        //Note: The old code wrote a message to standard output,
        //      if the matrix was not symmetric.  This code doesn't.
    }
    void calculateScaledRowTotals() const {
        NJFloat nless2      = ( n - 2 );
        NJFloat tMultiplier = ( n <= 2 ) ? 0 : (1 / nless2);
        #pragma omp parallel for
        for (size_t r=0; r<n; ++r) {
            scaledRowTotals[r] = rowTotals[r] * tMultiplier;
        }
    }
    virtual void getRowMinima() const {
        //
        //Note: Rather than multiplying distances by (n-2)
        //      repeatedly, it is cheaper to work with row
        //      totals multiplied by (1/(NJFloat)(n-2)).
        //      Better n multiplications than n*(n-1)/2.
        //
        NJFloat nless2      = ( n - 2 );
        NJFloat tMultiplier = ( n <= 2 ) ? 0 : (1 / nless2);
        calculateScaledRowTotals();
        auto tot = scaledRowTotals.data();
        for (size_t r=0; r<n; ++r) {
            tot[r] = rowTotals[r] * tMultiplier;
        }
        rowMinima.resize(n);
        rowMinima[0].value = infiniteDistance;
        #pragma omp parallel for
        for (size_t row=1; row<n; ++row) {
            Position pos(row, 0, infiniteDistance);
            const NJFloat* rowData = rows[row];
            for (size_t col=0; col<row; ++col) {
                NJFloat v = rowData[col] - tot[col];
                if (v < pos.value) {
                    pos.column = col;
                    pos.value  = v;
                }
            }
            pos.value -= tot [row];
            rowMinima[row] = pos;
        }
    }
    void getMinimumEntry(Position &best) {
        getRowMinima();
        best.value = infiniteDistance;
        for (size_t r=0; r<n; ++r) {
            if (rowMinima[r].value < best.value) {
                best = rowMinima[r];
            }
        }
    }
    virtual void cluster(size_t a, size_t b) {
        //Assumed 0<=a<b<n
        NJFloat nless2        = n-2;
        NJFloat tMultiplier   = (n<3) ? 0 : (0.5 / nless2);
        NJFloat medianLength  = 0.5 * rows[a][b];
        NJFloat fudge         = (rowTotals[a] - rowTotals[b]) * tMultiplier;
        NJFloat aLength       = medianLength + fudge;
        NJFloat bLength       = medianLength - fudge;
        NJFloat lambda        = 0.5;
        NJFloat mu            = 1.0 - lambda;
        NJFloat dCorrection   = - lambda * aLength - mu * bLength;
        for (int i=0; i<n; ++i) {
            if (i!=a && i!=b) {
                size_t  x     = (i<a) ? i : a;
                size_t  y     = (a<i) ? a : i;
                
                NJFloat Dai   = rows[a][i];
                NJFloat Dbi   = rows[b][i];
                NJFloat Dci   = lambda * Dai + mu * Dbi + dCorrection;
                rows[a][i]    = Dci;
                rows[i][a]    = Dci;
                rowTotals[i] += Dci - Dai - Dbi; //JB2020-06-18 Adjust row totals
                rowTotals[a] += Dci - Dai;       //on the fly.
            }
        }
        rowTotals[a] -= rows[a][b];
        clusters.emplace_back ( rowToCluster[a], aLength,
                                rowToCluster[b], bLength);
        rowToCluster[a] = clusters.size()-1;
        rowToCluster[b] = rowToCluster[n-1];
        removeRowAndColumn(b);
    }
    void finishClustering() {
        //Assumes that n is 3
        NJFloat halfD01 = 0.5 * rows[0][1];
        NJFloat halfD02 = 0.5 * rows[0][2];
        NJFloat halfD12 = 0.5 * rows[1][2];
        clusters.emplace_back
            ( rowToCluster[0], halfD01 + halfD02 - halfD12
            , rowToCluster[1], halfD01 + halfD12 - halfD02
            , rowToCluster[2], halfD02 + halfD12 - halfD01);
        n = 0;
    }
    virtual void doClustering() {
        while (3<n) {
            Position best;
            getMinimumEntry(best);
            cluster(best.column, best.row);
        }
        finishClustering();
    }
    void writeTreeFile(const std::string &treeFilePath) const {
        std::vector<Place> stack;
        std::fstream out;
        out.open(treeFilePath, std::ios_base::out);
        out.precision(8);
        bool failed = false; //Becomes true if clusters
                             //defines cycles (should never happen)
                             //Indicates a fatal logic error
        int maxLoop = 3 * clusters.size();
                             //More than this, and there must be
                             //a cycle.  Or something.

        stack.emplace_back(clusters.size()-1, 0);
        do {
            --maxLoop;
            if (maxLoop==0) {
                failed = true;
                break;
            }
            Place here = stack.back();
            const Cluster& cluster = clusters[here.clusterIndex];
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
            const Link & linkPrev = cluster.links[nextChildNum-1];
            out << ":" << linkPrev.linkDistance;
            if (nextChildNum<cluster.links.size()) {
                out << ",";
                const Link & linkNext = cluster.links[nextChildNum];
                stack.emplace_back(here.clusterIndex, nextChildNum+1);
                stack.emplace_back(linkNext.clusterIndex, 0);
            } else {
                out << ")";
            }
        } while (0 < stack.size());
        out << ";" << std::endl;
        out.close();
    }
};

class BIONJMatrix : public NJMatrix {
protected:
    Matrix  variance;
    typedef NJMatrix super;
public:
    explicit BIONJMatrix(const std::string &distanceMatrixFilePath)
        : super(distanceMatrixFilePath){
        variance = *this;
    }
    inline NJFloat chooseLambda(size_t a, size_t b, NJFloat Vab) {
        //Assumed 0<=a<b<n
        NJFloat lambda = 0;
        if (Vab==0.0) {
            return 0.5;
        }
        for (int i=0; i<a; ++i) {
            lambda += variance.rows[b][i] - variance.rows[a][i];
        }
        for (int i=a+1; i<b; ++i) {
            lambda += variance.rows[b][i] - variance.rows[a][i];
        }
        for (int i=b+1; i<n; ++i) {
            lambda += variance.rows[b][i] - variance.rows[a][i];
        }
        lambda = 0.5 + lambda / (2.0*((NJFloat)n-2)*Vab);
        if (1.0<lambda) lambda=1.0;
        if (lambda<0.0) lambda=0.0;
        return lambda;
    }
    virtual void cluster(size_t a, size_t b) {
        //Assumed 0<=a<b<n
        //Bits that differ from super::cluster tagged BIO
        NJFloat nless2        = n - 2 ;
        NJFloat tMultiplier   = ( n < 3 ) ? 0 : ( 0.5 / nless2 );
        NJFloat medianLength  = 0.5 * rows[b][a];
        NJFloat fudge         = (rowTotals[a] - rowTotals[b]) * tMultiplier;
        NJFloat aLength       = medianLength + fudge;
        NJFloat bLength       = medianLength - fudge;
        NJFloat Vab           = variance.rows[b][a];     //BIO
        NJFloat lambda        = chooseLambda(a, b, Vab); //BIO
        NJFloat mu            = 1.0 - lambda;
        NJFloat dCorrection   = - lambda * aLength - mu * bLength;
        NJFloat vCorrection   = - lambda * mu * Vab;
        NJFloat replacementRowTotal = 0;
        #pragma omp parallel for
        for (int i=0; i<n; ++i) {
            if (i!=a && i!=b) {
                size_t  x     = (i<a) ? i : a;
                size_t  y     = (a<i) ? a : i;
              
                //Dci as per reduction 4 in [Gascuel]
                NJFloat Dai   = rows[a][i];
                NJFloat Dbi   = rows[b][i];
                NJFloat Dci   = lambda * Dai + mu * Dbi + dCorrection;
                rows[a][i]    = Dci;
                rows[i][a]    = Dci;
                rowTotals[i] += Dci - Dai - Dbi; //JB2020-06-18 Adjust row totals
                
                //BIO begin (Reduction 10 on variance estimates)
                NJFloat Vci   = lambda * variance.rows[a][i]
                              + mu * variance.rows[b][i]
                              + vCorrection;
                variance.rows[a][i] = Vci;
                variance.rows[i][a] = Vci;
                //BIO finish
            }
        }
        replacementRowTotal = 0;
        for (int i=0; i<a; ++i) {
            replacementRowTotal += rows[a][i];
        }
        for (int i=a+1; i<b; ++i) {
            replacementRowTotal += rows[a][i];
        }
        for (int i=b+1; i<n; ++i) {
            replacementRowTotal += rows[a][i];
        }
        rowTotals[a] = replacementRowTotal;
        clusters.emplace_back ( rowToCluster[a], aLength,
                                rowToCluster[b], bLength);
        rowToCluster[a] = clusters.size()-1;
        rowToCluster[b] = rowToCluster[n-1];
        removeRowAndColumn(b);
        variance.removeRowAndColumn(b); //BIO
    }
};

template <class super=BIONJMatrix> class BoundingMatrix: public super
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
    std::vector<int>     clusterToRow;   //Maps clusters to their rows
    std::vector<NJFloat> clusterTotals;  //"Row" totals indexed by cluster

    mutable std::vector<NJFloat> scaledClusterTotals;   //The same, multiplied by
                                                        //(1.0 / (n-2)).
    mutable std::vector<NJFloat> scaledMaxEarlierClusterTotal;
    mutable std::vector<bool>    rowOrderChosen; //Indicates if row order chosen
    mutable std::vector<size_t>  rowScanOrder;   //Order in which rows are to be scanned
                                                 //Only used in... getRowMinima().
    
    Matrix<NJFloat>      entriesSorted; //Entries in distance matrix
                                        //(each row sorted by ascending value)
    Matrix<int>          entryToCluster;//
    
public:
    BoundingMatrix(const std::string &distanceFilePath)
    : super(distanceFilePath) {
    }
    virtual void doClustering() {
        //1. Set up vectors indexed by cluster number,
        clusterToRow.resize(n);
        clusterTotals.resize(n);
        for (size_t r=0; r<n; ++r) {
            clusterToRow[r]  = r;
            clusterTotals[r] = rowTotals[r];
        }
        
        //2. Set up "scratch" vectors used in getRowMinima
        //   so that it won't be necessary to reallocate them
        //   for each call.
        scaledClusterTotals.resize(n);
        scaledMaxEarlierClusterTotal.resize(n);
        rowOrderChosen.resize(n);
        rowScanOrder.resize(n);

        //2. Set up the matrix with row sorted by distance
        //   And the matrix that tracks which distance is
        //   to which cluster (the S and I matrices, in the
        //   RapidNJ papers).
        entriesSorted.setSize(n);
        entryToCluster.setSize(n);
        #pragma omp parallel for
        for (size_t r=0; r<n; ++r) {
            sortRow(r,r);
            //copies the "left of the diagonal" portion of
            //row r from the D matrix and sorts it
            //into ascending order.
        }
        size_t nextPurge = 0; //n*2/3;
        while (3<n) {
            Position best;
            super::getMinimumEntry(best);
            cluster(best.column, best.row);
            if ( n == nextPurge ) {
                #pragma omp parallel for
                for (int r=0; r<n; ++r) {
                    purgeRow(r);
                }
                nextPurge = n*2/3;
            }
        }
        super::finishClustering();
    }
    void sortRow(size_t r /*row index*/, size_t c /*upper bound on cluster index*/) {
        //1. copy data from a row of the D matrix into the S matrix
        //   (and write the cluster identifiers that correspond to
        //    the values in the D row into the same-numbered
        //    row in the I matrix), for distances between the cluster
        //    in that row, and other live clusters (up to, but not including c).
        NJFloat* sourceRow      = rows[r];
        NJFloat* values         = entriesSorted.rows[r];
        int*     clusterIndices = entryToCluster.rows[r];
        size_t   w = 0;
        for (size_t i=0; i<n; ++i) {
            values[w]         = sourceRow[i];
            clusterIndices[w] = rowToCluster[i];
            if ( i!=r && clusterIndices[w] < c ) {
                ++w;
            }
        }
        values[w]         = infiniteDistance; //sentinel value, to stop row search
        clusterIndices[w] = 0; //Always room for this, because distance to self
                               //was excluded via the i!=r check above.
        
        //2. Sort the row in the S matrix and mirror the sort
        //   on the same row of the I matrix.
        mirroredHeapsort(values, 0, w, clusterIndices);
    }
    void purgeRow(size_t r /*row index*/) {
        //Scan a row of the I matrix, so as to remove
        //entries that refer to clusters that are no longer
        //being processed. Remove the corresponding values
        //in the same row of the S matrix.
        NJFloat* values         = entriesSorted.rows[r];
        int*     clusterIndices = entryToCluster.rows[r];
        size_t w = 0;
        size_t i = 0;
        for (; i<n ; ++i ) {
            values[w]         = values[i];
            clusterIndices[w] = clusterIndices[i];
            if ( infiniteDistance <= values[i] ) {
                break;
            }
            if ( 0 <= clusterToRow[clusterIndices[i]] ) {
                ++w;
            }
        }
    }
    virtual void cluster(size_t a, size_t b) {
        size_t clusterA = rowToCluster[a];
        size_t clusterB = rowToCluster[b];
        size_t clusterMoved = rowToCluster[n-1];
        clusterToRow[clusterA] = -1;
        clusterToRow[clusterB] = -1;
        size_t c = clusters.size(); //cluster # of new cluster
        super::cluster(a,b);
        clusterToRow.emplace_back(a);
        clusterTotals.emplace_back(rowTotals[a]);
        scaledClusterTotals.emplace_back(rowTotals[a] / (n-1));
        scaledMaxEarlierClusterTotal.emplace_back(0.0);
        if (b<n) {
            clusterToRow[clusterMoved] = b;
        }
        //Mirror row rearrangement done on the D (distance) matrix
        //(and possibly also on the V (variance estimate) matrix),
        //onto the S and I matrices.
        entriesSorted.removeRowOnly(b);
        entryToCluster.removeRowOnly(b);
        
        //Recalculate cluster totals.
        for (size_t wipe = 0; wipe<c; ++wipe) {
            clusterTotals[wipe] = -infiniteDistance;
            //A trick.  This way we don't need to check if clusters
            //are still "live" in the inner loop of getRowMinimum().
        }
        for (size_t r = 0; r<n; ++r) {
            size_t cluster = rowToCluster[r];
            clusterTotals[cluster] = rowTotals[r];
            //std::cout << "row " << r << ", cluster " << cluster
            //    << " have total " << rowTotals[r] << std::endl;
        }
        sortRow(a, clusters.size());
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
        //Since we always have to check these entries when we process
        //the row, why not process them up front, hoping to get a
        //better bound on min(V) (and "rule out" entire rows with that
        //bound). -James B).
        //
        std::sort(rowMinima.begin(), rowMinima.end());
        for (size_t i=0; i<n; ++i) {
            rowOrderChosen[i]=false;
        }
        size_t w = 0;
        for (size_t r=0; r < rowMinima.size()
             && rowMinima[r].value < infiniteDistance; ++r) {
            size_t rowA = rowMinima[r].row;
            size_t rowB = rowMinima[r].column;
            size_t clusterA = rowToCluster[rowA];
            size_t clusterB = rowToCluster[rowB];
            size_t row = (clusterA<clusterB) ? rowB : rowA;
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
        //      totals multiplied by (1.0/(NJFloat)(n-2)).
        //      Better n multiplications than 0.5*n*(n-1).
        //Note 2: Note that these are indexed by cluster number,
        //      and *not* by row number.
        //
        size_t  c           = clusters.size();
        NJFloat nless2      = ( n - 2 );
        NJFloat tMultiplier = ( n <= 2 ) ? 0 : (1 / nless2);
        NJFloat maxTot      = -infiniteDistance; //maximum row total divided by (n-2)
        for (size_t i=0; i<c; ++i) {
            scaledClusterTotals[i] = clusterTotals[i] * tMultiplier;
            scaledMaxEarlierClusterTotal[i] = maxTot;
            if ( 0 <= clusterToRow[i] && maxTot < scaledClusterTotals[i] ) {
                maxTot=scaledClusterTotals[i];
            }
        }
        NJFloat qBest = infiniteDistance;
            //upper bound on minimum Q[row,col]
            //  = D[row,col] - R[row]*tMultipler - R[col]*tMultiplier
            //

        decideOnRowScanningOrder();
        rowMinima.resize(n);
        #pragma omp parallel for
        for (size_t r=0; r<n; ++r) {
            size_t row             = rowScanOrder[r];
            size_t cluster         = rowToCluster[row];
            size_t maxEarlierTotal = scaledMaxEarlierClusterTotal[cluster];
            //Note: Older versions of RapidNJ used maxTot rather than
            //      maxEarlierTotal here...
            rowMinima[row]         = getRowMinimum(row, maxEarlierTotal, qBest);
            NJFloat v              = rowMinima[row].value;
            //#pragma omp critical
            {
                if ( v < qBest ) {
                    qBest = v;
                }
            }
        }
    }
    Position getRowMinimum(size_t row, NJFloat maxTot, NJFloat qBest) const {
        NJFloat nless2      = ( n - 2 );
        NJFloat tMultiplier = ( n <= 2 ) ? 0 : (1 / nless2);
        auto    tot         = scaledClusterTotals.data();
        NJFloat rowTotal    = rowTotals[row] * tMultiplier;
        NJFloat rowBound    = qBest + maxTot + rowTotal;
                //Upper bound for distance, in this row, that
                //could (after row totals subtracted) provide a
                //better min(Q).

        Position pos(row, 0, infiniteDistance);
        const NJFloat* rowData   = entriesSorted.rows[row];
        const int*     toCluster = entryToCluster.rows[row];
        size_t i = 0;
        NJFloat Drc;
        for (i=0; (Drc=rowData[i])<rowBound; ++i) {
            size_t  cluster = toCluster[i];
                //The cluster associated with this distance
                //The c in Qrc and Drc.
            NJFloat Qrc = Drc - tot[cluster] - rowTotal;
            if (Qrc < pos.value) {
                int otherRow = clusterToRow[cluster];
                if ( 0 <= otherRow ) { //I *think* this check is still necessary,
                                       //despite setting "out of matrix" cluster totals
                                       //to (0-infiniteDistance).
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

class VectorizedBIONJMatrix: public BIONJMatrix
{
    //
    //Note: this is a first attempt at hand-vectorizing
    //      BIONJMatrix::getRowMinima (via Agner Fog's
    //      vectorclass library).
    //
public:
    typedef BIONJMatrix super;
    VectorizedBIONJMatrix(const std::string &distanceMatrixFilePath)
    : super(distanceMatrixFilePath)
    {
    }
    virtual void getRowMinima(Positions& rowMinima) const {
        NJFloat nless2      = ( n - 2 );
        NJFloat tMultiplier = ( n <= 2 ) ? 0 : (1 / nless2);
        boost::scoped_array<NJFloat> scratchTotals(new NJFloat[n]);
        boost::scoped_array<NJFloat> scratchColumnNumbers(new NJFloat[n]);
        auto tot  = scratchTotals.get();
        auto nums = scratchColumnNumbers.get();
        for (size_t r=0; r<n; ++r) {
            tot[r]  = rowTotals[r] * tMultiplier;
            nums[r] = r;
        }
        Vec4d v;
        size_t blockSize = v.size();
        rowMinima[0].value = infiniteDistance;
        #pragma omp parallel for
        for (size_t row=1; row<n; ++row) {
            Position pos(row, 0, infiniteDistance);
            const NJFloat* rowData = rows[row];
            size_t col;
            Vec4d minVector(infiniteDistance, infiniteDistance, infiniteDistance, infiniteDistance);
                //The minima of columns with indices "congruent modulo 4"
                //For example minVector[1] holds the minimum of
                //columns 1,5,9,13,17,...
            Vec4d ixVector(-1,     -1,     -1,     -1);
                //For each entry in minVector, the column from which
                //that value came.
            
            //Take#2: Examine four columns at a time
            for (col=0; col+blockSize<row; col+=blockSize) {
                Vec4d  rowVector; rowVector.load(rowData+col);
                Vec4d  totVector; totVector.load(tot+col);
                Vec4d  adjVector = rowVector - totVector;
                Vec4db less      = adjVector < minVector;
                Vec4d  numVector; numVector.load(nums+col);
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
            
            /*Take#1 was worse.  Or-ing each 4 right away.  Ee-uw.
            for (col=0; col+blockSize<row; col+=blockSize) {
                Vec4d rowVector;rowVector.load(rowData+col);
                Vec4d totVector;totVector.load(tot+col);
                Vec4d adjVector = rowVector - totVector;
                if (horizontal_or(adjVector < minVector)) {
                    for (int c=0; c<blockSize; ++c) {
                        if (minVector[c] < pos.value) {
                            pos.column = col+c;
                            pos.value  = minVector[c];
                        }
                    }
                    minVector = pos.value;
                }
            }
            */
            
            for (; col<row; ++col) {
                NJFloat v = rowData[col] - tot[col];
                if (v < pos.value) {
                    pos.column = col;
                    pos.value  = v;
                }
            }
            pos.value -= tot [row];
            rowMinima[row] = pos;
        }
    }
};//end of class

} //end of anonymous namespace

void BIONJ2::constructTree
    ( const std::string & distanceMatrixFilePath
    , const std::string & newickTreeFilePath)
{
    /*Vectorized*/BIONJMatrix d(distanceMatrixFilePath);
    double joinStart = getRealTime();
    d.doClustering();
    double joinElapsed = getRealTime() - joinStart;
    std::cout.precision(6);
    std::cout << "Elapsed time for neighbour joining proper (in BIONJ2), " << joinElapsed << std::endl;
    std::cout.precision(3);
    d.writeTreeFile(newickTreeFilePath);
}

void BIONJ2::constructTreeRapid
    ( const std::string & distanceMatrixFilePath
    , const std::string & newickTreeFilePath)
{
    BoundingMatrix<BIONJMatrix> d(distanceMatrixFilePath);
    double joinStart = getRealTime();
    d.doClustering();
    double joinElapsed = getRealTime() - joinStart;
    std::cout.precision(6);
    std::cout << "Elapsed time for neighbour joining proper (in BIONJ2/rapidNJ), "
        << joinElapsed << std::endl;
    std::cout.precision(3);
    d.writeTreeFile(newickTreeFilePath);

    VectorizedBIONJMatrix d2(distanceMatrixFilePath);
    joinStart = getRealTime();
    d2.doClustering();
    joinElapsed = getRealTime() - joinStart;
    std::cout.precision(6);
    std::cout << "Elapsed time for neighbour joining proper (in BIONJ2/Hand-Vectorized), " << joinElapsed << std::endl;
    std::cout.precision(3);
}
