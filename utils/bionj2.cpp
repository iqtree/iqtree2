//
//  bionj2.cpp - Implementations of NJ and BIONJ algorithms
//               (that work in terms of .mldist inputs and
//                NEWICK outputs).
//
//  BIONJ implementation based on http://www.lirmm.fr/~w3ifa/MAAS/BIONJ/BIONJ.html
//  NJ    implementation reverse engineered from same,
//
//  Created by James Barbetti on 18/6/20.
//

#include "bionj2.h"
#include <vector>
#include <string>
#include <fstream>
#include <iostream>               //for std::istream
#include <boost/scoped_array.hpp> //for boost::scoped_array
#include "utils/timeutil.h"       //for getRealTime()

typedef double NJFloat;

namespace
{

struct Position
{
    //A position (row, column) in an NJ matrix
    //Note that column is always less than row.
    //(Because that's the convention in RapidNJ).
public:
    size_t  column;
    size_t  row;
    NJFloat value;
    Position() : row(0), column(0), value(0) {}
    Position& operator = (const Position &rhs) {
        row    = rhs.row;
        column = rhs.column;
        value  = rhs.value;
        return *this;
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

struct Cluster: Position
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

class Matrix
{
    //Note: This is a separate class so that it can be
    //used for variance as well as distance matrices.
    //Lines that access the upper-right triangle
    //of the matrix are tagged with U-R.
    friend class NJMatrix;
    friend class BIONJMatrix;
protected:
    size_t n;
    NJFloat *data;
    NJFloat **rows;
    NJFloat *rowTotals;
    void setSize(size_t rank) {
        n         = rank;
        data      = (rank==0) ? nullptr : new NJFloat[rank*rank];
        rows      = (rank==0) ? nullptr : new NJFloat*[rank];
        rowTotals = (rank==0) ? nullptr : new NJFloat[rank];
        NJFloat *rowStart = data;
        for (int r=0; r<n; ++r) {
            rows[r] = rowStart;
            rowStart += n;
            rowTotals[r] = 0.0;
        }
    }
    void assign(const Matrix& rhs) {
        setSize(rhs.n);
        #pragma omp for
        for (size_t r=0; r<n; ++r) {
            NJFloat * destRow = rows[r];
            NJFloat const * sourceRow = rhs.rows[r];
            NJFloat const * const endSourceRow = sourceRow + n;
            for (; sourceRow<endSourceRow; ++destRow, ++sourceRow) {
                *destRow = *sourceRow;
            }
            rowTotals[r] = rhs.rowTotals[r];
        }
    }
public:
    Matrix() {
        setSize(0);
    }
    Matrix(const Matrix& rhs) {
        assign(rhs);
    }
    virtual ~Matrix() {
        clear();
    }
    void clear() {
        delete [] data;
        delete [] rows;
        delete [] rowTotals;
        data = nullptr;
        rows = nullptr;
        rowTotals = nullptr;
    }
    Matrix& operator=(const Matrix& rhs) {
        if (&rhs!=this) {
            clear();
            assign(rhs);
        }
        return *this;
    }
    size_t size() {
        return n;
    }
    void calculateRowTotals() const {
        #pragma omp for
        for (size_t r=0; r<n; ++r) {
            NJFloat total = 0;
            const NJFloat* rowData = rows[r];
            for (size_t c=0; c<r; ++c) {
                total += rowData[c];
            }
            for (size_t c=r+1; c<n; ++c) {
                total += rowData[c]; //U-R
            }
            
            rowTotals[r] = total;
        }
    }
    void removeRow(size_t rowNum)  {
        #pragma omp for
        for (size_t r=0; r<n; ++r) {
            NJFloat* rowData = rows[r];
            rowData[rowNum] = rowData[n-1]; //U-R
        }
        rows[rowNum] = rows[n-1];
        rows[n-1] = nullptr;
        --n;
    }
};

class NJMatrix: public Matrix
{
protected:
    std::vector<size_t>      rowToCluster;
    std::vector<Cluster>     clusters;
    Positions                rowMinima;
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
        //Note: The old code wrote a message to standard output,
        //      if the matrix was not symmetric.  This code doesn't.
    }
    virtual void getRowMinima(Positions& rowMinima) const {
        //
        //Note: Rather than multiplying distances by (n-2)
        //      repeatedly, it is cheaper to work with row
        //      totals multiplied by (1/(NJFloat)(n-2)).
        //      Better n multiplications than n*(n-1)/2.
        //
        NJFloat nless2      = ( n - 2 );
        NJFloat tMultiplier = ( n <= 2 ) ? 0 : (1 / nless2);
        boost::scoped_array<NJFloat> scratchTotals(new NJFloat[n]);
        auto tot = scratchTotals.get();
        for (size_t r=0; r<n; ++r) {
            tot[r] = rowTotals[r] * tMultiplier;
        }
        rowMinima.reserve(n);
        #pragma omp for
        for (size_t row=1; row<n; ++row) {
            Position pos;
            pos.row    = row;
            pos.column = 0;
            pos.value  = 1e+300;
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
        getRowMinima(rowMinima);
        for (size_t r=1; r<n; ++r) {
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
        removeRow(b);
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
    void doClustering() {
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
        lambda = 0.5 + lambda / (2.0*(n-2)*Vab);
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
                rowTotals[a] += Dci - Dai;       //on the fly.
                
                //BIO begin (Reduction 10 on variance estimates)
                NJFloat Vci   = lambda * variance.rows[a][i]
                              + mu * variance.rows[b][i]
                              + vCorrection;
                variance.rows[a][i] = Vci;
                variance.rows[i][a] = Vci;
                //BIO finish
            }
        }
        rowTotals[a] -= rows[a][b];
        clusters.emplace_back ( rowToCluster[a], aLength,
                                rowToCluster[b], bLength);
        rowToCluster[a] = clusters.size()-1;
        rowToCluster[b] = rowToCluster[n-1];
        removeRow(b);
        variance.removeRow(b); //BIO
    }
}; //end of class
} //end of anonymous namespace

void BIONJ2::constructTree
    ( const std::string &distanceMatrixFilePath
    , const std::string & newickTreeFilePath)
{
    BIONJMatrix d(distanceMatrixFilePath);
    double joinStart = getRealTime();
    d.doClustering();
    double joinElapsed = getRealTime() - joinStart;
    std::cout.precision(6);
    std::cout << "Elapsed time for neighbour joining proper (in BIONJ2), " << joinElapsed << std::endl;
    std::cout.precision(3);
    d.writeTreeFile(newickTreeFilePath);
}
