//  stitchup.cpp - Implements the "Family Stitch-up" (distance matrix)
//                 tree construction algorithm, which works by "stitching up"
//                 a graph, based on probable familial relationships
//                 (steps 1 through 3), and then "removing the excess
//                 stitches" (in step 4).
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
//  1. For each leaf node, a "caterpillar chain" of nodes is maintained,
//     (initially, each leaf node is the only node in its chain). Interior
//      nodes are added only on the ends of chains.
//  2. The pair of leaf nodes, A, B, with the lowest observed distance, d(A,B)
//     that are not already connected are selected, and new nodes Ai, Bi,
//     are added to the end their caterpillar chains, with Ai and Bi
//     connected to the previous ends of the "caterpillar chains",
//     and connected to each other by a edge with length d(A,B)*STAPLE_ARCH
//     (The length of the edge that connects Ai to caterpillar chain for A
//      is STAPLE_LEG*(d(A,B)-d(A,Ap)) where Ap is the node that was previously
//      at the end of the chain);
//     (so: STAPLE_LEG is (0.5*(1.0-STAPLE_ARCH)).
//
//     (Implementation details: A min-heap is used to find short edges,
//      and connectedness is tracked via a union-find structure).
//  3. Step 2 is repeated until all of the leaf nodes are connected.
//     At this point, each leaf node will be degree 1, and there will be 2*(n-1)
//     interior nodes, about half of which (the nodes at the "top" of caterpillar
//     chains) will be degree 2.
//  4. Remove all the nodes of degree 2 (by directly linking the two nodes
//     they formerly connected, with an edge of length equal to the sum of those
//     of the edges, incident to the degree 2 node, just removed).
//
//  Neighbour Joining (NJ) and BIONJ work by trying to guess where interior nodes
//  will be, relative to the nodes that they are linked to (and later joins depend
//  on those positional guesses) (the strategy is "guess the geometry to use to
//  decide on the structure") (but then: the tree topology depends on *guesses*).
//
//  Family stitch-up places an each way bet, inserting *two* internal nodes (each
//  close to one of the leaf nodes getting linked) (Step 2), and only later removes
//  the (degree 2) nodes that correspond to "possibilities that didn't pan out"
//  (in Step 4) (in a nutshell: the strategy is:
//  "let the leaf-distances, alone, decide the topology", and then, only later,
//  "let the topology decide the geometry").
//
//  It is up to some later algorithm to choose *better* lengths for the
//  zero-length edges (some of which might exist, along the caterpillar chains).
//  In practice the trees that come out of STITCHUP are fed into a Maximum Likelihood
//  framework that will choose better lengths for those ("caterpillar") edges, so:
//  it's not something I felt that I needed to worry about.
//
//  Running time: O((n^2).ln(n)) in the worst case (dominated by: heap extraction).
//                A little worse than O(n^2) in practice
//                (dominated by: heap construction).
//  Notes:     1. The union-find structure used here has a ~ ln(n)
//                worst case. And the time to remove the degree-2 nodes
//                is proportional to n.ln(n) (it is dominated by the time
//                taken by the... stitches.insert()... line).
//             2. Both could be more efficient (union-find's running time
//                could easily be lowered to n*A(n) where A is the inverse 
//                of Ackermann's function), and the removal of the degree-2
//                nodes can be done in time linear in n.  But given that other 
//                steps take so much longer there is no need to make these steps 
//                more efficient.
//             3. The union-find structure used in NearestTaxonClusterJoiningMatrix
//                is even less efficient (it is quadratic; ~ n*n/2 operations).
//                Again, this doesn't cost much: there are other parts of the
//                algorithm that take far longer.
//
//  Created by James Barbetti on 12-Aug-2020 (tree construction)
//  and 24-Aug-20 (generating the newick file).
//

#include "starttree.h"
#include "distancematrix.h"
#include "nj.h"
#include "clustertree.h"
#include "progress.h"
#include "heapsort.h" //for MinHeap template class
#include <set>
#include <math.h> //for floor()
namespace StartTree
{

template <class T=double> class Stitch { //an Edge in a stitch-up graph
public:
    size_t source;      //
    size_t destination; //
    T      length;      //
    Stitch() : source(0), destination(0), length(0) { }
    Stitch(size_t sourceIndex, size_t destinationIndex, T edgeLength):
        source(sourceIndex), destination(destinationIndex), length(edgeLength) {
    }
    Stitch& operator= (const Stitch& rhs) {
        source = rhs.source;
        destination = rhs.destination;
        length = rhs.length;
        return *this;
    }
    bool operator < (const Stitch<T>& rhs) const {
        if (source<rhs.source) return true;
        if (rhs.source<source) return false;
        return destination<rhs.destination;
    }
    bool operator <= (const Stitch<T>& rhs) const {
        if (source<rhs.source) return true;
        if (rhs.source<source) return false;
        return destination<=rhs.destination;
    }
    Stitch converse() const {
        return Stitch(destination, source, length);
    }
};

namespace {
    size_t lastHack = 1;
}

template <class T=double> struct LengthSortedStitch: public Stitch<T> {
public:
    typedef Stitch<T> super;
    using   super::length;
    size_t  hack; //Used to impose a pseudo-random ordering on equal-length edges
    LengthSortedStitch() : super(0,0,0.0) {}
    LengthSortedStitch(size_t sourceIndex, size_t destinationIndex, T edgeLength):
        super(sourceIndex, destinationIndex, edgeLength ) {
        lastHack = lastHack * 2862933555777941757UL + 3037000493UL;
        hack     = lastHack;
    }
    bool operator < (const LengthSortedStitch<T>& rhs) const {
        if (length<rhs.length) return true;
        if (rhs.length<length) return false;
        return (hack<rhs.hack);
    }
    bool operator <= (const LengthSortedStitch<T>& rhs) const {
        if (length<rhs.length) return true;
        if (rhs.length<length) return false;
        return (hack<=rhs.hack);
    }
};

#define STAPLE_ARCH (1.0/3.0)
#define STAPLE_LEG  (0.5*(1.0-STAPLE_ARCH))

template <class T=double> struct StitchupGraph {
    StrVector                leafNames;
    std::set< Stitch<T> >    stitches;
    std::vector< int >       taxonToSetNumber;
    std::vector< int >       taxonToNodeNumber;
    std::vector< T   >       taxonToDistance;
    std::vector< IntVector > setMembers;
    int                      nodeCount;
    bool                     silent;
    StitchupGraph() : nodeCount(0) {
    }
    void clear() {
        StitchupGraph temp;
        std::swap(*this, temp);
        nodeCount = 0;
    }
    const std::string& operator[] (size_t index) const {
        return leafNames[index];
    }
    void addLeaf(const std::string& name) {
        leafNames.emplace_back(name);
        taxonToSetNumber.emplace_back(nodeCount);
        taxonToNodeNumber.emplace_back(nodeCount);
        taxonToDistance.emplace_back(0);
        std::vector<int> singletonSet;
        singletonSet.push_back(nodeCount);
        setMembers.push_back(singletonSet);
        ++nodeCount;
    }
    bool areLeavesInSameSet(size_t leafA, size_t leafB) {
        return taxonToSetNumber[leafA]
                == taxonToSetNumber[leafB];
    }
    int staple(size_t leafA, size_t leafB, T length) {
        int interiorA = nodeCount;
        T legLengthA = (length - taxonToDistance[leafA]) * STAPLE_LEG;
        stitchLink(taxonToNodeNumber[leafA], interiorA, legLengthA);
        taxonToNodeNumber[leafA] = interiorA;
        taxonToDistance[leafA] = legLengthA;
        ++nodeCount;
        
        int interiorB = nodeCount;
        T legLengthB = (length - taxonToDistance[leafB]) * STAPLE_LEG;
        stitchLink(taxonToNodeNumber[leafB], interiorB, legLengthB);
        taxonToNodeNumber[leafB] = interiorB;
        taxonToDistance[leafB] = legLengthB;
        ++nodeCount;
        
        stitchLink(interiorA, interiorB, length * STAPLE_ARCH);
        
        int setA = taxonToSetNumber[leafA];
        int setB = taxonToSetNumber[leafB];
        #if (0)
            int setASize = setMembers[setA].size();
            int setBSize = setMembers[setB].size();
        #endif
        int setC = mergeSets(setA, setB);
        
        #if (0)
            std::cout << "Staple " << leafA << ":" << interiorA << "-" << length << "-"
                << interiorB << ":" << leafB
                << " (sets " << setA << " (size " << setASize << ")"
                << " and " << setB << " (size " << setBSize << ") to " << setC << ")"
                << " " << leafNames[leafA] << " to " << leafNames[leafB]
                << std::endl;
        #endif
        return setC;
    }
    void stitchLink(int nodeA, int nodeB, T length) {
        stitches.insert(Stitch<T>(nodeA, nodeB, length));
        stitches.insert(Stitch<T>(nodeB, nodeA, length));
    }
    int mergeSets(int setA, int setB) {
        if (setA == setB) {
            return setA;
        }
        std::vector<int>& membersA = setMembers[setA];
        std::vector<int>& membersB = setMembers[setB];
        if (membersA.size() < membersB.size()) {
            for (auto it = membersA.begin(); it != membersA.end(); ++it) {
                int a = *it;
                taxonToSetNumber[a] = setB;
                membersB.push_back(a);
            }
            membersA.clear();
            return setB;
        } else {
            for (auto it = membersB.begin(); it != membersB.end(); ++it) {
                int b = *it;
                taxonToSetNumber[b] = setA;
                membersA.push_back(b);
            }
            membersB.clear();
            return setA;
        }
    }
    void removeThroughThroughNodes() {
        //Removes any "through-through" interior nodes of degree 2.
        #if USE_PROGRESS_DISPLAY
        const char* taskDescription = silent
            ? "" : "Removing degree-2 nodes from stitchup graph";
        progress_display progress ( stitches.size()*2,
                                    taskDescription, "", "");
        #else
        double progress = 0.0;
        #endif
        std::vector<intptr_t> replacements;
        std::vector<T>        replacementLengths;
        replacements.reserve(nodeCount);
        replacementLengths.resize(nodeCount, 0);
        for (int i=0; i<nodeCount; ++i) {
            replacements.push_back(i);
        }
        intptr_t node      = -1; //Source node of last edge
        size_t   degree    = 0;  //Degree of that node
        for (auto it=stitches.begin(); it!=stitches.end(); ++it) {
            if (it->source != node) {
                if (node!=-1) {
                    if (degree!=2) {
                        replacements[node] = node;
                        replacementLengths[node] = 0;
                    } else {
                        //std::cout << "replacing " << node
                        //  << " with " << replacements[node] << std::endl;
                    }
                }
                node   = it->source;
                degree = 1;
                if ((intptr_t)it->destination < node) {
                    replacements[node]       = it->destination;
                    replacementLengths[node] = it->length;
                }
            } else {
                ++degree;
            }
            ++progress;
        }
        if (degree!=2 && node!=-1) {
            replacements[node] = node;
            replacementLengths[node] = 0;
        }
        //Remove them (adjusting the lengths of later edges
        //that have to take over from them, at the same time).
        std::set< Stitch<T> > oldStitches;
        std::swap(stitches, oldStitches);
        for (auto it=oldStitches.begin(); it!=oldStitches.end(); ++it) {
            T        length = it->length;
            intptr_t source = replacements[it->source];
            intptr_t dest   = replacements[it->destination];
            if (source!=dest) {
                length += replacementLengths[it->source];
                length += replacementLengths[it->destination];
                stitches.insert(Stitch<T>(source, dest, length));
            }
            ++progress;
        }
        #if USE_PROGRESS_DISPLAY
        progress.done();
        #endif
    }
    template <class F>
    void dumpTreeToFile ( const std::string &treeFilePath, F& out ) const {
        int cols = 0;
        for (auto it=stitches.begin(); it!=stitches.end(); ++it) {
            if (it->source<it->destination) {
                ++cols;
                std::cout << it->source << ":" << it->destination 
                    << " " << it->length << "\t";
                if (cols==4) {
                    std::cout << std::endl;
                    cols = 0;
                }
            }
        }
        std::cout << std::endl;
    }
    template <class F>
    bool writeTreeToOpenFile (int precision, 
                              progress_display_ptr progress,
                              F& out) const {
        out.precision(precision);
        auto lastEdge = stitches.end();
        --lastEdge;
        size_t lastNodeIndex = lastEdge->source;
        size_t edgeCount = stitches.size();
        std::vector<Stitch<T>> stitchVector;
        std::vector<size_t>    nodeToEdge;
        nodeToEdge.resize(lastNodeIndex+1, edgeCount);
        int j = 0;
        for (auto it=stitches.begin(); it!=stitches.end(); ++it, ++j) {
            stitchVector.push_back(*it);
            size_t i = it->source;
            if (nodeToEdge[i]==edgeCount) {
                nodeToEdge[i] = j;
            }
        }
        writeSubtree(stitchVector, nodeToEdge, nullptr, 
                     lastNodeIndex, progress, out);
        out << ";" << std::endl;
        return true;
    }
    template <class F>
    bool writeTreeToFile(int precision, 
                         const std::string &treeFilePath, 
                         F& out ) const {
        bool success = false;
        std::string desc = "Writing STITCH tree to ";
        desc+=treeFilePath;
        #if USE_PROGRESS_DISPLAY
        progress_display progress(stitches.size(), desc.c_str(), "wrote", "edge");
        #else
        double progress = 0.0;
        #endif
        out.exceptions(std::ios::failbit | std::ios::badbit);
        try {
            out.open(treeFilePath.c_str(), std::ios_base::out);
            success = writeTreeToOpenFile(precision, &progress, out);
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
        out.close();
        #if USE_PROGRESS_DISPLAY
        progress.done();
        #endif
        return success;
    }
    template <class F>
    void writeSubtree ( const std::vector<Stitch<T>> stitchVector, 
                        std::vector<size_t>  nodeToEdge,
                        const Stitch<T>* backstop, size_t nodeIndex,
                        progress_display_ptr progress, F& out) const {
        bool isLeaf = ( nodeIndex < leafNames.size() );
        if (isLeaf) {
            out << leafNames [ nodeIndex ] ;
        } else {
            out << "(";
            const char* sep = "";
            size_t x = nodeToEdge[nodeIndex];
            size_t y = stitchVector.size();
            nodeToEdge[nodeIndex] = y;
            for (; x<y && stitchVector[x].source == nodeIndex; ++x) {
                size_t child = stitchVector[x].destination;
                if ( nodeToEdge[child] != y /*no backsies*/ ) {
                    out << sep;
                    sep = ",";
                    writeSubtree(stitchVector, nodeToEdge, &stitchVector[x], 
                                 child, progress, out);
                }
                ++progress;
            }
            out << ")";
        }
        if (backstop!=nullptr) {
            out << ":" << backstop->length;
        }
    }
    bool writeTreeFile(bool zipIt, int precision, 
                       const std::string &treeFilePath) const {
        if (treeFilePath == "STDOUT") {
            #if USE_PROGRESS_DISPLAY
            progress_display progress(stitches.size(), "", "", "");
            #else
            double progress=0;
            #endif
            return writeTreeToOpenFile(precision, &progress, std::cout);
        } else if (zipIt) {
            #if USE_GZSTREAM
            ogzstream     out;
            #else
            std::ofstream out;
            #endif
            return writeTreeToFile(precision, treeFilePath, out);
        } else {
            std::fstream out;
            return writeTreeToFile(precision, treeFilePath, out);
        }
    }
};

template < class T=double> class StitchupMatrix: public SquareMatrix<T> {
public:
    typedef SquareMatrix<T> super;
    using super::rows;
    using super::setSize;
    using super::row_count;
    using super::column_count;
    using super::loadDistancesFromFlatArray;
    bool silent;
    StitchupMatrix(): isOutputToBeZipped(false) {
    }
    virtual std::string getAlgorithmName() const {
        return "STITCHUP";
    }
    virtual void addCluster(const std::string &name) {
        graph.addLeaf(name);
    }
    virtual bool loadMatrixFromFile(const std::string &distanceMatrixFilePath) {
        graph.clear();
        return loadDistanceMatrixInto(distanceMatrixFilePath, true, *this);
    }
    virtual bool loadMatrix
        ( const StrVector& names, const double* matrix ) {
        //Assumptions: 2 < names.size(), all names distinct
        //  matrix is symmetric, with matrix[row*names.size()+col]
        //  containing the distance between taxon row and taxon col.
        setSize(static_cast<intptr_t>(names.size()));
        graph.clear();
        for (auto it = names.begin(); it != names.end(); ++it) {
            addCluster(*it);
        }
        loadDistancesFromFlatArray(matrix);
        return true;
    }
    bool writeTreeFile(int precision, const std::string &treeFilePath) const {
        return graph.writeTreeFile(isOutputToBeZipped, precision, treeFilePath);
    }
    virtual void setZippedOutput(bool zipIt) {
        isOutputToBeZipped = zipIt;
    }
    virtual void beSilent() {
        silent = true;
    }
    virtual bool constructTree() {
        if (row_count<3) {
            return false;
        }
        std::vector<LengthSortedStitch<T>> stitches;
        stitches.reserve(row_count * (row_count-1) / 2);
        for (intptr_t row=0; row<row_count; ++row) {
            const T* rowData = rows[row];
            for (intptr_t col=0; col<row; ++col) {
                stitches.emplace_back(row, col, rowData[col]);
            }
        }
        MinHeapOnArray< LengthSortedStitch<T> >
            heap ( stitches.data(), stitches.size()
                 , silent ? "" : "Constructing min-heap of possible edges" );
        size_t iterations = 0;
        #if USE_PROGRESS_DISPLAY
        double row_count_triangle = 0.5*(double)row_count*(double)(row_count+1);
        const char* task_name = silent ? "" : "Assembling Stitch-up Graph";
        progress_display progress(row_count_triangle, task_name );
        #else
        double progress = 0.0;
        #endif
        std::cout.precision(12);
        for (intptr_t join = 0; join + 1 < row_count; ++join) {
            LengthSortedStitch<T> shortest;
            size_t source;
            size_t dest;
            do {
                shortest = heap.pop_min();
                source   = shortest.source;
                dest     = shortest.destination;
                ++iterations;
            } while ( graph.areLeavesInSameSet(source,dest) );
            graph.staple(source, dest, shortest.length);
            progress += (join+1);
        }
        #if USE_PROGRESS_DISPLAY
        progress.done();
        #endif
        graph.removeThroughThroughNodes();
        return true;
    }
protected:
    StitchupGraph<T>  graph;
    bool            isOutputToBeZipped;
};

template < class T=double, class S=NJMatrix<T> > 
class NearestTaxonClusterJoiningMatrix: public S {
    //
    //This is a mash-up of StitchupMatrix and Neighbor-Joining
    //(It works by considering the initial taxa distances, as
    //per the NJ distance metric, between points in different
    //clusters, joining the two clusters containing the taxa
    //that are closest, according to that metric).
    //It is somewhat faster than NJ but it gets worse answers.
    //For now, this seems to be a failed experiment.
    //
public:
    typedef NJMatrix<T> super;
    //member variables from super-class
    using super::row_count;
    using super::rows;
    using super::silent;
    using super::clusters;
    using super::rowTotals;
    //member functions from super-class
    using super::cluster;
    using super::finishClustering;
    
    virtual std::string getAlgorithmName() const {
        return "NTCJ";
    }
    virtual bool constructTree() {
        if (row_count<3) {
            return false;
        }
        class TaxonEdge {
        public:
            size_t taxon1, taxon2;
            T      length;
            bool operator <  ( const TaxonEdge &r) const {
                return length < r.length;
            }
            bool operator <= ( const TaxonEdge &r) const {
                return length < r.length;
            }
            TaxonEdge(): taxon1(0), taxon2(0), length(0) {}
            TaxonEdge(size_t t1, size_t t2, T dist) :
                taxon1(t1), taxon2(t2), length(dist) {}
            TaxonEdge& operator=(const TaxonEdge& rhs) = default;
        };
        T multiplier = 1.0 / (T)row_count;
        size_t row_count_triangle = row_count*(row_count-1)/2;
        std::vector<TaxonEdge> edges;
        edges.reserve(row_count_triangle);
        for (intptr_t row=0; row<row_count; ++row) {
            const T* rowData = rows[row];
            for (intptr_t col=0; col<row; ++col) {
                double d = rowData[col] - (rowTotals[row] + rowTotals[col]) *multiplier;
                edges.emplace_back(col, row, d);
            }
        }
        MinHeapOnArray< TaxonEdge >
            heap ( edges.data(), edges.size()
                 , silent ? "" : "Constructing min-heap of possible edges" );
        size_t iterations = 0;
        
        #if USE_PROGRESS_DISPLAY
        const char* task_name = silent ? "" : "Assembling NTCJ Tree";
        progress_display progress(row_count_triangle, task_name );
        #else
        double progress = 0.0;
        #endif

        intptr_t taxon_count = row_count;
        std::vector<size_t> taxonToRow;
        taxonToRow.resize(taxon_count);
        size_t* tr = taxonToRow.data();
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (intptr_t t=0; t<taxon_count; ++t) {
            tr[t] = t;
        }
        while (3<row_count) {
            TaxonEdge shortest;
            do {
                shortest = heap.pop_min();
                ++iterations;
            } while ( tr[shortest.taxon1] == tr[shortest.taxon2] );
            size_t rA = tr[shortest.taxon1];
            size_t rB = tr[shortest.taxon2];
            size_t r1 = (rA<rB) ? rA : rB;
            size_t r2 = (rB<rA) ? rA : rB;
            cluster( r1, r2);
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for (intptr_t t=0; t<taxon_count; ++t) {
                if (tr[t] == r2) {
                    tr[t] = r1;
                }
                else if (tr[t] == row_count) {
                    tr[t] = r2;
                }
            }
            progress += (taxon_count-row_count);
        }
        finishClustering();
        #if USE_PROGRESS_DISPLAY
        progress.done();
        #endif
        return true;
    }
};

void addStitchupTreeBuilders(Factory& f) {
    f.advertiseTreeBuilder( new Builder<StitchupMatrix<double>>
        ("STITCH",      "Family Stitch-up (Lowest Cost)"));
    f.advertiseTreeBuilder( new Builder<NearestTaxonClusterJoiningMatrix<double>>
        ("NTCJ",        "Cluster joining by nearest (NJ) taxon distance"));
}

} //end of StartTree namespace.
