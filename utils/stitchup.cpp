//
//  stitchup.cpp - Implements the "Family Stitch-up" (distance matrix)
//                 tree construction algorithm, which works by "stitching up"
//                 a graph, based on probable familial relationships
//                 (steps 1 through 3), and then "removing the excess
//                 stitches" (in step 4).
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
//  the (degree 2) nodes that correspond to "possibilities that didn't pan out" (in Step 4)
//  (in a nutshell: the strategy is
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
//  Notes:        The union-find structure used here has a ~ n.ln(n)/ln(2)
//                worst case. And the time to remove the degree-2 nodes
//                is linear in n.
//
//  Created by James Barbetti on 12-Aug-2020 (tree construction)
//  and 24-Aug-20 (generating the newick file).
//

#include "starttree.h"
#include "distancematrix.h"
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
    std::vector<std::string>        leafNames;
    std::set< Stitch<T> >           stitches;
    std::vector< int >              taxonToSetNumber;
    std::vector< int >              taxonToNodeNumber;
    std::vector< T   >              taxonToDistance;
    std::vector< std::vector<int> > setMembers;
    int                             nodeCount;
    bool                            silent;
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
    bool areLeavesInSameSet(int leafA, int leafB) {
        return taxonToSetNumber[leafA]
                == taxonToSetNumber[leafB];
    }
    int staple(int leafA, int leafB, T length) {
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
        const char* taskDescription = silent ? "" : "Removing degree-2 nodes from stitchup graph";
        progress_display progress ( stitches.size()*2,
                                    taskDescription, "", "");
        std::vector<int> replacements;
        std::vector<T>   replacementLengths;
        replacements.reserve(nodeCount);
        replacementLengths.resize(nodeCount, 0);
        for (int i=0; i<nodeCount; ++i) {
            replacements.push_back(i);
        }
        int    node      = -1; //Source node of last edge
        size_t degree    = 0;  //Degree of that node
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
                if (it->destination < node) {
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
            T   length = it->length;
            int source = replacements[it->source];
            int dest   = replacements[it->destination];
            if (source!=dest) {
                length += replacementLengths[it->source];
                length += replacementLengths[it->destination];
                stitches.insert(Stitch<T>(source, dest, length));
            }
            ++progress;
        }
        progress.done();
    }
    template <class F>
    void dumpTreeToFile ( const std::string &treeFilePath, F& out ) const {
        int cols = 0;
        for (auto it=stitches.begin(); it!=stitches.end(); ++it) {
            if (it->source<it->destination) {
                ++cols;
                std::cout << it->source << ":" << it->destination << " " << it->length << "\t";
                if (cols==4) {
                    std::cout << std::endl;
                    cols = 0;
                }
            }
        }
        std::cout << std::endl;
    }
    template <class F>
    bool writeTreeToFile ( const std::string &treeFilePath, F& out ) const {
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
            int i = it->source;
            if (nodeToEdge[i]==edgeCount) {
                nodeToEdge[i] = j;
            }
        }
        std::string desc = "Writing STITCH tree to ";
        desc+=treeFilePath;
        progress_display progress(edgeCount, desc.c_str(), "wrote", "edge");
        out.exceptions(std::ios::failbit | std::ios::badbit);
        try {
            out.open(treeFilePath.c_str(), std::ios_base::out);
            out.precision(8);
            writeSubtree(stitchVector, nodeToEdge, nullptr, lastNodeIndex, progress, out);
            out << ";" << std::endl;
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
        progress.done();
        return true;
    }
    template <class F>
    void writeSubtree ( const std::vector<Stitch<T>> stitchVector, std::vector<size_t>  nodeToEdge,
                        const Stitch<T>* backstop, int nodeIndex,
                        progress_display &progress, F& out) const {
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
                int child = stitchVector[x].destination;
                if ( nodeToEdge[child] != y /*no backsies*/ ) {
                    out << sep;
                    sep = ",";
                    writeSubtree( stitchVector, nodeToEdge, &stitchVector[x], child, progress, out);
                }
                ++progress;
            }
            out << ")";
        }
        if (backstop!=nullptr) {
            out << ":" << backstop->length;
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

template < class T=double> class StitchupMatrix: public Matrix<T> {
public:
    typedef Matrix<T> super;
    using super::rows;
    using super::setSize;
    using super::n;
    bool silent;
    StitchupMatrix(): isOutputToBeZipped(false) {
    }
    virtual std::string getAlgorithmName() const {
        return "STITCHUP";
    }
    virtual void addCluster(const std::string &name) {
        graph.addLeaf(name);
    }
    bool loadMatrixFromFile(const std::string &distanceMatrixFilePath) {
        graph.clear();
        return loadDistanceMatrixInto(distanceMatrixFilePath, true, *this);
    }
    virtual bool loadMatrix
        ( const std::vector<std::string>& names, const double* matrix ) {
        //Assumptions: 2 < names.size(), all names distinct
        //  matrix is symmetric, with matrix[row*names.size()+col]
        //  containing the distance between taxon row and taxon col.
        setSize(names.size());
        graph.clear();
        for (auto it = names.begin(); it != names.end(); ++it) {
            addCluster(*it);
        }
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (size_t row=0; row<n; ++row) {
            const double* sourceStart = matrix + row * n;
            const double* sourceStop  = sourceStart + n;
            auto    dest        = rows[row];
            for (const double* source=sourceStart; source<sourceStop
                 ; ++source, ++dest ) {
                *dest = (T) *source;
            }
        }
        return true;
    }
    bool writeTreeFile(const std::string &treeFilePath) const {
        return graph.writeTreeFile(isOutputToBeZipped, treeFilePath);
    }
    virtual void setZippedOutput(bool zipIt) {
        isOutputToBeZipped = zipIt;
    }
    virtual void beSilent() {
        silent = true;
    }
    virtual bool constructTree() {
        if (n<3) {
            return false;
        }
        std::vector<LengthSortedStitch<T>> stitches;
        stitches.reserve(n * n);
        for (size_t row=0; row<n; ++row) {
            const T* rowData = rows[row];
            for (size_t col=0; col<row; ++col) {
                stitches.emplace_back(row, col, rowData[col]);
            }
        }
        MinHeapOnArray< LengthSortedStitch<T> >
            heap ( stitches.data(), stitches.size()
                 , silent ? "" : "Constructing min-heap of possible edges" );
        size_t iterations = 0;
        progress_display progress(0.5*n*(n+1), silent ? "" : "Assembling Stitch-up Graph");
        std::cout.precision(12);
        for (size_t join = 0; join + 1 < n; ++join) {
            LengthSortedStitch<T> shortest;
            int source;
            int dest;
            do {
                shortest = heap.pop_min();
                source   = shortest.source;
                dest     = shortest.destination;
                ++iterations;
            } while ( graph.areLeavesInSameSet(source,dest) );
            graph.staple(source, dest, shortest.length);
            progress += (join+1);
        }
        progress.done();
        graph.removeThroughThroughNodes();
        return true;
    }
protected:
    StitchupGraph<T>  graph;
    bool            isOutputToBeZipped;
};

void addStitchupTreeBuilders(Factory& f) {
    f.advertiseTreeBuilder( new Builder<StitchupMatrix<double>>
        ("STITCH",      "Family Stitch-up (Lowest Cost)"));
}

} //end of StartTree namespace.
