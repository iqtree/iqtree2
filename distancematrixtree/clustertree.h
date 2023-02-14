//
//  clustertree.h - Classes for tracking clusters during distance matrix
//                  tree construction, and for writing out a Newick format
//                  file, based on those clusters.
//
//  Created by James Barbetti on 12-Aug-2020 (but Link, Cluster, and ClusterTree
//  were originally, from 18-Jun-2020, in bionj2.cpp).
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

#ifndef clustertree_h
#define clustertree_h

#if      USE_GZSTREAM
#include <utils/gzstream.h>          //for igzstream
#endif
#include <fstream>
#include <iostream>                  //for std::istream
#include <sstream>                   //for std::stringstream
#include <utils/progress.h>          //for progress_display
#include <utils/parallel_mergesort.h>//for MergeSorter

/**
 * @brief  Describes a link between an interior node and
 *         a cluster (clusters are identified by index).
 * @tparam T - the floating point type used to represent a distance
 *             between two taxa (or two clusters of taxa).
 */
template <class T=double> struct Link {
public:
    intptr_t clusterIndex;
    T        linkDistance;
    Link(size_t index, T distance) 
        : clusterIndex(index), linkDistance(distance) {
    }
};

/**
 * @brief  Describes a cluster (either a single exterior
 *         node (for a taxon, identified in the input distance matrix), 
 *         with an index equal to the row number that taxon had, in the 
 *         input distance matrix) with no links out from it, or an inerior
 *         node, with links to clusters that were formed earlier.
 * @tparam T - the floating point type used to represent a distance
 *             between two taxa (or two clusters of taxa).
 * @note   Tracking of the number of exterior nodes may be used to prioritize
 *         the construction of better balanced phylogenetic trees
 *         (if there are multiple pairs of clusters, that are the same 
 *          adjusted distance apart, the pair of clusters for which the
 *          absolute difference, between the exterior node counts of the 
 *          clusters, is minimized, is the one that should be chosen).
 */
template <class T=double> struct Cluster
{
public:
    size_t               countOfExteriorNodes;
    std::string          name;
    std::vector<Link<T>> links;
    Cluster(): countOfExteriorNodes(0) {
    }
    explicit Cluster(const std::string &taxon_name)
        : countOfExteriorNodes(1), name(taxon_name) {
    }
};

/**
 * @brief  A representation of a (partially or entirely) constructed
 *           phylogenetic tree. It is a vector of Cluster<T>, with 
 *           the links, between clusters (in the links member variables
 *           of the Cluster<T> instances in the vector), recording
 *           links, and distances, between subtrees.
 * @tparam T - the floating point type used to represent a distance
 *             between two taxa (or two clusters of taxa).
 */
template <class T> class ClusterTree: public std::vector<Cluster<T>>
{
    /**
     * @brief Tracks where we are up to, when we are writing out the description
     *        of a cluster, in a phylogenetic tree.
     */
    struct Place
    {
    public:
        size_t clusterIndex; 
        size_t linkNumber;   
        Place(size_t ix, size_t num) {
            clusterIndex = ix;
            linkNumber = num;
        }
    };

public:
    typedef std::vector<Cluster<T>> super;
    using super::at;
    using super::back;
    using super::emplace_back;
    using super::push_back;
    using super::size;

    /**
     * @brief  Add a leaf node (for an input taxon) to a phylogenetic tree
     * @param  taxon_name The taxon number
     * @return Cluster<T>& a (non-const) reference to the newly added Cluster<T>
     */
    Cluster<T>& addCluster(const std::string& taxon_name) {
        emplace_back(taxon_name);
        return back();
    }
    /**
     * @brief  Add a degree-2 node (or rather, a cluster, representing the
     *         subtree beneath that node), to an existing (partially-assembled)
     *         phylogenetic tree; by joining two existing clusters.
     * @param  a       - the index of the first cluster to be joined.
     * @param  aLength - the length, from the new degree-2 node, to the
     *                   first cluster.
     * @param  b       - the index of the second cluster to be joined.
     * @param  bLength - the length, from the new degree-2 node, to the
     *                   second cluster.
     * @return Cluster<T>& 
     * @note   It is assumed that 0<=a<b<size(), that the cluster indices
     *         are both valid, and are distinct.
     * @note   There is no check that a<b.
     * @note   It is assumed that there aren't any other existing clusters
     *         that have links "down" to the subtrees for clusters a or b.
     */
    Cluster<T>& addCluster(size_t a, T aLength, size_t b, T bLength) {
        push_back(Cluster<T>());
        Cluster<T>& cluster = back();
        cluster.links.emplace_back(a, aLength);
        cluster.links.emplace_back(b, bLength);
        cluster.countOfExteriorNodes = at(a).countOfExteriorNodes 
                                     + at(b).countOfExteriorNodes;
        //std::cout << "Cluster " << (size()-1)
        //<< " is (" << a << ", " << b << ")" << std::endl;
        return cluster;
    }
    /**
     * @brief  Add a degree-3 node (or rather, a cluster, representing the
     *         subtree beneath that node), to an existing (partially-assembled)
     *         phylogenetic tree; by joining THREE existing clusters.
     *         In practice this is only called to join the last three clusters,
     *         in an unrooted phylogenetic tree.
     * @param  a       - the index of the first cluster to be joined.
     * @param  aLength - the length, from the new degree-3 node, to the
     *                   first cluster.
     * @param  b       - the index of the second cluster to be joined.
     * @param  bLength - the length, from the new degree-3 node, to the
     *                   second cluster.
     * @param  c       - the index of the third cluster to be joined.
     * @param  cLength - the length, from the new degree-3 node, to the
     *                   third cluster.
     * @return Cluster<T>& 
     * @note   It is assumed that 0<=a<b<c<size(), that the cluster indices
     *         are all valid, and are distinct.
     * @note   There is no check that a<b<c.
     * @note   It is assumed that there aren't any other existing clusters
     *         that have links "down" to the subtrees for clusters a, b, or c.
     */    
    Cluster<T>& addCluster
    ( size_t a, T aLength, size_t b, T bLength
     , size_t c, T cLength)  {
        Cluster<T>& cluster = addCluster(a, aLength, b, bLength);
        cluster.links.emplace_back(c, cLength);
        cluster.countOfExteriorNodes += at(c).countOfExteriorNodes;
        //std::cout << "Final cluster " << (size()-1)
        //<< " is (" << a << ", " << b << ", " << c << ")" << std::endl;
        return cluster;
    }
    /**
     * @brief  "Append" a cluster, c, to the most recently created cluster, x
     *         (add a link, from cluster x, to cluster c, with the specified
     *         length).
     * @param  c      - the cluster index of the cluster to "append".
     * @param  length - the length of the link from the last cluster, x, 
     *                  to cluster c.
     * @return Cluster<T>& 
     * @note   It is assumed that 0 <= c < (size()-1).  But this is not checked.
     * @note   It is assumed that there are no downward links to c.  But this
     *         is not checked, either.
     */
    Cluster<T>& appendToLastCluster
    ( size_t c, T length) {
        back().links.emplace_back(c, length);
        back().countOfExteriorNodes += at(c).countOfExteriorNodes;
        return back();
    }
    typedef std::pair<intptr_t /*taxon*/, T /*distance*/ >     LeafDistance;
    typedef std::vector<LeafDistance> LeafDistanceVector;
    
    /**
     * @brief Calculate, for a given cluster, the sums of the lengths of
     *        the path, leading to each of the exterior clusters (the 
     *        taxon nodes), and return a vector of leaf distances.
     * @param top 
     * @param branch_length 
     * @param ldv a vector of LeafDistance (output)
     * @param ldv will be sorted by leaf cluster index.
     * @note  It is assumed (but not checked) that ldv is empty on input.
     * @note  The run-time complexity of this implementation is
     *        (guaranteed to be) O(n.log_2(n)), because a mergesort
     *        is used to sort by distances.
     */
    void calculateDistancesToLeaves(intptr_t top, double branch_length,
                                    LeafDistanceVector& ldv) {
        //Running time proportional to subtree size
        std::vector<LeafDistance> stack;
        LeafDistance top_node;
        top_node.first  = top;
        top_node.second = branch_length;
        stack.emplace_back(top_node);
        while (!stack.empty()) {
            LeafDistance b = stack.back();
            stack.pop_back();
            Cluster<T>& node = at(b.first);
            if (node.links.empty()) {
                ldv.emplace_back(b);
                continue;
            } 
            for (Link<T>& link : node.links ) {
                stack.emplace_back(link.clusterIndex, 
                                   b.second + link.linkDistance );
            }
        }
        //
        //Hack: Reordering makes memory accesses "cache-friendlier"
        //      for large clusters, in calculateRMSOfTMinusD().
        //Note: Single-threaded sorting is used, here, because 
        //      calculateRMSOfTMinusD() is already executing 
        //      calls to this function in parallel.
        //
        MergeSorter<LeafDistance> sorter;
        sorter.single_thread_sort(ldv.data(), ldv.size());
    }
    /**
     * @brief  Calculate the root mean square of the difference between the
     *         inter-taxon distances, implicit in the tree, and those made
     *         explicit in an input distance matrix.
     * @param  matrix - a pointer to the first element of a flat distance matrix,
     *                  recording (input) distances between taxa, in row-major 
     *                  order: matrix[r*rank+c] is the distance between taxa r
     *                  and c.
     * @param  rank   - the rank of the distance matrix
     * @param  rms    - the root mean square, of the differences, between the
     *                  distances, between taxa, according to (i) the tree,
     *                  and (ii) the distance matrix provided in (matrix).
     * @return true - because it always works (unless it throws?! Out of memory?!)
     * @note   Total running time: asymptotically proportional to 
     *         rank*(rank-1)/2 (that's the number of additions 
     *         to sum_of_squares)
     * @note   Total memory consumption:  depends on thread count, X.
     *         In theory, on the order of 
     *         (X+1)*rank*(sizeof(intptr_t)+sizeof(T)):
     *         2 * X * (maximum sum of sizes of leaf distance vectors).
     *         The 2 is because of the use of a MergeSorter in 
     *         calculateDistancesToLeaves().  An in-place sort would 
     *         lower that to a 1.
     * @note   It is assumed that 3<=rank.
     */
    bool calculateRMSOfTMinusD(const double* matrix, intptr_t rank, double& rms) {
        intptr_t cluster_count  = size();
        double   sum_of_squares = 0.0;
        #ifdef _OPENMP
        #pragma omp parallel for reduction(+:sum_of_squares)
        #endif
        for (intptr_t h=rank; h<cluster_count; ++h) {
            //For each (non-leaf) cluster...
            //1. calculate distances to all the 
            //   leaves, for each contributing cluster.
            const Cluster<T>& node = at(h);
            std::vector<LeafDistanceVector> subtrees;
            for (const Link<T>& link : node.links ) {
                LeafDistanceVector ldv;
                calculateDistancesToLeaves(link.clusterIndex, 
                                           link.linkDistance, ldv);
                subtrees.emplace_back(ldv);
            }
            //2. for each pair of LeafDistanceVectors, A, B,
            //   (for separate contributing clusters)...
            size_t subtree_count = subtrees.size();
            for (size_t i=0; i+1<subtree_count; ++i) {
                for (size_t j=i+1; j<subtree_count; ++j) {
                    //2b. for each leaf (a.first) from A
                    for (const LeafDistance& a: subtrees[i]) {
                        auto row = matrix + a.first * rank;
                        //2c. for each leaf (b.first) from B
                        //    calculate error; difference between
                        //    distance(a.first) + distance(b.first)
                        //    and D[a.first * rank + b.first].
                        for (const LeafDistance& b: subtrees[j]) {
                            double diff     = a.second + b.second
                                            - row[b.first];
                            sum_of_squares += (diff * diff);
                            #if (0)
                            std::cout << a.first << " " << b.first << " "
                                      << a.second + b.second << " "
                                      << row[b.first] << " "
                                      << diff << "\n";
                            #endif
                        }
                    }
                }
            }
            //std::cout << "\n";
        }
        double double_rank = static_cast<double>(rank);
        rms = sqrt( sum_of_squares * 2.0 / double_rank / (double_rank - 1.0) );
        //std::cout << "rank " << (double_rank) << std::endl;
        //std::cout << "sum " << sum_of_squares << ","
        //          << " divisor " << double_rank*(double_rank-1)*0.5 << std::endl;
        return true;
    }

    /**
     * @brief  write the phylogenetic tree (implicit in the Cluster<T>
     *         instances in the vector), to a file, with the specified 
     *         file path, using a stream of the specified type, F.
     * @tparam F - the stream type to use
     * @param  precision    - the number of digits of precision to
     *                        use when outputting genetic distances
     * @param  treeFilePath - the path of the file to write to
     * @param  isOutputToBeAppended - true if the file is to be appended
     * @param  isSubtreeOnly        - true if the leading ( and trailing );
     *                                are to be omitted, false otherwise
     * @param  out   - reference to a stream instance to use
     * @return true  - if the file can be written to
     * @return false - if it can't (error messages saying why not will
     *                 be written to std::err before false is returned).
     */
    template <class F> 
    bool writeTreeToFile(int precision, 
                         const std::string &treeFilePath,
                         bool isOutputToBeAppended,
                         bool isSubtreeOnly,
                         F& out) const {
        out.exceptions(std::ios::failbit | std::ios::badbit);
        try {
            std::ios_base::openmode openMode = isOutputToBeAppended
                          ? std::ios_base::app : std::ios_base::trunc;
            openMode |= std::ios_base::out;  
            out.open(treeFilePath.c_str(), openMode );
            out.precision(precision);
            writeTreeToOpenFile(isSubtreeOnly, out);
            out.close();
            return true;
        } catch (std::ios::failure &) {
            std::cerr << "IO error"
            << " opening/writing file: " << treeFilePath << std::endl;
            return false;
        } catch (const char *str) {
            std::cerr << "Writing newick file failed: " << str << std::endl;
            return false;
        } catch (const std::string &str) {
            std::cerr << "Writing newick file failed: " << str << std::endl;
            return false;
        }
    }
    bool isRootCluster(const Place& here) const {
        return here.clusterIndex + 1 == size();
    }
    /**
     * @brief  append the phylogenetic tree (implicit in the Cluster<T>
     *         instances in the vector), to a n output stream.
     * @tparam F - the type of the outputstream
     * @param  isSubtreeOnly - true if the leading ( and trailing );
     *                         are to be omitted, false otherwise
     * @param  out           - the output stream
     * @return true  - if the tree can be written
     * @return false - if there is an error (error messages will be
     *                 written to std::err, before false is returned).
     * @note   Cycles (where a cluster is an ancestor of itself),
     *         and reticulations, where the same cluster appears in
     *         more than one subtree) are not detected (as such).  But
     *         if the number of nodes that have been written to (out)
     *         indicates that there is (either) a cycle, or a 
     *         reticulation in the tree, this function will error out.
     */
    template <class F> bool writeTreeToOpenFile
        ( bool isSubtreeOnly, F& out) const {
        std::vector<Place> stack;
        bool failed = false; //Becomes true if clusters
        //defines cycles (should never happen)
        //Indicates a fatal logic error
        size_t maxLoop = 3 * size();
        //More than this, and there must be
        //a cycle.  Or something.
        
        ASSERT(0<size());
        stack.emplace_back(size()-1, 0);
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
                if (!isSubtreeOnly || !isRootCluster(here))
                {
                    out << "(";
                }
                stack.emplace_back(here.clusterIndex, 1);
                stack.emplace_back(cluster.links[0].clusterIndex, 0);
                continue;
            }
            size_t nextChildNum = here.linkNumber;
            const Link<T> & linkPrev = cluster.links[nextChildNum-1];
            out << ":" << linkPrev.linkDistance;
            if (nextChildNum<cluster.links.size()) {
                out << ",";
                const Link<T> & linkBelow = cluster.links[nextChildNum];
                stack.emplace_back(here.clusterIndex, nextChildNum+1);
                stack.emplace_back(linkBelow.clusterIndex, 0);
            } else if (!isSubtreeOnly || !isRootCluster(here)) {
                out << ")";
            }
        } while (0 < stack.size());
        if (!isSubtreeOnly) {
            out << ";" << std::endl;
        }
        return !failed;
    }
    /**
     * @brief write the phylogenetic tree (implicit in the Cluster<T>
     *        instances in the vector), to a file, with the specified 
     *        file path (or to standard output - see below).
     * @param zipIt        - indicates whether the output should be compressed
     *                       (with gzip compression). This is *only* honoured
     *                       if the USE_GZSTREAM symbol is defined and non-zero.
     * @param precision    - the number of digits of precision to use for each
     *                       genetic distance.
     * @param treeFilePath - either "STDOUT" (write to standard output), or
     *                       the file path, to which the tree is to be written
     * @param isOutputToBeAppended - true if the file is (if it exists) to be 
     *                       appended, false if the file is to be truncated.
     * @param subtreeOnly  - true, if the leading "(" and trailing ");" are to
     *                       be suppressed, in the newick output.
     * @return true  - on success
     * @return false - on failure (error messages will be written to std::cerr)
     */
    bool writeTreeFile(bool zipIt, int precision,
                       const std::string &treeFilePath,
                       bool isOutputToBeAppended,
                       bool subtreeOnly) const {
        if (treeFilePath == "STDOUT") {
            std::cout.precision(precision);
            return writeTreeToOpenFile(subtreeOnly, std::cout );
        } else if (zipIt) {
            #if USE_GZSTREAM
            ogzstream     out;
            #else
            std::ofstream out;
            #endif
            return writeTreeToFile(precision, treeFilePath, 
                                   isOutputToBeAppended, 
                                   subtreeOnly, out);
        } else {
            std::fstream out;
            return writeTreeToFile(precision, treeFilePath, 
                                   isOutputToBeAppended, 
                                   subtreeOnly, out);
        }
    }
}; //end of ClusterTree class

#endif /* clustertree_h */
