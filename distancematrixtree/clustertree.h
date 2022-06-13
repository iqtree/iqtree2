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

template <class T=double> struct Link {
    //
    //Describes a link between an interior node and
    //a cluster (clusters are identified by index).
    //
public:
    size_t  clusterIndex;
    T       linkDistance;
    Link(size_t index, T distance) 
        : clusterIndex(index), linkDistance(distance) {
    }
};

template <class T=double> struct Cluster
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
    explicit Cluster(const std::string &taxon_name)
        : countOfExteriorNodes(1), name(taxon_name) {
    }
};

template <class T> class ClusterTree: public std::vector<Cluster<T>>
{
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
        cluster.countOfExteriorNodes = at(a).countOfExteriorNodes 
                                     + at(b).countOfExteriorNodes;
        //std::cout << "Cluster " << (size()-1)
        //<< " is (" << a << ", " << b << ")" << std::endl;
        return cluster;
    }
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
    Cluster<T>& appendToLastCluster
    ( size_t c, T length) {
        back().links.emplace_back(c, length);
        back().countOfExteriorNodes += at(c).countOfExteriorNodes;
        return back();
    }
    typedef std::pair<intptr_t /*taxon*/, T /*distance*/ >     LeafDistance;
    typedef std::vector<LeafDistance> LeafDistanceVector;
    void calculateDistancesToLeaves(intptr_t top, double branch_length,
                                    LeafDistanceVector& ldv) {
        //Running time proportional to subtree size
        std::vector<LeafDistance> stack;
        stack.emplace_back(top, branch_length);
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
    bool calculateRMSOfTMinusD(const double* matrix, intptr_t rank, double& rms) {
        //Assumes: rank is at least 3.
        //
        //Total running time: proportional to rank*(rank-1)/2.
        //(that's the number of additions to sum_of_squares).
        //
        //Total memory consumption:  depends on thread count, X.
        //In theory, on the order of (X+1)*rank*(sizeof(intptr_t)+sizeof(T)).
        //2 * X * (maximum sum of sizes of leaf distance vectors).
        //The 2 is because of the use of a MergeSorter in 
        //calculateDistancesToLeaves().  An in-place sort would 
        //lower that to a 1.
        //
        intptr_t cluster_count  = size();
        double   sum_of_squares = 0.0;
        #ifdef _OPENMP
        #pragma omp parallel for reduction(+:sum_of_squares)
        #endif
        for (intptr_t h=rank; h<cluster_count; ++h) {
            //For each (non-leaf) cluster...
            //1. calculate distances to all the 
            //   leaves, for each contributing cluster.
            Cluster<T>& node = at(h);
            std::vector<LeafDistanceVector> subtrees;
            for (Link<T>& link : node.links ) {
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
                    for (LeafDistance& a: subtrees[i]) {
                        auto row = matrix + a.first * rank;
                        //2c. for each leaf (b.first) from B
                        //    calculate error; difference between
                        //    distance(a.first) + distance(b.first)
                        //    and D[a.first * rank + b.first].
                        for (LeafDistance& b: subtrees[j]) {
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

    template <class F> 
    bool writeTreeToFile(int precision, 
                         const std::string &treeFilePath,
                         bool isOutputToBeAppended,
                         bool isSubtreeOnly,
                         F& out) const {
        out.exceptions(std::ios::failbit | std::ios::badbit);
        try {
            auto openMode = isOutputToBeAppended
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
        } catch (std::string &str) {
            std::cerr << "Writing newick file failed: " << str << std::endl;
            return false;
        }
    }
    bool isRootCluster(const Place& here) const {
        return here.clusterIndex + 1 == size();
    }
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
