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
    template <class F> 
    bool writeTreeToFile(int precision, 
                         const std::string &treeFilePath, 
                         F& out) const {
        out.exceptions(std::ios::failbit | std::ios::badbit);
        try {
            out.open(treeFilePath.c_str(), std::ios_base::out);
            writeTreeToOpenFile(precision, out);
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
    template <class F> bool writeTreeToOpenFile(int precision, F& out) const {
        out.precision(precision);
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
                const Link<T> & linkBelow = cluster.links[nextChildNum];
                stack.emplace_back(here.clusterIndex, nextChildNum+1);
                stack.emplace_back(linkBelow.clusterIndex, 0);
            } else {
                out << ")";
            }
        } while (0 < stack.size());
        out << ";" << std::endl;
        return !failed;
    }
    bool writeTreeFile(bool zipIt, int precision,
                       const std::string &treeFilePath) const {
        if (treeFilePath == "STDOUT") {
            return writeTreeToOpenFile(precision, std::cout);
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
}; //end of ClusterTree class

#endif /* clustertree_h */
