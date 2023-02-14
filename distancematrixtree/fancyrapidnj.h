//
//  fancyrapidnj.h - Template classes for a more efficient version of 
//                   RapidNJ.  Unlike NJMatrix, which maintains D, S 
//                   and I matrices, it maintains per-cluster rows
//                   of (cluster number, distance pairs)
//                   (S+I) - distances to earlier clusters
//                           ordered by ascending distance.
//                   (D+I) - distances to earlier clusters
//                           ordereed by ascending cluster number.
//                   Clustering merges (S+I) *distance-sorted* information,
//                   for the clusters (x, and y) being merged, and 
//                   "looks up" the cluster-number-sorted information,
//                   from clusters numbered "higher" than x and y
//                   by interpolation searching their (D+I) rows.
//
//                   It avoids *writing* (S+I) and (D+I) memory,
//                   except for the next cluster's (S+I) and (D+I) rows.
//
//  See:  https://birc.au.dk/software/rapidnj/.
//        Paper: "Inference of Large Phylogenies using Neighbour-Joining."
//               Martin Simonsen, Thomas Mailund, Christian N. S. Pedersen.
//               Communications in Computer and Information Science
//               (Biomedical Engineering Systems and Technologies:
//               3rd International Joint Conference, BIOSTEC 2010,
//               Revised Selected Papers), volume 127, pages 334-344,
//               Springer Verlag, 2011. [Simonsen+Mailund+Pedersen]
//
//  FancyNJMatrix<T> differs from BoundingMatrix<T, NJMatrix<T>> in several
//  respects, because it implemented with space efficiency performance in mind 
//  (BoundingMatrix was NOT; it was implemented to be easy to read, easy
//  to "map" back to the paper it is based on).
//
//  Copyright James Barbetti (2021-22)
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

#pragma  once
#include "upgma.h"
#include "utils/gzstream.h"
#include "utils/timeutil.h"
#ifndef  fancy_rapid_nj_h
#define  fancy_rapid_nj_h

#include "nj.h"                       //for NJFloat (should be separate header for that)
#include "clustertree.h"              //for ClusterTree template class
#include "flatmatrix.h"               //for FlatMatrix
#include "hashrow.h"                  //for HashRow
#include "utils/parallel_mergesort.h" //for MergeSorter

#define  FNJ_TRACE(x) {}

namespace StartTree {
/**
 * @brief  An alternative implementation of the Rapid NJ algorithm
 *         (somewhat closer to the [Simonsen+Mailund+Pedersen] one,
 *         with { cluster, distance } entries in a matrix vector 
 *         of row vectors (rather than a single flattened out array).
 * @tparam T 
 */
template <class T=NJFloat> class FancyNJMatrix {
protected:
    bool           be_silent;      //true if log messages are to be suppressed
    bool           zip_it;         //true if output is to be compressed
    bool           append_file;    //true if the output file is to be appended
                                   //(false if it is to be truncated)
    bool           is_rooted;      //true if the "top" node is to be degree 2
                                   //(false if it to be dgree 3)
    bool           omit_semicolon; //true if the semi-colon at the end of the
                                   //output is to be omitted (e.g. when 
                                   //calculating a subtree rather than a tree).
    ClusterTree<T> clusters;
    #if USE_GZSTREAM
        #if USE_PROGRESS_DISPLAY
            typedef pigzstream InFile;
            #define INFILE(name) InFile name("matrix")
        #else
            typedef igzstream  InFile;
            #define INFILE(name) Infile name
        #endif
    #else
        typedef std::ifstream InFile;
        #define INFILE(name) InFile name
    #endif
    typedef std::vector<T>    DistanceVector;
    struct MatrixEntry {
        public:
            T   distance;
            int cluster_num;
            bool operator < (const MatrixEntry &rhs) const {
                return  distance    <  rhs.distance
                    || (distance    == rhs.distance && 
                        cluster_num <  rhs.cluster_num);
            }
            bool operator <= (const MatrixEntry& rhs) const {
                return  distance    <  rhs.distance
                    || (distance    == rhs.distance && 
                        cluster_num <= rhs.cluster_num);
            }
    };
    typedef std::vector<MatrixEntry>  EntryVector;
    typedef std::vector<MatrixEntry*> EntryPtrVector;

    intptr_t       original_rank;        //the rank of the matrix before any
                                         //clusters were joined
    intptr_t       next_cluster_number;  //the cluster number to use next
    IntVector      cluster_in_play;      //indicates which clusters are in use
                                         //(it's an IntVector, rather than a vector
                                         //of bool, because std::vector<bool> tends
                                         //*not* to be multi-thread-friendly).
                                         //For each cluster:
    DistanceVector cluster_total;        // - The total of distances to other clsuters
    DistanceVector cluster_total_scaled; // - The same, scaled (divided by (n-2))
    DistanceVector cluster_cutoff;
    EntryPtrVector cluster_sorted_start; 
    EntryPtrVector cluster_sorted_stop;
    EntryPtrVector cluster_unsorted_start;
    EntryPtrVector cluster_unsorted_stop;
    IntVector      cluster_row;            //maps a cluster number to a row number

                                  //when searching for clusters to merge:...
    IntVector      row_cluster;   //cluster number, y, of i(th) clusters still in play
    DistanceVector row_raw_dist;  //raw distance, Dxy, between cluster y and
                                  //the cluster x, closest to it (according to Dxy).
    DistanceVector row_best_dist; //Dxy - Rx - Ry distance for same
    IntVector      row_choice;    //Indicates which cluster, x<y, had 
                                  //the best Dxy - Rx - Ry distance, for 
                                  //cluster y.

                                  //when merging clusters:...
    DistanceVector x_distances;   //distances to cluster x
    DistanceVector y_distances;   //distances to cluster y

    int            threadCount;               //the number of threads
    typedef MergeSorter<MatrixEntry> Sorter;  //used for sorting matrix rows
    std::vector<Sorter> sorters;              //mergesorting contexts (one per thread)
                                              //(mergesorting requires an auxiliary array,
                                              // and allocation of those arrays is amortized)

    volatile T     global_best_dist;          //best adusted distance found in curent iteration

    /**
    * @brief  Get the thread number of the current thread (e.g. for looking)
    *         up the entry in sorters, that keeps track of the auxiliary arrays
    *         allocated for the current thread.
    * @return int - either the current thread number or (if not multithreading), 0.
    */
    int getThreadNumber() const {
        #ifdef _OPENMP
            return omp_get_thread_num();
        #else
            return 0;
        #endif
    }

    /**
     * @brief  return the number of threads of execution
     * @return int - the number of threads (always 1 if _OPENMP is not defined)
     */
    int getThreadCount() const {
        #ifdef _OPENMP
            return omp_get_num_threads();
        #else
            return 1;
        #endif
    }

    DuplicateTaxa duplicate_taxa;

public:
    FancyNJMatrix() : be_silent(false), zip_it(false), 
                      append_file(0), is_rooted(false), 
                      omit_semicolon(false), original_rank(0), 
                      next_cluster_number(0),  global_best_dist(0) {
        #ifdef _OPENMP
            threadCount = omp_get_max_threads();
        #else
            threadCount = 1;
        #endif
        sorters.resize(threadCount);
    }

    ~FancyNJMatrix() {
        for (int c=0; c<original_rank; ++c) {
            deallocateCluster(c);
        }
    }

    std::string getAlgorithmName() const {
        return "FancyNJ";
    }

    /**
     * @brief Turn off logging messages.
     */
    void beSilent() { 
        be_silent = true; 
    }

    /**
     * @brief  Control whether output files are to be appended (true)
     *         or overwritten (false).
     * @param  appendIt - whether output files are to be appended
     * @return true - always succeeds (this implementation appends when asked to)
     */
    virtual bool setAppendFile(bool appendIt) {
        append_file = appendIt;
        return true;
    }

    /**
     * @brief  Control whether output files are to be compressed (true)
     *         or not (false).
     * @param  appendIt - whether output files are to be compressed
     * @return true - always succeeds (this implementation supports 
     *         compression)
     */
    virtual bool setZippedOutput(bool zipIt) { 
        zip_it = zipIt;
        return true;
    }

    /**
     * @brief  load sequence names and a distance matrix, from the file
     *         matching the specified file path
     * @param  distanceMatrixFilePath - the file path
     * @return true  - on success
     * @return false - on failure (error message will be logged to std::cerr)
     */
    virtual bool loadMatrixFromFile(const std::string &distanceMatrixFilePath) {
        INFILE(in);
        try {
            in.exceptions(std::ios::failbit | std::ios::badbit);
            in.open(distanceMatrixFilePath.c_str(), std::ios_base::in);
            bool rc = loadMatrixFromOpenFile(in);
            in.close();
            return rc;
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
    }

    /**
     * @brief  Load sequence names and a distance matrix from an open file
     * @param  in the file to load it from 
     * @return true  - if the load was successful
     * @return false - if there was an error (error message logged to std::cerr)
     */
    virtual bool loadMatrixFromOpenFile(InFile& in) {
        FlatMatrix dummy;
        loadDistanceMatrixFromOpenFile(in, !be_silent, dummy);
        return this->loadMatrix(dummy.getSequenceNames(), dummy.getDistanceMatrix());
    }

    /**
     * @brief  Load a vector of sequence names, and a matrix of distances
     * @param  names   a vector of n (unique!) sequence names
     * @param  matrix  a pointer to the first element of a flat matrix of 
     *                 size n*n (the distances between the sequences,
     *                 stored in row-major order).
     * @return true  - on success
     * @return false - on failure
     * @note   loading of data from (matrix) is parallelized, over rows,
     *         if _OPENMP is defined.
     */
    virtual bool loadMatrix(const StrVector& names,
                            const double* matrix) {
        #if USE_PROGRESS_DISPLAY
        double row_count  = names.size();
        double work_to_do = 0.5 * row_count * (row_count - 1.0);
        const char* taskName = be_silent
            ? "" :  "Setting up triangular I+S and I+D matrices";
        progress_display setupProgress(work_to_do, taskName, "", "");
        #endif
        //Assumes: matrix is symmetrical
        clusters.clear();
        for (const std::string& name : names) {
            clusters.addCluster(name);
        }
        setRank(names.size()); //sets original_rank

        //
        //A weird thing here. Visual Studio C++ doesn't like 0<=r.
        //But it is perfectly happy with r>=0.  That's a bug,
        //but it can't be helped.
        //
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) if(1<threadCount)
        #endif
        for (intptr_t r=original_rank-1; r>=0; --r ) {
            const double* row       = matrix + r * original_rank;
            MatrixEntry*  data      = cluster_sorted_start[r];
            MatrixEntry*  unsorted  = cluster_unsorted_start[r];
            T             total     = 0;
            FNJ_TRACE("Z" << r << "...");
            for (int c=0; c<r; ++c) {
                data[c].distance    = row[c];
                data[c].cluster_num = c;
                unsorted[c]         = data[c];
                total              += row[c];
                FNJ_TRACE(" " << row[c]);
            }
            FNJ_TRACE("X");
            for (int c=r+1; c<original_rank; ++c) {
                total              += row[c];
                FNJ_TRACE(" " << row[c]);
            }
            FNJ_TRACE(" total " << total << "\n");
            cluster_total[r] = total;

            //formerly: std::sort(data, data+r);
            sorters[getThreadNumber()].single_thread_sort(data, r);

            #if USE_PROGRESS_DISPLAY
                setupProgress += ((double)r-1);
            #endif
        }        
        #if USE_PROGRESS_DISPLAY
            setupProgress.done();
        #endif

        identifyDuplicateTaxa(matrix);
        return true;
    }
    /**
     * @brief Hook for code that must execute before a tree is constructed
     *        (there is any, so this does nothing at all)
     */
    virtual void prepareToConstructTree() {
    }

    /**
     * @brief  Indicate whether the tree, to be inferred, is to have a top node
     *         with a degree of 2 (if true), or of 3 (if false)
     * @param  rootIt true, for a subtree (top node degree 2), false otherwise
     * @return true   - if last/top "last-joined" node's degree is to be 2
     * @return false  - if last node's degree is to be 3
     */
    virtual bool setIsRooted(bool rootIt) {
        is_rooted = rootIt;
        return true;
    }

    /**
     * @brief  Indicate whether output is to include ( and ) characters
     *         around the output for the tree or subtree (true if no '('
     *         or ')' characters are wanted, false otherwise).
     * @param  wantSubtree true if '(' and ')' characters to be dropped
     * @return true  - this implementation supports outputting subtrees.
     */
    virtual bool setSubtreeOnly(bool wantSubtree) {
        omit_semicolon = wantSubtree;
        return true;
    }

    /**
     * @brief  construct a tree via phylogenetic inference, using the
     *         Rapid NJ algorithm.
     * @return true  - on success (always!)
     * @return false - on failure (never!)
     * @note   periodically, items in per-cluster rows, that refer to
     *         clusters that are "no longer in play" (that have already
     *         been joined into larger clusters), are stripped out.
     *         the code for this is tagged PURGE.
     */
    virtual bool constructTree() {
        prepareToConstructTree();
        if (original_rank<3) {
            return false;
        }

        int  n = original_rank;
        int  q = n + n - 2;

        double purgeTime   = 0.0;
        double recalcTime  = 0.0;
        double previewTime = 0.0;
        double mergeTime   = 0.0;
        int    next_purge  = n * 7 / 8;

        intptr_t duplicate_merges = clusterDuplicateTaxa 
                                    (n, recalcTime, mergeTime);

        #if USE_PROGRESS_DISPLAY
        std::string taskName;
        if (!be_silent) {
            taskName = "Constructing " + getAlgorithmName() + " tree";
        }
        double triangle  = (double)n * ((double)n + 1.0) * 0.5;
        progress_display show_progress(triangle, taskName.c_str(), "", "");
        #endif

        for ( ; 1 < n ; --n) {
            if (n<=next_purge) {
                //PURGE
                purgeTime -= getRealTime();
                removeOutOfPlayClusters(next_cluster_number);
                next_purge = n*7/8;
                if (next_purge<100) {
                    next_purge = 0;
                }
                purgeTime += getRealTime();
            }
            recalcTime  -= getRealTime();
            recalculateTotals    (n, next_cluster_number);
            recalcTime  += getRealTime();
            previewTime -= getRealTime();
            previewRows          (n, q);
            previewTime += getRealTime();
            findPreferredPartners(n, q);
            int best_row      = chooseBestRow(n, q);
            int best_cluster  = row_cluster[best_row];
            int other_cluster = row_choice[best_row];
            int other_row     = cluster_row[other_cluster];
            int low_row       = (best_row < other_row) ? best_row  : other_row;
            int high_row      = (best_row < other_row) ? other_row : best_row;
            T   raw_Dxy       = clusterDistance(other_cluster, best_cluster);
            row_cluster[low_row]  = next_cluster_number;
            int n_less_1          = n - 1;
            row_cluster[high_row] = row_cluster[n_less_1];
            mergeTime -= getRealTime();
            mergeClusters(best_cluster, other_cluster, next_cluster_number,
                          raw_Dxy, n, is_rooted);
            mergeTime += getRealTime();
            ++next_cluster_number;
            #if USE_PROGRESS_DISPLAY
            show_progress+=(double)n;
            #endif
        }

        #if USE_PROGRESS_DISPLAY
        show_progress.done();
        #endif

        reportConstructionDone(duplicate_merges, purgeTime, recalcTime, previewTime, mergeTime);
        return true;
    }

    /**
     * @brief  Calculate the root-mean-square of the matrix of 
     *         differences between the tree-distances (between taxa) 
     *         implied by the edge lengths in (clusters), and the
     *         input distances (in a distance matrix, which should be
     *         a copy of the input that was provided to the phylogenetic
     *         inference algorithm).
     * @param  matrix a flat matrix (rank*rank) of doubles in row-major 
     *                order.
     * @param  rank   the rank of the flat matrix.
     * @param  rms    the root mean square will be written here.
     * @return true if the calculation succeeded
     * @note   it is assumed that rank is equal to the number of leaf
     *         taxa in clusters, and that the rows and columns in
     *         the distance matrix (pointed to by matrix) correspond
     *         one-to-one and in order to the first (rank) clusters.
     * @note   it is assumed that the distance matrix is symmetric.
     *         it is assumed that the diagnal entries of the distance 
     *         matrix are all zeroes.
     */
    virtual bool calculateRMSOfTMinusD(const double* matrix, 
                                       intptr_t rank, double& rms) {
        return clusters.calculateRMSOfTMinusD(matrix, rank, rms);
    }

    /**
     * @brief  Write the current tree, in newick format to the file
     *         matching the supplied file path.
     * @param  precision the number of digits after the decimal point
     *                   in distances, reported between nodes in the tree.
     * @param  treeFilePath the file path to write to
     * @return true  - if the write succeeds
     * @return false - if the write fails
     * @note   zip_it (compression on/off), append_file (append rather than overwrite,
     *         and omit_semicolon (subtree only, no leading '(', trailing ");") 
     *         are all honoured).
     */
    virtual bool writeTreeFile     (int precision,
                                    const std::string &treeFilePath) const { 
        return clusters.writeTreeFile
               ( zip_it, precision, treeFilePath
               , append_file, omit_semicolon );
    }
    /**
     * @brief  Append the current tree, in newick format to an open I/O stream.
     * @param  stream the output stream (std::ostream) to write to.
     * @return true  - if the write succeeds
     * @return false - if the write fails
     * @note   zip_it (compression on/off), append_file (append rather than overwrite,
     *         and omit_semicolon (subtree only, no leading '(', trailing ");") 
     *         are all honoured).
     * @note   it is expected that the caller will have set the precision already.
     */
    virtual bool writeTreeToOpenFile(std::iostream &stream) const { 
        return clusters.writeTreeToOpenFile
               ( omit_semicolon, stream );
    }

protected:
    /**
     * @brief Set the number of leaf taxa (or sequences)
     * @param n the number of leaf taxa
     */
    virtual void setRank(size_t n) {
        original_rank       = n;
        next_cluster_number = n;
        size_t q = n+n-1; //number of clusters that will be needed
                          //in total, during the course of the tree
                          //construction: n leaves, and n-1 interiors
                          //(at the very end of processing, for an 
                          //unrooted tree, the last two clusters 
                          //will be merged; so n-1 not n-2), see
                          //the appendToLastCluster call, near the
                          //end of mergeClusters().
                          //
        cluster_in_play.resize       (q, 1);
        cluster_total.resize         (q, 0.0);
        cluster_total_scaled.resize  (q, 0.0);
        cluster_sorted_start.resize  (q, nullptr);
        cluster_sorted_stop.resize   (q, nullptr);
        cluster_unsorted_start.resize(q, nullptr);
        cluster_unsorted_stop.resize (q, nullptr);
        cluster_row.resize           (q, -1);
        cluster_cutoff.resize        (q, 0.0);
        for (size_t c=0; c<n; ++c) {
            allocateCluster(c);
        }
        row_cluster.resize(original_rank);
        row_raw_dist.resize(original_rank);
        row_best_dist.resize(original_rank);
        row_choice.resize(original_rank);

        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (int i=0; i<original_rank; ++i) {
            row_cluster[i] = i;
        }
    }

    /**
     * @brief Allocate resources (in particular a MatrixEntry row)
     *        for a cluster.  The matrix entry row will be large 
     *        enough for TWO copies of the entries needed for the
     *        cluster with this index ("just big enough").
     * @param c the cluster index
     */
    void allocateCluster(int c) {
        ASSERT(cluster_sorted_start[c] == nullptr);
        size_t       q             = original_rank + original_rank - 2;
        //Cluster, c, 0 through n-1, has c previous clusters
        //When     c, n through q-1, is created, it'll have q-c
        size_t       entries       = (c<original_rank) ? c : ( q - c );
        MatrixEntry* data          = new MatrixEntry[entries*2];
        cluster_sorted_start[c]    = data;
        data                      += entries;
        cluster_sorted_stop[c]     = data;
        cluster_unsorted_start[c]  = data;
        data                      += entries;
        cluster_unsorted_stop[c]   = data;
    }

    /**
     * @brief  identify duplicate taxa.  Record equivalence classes
     *         in the duplicate_taxa member.
     * @param  matrix the input distance matrix (an n*n matrix, 
     *         of inter-taxon distances, stored in row-major order).
     */
    void identifyDuplicateTaxa(const double* matrix) {
        //1. Calculate row hashes
        #if USE_PROGRESS_DISPLAY
        const char* task_name = "";
        if (!be_silent) {
            task_name = "Identifying duplicate clusters";
        }
        progress_display show_progress((double)original_rank, task_name, "", "");
        #endif

        std::vector< HashRow<double> > hashed_rows(original_rank);
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (intptr_t i=0; i<original_rank; ++i) {
            hashed_rows[i] = HashRow<double>(i, matrix+original_rank*i, 
                                             original_rank);
            #if USE_PROGRESS_DISPLAY
                if ((i%1000)==999) {
                    show_progress += (double)(1000.0);
                }
            #endif
        }
        //2. Sort rows by hash (and tiebreak on row content)
        std::sort(hashed_rows.begin(), hashed_rows.end());

        //3. Now, identify any identical rows.
        duplicate_taxa.clear();
        HashRow<double>::identifyDuplicateClusters(hashed_rows, duplicate_taxa);

        #if USE_PROGRESS_DISPLAY
            show_progress += (double)(original_rank%1000);
            show_progress.done();
        #endif
    }

    /**
     * @brief  Identify, and cluster, taxa that are identical (or, at least,
     *         have identical rows in the supplied distance matrix).
     * @param  n          - the number of taxa
     * @param  recalcTime - used to return how long it took to identify
     *                      clusters of duplciate taxa (elapsed seconds).
     * @param  mergeTime  - used to return how long it took to join the
     *                      clusters of duplicate tax (elapsed seconds).
     * @return intptr_t   - how many cluster joins were done (the number
     *                      of taxa that were duplicates, minus the number
     *                      of taxon equivalence classes).
     * @note on entry, duplciate_taxa has already been set.
     *       usually, in identifyDuplicateTaxa().
     *       (yes, this member function is "episodic").
     */
    intptr_t clusterDuplicateTaxa(int& n, double& recalcTime, 
                                  double& mergeTime) {
        if (duplicate_taxa.empty()) {
            return 0;
        }

        #if USE_PROGRESS_DISPLAY
        double work_estimate = 0;
        for (std::vector<intptr_t>& cluster_members: duplicate_taxa) {
            work_estimate += cluster_members.size() - 1.0;
        }
        const char* task_name = "";
        if (!be_silent) {
            task_name = "Merging duplicate clusters";
        }
        progress_display show_progress(work_estimate, task_name, "", "");
        #endif

        intptr_t duplicate_merges = 0;
        for (std::vector<intptr_t>& cluster_members: duplicate_taxa) {
            //cluster_members: cluster numbers of the next set 
            //of duplicate taxa to merge...
            while (1<cluster_members.size()) {
                std::vector<intptr_t> next_level;
                for (size_t i=0; i<cluster_members.size(); i+=2) {
                    if (i+1==cluster_members.size()) {
                        next_level.push_back(cluster_members[i]);
                    } else {    
                        int x = cluster_members[i];
                        int y = cluster_members[i+1];

                        recalcTime  -= getRealTime();
                        recalculateTotals (n, next_cluster_number);
                        recalcTime  += getRealTime();

                        int low_row  = cluster_row[x];
                        int high_row = cluster_row[y];
                        int n_less_1 = n - 1;
                        if (high_row<low_row) {
                            std::swap(low_row, high_row);
                        }
                        row_cluster[low_row]  = next_cluster_number;
                        row_cluster[high_row] = row_cluster[n_less_1];   
                        mergeTime -= getRealTime();
                        mergeClusters(x, y, next_cluster_number,
                                        (T)0, n, false);
                        mergeTime += getRealTime();
                        next_level.push_back(next_cluster_number);
                        ++next_cluster_number;
                        #if USE_PROGRESS_DISPLAY
                        show_progress+=1.0;
                        #endif
                        --n;
                        ++duplicate_merges;
                    }
                }
                std::swap(cluster_members, next_level);
            }
        }
        #if USE_PROGRESS_DISPLAY
            show_progress.done();
        #endif
        return duplicate_merges;
    }

    /**
     * @brief Recalculate cluster totals, for clusters 0..next_cluster_num
     *        and map current rows to in-play clusters (an in-play
     *        cluster is one that has not yet been joined to make a
     *        larger cluster).
     * @param n - the number of clusters in play
     * @param next_cluster_num - the next unused cluster number 
     */
    void recalculateTotals(int n, int next_cluster_num) {
        //writes: row_cluster and cluster_row
        double cutoff           = -infiniteDistance;
        double one_on_n_minus_2 = (n<3) ? 0.0 : (1.0 / ((double)n - 2.0));
        int    r                = 0;
        FNJ_TRACE("\nRow Totals: ");
        for (int c=0; c<next_cluster_num; ++c) {
            if (0<cluster_in_play[c]) {
                cluster_total_scaled[c] = cluster_total[c] * one_on_n_minus_2;
                cluster_cutoff[c]       = cutoff;
                if (cutoff < cluster_total_scaled[c] ) {
                    cutoff = cluster_total_scaled[c];
                }
                cluster_row[c] = r;
                row_cluster[r] = c;
                ++r;
                FNJ_TRACE("c=" << c << " t=" << cluster_total_scaled[c]
                          << " s=" << cluster_total[c] << "\n");
            }
        }
    }
    /**
     * @brief For each cluster, y, that is currently in play, find the
     *        first MatrixEntry<T>, in that cluster's block of 
     *        matrix entries for LOWER numbered clusters (which is
     *        sorted by raw distance from y), that refers to a
     *        cluster that is till in play (if any!).
     * @param n - the number of clusters currently in play
     * @param q - the number of nodes there will be, in a
     *            complete unrooted tree:(n0-1)*2, where n0 is
     *            the number of rows in the input distance matrix
     * @note  reads:  row_cluster
     * @note  writes: row_raw_dist, row_best_distance, row_choice
     * @note  doesn't usually update cluster_sorted_start[y],
     *        but will deallocate the memory allocated for the 
     *        cluster, if its block of matrix entries no longer
     *        has any distances to (lower-numbered) in-play clusters
     *        in it. 
     */
    void previewRows(int n, int q) {
        global_best_dist = infiniteDistance;
        FNJ_TRACE("\nPreview for n=" << n << "\n");

        //Problem here: reduction(min) not supported in Visual Studio C++ 19 on Windows
        //because Visual Studio C++ only supports OpenMP 2.0, and reduction(min) is OpenMP 3.1
        #ifdef _OPENMP
        #ifndef _MSC_VER
        #pragma omp parallel for reduction(min:global_best_dist)
        #else
        #pragma omp parallel for
        #endif
        #endif
        for (int r = 0; r < n; ++r ) /* r is current row number */ {
            int  y     = row_cluster[r]; //y is cluster number
            auto scan  = cluster_sorted_start[y]; 
            auto stop  = cluster_sorted_stop[y];
            bool found = false;
            for (; scan<stop; ++scan) {
                int x = scan->cluster_num;
                if (0<cluster_in_play[x]) {
                    row_raw_dist[r]  = scan->distance;
                    row_best_dist[r] = scan->distance
                                     - cluster_total_scaled[x]
                                     - cluster_total_scaled[y];
                    row_choice[r]    = x;
                    found = true;
                    FNJ_TRACE("x=" << x << ", y=" << y
                              << ", Dxy=" << row_raw_dist[r]
                              << ", Rx=" << cluster_total_scaled[x]
                              << ", Ry=" << cluster_total_scaled[y]
                              << ", Dxy-Rx-Ry=" << row_best_dist[r]
                              << "\n");
                    #ifdef _MSC_VER
                    #ifdef _OPENMP
                    #pragma omp critical
                    #endif
                    #endif              
                    if (row_best_dist[r]<global_best_dist) {
                        global_best_dist = row_best_dist[r];
                    }                    
                    break;
                } else {
                    FNJ_TRACE("For y=" << y << ", cluster x=" 
                              << x << " is no longer in play.");
                }
            }
            if (!found) {
                row_raw_dist[r]  = infiniteDistance;
                row_best_dist[r] = infiniteDistance;
                row_choice[r]    = q;
                deallocateCluster(y);
            }
        }
    }

    /**
     * @brief For each row, r, and corresponding in-play cluster, y,
     *        find the best choice (of MatrixEntry) for cluster y 
     *        (and so, for row r), and also the raw, and adjusted 
     *        distance, of that best choice.
     * @param n - the number of clusters currently in play
     * @param q - the number of nodes there will be, in a
     *            complete unrooted tree:(n0-1)*2, where n0 is
     *            the number of rows in the input distance matrix
     * @note  reads:   rows_by_dist (ordered by increasing preview 
     *                               distance, to get r)
     *                 row_cluster  (indexed by r to get y)
     * @note  updates: row_raw_dist, row_best_dist, row_choice 
     *                 (all of these vectors are indexed by r)
     * @note  if _OPENMP is defined, parallelizes over rows.
     *        each thread, (i) of the (step) running threads looks
     *        at those rows that correspond to entries that are
     *        (i) modulo (step) in rows_by_dist.
     */
    void findPreferredPartners(int n, int q) {
        //reads:   row_cluster
        //updates: row_raw_dist, row_best_dist, row_choice
        FNJ_TRACE("\nSearch, n=" << n << "\n");

        std::vector< std::pair<T, int> > rows_by_dist;
        chooseRowSearchOrder(n, rows_by_dist);

        //For each cluster, find preferred partner
        #ifdef _OPENMP
        #pragma omp parallel 
        #endif
        {
            int step = getThreadCount();
            double per_thread_best_dist = global_best_dist;
            for (int i = getThreadNumber(); i < n; i+=step ) {
                int r = rows_by_dist[i].second;
                int y = row_cluster[r];
                if (cluster_in_play[y]==0) {
                    continue;                    
                }
                findPartnerForOneCluster(r,y);
                if (row_best_dist[r] < per_thread_best_dist) {
                    per_thread_best_dist = row_best_dist[r];
                    if (per_thread_best_dist < global_best_dist) {
                        #pragma omp critical
                        if (per_thread_best_dist < global_best_dist) {
                            global_best_dist = per_thread_best_dist;
                        }
                        else {
                            per_thread_best_dist = global_best_dist;
                        }
                    }
                }
            }
        }
    }

    /**
     * @brief Find an ordering of rows, ordering each row (r), according to 
     *        the (Dcj - Rc - Rj) "cooked" distance, to the first entry 
     *        (for cluster x) from the cluster (y) mapped to row r.
     *        (Here: Dcj is the distance between clusters c and j, Rc
     *         is the scaled row total, sum divided by (n-2), for cluster c,
     *         and Rj the same for cluster j).
     * @param n            - the number of clsuters currently inplay
     * @param rows_by_dist - a vector of pair<T,int>, where first is the
     *                       adjusted distance, and second the row number
     * @note  reads: row_best_dist.
     * @note  Sorts using std::sort().  Possibly for VERY large inputs, 
     *        it might be worth sorting with a parallel sorting 
     *        implementation.  But I haven't done that. -James B.
     */
    void chooseRowSearchOrder(int n, std::vector< std::pair<T, int> >& rows_by_dist) {
        for (int r = 0; r < n; ++r ) {
            rows_by_dist.emplace_back(row_best_dist[r], r);
        }
        std::sort(rows_by_dist.begin(), rows_by_dist.end());
    }

    /**
    * @brief For a given cluster, y, corresponding to a specified row, r,
    *        of the working distance matrix, find the MatrixEntry,
    *        and corresponding distance, for the (earlier-numbered)
    *        cluster, x, still in play, that has the lowest adjusted 
    *        distance to cluster y.
    * @param r - row number (of the cluster, in the working distance matrix)
    * @param y - cluster number
    * @note  reads: row_choice, row_best_dist.
    * @note  writes:row_best_dist and row_choice.
    * @note  Rather than searching (and comparing distances, until it
    *        is no longer possible, for a MatrixEntry<T> in y's block of
    *        lower-numbered MatrixEntry<T>, to be for a cluster, x, 
    *        such that the adjusted distance to x, from y, will be the 
    *        lowest), using the current cutoff, this code performs a 
    *        binary search, for an entry beyond the initial cutoff.
    *        The idea is to spend log_2(CX(y)) unpredictable dstiance 
    *        comparisons before the scanning loop, to save CX(y) in it.
    */
    virtual bool findPartnerForOneCluster(int r /*row*/, int y /*cluster*/) {
        int best_x         = row_choice[r];    //other cluster
        T   best_hc_dist   = row_best_dist[r]  + cluster_total_scaled[y];
        T   cutoff         = best_hc_dist      + cluster_cutoff[y];
        T   earlier_cutoff = global_best_dist 
                           + cluster_cutoff[y] + cluster_total_scaled[y];
        if (earlier_cutoff < cutoff) {
            cutoff = earlier_cutoff;
        }
        //Any MatrixEntry in the data for this cluster that has a
        //distance, Dxy greater tan cutoff will have Dxy - Rx - Ry
        //greater than max(global_best_dist, row_best_dist),
        //where global_best_dist is the best Dij - Ri - Rj found
        //for any i,j, so far, and row_best_dist is the best 
        //Dyi - Ri - Ry found for cluster y.

        //Can skip the first MatrixEntry at cluster_sorted_start[c].
        //As it has definitely been looked at, already, by 
        //previewRows().
        auto dataStart = cluster_sorted_start[y] + 1; 
        auto dataStop  = cluster_sorted_stop[y];
        bool found     = false;

        dataStop = findFirstGreaterDistance
                    (dataStart, dataStop, cutoff);
        for (auto scan=dataStart; scan<dataStop; ++scan) {
            int x                = scan->cluster_num;
            T   dist_half_cooked = scan->distance
                                 - cluster_total_scaled[x];
            if (best_hc_dist<=dist_half_cooked) {
                continue;
            }
            found         = true;
            best_hc_dist  = dist_half_cooked;
            best_x        = x; //best cluster found
        }
        row_best_dist[r] = best_hc_dist - cluster_total_scaled[y];
        row_choice[r]    = best_x;
        return found;
    }

    /**
     * @brief  Find, in a block of contiguous MatrixEntry<T> rows,
     *         indicated by (start, stop), sorted by distance, 
     *         the first entry that has a distance greater than (dist),
     *         and return a pointer to it. If there ISN'T one, 
     *         return stop.
     * @param  start - first MatrixEntry<T> to look at
     * @param  stop  - one more than the last MatrixEntry<T> to look at
     *                 (or, if you prefer, the first one, after than,
     *                 NOT to look at)
     * @param  dist  - the threshold distance, we are looking to find
     *                 a "boundary" element for.
     * @return MatrixEntry* the first such element. Or (stop) if there are
     *         no elements with a distance greater than (dist).
     */
    MatrixEntry* findFirstGreaterDistance(MatrixEntry* start, 
                                          MatrixEntry* stop, T dist) {
        while (start<stop) {
            MatrixEntry* middle  = start + (stop-start) / 2;
            bool         greater = dist < middle->distance;
            if (greater) {
                stop=middle; 
            } else {
                start = middle+1;
            }
        }
        return start;
    }

    /**
     * @brief  Given, row_choice and row_best_dist have been set,
     *         (indicating for each row, r, and corresponding cluster, y,
     *          which lower-numbered cluster x, if any, has the minimal
     *          adjusted distance to y), find the best row.
     * @param  n - the number of clusters currently in play
     * @param  q - a dummy cluster number (indicating no cluster!) 
     *           - in practice, (n0-1)*2 where n0 is the number of
     *             rows that there were in the initial distance matrix.
     * @return int - the row, r, which corresponds to the cluster, y, 
     *               which has the mimimum minimum distance to a 
     *               lower-numbered in-play cluster, x.
     * @note   reads:  row_choice, row_best_dist
     */
    int chooseBestRow(int n, int q) {
        //reads: row_best_dist, row_choice
        int best_row  = q;
        T   best_dist = infiniteDistance;
        //std::cout << "Choose Row:";
        int r;
        for (r=0; r<n; ++r) {
            if (row_choice[r] < q) {
                best_dist = row_best_dist[r];
                best_row  = r;
                //std::cout << " R" << best_row << "=D" << best_dist; 
                break;
            }
        }
        for (++r; r<n; ++r) {
            if (row_best_dist[r] < best_dist) {
                if (row_choice[r] < q) {
                    best_dist = row_best_dist[r];
                    best_row  = r;
                    //std::cout << " R" << best_row << "=D" << best_dist; 
                }    
            }
        }
        //std::cout << "\n";
        return best_row;
    }

    /**
     * @brief merge two existing clusters, into a new cluster
     * @param cluster_X - the first existing cluster
     * @param cluster_Y - the second
     * @param cluster_U - the new cluster to be built by joining
     *                    clusters cluster_X and cluster_Y
     * @param Dxy       - the raw (not the adjusted) distance between
     *                    the two existing clusters
     * @param n         - the number of clusters still in play
     * @param is_rooted - indicates whether the last cluster is to have
     *                    degree 2 (true), or degree 3 (false).
     * @note  it is assumed 0 <= cluster_X < cluster_Y < cluster_U.
     */
    void mergeClusters(int cluster_X, int cluster_Y, 
                       int cluster_U, T Dxy, int n,
                       bool is_rooted) {
        allocateCluster(cluster_U);
        ASSERT(cluster_Y != cluster_X);
        if (cluster_Y<cluster_X) {
            std::swap(cluster_X, cluster_Y);
        }
        T nless2        = (T)(n-2);
        T tMultiplier   = (n<3) ? (T)0.0 : ((T)0.5 / nless2);
        T lambda        = (T)0.5;
        T medianLength  = lambda * Dxy;
        T fudge         = (cluster_total[cluster_X] - 
                           cluster_total[cluster_Y]) * tMultiplier;
        T length_to_X   = medianLength + fudge;
        T length_to_Y   = medianLength - fudge;
        T mu            = (T)1.0 - lambda;
        T dCorrection   = - lambda * length_to_X - mu * length_to_Y;

        x_distances.resize(cluster_U, infiniteDistance);
        getDistances(cluster_X, cluster_U, x_distances);
        y_distances.resize(cluster_U, infiniteDistance);
        getDistances(cluster_Y, cluster_U, y_distances);
        auto start      = cluster_sorted_start[cluster_U];
        auto entry      = start;
        auto unsorted   = cluster_unsorted_start[cluster_U];
        T    cTotal     = 0.0;
        for (int c=0; c<cluster_U; ++c) {
            if ( c!=cluster_X && c!=cluster_Y && cluster_in_play[c] ) {
                T Dcx = x_distances[c];
                T Dcy = y_distances[c];
                T Dcu              = lambda * Dcx + mu * Dcy + dCorrection;
                entry->distance    = Dcu;
                entry->cluster_num = c;
                cluster_total[c]  += Dcu - Dcx - Dcy;
                cTotal            += Dcu;
                *unsorted          = *entry;
                ++entry;
                ++unsorted;
            }
        }
        cluster_total[cluster_U]   = cTotal;
        auto stop                  = cluster_sorted_stop[cluster_U];
        if (2<n || is_rooted) {
            ASSERT ( stop == entry );
            sorters[0].parallel_sort(start, stop-start);
            clusters.addCluster(cluster_X, length_to_X, 
                                cluster_Y, length_to_Y);
        } else {
            //we're not connecting up X and Y to U,
            //we're connecting X and Y with each other
            //  Not this:       But rather this:
            //  X         Y          
            //    \      /            
            //    dx   dy       X--(dx+dy)--Y
            //      \ /
            //       U
            //       |
            //       ?
            //
            clusters.appendToLastCluster(cluster_X, length_to_X + length_to_Y);
        }
        markClusterAsUsedUp(cluster_X);
        markClusterAsUsedUp(cluster_Y);
    }

    /**
    * @brief mark a cluster as no longer being in play (because it has
    *        been merged into another cluster).
    * @param c - the cluster number
    */
    void markClusterAsUsedUp(int c) {
        cluster_in_play[c]       = 0;
        cluster_total[c]         = -infiniteDistance;
        cluster_total_scaled[c]  = -infiniteDistance;
        cluster_cutoff[c]        = -infiniteDistance;
        deallocateCluster(c);
    }

    /**
    * @brief Free up the memory allocated on behalf of a cluser that is
    *        no longer in play.
    * @param c - the cluster number.
    */
    void deallocateCluster(int c) {
        delete [] cluster_sorted_start[c];
        cluster_sorted_start[c]   = nullptr;
        cluster_sorted_stop[c]    = nullptr;
        cluster_unsorted_start[c] = nullptr;
        cluster_unsorted_stop[c]  = nullptr;
    }

    /**
    * @brief Get the raw distances, from cluster c, to lower-numbered clusters
    * @param c - cluster number
    * @param u -
    * @param distances - a vector, of distances (i.e. a vector of T),
    *                    initialized on entry, to infiniteDistance, 
    *                    of size at least n (where n is the number of
    *                    clusters currently in play)).  
    *                    On exit, distance[b] will
    *                    be set to the raw distance between b and c, if
    *                    and only if b is currently in play.
    * @note  The sorted distances, from cluster cluster_X, and from 
    *        cluster cluster_Y, are likely to be in different orders.
    *        That is why unsorted distances are used, by getDistances();
    *        (despite the name, those are already in cluster-number order).
    * @note  For clusters b<c, distances can be read from cluster c's
    *        unsorted distances row. But for c<b clusters, distances have
    *        to be read from the index c entry in the unsorted distances
    *        row for cluster b.
    *        This is the whole reason that unsorted distances have to be 
    *        tracked in the first place!  So they can be looked up, here,
    *        like this, via clusterDistance() - see below.  -James B.
    */
    void getDistances(int c, int u, DistanceVector& distances) {
        //It's better to read from cluster_unsorted_start, 
        //since reading (and, more to the point, writing) 
        //in cluster number order is cache-friendlier.
        MatrixEntry* start       = cluster_unsorted_start[c];
        MatrixEntry* stop        = cluster_unsorted_stop[c];
        intptr_t     entry_count = stop - start;

        #ifdef _OPENMP
        #pragma omp parallel for if(1<threadCount)
        #endif
        for (intptr_t i = 0; i<entry_count; ++i) {
            MatrixEntry* scan = start + i;
            int       p  = scan->cluster_num;
            distances[p] = (0<cluster_in_play[p]) 
                         ? scan->distance : distances[p];
        }
        //Now, for clusters *after* c, for which 
        //distances to c were recorded.
        #ifdef _OPENMP
        #pragma omp parallel for if(1<threadCount)
        #endif
        for (int p=c+1; p<u; ++p) {
            if (0<cluster_in_play[p]) {
                distances[p] = clusterDistance(c, p);
            }
        }
    }

    /**
     * @brief Finds:  record of raw distance, to a, from b, using an 
     *                interpolation search in the "unsorted" entries 
     *                for cluster b (which, it so happens, are written 
     *                in cluster order).
     * @param  a lower-numbered cluster's cluster index
     * @param  b higher-numberd cluster's cluster index
     * @return T distance type
     * @note   Assumes that a is less than b (but doesn't check that is!)
     * @note   Why:   Theoretically, an interpolation search over x entries... 
     *                has a running time proportional to 
     *                the log of the log of x, where x = max(b-1, n0+n0-2-b),
     *                where n0 is the number of rows there were in the input
     *                distance matrix.
     * @note   Since: An overhead of n*n*log(log(n)) reads isn't serious 
     *                (since, we can hope that each of the ~ order n*n ~ 
     *                interpolation searches only results in at most 
     *                two or three cache misses).
     * @note   A binary search would take time proportional to log(x)
     *         which is more by a factor proportional to log(x).
     */
    T clusterDistance(int a, int b) {
        auto start = cluster_unsorted_start[b];
        int  count = cluster_unsorted_stop[b] - start;
        int  low   = 0;
        int  hi    = b;

        while (0<count) {
            double dGuess = (double)(a-low) / (double)(hi-low) 
                          * (double)(count);
            int    iGuess = static_cast<int>(floor(dGuess+.5));
            if (count<=iGuess) {
                iGuess = count-1;
            }
            MatrixEntry* guess = start + iGuess;
            if (guess->cluster_num<a) {
                low    = guess->cluster_num + 1;
                count -= (iGuess+1);
                start  = guess + 1;
            } else if (a<guess->cluster_num) {
                count  = iGuess;
                hi     = guess->cluster_num;
            } else {
                //std::cout << "Found " << a << " in " << b << "'s D+I row\n";
                return guess->distance;
            }
            //This will always terminate, because the entry 
            //pointed to by guess (at least) is always removed 
            //from consideration.
        }

        //std::cout << "Did not find "<< a << " in " << b << "'s D+I row\n";
        return infiniteDistance;
    }

    /**
     * @brief For each cluster that is currently in play,
     *        examine the block of distance-sorted MatrixEntry<T>,
     *        and "shrink" the block (by adjusting the boundary
     *        from the left, to point to the first entry for an
     *        "in-play" cluster, shuffling "live" entries back 
     *        toward the left of the block, and then adjusting
     *        the boundary on the right, to point after the last
     *        live entry.
     * @param used_cluster_count 
     * @note  reads:  cluster_in_play
     *        reads:  cluster_sorted_start, cluster_sorted_stop
     * @note  writes: cluster_sorted_start, cluster_sorted_stop
     */
    void removeOutOfPlayClusters(int used_cluster_count) {
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (int c=0; c<used_cluster_count; ++c) {
            if (cluster_in_play[c]==0) {
                continue;
            }
            MatrixEntry* w          = cluster_sorted_start[c];
            MatrixEntry* scan       = w;
            for (; scan<cluster_sorted_stop[c]; ++scan) {
                *w = *scan;
                w += cluster_in_play[scan->cluster_num];
            }
            cluster_sorted_stop[c]    = w;
            scan                      = cluster_unsorted_start[c];
            cluster_unsorted_start[c] = w;
            for (; scan<cluster_unsorted_stop[c]; ++scan) {
                *w = *scan;
                w += cluster_in_play[scan->cluster_num];
            }
            cluster_unsorted_stop[c] = w;
        }
    }

    /**
     * @brief Log a summary of what has been done, neighbour joining 
     *        a phylogenetic tree, for the taxa, given the input distance
     *        matrix.  Most of the information here is about how much
     *        time was spent on various activities.
     * @param duplicate_merges - how many merges there were of clusters
     *                           containing duplicate taxa
     * @param purgeTime        - total time spent removing "out-of-play"
     *                           MatrixEntry<T> elements from the sorted
     *                           by distance blocks, for "in-play" clusters.
     * @param recalcTime       - time spent recalculating distances
     * @param previewTime      - time spent getting "previews" (first live
     *                           cluster in each cluster's sorted-by-distance
     *                           MatrixEntry<T> block).
     * @param mergeTime        - time spent merging clusters
     * @note  All times are in seconds.
     */
    void reportConstructionDone(intptr_t duplicate_merges, 
                                double purgeTime,
                                double recalcTime, 
                                double previewTime, 
                                double mergeTime) {
        if (!be_silent) {
            if (0<duplicate_merges) {
                std::cout << "Did " << duplicate_merges 
                          << " merges of duplicate clusters.\n";
            }
            std::cout << "Purging time was " << purgeTime << ","
                      << " Recalc time was " << recalcTime << ".\n";
            std::cout << "Preview time was " << previewTime << ","
                      << " Cluster merge time was " << mergeTime << ".\n";
        }
    }
}; //FancyNJMatrix template class

/**
 * @brief Vectorized version of FancyNJMatrix
 * @note  The only method that is vectorized is findPartnerForOneCluster().
 *        (Since that's the one that "matters" most, by far). Between
 *        n0*n0 and n0^3 (~ probably about no~2.5) operations on T.
 * @note  I didn't judge it worth vectorizing getDistances(), as that
 *        only does O(n0*n0) operations on T. -James B.
 * @note  Nor did I judge it worth vectorizing chooseBestRow(). 
 *        For the same reason. -James B.
 */
#if USE_VECTORCLASS_LIBRARY
template <class T=NJFloat, class V=FloatVector, class VB=FloatBoolVector> 
class VectorizedFancyNJMatrix: public FancyNJMatrix<T> {
public:
    typedef FancyNJMatrix<T> super;
    typedef typename super::DistanceVector DistanceVector;

    using super::cluster_total_scaled;
    using super::cluster_cutoff;
    using super::cluster_sorted_start;
    using super::cluster_sorted_stop;

    using super::row_choice;
    using super::row_raw_dist;
    using super::row_best_dist;

    using super::global_best_dist;

    using super::getThreadCount;
    using super::getThreadNumber;
    using super::clusterDistance;
    using super::findFirstGreaterDistance;

protected:
    const intptr_t  block_size;
    int             thread_count;
    std::vector<T>  vector_storage;
    std::vector<T*> per_thread_vector_storage;

public:
    VectorizedFancyNJMatrix(): super(), block_size(VB().size()),
        thread_count(1), 
        vector_storage  (block_size*4) {
        #ifdef _OPENMP
        #pragma omp parallel
        {
            if (getThreadNumber()==0) {
                thread_count = getThreadCount();
            }
        }
        vector_storage.resize(thread_count*block_size*4);
        #endif
        for (int i=0; i<thread_count; ++i) {
            T* ptr = vector_storage.data() + block_size * i * 4;
            per_thread_vector_storage.push_back( ptr ) ;
        }
    }

    virtual bool findPartnerForOneCluster(int r /*row*/, int y /*cluster*/) override {
        auto thread_num     = getThreadNumber();
        int  best_x         = row_choice[r];    //other cluster
        T    best_hc_dist   = row_best_dist[r]  + cluster_total_scaled[y];
        T    cutoff         = best_hc_dist      + cluster_cutoff[y];
        T    earlier_cutoff = global_best_dist 
                            + cluster_cutoff[y] + cluster_total_scaled[y];

        if (earlier_cutoff < cutoff) {
            cutoff = earlier_cutoff;
        }
        //Any MatrixEntry in the data for this cluster that has a
        //distance, Dxy greater tan cutoff will have Dxy - Rx - Ry
        //greater than max(global_best_dist, row_best_dist),
        //where global_best_dist is the best Dij - Ri - Rj found
        //for any i,j, so far, and row_best_dist is the best 
        //Dyi - Ri - Ry found for cluster y.

        //Can skip the first MatrixEntry at cluster_sorted_start[c].
        //As it has already been looked at by previewRows().
        auto dataStart = cluster_sorted_start[y] + 1; 
        auto dataStop  = cluster_sorted_stop[y];

        dataStop = findFirstGreaterDistance
                   (dataStart, dataStop, cutoff);
        intptr_t blockCount = ( dataStop - dataStart ) / block_size;
        auto blockStop      = dataStart + ( blockCount * block_size );
        T*   blockRawDist   = per_thread_vector_storage[thread_num];
        T*   blockCluster   = blockRawDist + block_size;
        T*   blockIndex     = blockCluster + block_size;
        T*   blockHCDist    = blockIndex   + block_size;

        V    best_hc_vector = best_hc_dist;
        V    best_ix_vector = (T)best_x;
        V    raw(0);
        V    tot(0);
        V    ix(0);
        bool found = false;
        for (auto scan=dataStart; scan<blockStop; scan+=block_size) {
            for (int i=0; i<block_size; ++i) {
                blockRawDist[i] = scan[i].distance;
                blockCluster[i] = cluster_total_scaled[scan[i].cluster_num];
                blockIndex[i]   = scan[i].cluster_num;
            }

            raw.load(blockRawDist);
            tot.load(blockCluster);
            ix.load(blockIndex);
            V  hc(raw - tot); //subtract cluster totals to get half-cooked distances?
            VB less = hc < best_hc_vector; //which are improvements?
            best_hc_vector = select(less, hc, best_hc_vector);
            best_ix_vector = select(less, ix, best_ix_vector);
        }

        if (dataStart<blockStop) {
            best_hc_vector.store(blockHCDist);
            best_ix_vector.store(blockIndex);
            for (int i=0; i<block_size; ++i) {
                if (blockHCDist[i] < best_hc_dist ) {
                    best_hc_dist = blockHCDist[i];
                    best_x       = static_cast<int>(blockIndex[i]);
                    found        = true;
                }
            }
        }

        for (auto scan=blockStop; scan<dataStop; ++scan) {
            int x                = scan->cluster_num;
            T   dist_raw         = scan->distance;
            T   dist_half_cooked = dist_raw 
                                 - cluster_total_scaled[x];
            if (best_hc_dist<=dist_half_cooked) {                
                continue;
            }
            found         = true;
            best_hc_dist  = dist_half_cooked;
            best_x        = x; //best cluster found
        }

        row_best_dist[r] = best_hc_dist - cluster_total_scaled[y];
        row_choice[r]    = best_x;
        return found;
    } //findPartnerForOneCluster
}; //VectorizedFancyNJMatrix
#endif //USE_VECTORCLASS_LIBRARY


}; //StartTree Namespace

#endif //fancy_rapid_nj_h