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
//               Springer Verlag, 2011.
//
//  FancyNJMatrix<T> differs from BoundingMatrix<T, NJMatrix<T>> in several
//  respects, because it implemented with space efficiency performance in mind 
//  (BoundingMatrix was NOT; it was implemented to be easy to read, easy
//  to "map" back to the paper it is based on).
//
//  Copyright James Barbetti (2021)
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
#include "distancematrixtree/upgma.h"
#include "utils/gzstream.h"
#include "utils/timeutil.h"
#ifndef  fancy_rapid_nj_h
#define  fancy_rapid_nj_h

#include "nj.h"                       //for NJFloat (should be separate header for that)
#include "clustertree.h"              //for ClusterTree template class
#include "flatmatrix.h"               //for FlatMatrix
#include "hashrow.h"                  //for HashRow
#include "utils/parallel_mergesort.h" //for MergeSorter

#define  FNJ_TRACE(x) (0)

namespace StartTree {
template <class T=NJFloat> class FancyNJMatrix {
protected:
    bool           be_silent;
    bool           zip_it;
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
        typedef ifstream InFile;
        #define INFILE(name) Infile name
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

    size_t         original_rank;
    size_t         next_cluster_number;
    EntryVector    storage;
    IntVector      cluster_in_play;
    DistanceVector cluster_total;
    DistanceVector cluster_total_scaled;
    DistanceVector cluster_cutoff;
    EntryPtrVector cluster_sorted_start;
    EntryPtrVector cluster_sorted_stop;
    EntryPtrVector cluster_unsorted_start;
    EntryPtrVector cluster_unsorted_stop;
    IntVector      cluster_row;

    //uint64_t       search_iterations;

    int            threadCount;
    typedef MergeSorter<MatrixEntry> Sorter;
    std::vector<Sorter> sorters;
    int getThreadNumber() {
        #ifdef _OPENMP
            return omp_get_thread_num();
        #else
            return 0;
        #endif
    }
    int getThreadCount() {
        #ifdef _OPENMP
            return omp_get_num_threads();
        #else
            return 1;
        #endif
    }
    DuplicateTaxa duplicate_taxa;

public:
    FancyNJMatrix() : be_silent(false), zip_it(false), original_rank(0), 
                      next_cluster_number(0) /*, search_iterations(0)*/ {
        #ifdef _OPENMP
            threadCount = omp_get_max_threads();
        #else
            threadCount = 1;
        #endif
        sorters.resize(threadCount);
    }
    std::string getAlgorithmName() const {
        return "FancyNJ";
    }

    void beSilent() { 
        be_silent = true; 
    } 
    virtual void setZippedOutput(bool zipIt) { 
        zip_it = zipIt; 
    }
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
    virtual bool loadMatrixFromOpenFile(InFile& in) {
        FlatMatrix dummy;
        loadDistanceMatrixFromOpenFile(in, !be_silent, dummy);
        return this->loadMatrix(dummy.getSequenceNames(), dummy.getDistanceMatrix());
    }
    virtual bool loadMatrix(const StrVector& names,
                            const double* matrix) {
        #if USE_PROGRESS_DISPLAY
        double row_count  = names.size();
        double work_to_do = 0.5 * row_count * (row_count - 1.0);
        const char* taskName = be_silent
            ? "" :  "Setting up triangular I+S and I+D matrices";
        progress_display setupProgress(work_to_do, taskName, "loading", "row");
        #endif
        //Assumes: matrix is symmetrical
        clusters.clear();
        for (const std::string& name : names) {
            clusters.addCluster(name);
        }
        setRank(names.size()); //sets original_rank
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) if(1<threadCount)
        #endif
        for (int r=original_rank-1; 0<=r; --r ) {
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
    virtual bool constructTree() {
        //search_iterations = 0;
        if (original_rank<3) {
            return false;
        }

        int  n = original_rank;
        int  q = n + n - 2;
        IntVector      row_cluster(original_rank);
        DistanceVector row_raw_dist(original_rank);
        DistanceVector row_best_dist(original_rank);
        IntVector      row_choice(original_rank);
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (int i=0; i<original_rank; ++i) {
            row_cluster[i] = i;
        }

        double recalcTime  = 0;
        double previewTime = 0;
        double mergeTime   = 0;

        intptr_t duplicate_merges = clusterDuplicateTaxa(row_cluster, n, 
                                                         recalcTime, mergeTime);

        #if USE_PROGRESS_DISPLAY
        std::string taskName;
        if (!be_silent) {
            taskName = "Constructing " + getAlgorithmName() + " tree";
        }
        double triangle  = (double)n * ((double)n + 1.0) * 0.5;
        progress_display show_progress(triangle, taskName.c_str(), "", "");
        #endif

        for ( ; 1 < n ; --n) {
            T best_dist;
            recalcTime  -= getRealTime();
            recalculateTotals    (n, next_cluster_number, 
                                  row_cluster, cluster_row);
            recalcTime  += getRealTime();
            previewTime -= getRealTime();
            previewRows          (n, q, row_cluster,
                                  row_raw_dist, row_best_dist,
                                  row_choice, best_dist);
            previewTime += getRealTime();
            findPreferredPartners(n, q, best_dist, row_cluster, 
                                  row_raw_dist, row_best_dist, 
                                  row_choice);
            int best_row      = chooseBestRow(n, q, row_best_dist, row_choice);
            int best_cluster  = row_cluster[best_row];
            int other_cluster = row_choice[best_row];
            int other_row     = cluster_row[other_cluster];
            int low_row       = (best_row < other_row) ? best_row  : other_row;
            int high_row      = (best_row < other_row) ? other_row : best_row;
            row_cluster[low_row]  = next_cluster_number;
            row_cluster[high_row] = row_cluster[n-1];   
            mergeTime -= getRealTime();
            mergeClusters(best_cluster, other_cluster, next_cluster_number,
                          row_raw_dist[best_row], n);
            mergeTime += getRealTime();
            ++next_cluster_number;
            #if USE_PROGRESS_DISPLAY
            show_progress+=(double)n;
            #endif
        }

        #if USE_PROGRESS_DISPLAY
        show_progress.done();
        #endif
        if (!be_silent) {
            if (0<duplicate_merges) {
                std::cout << "Did " << duplicate_merges 
                          << " merges of duplicate clusters.\n";
            }
            std::cout << "Recalc time " << recalcTime << ","
                      << " Preview time " << previewTime << " and"
                      << " Cluster merge time was " << mergeTime << ".\n";
            #if (0)
            double coefficient = (double)search_iterations
                               / (double)original_rank 
                               / (double)original_rank;
            std::cout << "Number of search iterations was " 
                      << search_iterations
                      << " (coefficient " << coefficient << ").\n";
            #endif
        }
        return true;
    }
    virtual bool calculateRMSOfTMinusD(const double* matrix, 
                                       intptr_t rank, double& rms) {
        return clusters.calculateRMSOfTMinusD(matrix, rank, rms);
    }
    virtual bool writeTreeFile     (int precision,
                                    const std::string &treeFilePath) const { 
        return clusters.writeTreeFile(zip_it, precision, treeFilePath);
    }
protected:
    virtual void setRank(size_t n) {
        original_rank       = n;
        next_cluster_number = n;
        storage.resize(original_rank*original_rank*size_t(2));
        size_t q = n+n-2; //number of clusters that will be needed
                          //in total, during the course of the tree
                          //construction: n leaves, n-2 interiors.
        cluster_in_play.resize       (q, 1);
        cluster_total.resize         (q, 0.0);
        cluster_total_scaled.resize  (q, 0.0);
        cluster_sorted_start.resize  (q, nullptr);
        cluster_sorted_stop.resize   (q, nullptr);
        cluster_unsorted_start.resize(q, nullptr);
        cluster_unsorted_stop.resize (q, nullptr);
        cluster_row.resize           (q, -1);
        cluster_cutoff.resize        (q, 0.0);
        MatrixEntry* data   = storage.data();
        for (int c=0; c<q; ++c) {
            //Cluster, c, 0 through n-1, has c previous clusters
            //When     c, n through q-1, is created, it'll have q-c
            size_t entries = (c<n) ? c : q-c;
            cluster_sorted_start[c]    = data;
            data                      += entries;
            cluster_sorted_stop[c]     = data;
            cluster_unsorted_start[c]  = data;
            data                      += entries;
            cluster_unsorted_stop[c]   = data;
        }
    }

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

    intptr_t clusterDuplicateTaxa(IntVector& row_cluster, int& n,
                                  double& recalcTime, double& mergeTime) {
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
                for (int i=0; i<cluster_members.size(); i+=2) {
                    if (i+1==cluster_members.size()) {
                        next_level.push_back(cluster_members[i]);
                    } else {    
                        int x = cluster_members[i];
                        int y = cluster_members[i+1];

                        recalcTime  -= getRealTime();
                        recalculateTotals    (n, next_cluster_number, 
                                            row_cluster, cluster_row);
                        recalcTime  += getRealTime();

                        int low_row  = cluster_row[x];
                        int high_row = cluster_row[y];
                        if (high_row<low_row) {
                            std::swap(low_row, high_row);
                        }
                        row_cluster[low_row]  = next_cluster_number;
                        row_cluster[high_row] = row_cluster[n-1];   
                        mergeTime -= getRealTime();
                        mergeClusters(x, y, next_cluster_number,
                                        (T)0, n);
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

    void recalculateTotals(int n, int next_cluster_num, 
                           IntVector& row_cluster,
                           IntVector& cluster_row) {
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
    void previewRows(int n, int q,
                     const IntVector& row_cluster,
                     DistanceVector&  row_raw_dist,
                     DistanceVector&  row_best_dist,
                     IntVector&       row_choice,
                     T&               best_dist) {

        FNJ_TRACE("\nPreview for n=" << n << "\n");
        intptr_t iterations = 0;
        #ifdef _OPENMP
        #pragma omp parallel for if(1<threadCount) reduction(+:iterations)
        #endif
        for (int r = 0; r < n; ++r ) {
            int  y     = row_cluster[r];
            auto scan  = cluster_sorted_start[y]; 
            auto stop  = cluster_sorted_stop[y];
            bool found = false;
            for (; scan<stop; ++scan) {
                ++iterations;
                int x = scan->cluster_num;
                if (0<cluster_in_play[x]) {
                    row_raw_dist[r]  = scan->distance;
                    row_best_dist[r] = scan->distance
                                     - cluster_total_scaled[x]
                                     - cluster_total_scaled[y];
                    row_choice[r]    = x;
                    cluster_sorted_start[y] = scan;
                    found = true;
                    FNJ_TRACE("x=" << x << ", y=" << y
                              << ", Dxy=" << row_raw_dist[r]
                              << ", Rx=" << cluster_total_scaled[x]
                              << ", Ry=" << cluster_total_scaled[y]
                              << ", Dxy-Rx-Ry=" << row_best_dist[r]
                              << "\n");
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
                cluster_sorted_start[y] = cluster_sorted_stop[y];
            }
        }

        best_dist = row_best_dist[0];
        for (int r=1; r<n; ++r) {
            if (row_best_dist[r]<best_dist) {
                best_dist = row_best_dist[r];
            }
        }
        //search_iterations += iterations;
    }
    void findPreferredPartners(int n, int q, T global_best_dist,
                               const IntVector& row_cluster,
                               DistanceVector&  row_raw_dist,
                               DistanceVector&  row_best_dist,
                               IntVector&       row_choice) {

        FNJ_TRACE("\nSearch, n=" << n << "\n");
        //
        //Find an ordering of rows, ordering each row (r),
        //according to the (Dcj - Rc - Rj) "cooked" distance, 
        //to the first entry (for cluster x) from the cluster (y)
        //mapped to row r.
        //
        std::vector< std::pair<T, int> > rows_by_dist;
        for (int r = 0; r < n; ++r ) {
            rows_by_dist.emplace_back(row_best_dist[r], r);
        }
        std::sort(rows_by_dist.begin(), rows_by_dist.end());

        //intptr_t iterations = 0;
        //For each cluster, find preferred partner
        #ifdef _OPENMP
        #pragma omp parallel 
                //for schedule(static) if(1<threadCount) 
                //reduction(+:iterations)
        #endif
        {
            int step = getThreadCount();
            for (int i = getThreadNumber(); i < n; i+=step ) {
                int r = rows_by_dist[i].second;
                int y = row_cluster[r];
                if (cluster_in_play[y]==0) {
                    continue;
                }
                int best_x         = row_choice[r];    //other cluster
                T   best_raw_dist  = row_raw_dist[r];  //set by previewRows().
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
                //As it has already been looked at by previewRows().
                auto dataStart = cluster_sorted_start[y] + 1; 
                auto dataStop  = cluster_sorted_stop[y];

    #if (0)
                if (r+r < dataStop-dataStart) {
                    //At most half of the information recorded for this cluster
                    //is still of any use.  Purge it of references to any cluster x
                    //that is no longer in play (for which cluster_in_play[x]==0).
                    auto oldStop = dataStop;
                    for (auto scan=dataStart; scan<oldStop; ++scan) {
                        int x = scan->cluster_num;
                        *dataStop = *scan;
                        dataStop += (cluster_in_play[x]==0) ? 0 : 1;
                    }
                    cluster_sorted_stop[y] = dataStop;
                }
    #endif
                dataStop = findFirstGreaterDistance
                           (dataStart, dataStop, cutoff);
                for (auto scan=dataStart; scan<dataStop; ++scan) {
                    //++iterations;
                    int x                = scan->cluster_num;
                    T   dist_raw         = scan->distance;
                    T   dist_half_cooked = dist_raw 
                                         - cluster_total_scaled[x];
                    FNJ_TRACE("x=" << x << ",y=" << y
                            << ", Dxy=" << dist_raw 
                            << " versus cutoff " << cutoff
                            << ", Dxy-Rx=" << dist_half_cooked
                            << " versus best " << best_hc_dist << "\n");
                    if (best_hc_dist<=dist_half_cooked) {
                        continue;
                    }
                    best_raw_dist = dist_raw;
                    best_hc_dist  = dist_half_cooked;
                    best_x        = x; //best cluster found
                    cutoff        = dist_half_cooked
                                  + cluster_cutoff[y];
                                
                    if (earlier_cutoff < cutoff) {                    
                        cutoff = earlier_cutoff;
                    }
                    FNJ_TRACE("New cutoff " << cutoff << "\n");
                }
                row_raw_dist[r]  = best_raw_dist;
                row_best_dist[r] = best_hc_dist - cluster_total_scaled[y];
                row_choice[r]    = best_x;
                //#pragma omp critical
                if (row_best_dist[r] < global_best_dist) {
                    global_best_dist = row_best_dist[r];
                }
            }
        }
        #if (0)
        for (int r = 0; r < n; ++r ) {
            int c             = row_cluster[r];
            FNJ_TRACE("c=" << c << " raw=" << row_raw_dist[r]
                      << " d=" << row_best_dist[r] 
                      << " o=" << row_choice[r] << "\n");
        }
        #endif
        //search_iterations += iterations;
    }
    MatrixEntry* findFirstGreaterDistance(MatrixEntry* start, 
                                          MatrixEntry* stop, T dist) {
        while (start<stop) {
            MatrixEntry* middle  = start + (stop-start) / 2;
            bool         greater = dist < middle->distance;
#if (0)
            //If two conditional moves 
            //cheaper than one if statement
            stop                 = greater ? middle : stop;
            ++middle;
            start                = greater ? start : middle;
#else
            if (greater) stop=middle; else start = middle+1;
#endif
        }
        return start;
    }
    int chooseBestRow(int n, int q,
                      const DistanceVector& row_best_dist,
                      const IntVector&      row_choice) {
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
    void mergeClusters(int cluster_X, int cluster_Y, 
                       int cluster_U, T Dxy, int n) {
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

        std::vector<T> x_distances(cluster_U, infiniteDistance);
        getDistances(cluster_X, cluster_U, x_distances);
        std::vector<T> y_distances(cluster_U, infiniteDistance);
        getDistances(cluster_Y, cluster_U, y_distances);
        auto start      = cluster_sorted_start[cluster_U];
        auto entry      = start;
        auto unsorted   = cluster_unsorted_start[cluster_U];
        T    cTotal     = 0.0;
        #if (0)
        std::cout << "Merging clusters x=" << cluster_X 
                  << " and y=" << cluster_Y << " with dXY=" << Dxy << "\n";
        std::cout << "D+I for cluster " << cluster_U << "is:\n";
        #endif
        for (int c=0; c<cluster_U; ++c) {
            T Dcx = x_distances[c];
            T Dcy = y_distances[c];
            if (Dcx<infiniteDistance && Dcy<infiniteDistance &&
                c!=cluster_X && c!=cluster_Y) {
                T Dcu              = lambda * Dcx + mu * Dcy + dCorrection;
                //std::cout << c << " " << Dcu << " \t";
                entry->distance    = Dcu;
                entry->cluster_num = c;
                cluster_total[c]  += Dcu - Dcx - Dcy;
                cTotal            += Dcu;
                *unsorted          = *entry;
                ++entry;
                ++unsorted;
            }
        }
        //std::cout << "\n";
        cluster_total[cluster_U]   = cTotal;
        auto stop                  = cluster_sorted_stop[cluster_U];
        if (2<n) {
            ASSERT ( stop == entry );
            //formerly: std::sort(start, stop);
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
    void markClusterAsUsedUp(int c) {
        cluster_in_play[c]       = 0;
        cluster_total[c]         = -infiniteDistance;
        cluster_sorted_stop[c]   = cluster_sorted_start[c];
        cluster_total_scaled[c]  = -infiniteDistance;
        cluster_cutoff[c]        = -infiniteDistance;
        cluster_unsorted_stop[c] = cluster_unsorted_start[c];
    }
    void getDistances(int c, int u, DistanceVector& distances) {
        //Have to read from cluster_unsorted_start, because 
        //cluster_sorted_start might have "skipped" over entries
        //that have been removed from consideration; see 
        //the implementation of previewRows().  Besides, reading
        //in order is cache-friendlier.
        auto start = cluster_unsorted_start[c];
        auto stop  = cluster_unsorted_stop[c];
        #ifdef _OPENMP
        #pragma omp parallel for if(1<threadCount)
        #endif
        for (auto scan = start; scan<stop; ++scan) {
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
    T clusterDistance(int a, int b) {
        //Assumes:  a is less than b
        //Finds:    record of raw distance, to a, from b,
        //          using an interpolation search in the
        //          "unsorted" entries for cluster b
        //          (which, it so happens, are written 
        //           in cluster order).
        //Why?:     Theoretically, an interpolation search 
        //          over x entries... has a running time 
        //          proportional to the log of the log of x,
        //          where x = max(b-1, n+n-2-b).
        //Since:    An overhead of n*n*log(log(n)) reads 
        //          isn't serious (since, we can hope that
        //          each of the ~ order n*n ~ interpolation 
        //          searches only results in one or two cache 
        //          misses).
        //         
        auto start = cluster_unsorted_start[b];
        int  count = cluster_unsorted_stop[b] - start;
        int  low   = 0;
        int  hi    = b;

#if (0)
        std::stringstream s;
        s << "Find " << a << " in " << b << "'s D+I row, which is: \n";
        for (int i=0; i<count; ++i) {
            s << "/" << start[i].cluster_num << "/" << start[i].distance << " ";
        }
        s << "\n";
        std::cout << s.str();
#endif

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
}; //FancyNJMatrix template class
} //StartTree Namespace

#endif