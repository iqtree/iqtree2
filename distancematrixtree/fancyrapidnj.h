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

#define  FNJ_TRACE(x) (0)

namespace StartTree {
template <class T=NJFloat> class FancyNJMatrix {
protected:
    bool           be_silent;
    bool           zip_it;
    bool           append_file;          
    bool           is_rooted;
    bool           omit_semicolon;
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
    IntVector      cluster_in_play;
    DistanceVector cluster_total;
    DistanceVector cluster_total_scaled;
    DistanceVector cluster_cutoff;
    EntryPtrVector cluster_sorted_start;
    EntryPtrVector cluster_sorted_stop;
    EntryPtrVector cluster_unsorted_start;
    EntryPtrVector cluster_unsorted_stop;
    IntVector      cluster_row;

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

    int            threadCount;
    typedef MergeSorter<MatrixEntry> Sorter;
    std::vector<Sorter> sorters;

    volatile T     global_best_dist;

    int getThreadNumber() const {
        #ifdef _OPENMP
            return omp_get_thread_num();
        #else
            return 0;
        #endif
    }
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
    void beSilent() { 
        be_silent = true; 
    } 
    virtual bool setAppendFile(bool appendIt) {
        append_file = appendIt;
        return true;
    }
    virtual bool setZippedOutput(bool zipIt) { 
        zip_it = zipIt;
        return true;
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
        progress_display setupProgress(work_to_do, taskName, "", "");
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
    virtual void prepareToConstructTree() {
    }
    virtual bool setIsRooted(bool rootIt) {
        is_rooted = rootIt;
        return true;
    }
    virtual bool setSubtreeOnly(bool wantSubtree) {
        omit_semicolon = wantSubtree;
        return true;
    }
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

        int  next_purge = n * 7 / 8;

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
            row_cluster[high_row] = row_cluster[n-1];   
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
    virtual bool calculateRMSOfTMinusD(const double* matrix, 
                                       intptr_t rank, double& rms) {
        return clusters.calculateRMSOfTMinusD(matrix, rank, rms);
    }
    virtual bool writeTreeFile     (int precision,
                                    const std::string &treeFilePath) const { 
        return clusters.writeTreeFile
               ( zip_it, precision, treeFilePath
               , append_file, omit_semicolon );
    }
    virtual bool writeTreeToOpenFile(std::iostream &stream) const { 
        return clusters.writeTreeToOpenFile
               ( omit_semicolon, stream );
    }

protected:
    virtual void setRank(size_t n) {
        original_rank       = n;
        next_cluster_number = n;
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
        for (int c=0; c<n; ++c) {
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

    void allocateCluster(int c) {
        ASSERT(cluster_sorted_start[c] == nullptr);
        size_t       q       = original_rank + original_rank-2;
        //Cluster, c, 0 through n-1, has c previous clusters
        //When     c, n through q-1, is created, it'll have q-c
        size_t       entries = (c<original_rank) ? c : q-c;
        MatrixEntry* data    = new MatrixEntry[entries*2];
        cluster_sorted_start[c]    = data;
        data                      += entries;
        cluster_sorted_stop[c]     = data;
        cluster_unsorted_start[c]  = data;
        data                      += entries;
        cluster_unsorted_stop[c]   = data;
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

    intptr_t clusterDuplicateTaxa(int& n, double& recalcTime, 
                                  double& mergeTime) {
        //writes: row_cluster                                      
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
                        recalculateTotals (n, next_cluster_number);
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
    void previewRows(int n, int q) {
        //reads: row_cluster
        //writes: row_rw_dist, row_best_distance, row_choice

        global_best_dist = infiniteDistance;
        FNJ_TRACE("\nPreview for n=" << n << "\n");
        #ifdef _OPENMP
        #pragma omp parallel for reduction(min:global_best_dist)
        #endif
        for (int r = 0; r < n; ++r ) {
            int  y     = row_cluster[r];
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
            for (int i = getThreadNumber(); i < n; i+=step ) {
                int r = rows_by_dist[i].second;
                int y = row_cluster[r];
                if (cluster_in_play[y]==0) {
                    continue;                    
                }
                findPartnerForOneCluster(r,y);
                #pragma omp critical
                if (row_best_dist[r] < global_best_dist) {
                    global_best_dist = row_best_dist[r];
                }
            }
        }
    }

    void chooseRowSearchOrder(int n, std::vector< std::pair<T, int> >& rows_by_dist) {
        //
        //Find an ordering of rows, ordering each row (r),
        //according to the (Dcj - Rc - Rj) "cooked" distance, 
        //to the first entry (for cluster x) from the cluster (y)
        //mapped to row r.
        //
        for (int r = 0; r < n; ++r ) {
            rows_by_dist.emplace_back(row_best_dist[r], r);
        }
        std::sort(rows_by_dist.begin(), rows_by_dist.end());
    }

    virtual void findPartnerForOneCluster(int r /*row*/, int y /*cluster*/) {
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

        dataStop = findFirstGreaterDistance
                    (dataStart, dataStop, cutoff);
        for (auto scan=dataStart; scan<dataStop; ++scan) {
            int x                = scan->cluster_num;
            T   dist_half_cooked = scan->distance
                                 - cluster_total_scaled[x];
            if (best_hc_dist<=dist_half_cooked) {
                continue;
            }
            best_hc_dist  = dist_half_cooked;
            best_x        = x; //best cluster found
        }
        row_best_dist[r] = best_hc_dist - cluster_total_scaled[y];
        row_choice[r]    = best_x;
    }

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

    void markClusterAsUsedUp(int c) {
        cluster_in_play[c]       = 0;
        cluster_total[c]         = -infiniteDistance;
        cluster_total_scaled[c]  = -infiniteDistance;
        cluster_cutoff[c]        = -infiniteDistance;
        deallocateCluster(c);
    }

    void deallocateCluster(int c) {
        delete [] cluster_sorted_start[c];
        cluster_sorted_start[c]   = nullptr;
        cluster_sorted_stop[c]    = nullptr;
        cluster_unsorted_start[c] = nullptr;
        cluster_unsorted_stop[c]  = nullptr;
    }

    void getDistances(int c, int u, DistanceVector& distances) {
        //It's better to read from cluster_unsorted_start, 
        //since reading (and, more to the point, writing) 
        //in cluster number order is cache-friendlier.
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

#ifdef USE_VECTORCLASS_LIBRARY
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

    virtual void findPartnerForOneCluster(int r /*row*/, int y /*cluster*/) override {
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

        for (auto scan=dataStart; scan<blockStop; scan+=block_size) {
            for (int i=0; i<block_size; ++i) {
                blockRawDist[i] = scan[i].distance;
                blockCluster[i] = cluster_total_scaled[scan[i].cluster_num];
                blockIndex[i]   = scan[i].cluster_num;
            }

            V  raw;  raw.load(blockRawDist);
            V  tot;  tot.load(blockCluster);
            V  ix;   ix.load(blockIndex);
            V  hc   = raw - tot; //subtract cluster totals to get half-cooked distances?
            VB less = hc < best_hc_vector; //which are improvements?
            best_hc_vector = select(less, hc, best_hc_vector);
            best_ix_vector = select(less, ix, best_ix_vector);
        }

        if (dataStart<blockStop) {
            best_hc_vector.store(blockHCDist);
            best_ix_vector.store(blockIndex);
            bool found = false;
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
            best_hc_dist  = dist_half_cooked;
            best_x        = x; //best cluster found
        }

        row_best_dist[r] = best_hc_dist - cluster_total_scaled[y];
        row_choice[r]    = best_x;
    } //findPartnerForOneCluster
}; //VectorizedFancyNJMatrix
#endif //USE_VECTORCLASS_LIBRARY


}; //StartTree Namespace

#endif //fancy_rapid_nj_h