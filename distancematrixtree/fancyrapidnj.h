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

#include "nj.h"                       //for NJFloat
#include "clustertree.h"              //for ClusterTree template class
#include "utils/parallel_mergesort.h" //for MergeSorter

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
    IntVector      cluster_alive;
    DistanceVector cluster_total;
    DistanceVector cluster_total_scaled;
    DistanceVector cluster_cutoff;
    EntryPtrVector cluster_sorted_start;
    EntryPtrVector cluster_sorted_stop;
    EntryPtrVector cluster_unsorted_start;
    EntryPtrVector cluster_unsorted_stop;
    IntVector      cluster_row;

    int    threadCount;
    typedef MergeSorter<MatrixEntry> Sorter;
    std::vector<Sorter> sorters;
    int getThreadNumber() {
        #ifdef _OPENMP
            return omp_get_thread_num();
        #else
            return 0;
        #endif
    }

public:
    FancyNJMatrix() : be_silent(false), zip_it(false), 
                      original_rank(0), next_cluster_number(0) {
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
            //std::cout << "Z" << r << "...";
            for (int c=0; c<r; ++c) {
                data[c].distance    = row[c];
                data[c].cluster_num = c;
                unsorted[c]         = data[c];
                total              += row[c];
                //std::cout << " " << row[c];
            }
            for (int c=r+1; c<original_rank; ++c) {
                total              += row[c];
            }
            //std::cout << " ... total " << total << "\n";
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
        return true;
    }
    virtual bool constructTree() {
        if (original_rank<3) {
            return false;
        }

        #if USE_PROGRESS_DISPLAY
        std::string taskName;
        if (!be_silent) {
            taskName = "Constructing " + getAlgorithmName() + " tree";
        }
        double row_count = (double)original_rank;
        double triangle  = row_count * (row_count + 1.0) * 0.5;
        progress_display show_progress(triangle, taskName.c_str(), "", "");
        #endif

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

        double previewTime = 0;
        double mergeTime   = 0;
        for ( ; 1 < n ; --n) {
            T best_dist;
            recalculateTotals    (n, next_cluster_number, 
                                  row_cluster, cluster_row);
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
            //std::cout << "br=" << best_row << ", or=" << other_row << "\n";
            //std::cout << "bc=" << best_cluster << ", oc=" << other_cluster << "\n";
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
            std::cout << "Preview time " << previewTime << " and"
                      << " Cluster merge time was " << mergeTime << ".\n";
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
    virtual bool loadMatrixFromOpenFile(InFile& in) {
        return false;
    }
    virtual void setRank(size_t n) {
        original_rank       = n;
        next_cluster_number = n;
        storage.resize(original_rank*original_rank*size_t(2));
        size_t q = n+n-2; //number of clusters that will be needed
                          //in total, during the course of the tree
                          //construction: n leaves, n-2 interiors.
        cluster_alive.resize         (q, 1);
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
    void recalculateTotals(int n, int next_cluster_num, 
                           IntVector& row_cluster,
                           IntVector& cluster_row) {
        double cutoff           = -infiniteDistance;
        double one_on_n_minus_2 = (n<3) ? 0.0 : (1.0 / ((double)n - 2.0));
        int    r                = 0;
        //std::cout << "\nRow Totals: ";
        for (int c=0; c<next_cluster_num; ++c) {
            if (0<cluster_alive[c]) {
                cluster_total_scaled[c] = cluster_total[c] * one_on_n_minus_2;
                cluster_cutoff[c]       = cutoff;
                if (cutoff < cluster_total_scaled[c] ) {
                    cutoff = cluster_total_scaled[c];
                }
                cluster_row[c] = r;
                row_cluster[r] = c;
                ++r;
                //std::cout << "c=" << c << " t=" << cluster_total_scaled[c] << " ";
            }
        }
        //std::cout << "\n";
    }
    void previewRows(int n, int q,
                     const IntVector& row_cluster,
                     DistanceVector&  row_raw_dist,
                     DistanceVector&  row_best_dist,
                     IntVector&       row_choice,
                     T&               best_dist) {
        #ifdef _OPENMP
        #pragma omp parallel for if(1<threadCount)
        #endif
        for (int r = 0; r < n; ++r ) {
            int  c     = row_cluster[r];
            auto scan  = cluster_sorted_start[c]; 
            auto stop  = cluster_sorted_stop[c];
            bool found = false;
            for (; scan<stop; ++scan) {
                int b = scan->cluster_num;
                if (0<cluster_alive[b]) {
                    row_raw_dist[r]  = scan->distance;
                    row_best_dist[r] = scan->distance
                                     - cluster_total_scaled[c]
                                     - cluster_total_scaled[b];
                    row_choice[r]    = b;
                    cluster_sorted_start[c] = scan;
                    found = true;
                    break;
                }
            }
            if (!found) {
                    row_raw_dist[r]  = infiniteDistance;
                    row_best_dist[r] = infiniteDistance;
                    row_choice[r]    = q;
            }
        }

        best_dist = row_best_dist[0];
        for (int r=1; r<n; ++r) {
            if (row_best_dist[r]<best_dist) {
                best_dist = row_best_dist[r];
            }
        }
    }
    void findPreferredPartners(int n, int q, T global_best_dist,
                               const IntVector& row_cluster,
                               DistanceVector&  row_raw_dist,
                               DistanceVector&  row_best_dist,
                               IntVector&       row_choice) {
        //For each cluster, find preferred partner
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic) if(1<threadCount)
        #endif
        for (int r = n-1; 0 <= r; --r ) {
            int c             = row_cluster[r];
            if (cluster_alive[c]==0) {
                continue;
            }
            int best_cluster   = row_choice[r];    //set by previewRows().
            T   best_raw_dist  = row_raw_dist[r];  //set by previewRows().
            T   best_dist      = row_best_dist[r]; //set by previewRows().
            T   cutoff         = best_raw_dist + cluster_cutoff[c];
            T   earlier_cutoff = global_best_dist 
                                + cluster_cutoff[c] + cluster_total[c];
            if (earlier_cutoff < cutoff) {
                cutoff = earlier_cutoff;
            }
            //Can skip the first MatrixEntry at cluster_sorted_start[c].
            //It's already been looked at by previewRows().
            auto dataStart = cluster_sorted_start[c] + 1; 
            auto dataStop  = cluster_sorted_stop[c];

#if (0)
            if (r+r < dataStop-dataStart) {
                //At most half of the information recorded for this cluster
                //is still of any use.  Purge it.
                auto oldStop = dataStop;
                for (auto scan=dataStart; scan<oldStop; ++scan) {
                    *dataStop = *scan;
                    dataStop += cluster_alive[scan->cluster_num];
                }
                cluster_sorted_stop[c] = dataStop;
            }
#endif 

            for (auto scan=dataStart; scan<dataStop; ++scan) {
                int p                = scan->cluster_num;
                if (cluster_alive[p]==0) {
                    continue;
                }
                T   dist_raw         = scan->distance;
                if (cutoff < dist_raw) {
                    break;
                }
                T   dist_half_cooked = dist_raw - cluster_total_scaled[p];
                if (best_dist<=dist_half_cooked) {
                    continue;
                }
                best_raw_dist = dist_raw;
                best_dist     = dist_half_cooked;
                best_cluster  = p;
                cutoff        = best_raw_dist + cluster_cutoff[c];
                if (earlier_cutoff < cutoff) {
                    cutoff = earlier_cutoff;
                }
            }
            row_raw_dist[r]  = best_raw_dist;
            row_best_dist[r] = best_dist - cluster_total_scaled[c];
            row_choice[r]    = best_cluster;
            if (row_best_dist[r] < global_best_dist) {
                global_best_dist = row_best_dist[r];
            }
        }
        #if (0)
        for (int r = 0; r < n; ++r ) {
            int c             = row_cluster[r];
            std::cout << "c=" << c << " raw=" << row_raw_dist[r]
                      << " d=" << row_best_dist[r] << " o=" << row_choice[r] << "\n";
        }
        #endif
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
        cluster_alive[c]         = 0;
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
            int p = scan->cluster_num;
            distances[p] = (0<cluster_alive[p]) ? scan->distance : distances[p];
        }
        //Now, for clusters *after* c, for which 
        //distances to c were recorded.
        #ifdef _OPENMP
        #pragma omp parallel for if(1<threadCount)
        #endif
        for (int p=c+1; p<u; ++p) {
            if (cluster_alive[p]) {
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
                low    = guess[1].cluster_num;
                //looking at guess[1] would be naughty,
                //if we did not know that there *is* a
                //MatrixEntry with a cluster_num of a
                //somewhere to the right of guess.
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