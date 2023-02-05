//
//  HashRow<T> type for keeping track of rows of a distance matrix of T
//  (hashing them, to make it cheaper to determine which rows
//  are identical).
//
//  DuplicateTaxa type, for reporting the row numbers of rows
//  that indicate identical taxa (or at least, taxa whose inter-taxa
//  differences are the same, and that are all distance 0 from each 
//  other).
//
//  Copyright James Barbetti (2021-22)
//

#pragma once
#ifndef hashrow_h
#define hashrow_h

#include <inttypes.h> //for intptr_t
#include <vector>     //for std::vector template class
#include <algorithm>  //for std::sort template function

typedef std::vector< std::vector< intptr_t > > DuplicateTaxa;

/**
 * @brief  A pair, of a hash, and contiguous block of
 *         something (T) that is hashable and comparable
 *         via both operator!= and operator<.
 *         Ordering is hash primary, with lexicographic 
 *         block order to tie-break.
 * @tparam T something hashable (via std::hash) and 
 *         comparable.  
 * @note   if blocks are different sizes, and the contents
 *         of the blocks compare equal (to the size of the
 *         smaller block, the shorter block is treated as
 *         being less than the longer).
 */
template <class T> class HashRow {
public:
    intptr_t row_num;    //an identifying row number (not used for ordering)
    const T* row_data;   //the data block for the row
    size_t   row_length; //the size of the data block
    size_t   row_hash;   //the hash of the elements in the block

    HashRow(): row_num(-1), row_data(nullptr), row_length(0), row_hash(0) {}
    HashRow(intptr_t num, const T* row_start, size_t length)
        : row_num(num), row_data(row_start), row_length(length), row_hash(0) {
        for (size_t col=0; col<length;++col ) {
            row_hash ^= std::hash<double>()(static_cast<double>(row_data[col]))
                + 0x9e3779b9 + (row_hash<<6) + (row_hash>>2);
        }
    }
    HashRow(const HashRow& rhs) = default;
    HashRow& operator= (const HashRow& rhs) = default;
    
    /**
     * @brief  Compare *this with rhs
     * @param  rhs 
     * @return int -1 if *this compares less, +1 if rhs compares less,
     *         and 0, if *this and rhs compare equal
     */
    int compare(const HashRow& rhs) const {
        //First compare hashes
        if (row_hash<rhs.row_hash) { return -1; }
        if (rhs.row_hash<row_hash) { return 1;  }
        int    count_diff    = 0;
        size_t min_col_count = row_length; 
        if (row_length != rhs.row_length) {
            bool left_less = (row_length < rhs.row_length);
            count_diff     = left_less ? -1 : 1;
            min_col_count  = left_less ? row_length : rhs.row_length;
        }
        //compare elements in the blocks, up to the size of
        //the smaller block.
        for ( size_t col=0; col<min_col_count ; ++col) {
            if (row_data[col]!=rhs.row_data[col]) {
                return (row_data[col]<rhs.row_data[col]) ? -1 : 1;
            }
        }
        return count_diff;
    }
    bool operator< (const HashRow& rhs) const {
        return compare(rhs)<0;
    }
    bool operator<= (const HashRow& rhs) const {
        return compare(rhs)==0;
    }
    
    /**
     * @brief  Given a sorted vector of HashRow<T>, that represents the
     *         rows of a distance matrix, determine which groups of rows
     *         in the distance matrix are duplicates.
     * @param  hashed_rows the hashed distance matrix rows, in sorted order
     * @param  vvc reference to a DuplicateTaxa instance (a vector of 
     *             vectors of intptr_t).  Each of the vectors will 
     *             contains the row numbers of the (duplicate) rows
     *             in an equivalence class (in no particular order).
     * @note   it assumed (but not checked) that hashed_rows has been sorted
     * @note   row numbers, within an equivalence class, appear in vvc
     *         in their order of appearance in hashed_rows.
     */
    static void identifyDuplicateClusters(const std::vector< HashRow<T> >& hashed_rows,
                                          DuplicateTaxa& vvc) {
        std::vector< intptr_t> vc; //vector of cluster #s
        intptr_t row_count = hashed_rows.size();
        for (intptr_t i=1; i<row_count; ++i) {
            bool is_duplicate = hashed_rows[i].compare(hashed_rows[i-1])==0;
            if (is_duplicate) {
                intptr_t   h = hashed_rows[i-1].row_num;
                is_duplicate = hashed_rows[i].row_data[h] == 0 ;
            }
            if (!is_duplicate) {
                //Not a duplicate of the previous row.
                if (!vc.empty()) {
                    //Add vector of the clusters in previous group of
                    //duplicates, to vvc.  And start a new one (empty).
                    vvc.emplace_back(vc);
                    vc.clear();
                }
            } else { 
                //duplicate of previous row
                if (vc.empty()) {
                    //This row and the previous one are to be the first
                    //two rows whose clusters belong in the current
                    //vector of duplicate clusters, which is empty.
                    //Add the cluster that corresponds to the previous row.
                    vc.push_back(hashed_rows[i-1].row_num);
                }
                vc.push_back(hashed_rows[i].row_num);
            }
        }
        if (!vc.empty()) {
            vvc.emplace_back(vc);
        }
    }
};

#endif