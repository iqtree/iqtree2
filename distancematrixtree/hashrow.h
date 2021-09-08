#pragma once
#ifndef hashrow_h
#define hashrow_h

#include <inttypes.h> //for intptr_t
#include <vector>     //for std::vector template class
#include <algorithm>  //for std::sort template function

typedef std::vector< std::vector< intptr_t > > DuplicateTaxa;

template <class T> class HashRow {
public:
    intptr_t row_num;
    const T* row_data;
    size_t   row_length;
    size_t   row_hash;
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
    

    int compare(const HashRow& rhs) const {
        if (row_hash<rhs.row_hash) { return -1; }
        if (rhs.row_hash<row_hash) { return 1;  }
        for ( size_t col=0; col<row_length ; ++col) {
            if (row_data[col]!=rhs.row_data[col]) {
                return (row_data[col]<rhs.row_data[col]) ? -1 : 1;
            }
        }
        return 0;
    }
    bool operator< (const HashRow& rhs) const {
        return compare(rhs)<0;
    }
    bool operator<= (const HashRow& rhs) const {
        return compare(rhs)==0;
    }
    
    static void identifyDuplicateClusters(const std::vector< HashRow<T> >& hashed_rows,
                                          DuplicateTaxa& vvc) {
        std::vector< intptr_t> vc; //vector of cluster #s
        intptr_t row_count = hashed_rows.size();
        for (intptr_t i=1; i<row_count; ++i) {
            bool is_duplicate = hashed_rows[i].compare(hashed_rows[i-1])==0;
            if (is_duplicate) {
                intptr_t   h = hashed_rows[i-1].row_num;
                is_duplicate = hashed_rows[i].row_data[h] ==0;
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