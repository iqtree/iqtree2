//
// likelihoodbufferset.cpp
// Implementation of the LikelihoodBufferSet class
//
// Created by James Barbetti on 09-Oct-2020.
//

#include "likelihoodbufferset.h"
#include "alignedalloc.h"
    
LikelihoodBufferSet::LikelihoodBufferSet():
    theta_all(nullptr),         theta_block_size(0),
    theta_computed(false),      theta_borrowed(false),
    _pattern_lh(nullptr),       pattern_lh_block_size(0),     pattern_lh_borrowed(false),
    _pattern_lh_cat(nullptr),   pattern_lh_cat_block_size(0), pattern_lh_cat_borrowed(false),
    buffer_partial_lh(nullptr), partial_lh_block_size(0),     partial_lh_borrowed(false),
    buffer_scale_all(nullptr),  scale_all_block_size(0),      scale_all_borrowed(false) {
}

LikelihoodBufferSet::LikelihoodBufferSet(const LikelihoodBufferSet& copyMe) {
    theta_all = nullptr;
    ensureThetaAllocated(copyMe.theta_block_size);

    _pattern_lh = nullptr;
    ensurePatternLhAllocated(copyMe.pattern_lh_block_size);
    
    _pattern_lh_cat = nullptr;
    ensurePatternLhCatAllocated(copyMe.pattern_lh_cat_block_size);
    if (_pattern_lh_cat != nullptr ) {
        memcpy(_pattern_lh_cat, copyMe._pattern_lh_cat,
               copyMe.pattern_lh_cat_block_size * sizeof(double));
    }
    
    buffer_partial_lh = nullptr;
    ensurePartialLhAllocated(copyMe.partial_lh_block_size);
    
    buffer_scale_all = nullptr;
    ensureScaleAllAllocated(copyMe.scale_all_block_size);
}

void LikelihoodBufferSet::ensureThetaAllocated(size_t desired_block_size) {
    if (theta_all != nullptr) {
        return;
    }
    theta_all        = aligned_alloc<double>(desired_block_size);
    theta_block_size = desired_block_size;
    theta_borrowed   = false;
    for (int i=0; i<theta_block_size; ++i) theta_all[i]=0;
}

void LikelihoodBufferSet::borrowTheta(double* theta_start, size_t theta_size_in_doubles) {
    if (!theta_borrowed) {
        aligned_free(theta_all);
    }
    theta_all        = theta_start;
    theta_block_size = theta_size_in_doubles;
    theta_borrowed   = true;
}

void LikelihoodBufferSet::ensurePatternLhAllocated(size_t desired_block_size_in_doubles) {
    if (!_pattern_lh) {
        _pattern_lh = aligned_alloc<double>(desired_block_size_in_doubles);
        pattern_lh_block_size = desired_block_size_in_doubles;
        pattern_lh_borrowed = false;
        for (int i=0; i<pattern_lh_block_size; ++i) _pattern_lh[i]=0;
    }
}

void LikelihoodBufferSet::borrowPatternLh(double* borrowMe, size_t size_in_doubles) {
    if (_pattern_lh && !pattern_lh_borrowed) {
        aligned_free(_pattern_lh);
    }
    _pattern_lh = borrowMe;
    pattern_lh_block_size = size_in_doubles;
    pattern_lh_borrowed = true;
}
    
void LikelihoodBufferSet::ensurePatternLhCatAllocated(size_t desired_block_size_in_doubles) {
    if (!_pattern_lh_cat) {
        _pattern_lh_cat           = aligned_alloc<double>(desired_block_size_in_doubles);
        pattern_lh_cat_block_size = desired_block_size_in_doubles;
        pattern_lh_cat_borrowed   = false;
        for (int i=0; i<pattern_lh_cat_block_size; ++i) _pattern_lh_cat[i]=0;
    }
}

void LikelihoodBufferSet::borrowPatternLhCat(double* borrowMe, size_t size_in_doubles) {
    if (_pattern_lh_cat && !pattern_lh_cat_borrowed) {
        aligned_free(_pattern_lh_cat);
    }
    _pattern_lh_cat           = borrowMe;
    pattern_lh_cat_block_size = size_in_doubles;
    pattern_lh_cat_borrowed   = true;
}

void LikelihoodBufferSet::ensurePartialLhAllocated(size_t size_in_doubles) {
    if (!buffer_partial_lh) {
        buffer_partial_lh = aligned_alloc<double>(size_in_doubles);
        partial_lh_block_size   = size_in_doubles;
        pattern_lh_cat_borrowed = false;
        for (int i=0; i<partial_lh_block_size; ++i) buffer_partial_lh[i]=0;
    }
}

void LikelihoodBufferSet::borrowPartialLh(double* borrowMe, size_t size_in_doubles) {
    if (buffer_partial_lh && !partial_lh_borrowed) {
        aligned_free(buffer_partial_lh);
    }
    buffer_partial_lh     = borrowMe;
    partial_lh_block_size = (size_in_doubles==0) ? partial_lh_block_size : size_in_doubles;
    partial_lh_borrowed   = true;
}

void LikelihoodBufferSet::ensureScaleAllAllocated(size_t size_in_doubles) {
    if (!buffer_scale_all) {
        buffer_scale_all     = aligned_alloc<double>(size_in_doubles);
        scale_all_block_size = size_in_doubles;
        scale_all_borrowed   = false;
        for (int i=0; i<scale_all_block_size; ++i) buffer_scale_all[i]=0;
    }
}

void LikelihoodBufferSet::borrowScaleAll(double* borrowMe, size_t size_in_doubles) {
    if (buffer_scale_all && !scale_all_borrowed) {
        aligned_free(buffer_partial_lh);
    }
    buffer_partial_lh    = borrowMe;
    scale_all_block_size = size_in_doubles;
    scale_all_borrowed   = true;
}

void LikelihoodBufferSet::forget() {
    theta_all                 = nullptr;
    theta_borrowed            = false;
    _pattern_lh               = nullptr;
    pattern_lh_borrowed       = false;
    _pattern_lh_cat           = nullptr;
    pattern_lh_cat_borrowed   = false;
    buffer_partial_lh         = nullptr;
    partial_lh_borrowed       = false;
    buffer_scale_all          = nullptr;
    scale_all_block_size      = 0;
}

void LikelihoodBufferSet::freeBuffers() {
    if (!theta_borrowed) {
        aligned_free(theta_all);
    } else {
        theta_all = nullptr;
    }
    theta_borrowed   = false;
    theta_block_size = 0;
    
    if (!pattern_lh_borrowed) {
        aligned_free(_pattern_lh);
    } else {
        _pattern_lh = nullptr;
    }
    pattern_lh_borrowed   = false;
    pattern_lh_block_size = 0;
    
    if (!pattern_lh_cat_borrowed) {
        aligned_free(_pattern_lh_cat);
    } else {
        _pattern_lh_cat = nullptr;
    }
    pattern_lh_cat_borrowed   = false;
    pattern_lh_cat_block_size = 0;

    if (!partial_lh_borrowed) {
        aligned_free(buffer_partial_lh);
    } else {
        buffer_partial_lh = nullptr;
    }
    partial_lh_borrowed   = false;
    partial_lh_block_size = 0;

    if (!scale_all_borrowed) {
        aligned_free(buffer_scale_all);
    } else {
        buffer_scale_all  = nullptr;
    }
    scale_all_borrowed    = false;
    scale_all_block_size  = 0;
}

LikelihoodBufferSet::~LikelihoodBufferSet() {
    freeBuffers();
}
