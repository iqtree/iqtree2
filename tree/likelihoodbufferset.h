//
// likelihoodbufferset.h
// A set of buffers, used for likelihood calculations, that are external to
// the nodes in a tree, that cannot be shared between threads that are
// independently calculating likelihoods on mutually inconsistent *possible*
// subtrees.
//
// Each instance of PhyloTree has its own LikelihoodBufferSet (tree_buffers),
// (and that is what is used, most of the time).
//
// But, during likelihod placement, likelihoodfs for many possible candidate
// subtrees (that differ only in one leaf and one interior node), potentially
// connected to the tree at the same point, are calculated, possibly in parallel.
// Each such calculation (for a different potential subtree) needs its own
// set of buffers.
//
// See in particular, LikelihoodCostCalculator::assessPlacementCost,
// in placement/placementcostcalculator.cpp.
//
//
// Created by James Barbetti on 9/10/20.
//

#ifndef likelihoodbufferset_h
#define likelihoodbufferset_h

#include <stdlib.h> //for size_t
        
class LikelihoodBufferSet {
public:
    /**
     *    NSTATES x NUMCAT x (number of patterns) array
     *    Used to store precomputed values when optimizing branch length
     *    See Tung's report on 07.05.2012 for more information
     */
    double* theta_all;
    size_t  theta_block_size;   //How big it is
    bool    theta_computed;     //True if it's content is computed
    bool    theta_borrowed;     //True if it belongs to something else
                                //(other than this instance)

    /**
            internal pattern log-likelihoods, always stored after calling computeLikelihood()
            or related functions. Note that scaling factors are not incorporated here.
            If you want to get real pattern log-likelihoods, please use computePatternLikelihood()
     */
    double *_pattern_lh; //Todo: not distributed
    size_t pattern_lh_block_size; //
    bool   pattern_lh_borrowed;   //True if it belongs to something else

    /**
            internal pattern likelihoods per category,
    */
    double *_pattern_lh_cat; //Todo: not distributed
    size_t pattern_lh_cat_block_size;
    bool   pattern_lh_cat_borrowed;
    
    /** buffer used when computing partial_lh, to avoid repeated mem allocation */
    double *buffer_partial_lh; //Todo: not distributed
    size_t partial_lh_block_size;
    bool   partial_lh_borrowed;
    
    /** total scaling buffer */

    double *buffer_scale_all;
    size_t scale_all_block_size;
    bool   scale_all_borrowed;
    
    LikelihoodBufferSet();
    
    LikelihoodBufferSet(const LikelihoodBufferSet& copyMe);
    
    void ensureThetaAllocated(size_t desired_block_size);

    void borrowTheta(double* theta_start, size_t theta_size_in_doubles);
    
    void ensurePatternLhAllocated(size_t desired_block_size_in_doubles);
    
    void borrowPatternLh(double* borrowMe, size_t size_in_doubles);
        
    void ensurePatternLhCatAllocated(size_t desired_block_size_in_doubles);
    
    void borrowPatternLhCat(double* borrowMe, size_t size_in_doubles);
    
    void ensurePartialLhAllocated(size_t size_in_doubles);
    
    void borrowPartialLh(double* borrowMe, size_t size_in_doubles);
    
    void ensureScaleAllAllocated(size_t size_in_doubles);
    
    void borrowScaleAll(double* borrowMe, size_t size_in_doubles);

    void forget();
    
    void freeBuffers();
    
    ~LikelihoodBufferSet();
};

#endif /* likelihoodbufferset_h */
