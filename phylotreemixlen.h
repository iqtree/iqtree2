//
//  phylotreemixlen.h
//  iqtree
//
//  Created by Minh Bui on 24/08/15.
//
//

#ifndef __iqtree__phylotreemixlen__
#define __iqtree__phylotreemixlen__

#include <stdio.h>
#include "iqtree.h"


/**
    Phylogenetic tree with mixture of branch lengths
    Started within joint project with Stephen Crotty
*/
class PhyloTreeMixlen : public IQTree {

public:

    /**
            default constructor
     */
    PhyloTreeMixlen();

    PhyloTreeMixlen(Alignment *aln, int mixlen);

    /**
        @return true if this is a tree with mixture branch lengths, default: false
    */
    virtual bool isMixlen() { return true; }

    /**
        set number of mixture branch lengths
    */
    void setMixlen(int mixlen);

    int mixlen;

};

#endif /* defined(__iqtree__phylotreemixlen__) */
