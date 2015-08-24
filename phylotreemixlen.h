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


class PhyloTreeMixlen : public IQTree {

public:

    /**
            default constructor
     */
    PhyloTreeMixlen();

    PhyloTreeMixlen(Alignment *aln, int mixlen);

    /**
        set number of mixture branch lengths
    */
    void setMixlen(int mixlen);

    int mixlen;

};

#endif /* defined(__iqtree__phylotreemixlen__) */
