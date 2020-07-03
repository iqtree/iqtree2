//
//  superalignmentpairwiseplen.cpp
//  iqtree
//
//  Created by Olga on 04/05/17.
//
//

#include <stdio.h>
#include "superalignmentpairwiseplen.h"

/**********************************************************
 * class SuperAlignmentPairwisePlen
 **********************************************************/

SuperAlignmentPairwisePlen::SuperAlignmentPairwisePlen(PhyloSuperTreePlen *atree, int seq1, int seq2)
: SuperAlignmentPairwise((PhyloSuperTree*) atree, seq1, seq2)
{
    part_info = &(atree->part_info);
}

double SuperAlignmentPairwisePlen::computeFunction(double value) {
    int part = 0;
    double lh = 0.0;
    for (auto it = partitions.begin(); it != partitions.end(); it++, part++) {
        lh += it->computeFunction(part_info->at(part).part_rate*value);
    }
    return lh;
}

void SuperAlignmentPairwisePlen::computeFuncDerv(double value, double &df, double &ddf) {
    int part = 0;
    df = 0.0;
    ddf = 0.0;
    for (auto it = partitions.begin(); it != partitions.end(); it++, part++) {
        double d1, d2;
        it->computeFuncDerv(part_info->at(part).part_rate*value, d1, d2);
        df += part_info->at(part).part_rate*d1;
        ddf += part_info->at(part).part_rate*part_info->at(part).part_rate*d2;
    }
}

SuperAlignmentPairwisePlen::~SuperAlignmentPairwisePlen()
{}
