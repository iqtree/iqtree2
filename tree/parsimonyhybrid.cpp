//
//parsimonyhybrid.cpp
//Created by James Barbetti on 05-Mar-2021.
//

#include "parsimonyhybrid.h"
#include "parsimonysearch.h"

int PhyloTree::doParsimonyHybrid(VerboseMode how_loud) {
    if (leafNum<6) {
        return computeParsimony();
    }
    ParsimonySearchParameters s("Hybrid");
        
    s.iterations                 = params->parsimony_hybrid_iterations;
    s.lazy_mode                  = false;
    s.radius                     = params->spr_radius;
    s.calculate_connection_costs = false;
    s.be_quiet                   = verbose_mode < how_loud;
    
#if (0)
    //Adding TBR to the mix doesn't really seem to help.
    typedef ParsimonyHybridMove<ParsimonyNNIMove, ProperParsimonyTBRMove, false> NNI_Or_TBR;
    typedef ParsimonyHybridMove<ParsimonySPRMove, NNI_Or_TBR,             true>  PHM;
#endif
    
    typedef ParsimonyHybridMove<ParsimonySPRMove, ParsimonyNNIMove, false> PHM;

    return doParsimonySearch<PHM>(s);
}
