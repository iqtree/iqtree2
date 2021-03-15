//
//  parsimonysearchparameters.h
//  Created by James Barbetti on 19-Jan-2021.
//

#ifndef parsimonysearchparameters_h
#define parsimonysearchparameters_h

#include "phylotree.h"
#include <utils/timekeeper.h>

class ParsimonySearchParameters
{
public:
    std::string name;
    intptr_t    iterations;
    bool        lazy_mode;
    int         radius;
    bool        calculate_connection_costs;
    bool        be_quiet;

    TimeKeeper overall;
    TimeKeeper initializing;
    TimeKeeper rescoring;
    TimeKeeper evaluating;
    TimeKeeper sorting;
    TimeKeeper applying;

    ParsimonySearchParameters() = delete;
    ParsimonySearchParameters(const char* move_name):
        name(move_name), be_quiet(false), overall(move_name),
        initializing("initializing"), rescoring("rescoring parsimony"),
        evaluating(std::string("evaluating ") + name + " moves"),
        sorting(std::string("sorting ")       + name + " moves"),
        applying(std::string("applying ")     + name + " moves") {
    }
    
    void report() {
        if (VB_MED <= verbose_mode && !be_quiet) {
            std::cout.precision(4);
            if (!progress_display::getProgressDisplay()) {
                overall.report();
            }
            initializing.report();
            rescoring.report();
            evaluating.report();
            sorting.report();
            applying.report();
        }

    }

};

#endif /* parsimonysearchparameters_h */
