//
//  parsimonysearchparameters.h
//  Created by James Barbetti on 19-Jan-2021.
//

#ifndef parsimonysearchparameters_h
#define parsimonysearchparameters_h

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
    explicit ParsimonySearchParameters(const char* move_name);
    
    void report();
};

#endif /* parsimonysearchparameters_h */
