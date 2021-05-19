#include "parsimonysearchparameters.h"
#include "phylotree.h"

ParsimonySearchParameters::ParsimonySearchParameters(const char* move_name):
    name(move_name), iterations(0), lazy_mode(false), radius(0),
    calculate_connection_costs(false), be_quiet(false), overall(move_name),
    initializing("initializing"), rescoring("rescoring parsimony"),
    evaluating(std::string("evaluating ") + name + " moves"),
    sorting(std::string("sorting ")       + name + " moves"),
    applying(std::string("applying ")     + name + " moves") {
}

void ParsimonySearchParameters::report() {
    if (VerboseMode::VB_MED <= verbose_mode && !be_quiet) {
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
