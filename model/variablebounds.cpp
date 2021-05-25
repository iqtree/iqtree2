#include "variablebounds.h"

VariableBounds::VariableBounds(int dimension) {
    variables  = new double[dimension*4];
    variables2  = variables + dimension;
    upper_bound = variables + dimension*2;
    lower_bound = variables + dimension*3;
    bound_check = new bool[dimension];
    for (int i=0; i<dimension; ++i) {
        variables[i]   = 0.0;
        variables2[i]  = 0.0;
        upper_bound[i] = 0.0;
        lower_bound[i] = 0.0;
        bound_check[i] = false;
    }
}

VariableBounds::~VariableBounds() {
    delete [] variables;
    delete [] bound_check;
}
