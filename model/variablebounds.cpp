#include "variablebounds.h"

VariableBounds::VariableBounds(int dimension) : variable_count(dimension) {
    int doubles_needed = dimension * 4;
    variables   = new double[doubles_needed];
    variables2  = variables + dimension;
    upper_bound = variables + dimension*2;
    lower_bound = variables + dimension*3;
    for (int i=0; i < doubles_needed; ++i) {
        variables[i] = 0.0;
    }
    bound_check = new bool[dimension];
    for (int i=0; i<dimension; ++i) {
        bound_check[i] = false;
    }
}

VariableBounds::~VariableBounds() {
    delete [] variables;
    delete [] bound_check;
}
