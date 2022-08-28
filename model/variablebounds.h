#pragma once
#ifndef variable_bounds_h
#define variable_bounds_h

class VariableBounds {
public:
	int     variable_count;
	double* variables;
	double* variables2;
	double* upper_bound;
	double* lower_bound;
	bool*   bound_check;
	explicit VariableBounds(int dimension);
	~VariableBounds();
};

#endif //variable_bounds_h
