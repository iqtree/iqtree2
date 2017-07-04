//
// C++ Interface: optimization
//
// Description:
//
//
// Author: BUI Quang Minh, Steffen Klaere, Arndt von Haeseler <minh.bui@univie.ac.at>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

#include <iostream>

/**
Optimization class, implement some methods like Brent, Newton-Raphson (for 1 variable function), BFGS (for multi-dimensional function)

	@author BUI Quang Minh, Steffen Klaere, Arndt von Haeseler <minh.bui@univie.ac.at>
*/
class Optimization{
public:

    /** constructor */
    Optimization();

	/*****************************************************
		One dimensional optimization with Brent method
	*****************************************************/
	/**
		This function calculate f(value) of the f() function, used by other general optimization method to minimize it.
		Please always override this function to adapt to likelihood or parsimony score.
		The default is for function f(x)=x.
		@param value x-value of the function
		@return f(value) of function f you want to minimize
	*/
	virtual double computeFunction(double value) { return value; }

	/**
		the brent method to find the value that minimizes the computeFunction().
		@return the x-value that minimize the function
		@param xmin lower bound
		@param xmax upper bound
		@param xguess first guess
		@param tolerance tolerance
		@param fx (OUT) function value at the minimum x found
		@param ferror (OUT) Dont know
	*/
	double minimizeOneDimen(double xmin, double xguess, double xmax, double tolerance, double *fx, double *ferror);

	double minimizeOneDimenSafeMode(double xmin, double xguess, double xmax, double tolerance, double *fx);

	/*****************************************************
		One dimensional optimization with Newton Raphson
		only applicable if 1st and 2nd derivatives are easy to compute
	*****************************************************/

	/**
		This function calculate f(value), first derivative f'(value) and 2nd derivative f''(value).
		used by Newton raphson method to minimize the function.
		Please always override this function to adapt to likelihood or parsimony score.
		The default is for function f(x) = x^2.
		@param value x-value of the function
		@param df (OUT) first derivative
		@param ddf (OUT) second derivative
	*/
	virtual void computeFuncDerv(double value, double &df, double &ddf) {
		df = 2.0*value; ddf = 2.0;
	}

	/**
		Newton-Raphson method to minimize computeFuncDerv()
		@return the x-value that minimize the function
		@param xmin lower bound
		@param xmax upper bound
		@param xguess first guess
		@param tolerance tolerance of x-value to stop the iterations
		@param fx (OUT) function value at the minimum x found
		@param var (OUT) variance estimate of x
		@param maxNRStep max number of NR steps
	*/
	double minimizeNewton(double xmin, double xguess, double xmax, double tolerance, int maxNRStep = 100);

	double minimizeNewton(double xmin, double xguess, double xmax, double tolerance, double &d2l, int maxNRStep = 100);

	double minimizeNewtonSafeMode(double xmin, double xguess, double xmax, double tolerance, double &f);

	/*****************************************************
		MULTI dimensional optimization with Newton Raphson
		only applicable if 1st and 2nd derivatives are easy to compute
	*****************************************************/

	/**
		This function calculate f(value), first derivative f'(value) and 2nd derivative f''(value).
		used by Newton raphson method to minimize the function.
		Please always override this function to adapt to likelihood or parsimony score.
		The default is for function f(x) = x^2.
		@param value x-value of the function
		@param df (OUT) first derivative, with df[N] being function value!
		@param ddf (OUT) second derivative
	*/
	virtual void computeFuncDervMulti(double *value, double *df, double *ddf) {
        std::cerr << "Please override computeFuncDervMulti" << std::endl;
	}

	/**
		Newton-Raphson method to minimize computeFuncDervMulti()
		@param xmin lower bound
		@param xmax upper bound
		@param xguess[in/out] first guess
		@param tolerance tolerance of x-value to stop the iterations
        @param N number of variables
		@param maxNRStep max number of NR steps
        @return function value at optimum
	*/
	double minimizeNewtonMulti(double *xmin, double *xguess, double *xmax, double tolerance, int N, int maxNRStep = 100);


//	double rtsafe(double x1, double xguess, double x2, double xacc, double &f);

	/*****************************************************
		Multi dimensional optimization with BFGS method
	*****************************************************/

	/**
		return the number of dimensions
	*/
	virtual int getNDim() { return 0; }


	/**
		the target function which needs to be optimized
		@param x the input vector x
		@return the function value at x
	*/
	virtual double targetFunk(double x[]) { return 0.0; }

	/**
		the approximated derivative function
		@param x the input vector x
		@param dfx the derivative at x
		@return the function value at x
	*/
	virtual double derivativeFunk(double x[], double dfx[]);

	/**
	        Controls restarting of optimization if optimization gets
                stuck on the boundary. Models are free to override this
		to provide their own boundary avoidance strategies.
		@param guess the current optimized parameters, and
                       will get overwritten with new search start point
                       if required. (I.e. this is both input and output)
                @param ndim the number of dimensions
		@param lower lower bounds on parameters
		@param upper upper bounds on parameters
		@param bound_check whether to trigger a restart if
                       corresponding parameter is on boundary
		@param iteration how many times we've called restartParameters
		       for this model. Should start at 1, increase by one each 
		       call, until restartParameters returns false.
		@return true if guess has changed, false if not.
	 */
	virtual bool restartParameters(double guess[], int ndim, double lower[], double upper[], bool bound_check[], int iteration);

	/**
		multi dimensional optimization by BFGS method
		@param guess the initial starting point
		@param ndim number of dimension
		@param gtol tolerance
		@param lower the lower bound vector
		@param upper the upper bound vector
		@param bound_check bound checking vector
		@return the minimum function value obtained
	*/
	double minimizeMultiDimen(double guess[], int ndim, double lower[], double upper[],
		bool bound_check[], double gtol, double *hessian = NULL);

	/*****************************************************
		NEW 2015-08-19: Multi dimensional optimization with L-BFGS-B method
	*****************************************************/

    /**
     Function to access the L-BFGS-B function, taken from HAL_HAS software package
     
     1. int nvar : The number of the variables
     2. double* vars : initial values of the variables
     3. double* lower : lower bounds of the variables
     4. double* upper : upper bounds of the variables
     5. double pgtol: gradient tolerance
     5. int maxit : max # of iterations
     @return minimized function value
     After the function is invoked, the values of x will be updated
    */
    double L_BFGS_B(int nvar, double* vars, double* lower, double* upper, double pgtol = 1e-5, int maxit = 5); // changed maxit 1000 -> 5 by Thomas on Sept 11, 15

    /** internal function called by L_BFGS_B
        should return function value 
        @param nvar number of variables
        @param vars variables
    */
    virtual double optimFunc(int nvar, double *vars);
    
    /** internal function called by L_BFGS_B
        should return gradient value
        @param nvar number of variables
        @param vars variables
        @param gradient (OUT) function gradient
        @return function value
    */
    virtual double optimGradient(int nvar, double *vars, double *gradient);
    

    ~Optimization();

	/**
		original numerical recipes method
	*/
	double brent(double ax, double bx, double cx, double tol, double *xmin);

private:


	double brent_opt (double ax, double bx, double cx, double tol,
		double *foptx, double *f2optx, double fax, double fbx, double fcx);

	double dbrent(double ax, double bx, double cx, double tol, double *xmin);

	void dfpmin(double p[], int n, double lower[], double upper[], double gtol, int *iter, double *fret, double *hessian);

	void lnsrch(int n, double xold[], double fold, double g[], double p[], double x[],
		double *f, double stpmax, int *check, double lower[], double upper[]);

    void lbfgsb(int n, int m, double *x, double *l, double *u, int *nbd,
		double *Fmin, int *fail,
		double factr, double pgtol,
		int *fncount, int *grcount, int maxit, char *msg,
		int trace, int nREPORT);
    
};


void nrerror(const char *error_text);


#endif
