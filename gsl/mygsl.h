
/**
	stripped GSL (GNU Scientific library) used for IQ-TREE code
*/

#ifndef _MYGSL_H
#define _MYGSL_H

#include <stdio.h>

/*
    x power n (x^n)
    @return x^n
*/
double gsl_pow_uint(double x, unsigned int n);

/*
    binomial sampling
    @param p probability
    @param n sample size
    @return random value drawn from binominal distribution 
*/
unsigned int gsl_ran_binomial (double p, unsigned int n, int *rstream);

/*
    multinomial sampling
    @param K number of categories
    @param N sample size
    @param p probability vector of length K, will be normalized to 1 if not summing up to 1
    @param[out] n output vector of length K as drawn from multinomial distribution, sum to N
*/
void gsl_ran_multinomial (const size_t K, const unsigned int N, const double p[], unsigned int n[], int *rstream);


/*
    probability density function for standard normal distribution
    @param x x-value
    @return probability density p(x)
*/
double gsl_ran_ugaussian_pdf (const double x);

/*
    cumulative distribution function for standard normal distribution 
    @param x x-value
    @return CDF at x
*/
double gsl_cdf_ugaussian_P (const double x);

/*
 1.0 - cumulative distribution function for standard normal distribution
 @param x x-value
 @return 1.0-CDF at x
 */
double gsl_cdf_ugaussian_Q (const double x);

/*
    quantile function for standard normal distribution (or CDF-inverse function)
    @param P probability value
    @return x-value
*/
double gsl_cdf_ugaussian_Pinv (const double P);

#endif

