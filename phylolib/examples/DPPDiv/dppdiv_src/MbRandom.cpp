/*! 
 * \file
 * This file contains the implementation of MbRandom, a
 * class that deals with random variables. 
 *
 * \brief Implementation of MbRandom class 
 *
 * MrBayes version 4.0 beta
 *
 * (c) Copyright 2005.
 * \version 4.0 beta
 * \date Last modified: $Date: 2006/09/11 17:29:51 $
 * \author John Huelsenbeck (1)
 * \author Bret Larget (2)
 * \author Paul van der Mark (3)
 * \author Fredrik Ronquist (3)
 * \author Donald Simon (4)
 * \author (authors listed in alphabetical order)
 * (1) Department of Integrative Biology, University of California, Berkeley
 * (2) Departments of Botany and of Statistics, University of Wisconsin - Madison
 * (3) School of Computational Science, Florida State University
 * (4) Department of Mathematics/Computer Science, Duquesne University
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License (the file gpl.txt included with this
 * distribution or http://www.gnu.org/licenses/gpl.txt for more
 * details.
 *
 * $Id: MbRandom.cpp,v 1.3 2006/09/11 17:29:51 paulvdm Exp $
 */

#include <cmath>
#include <ctime>
#include <iostream>
#include <cassert>
#include <cstdio>

#include "MbRandom.h"
#include "MbVector.h"

using namespace std;



#pragma mark Constructors

/*!
 * Constructor for MbRandom class. This constructor does not take
 * any parameters and initializes the seed using the current 
 * system time.
 *
 * \brief Constructor for MbRandom, initializing seed with system time.
 * \param Takes no parameter.
 * \return Returns no value.
 * \throws Does not throw an error.
 */
MbRandom::MbRandom(void) {

	setSeed();
	initializedFacTable = false;
	availableNormalRv = false;
	
}

/*!
 * Constructor for MbRandom class. This constructor takes as a
 * parameter a long integer with the user-supplied seed for the
 * random number generator. 
 *
 * \brief Constructor for MbRandom, initializing seed with a user-supplied long integer.
 * \param x is a long integer with the user-supplied random number seed.
 * \return Returns no value.
 * \throws Does not throw an error.
 */
MbRandom::MbRandom(seedType x) {

	setSeed(x, 0);
	initializedFacTable = false;
	availableNormalRv = false;
	
}

#pragma mark Chi-Square Distribution

/*!
 * This function generates a chi square distributed random variable.
 *
 * \brief Chi-square random variable.
 * \param v is the degrees of freedom parameter of the chi-square distribution. 
 * \return Returns a chi-square distributed random variable.
 * \throws Does not throw an error.
 */
double MbRandom::chiSquareRv(double v) {

	/* Cast the degrees of freedom parameter as an integer. We will see
       if there is a decimal remainder later. */
	int n = (int)(v);
	
	double x2;
	if ( (double)(n) == v && n <= 100 )
		{
		/* If the degrees of freedom is an integer and less than 100, we
		   generate our chi-square random variable by generating v
		   standard normal random variables, squaring each, and taking the
		   sum of the squared random variables. */
		x2 = 0.0;
		for (int i=0; i<n; i++)
			{
			double x = normalRv();
			x2 += x * x;
			}
		}
	else
		{
		/* Otherwise, we use the relationship of the chi-square to a gamma
		   (it is a special case of the gamma) to generate the chi-square
		   random variable. */
		x2 = gammaRv(v/2.0, 0.5);
		}
	return x2;
	
}

/*!
 * This function calculates the probability density 
 * for a chi-square distributed random variable.
 *
 * \brief Chi-square probability density.
 * \param v is the degrees of freedom parameter of the chi-square. 
 * \param x is the chi-square random variable. 
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
double MbRandom::chiSquarePdf(double v, double x) {

	double pdf;
	if ( x < 0.0 )
		{
		pdf = 0.0;
		}
	else
		{
		double b = v / 2.0;
		pdf = exp ( -0.5 * x ) * pow ( x, ( b - 1.0 ) ) / ( pow ( 2.0, b ) * gamma ( b ) );
		}
	return pdf;
	
}

/*!
 * This function calculates the cumulative probability  
 * for a chi-square distributed random variable.
 *
 * \brief Chi-square cumulative probability.
 * \param v is the degrees of freedom parameter of the chi-square. 
 * \param x is the chi-square random variable. 
 * \return Returns the cumulative probability.
 * \throws Does not throw an error.
 */
double MbRandom::chiSquareCdf(double v, double x) {

	return gammaCdf( v / 2.0, 0.5, x );
	
}

/*!
 * This function calculates the natural log of the probability density 
 * for a chi-square distributed random variable.
 *
 * \brief Natural log of chi-square probability density.
 * \param v is the degrees of freedom parameter of the chi-square. 
 * \param x is the chi-square random variable. 
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
double MbRandom::lnChiSquarePdf(double v, double x) {

	double b = v / 2.0;
	return ( -(b * log(2.0) + lnGamma(b)) - b + (b - 1.0) * log(x) );
	
}

/*!
 * This function calculates the quantile of a chi square distribution with v
 * degrees of freedom.
 *
 * \brief Quantile of a chi square distribution.
 * \param v is the degrees of freedom of the chi square. 
 * \param prob is the probability up to the quantile. 
 * \return Returns quantile value (or -1 if in error). 
 * \throws Does not throw an error.
 */
double MbRandom::chiSquareQuantile(double prob, double v) {

	double 		e = 0.5e-6, aa = 0.6931471805, p = prob, g,
					xx, c, ch, a = 0.0, q = 0.0, p1 = 0.0, p2 = 0.0, t = 0.0, 
					x = 0.0, b = 0.0, s1, s2, s3, s4, s5, s6;

	if (p < 0.000002 || p > 0.999998 || v <= 0.0) 
		return (-1.0);
	g = lnGamma(v/2.0);
	xx = v/2.0;   
	c = xx - 1.0;
	if (v >= -1.24*log(p)) 
		goto l1;
	ch = pow((p*xx*exp(g+xx*aa)), 1.0/xx);
	if (ch-e < 0) 
		return (ch);
	goto l4;
	l1:
		if (v > 0.32) 
			goto l3;
		ch = 0.4;   
		a = log(1.0-p);
	l2:
		q = ch;  
		p1 = 1.0+ch*(4.67+ch);  
		p2 = ch*(6.73+ch*(6.66+ch));
		t = -0.5+(4.67+2.0*ch)/p1 - (6.73+ch*(13.32+3.0*ch))/p2;
		ch -= (1.0-exp(a+g+0.5*ch+c*aa)*p2/p1)/t;
		if (fabs(q/ch-1.0)-0.01 <= 0.0) 
			goto l4;
		else                       
			goto l2;
	l3: 
		x = pointNormal (p);
		p1 = 0.222222/v;   
		ch = v*pow((x*sqrt(p1)+1.0-p1), 3.0);
		if (ch > 2.2*v+6.0)  
			ch = -2.0*(log(1.0-p)-c*log(0.5*ch)+g);
	l4:
		q = ch;   
		p1 = 0.5*ch;
		if ((t = incompleteGamma (p1, xx, g)) < 0.0) 
			{
			printf ("\nerr IncompleteGamma");
			return (-1.0);
			}
		p2 = p-t;
		t = p2*exp(xx*aa+g+p1-c*log(ch));   
		b = t/ch;  
		a = 0.5*t-b*c;
		s1 = (210.0+a*(140.0+a*(105.0+a*(84.0+a*(70.0+60.0*a))))) / 420.0;
		s2 = (420.0+a*(735.0+a*(966.0+a*(1141.0+1278.0*a))))/2520.0;
		s3 = (210.0+a*(462.0+a*(707.0+932.0*a)))/2520.0;
		s4 = (252.0+a*(672.0+1182.0*a)+c*(294.0+a*(889.0+1740.0*a)))/5040.0;
		s5 = (84.0+264.0*a+c*(175.0+606.0*a))/2520.0;
		s6 = (120.0+c*(346.0+127.0*c))/5040.0;
		ch += t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
		if (fabs(q/ch-1.0) > e) 
			goto l4;
		return (ch);

}

#pragma mark Gamma Distribution

/*!
 * This function generates a gamma-distributed random variable.
 *
 * \brief Gamma random variable.
 * \param a is the shape parameter of the gamma. 
 * \param b is the scale parameter of the gamma. 
 * \return Returns a gamma-distributed random variable.
 * \throws Does not throw an error.
 */
double MbRandom::gammaRv(double a, double b) {

	return (rndGamma(a) / b);

}

/*!
 * This function calculates the probability density 
 * for a gamma-distributed random variable.
 *
 * \brief Gamma probability density.
 * \param a is the shape parameter of the gamma. 
 * \param b is the scale parameter of the gamma. 
 * \param x is the gamma random variable. 
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
double MbRandom::gammaPdf(double a, double b, double x) {

	return (pow(b, a) / gamma(a)) * pow(x, a - 1.0) * exp(-x * b);
	
}

/*!
 * This function calculates the cumulative probability  
 * for a gamma-distributed random variable.
 *
 * \brief Gamma cumulative probability.
 * \param a is the shape parameter of the gamma. 
 * \param b is the scale parameter of the gamma. 
 * \param x is the gamma random variable. 
 * \return Returns the cumulative probability.
 * \throws Does not throw an error.
 */
double MbRandom::gammaCdf(double a, double b, double x) {

	return incompleteGamma(b*x, a, lnGamma(a));
	
}

/*!
 * This function calculates the natural log of the probability density 
 * for a gamma-distributed random variable.
 *
 * \brief Natural log of gamma probability density.
 * \param a is the shape parameter of the gamma. 
 * \param b is the scale parameter of the gamma. 
 * \param x is the gamma random variable. 
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
double MbRandom::lnGammaPdf(double a, double b, double x) {

	return a * log(b) - lnGamma(a) + (a - 1.0) * log(x) - x * b;
	
}

/*!
 * This function calculates the average or median values of the  
 * categories for a discretized gamma distribution. The gamma is
 * broken into K categories, with each category having equal
 * probability. The mean or meadian value for each category is
 * then calculated.
 *
 * \brief Calculates the mean/median values for a discretized gamma.
 * \param catRate a pointer to a vector of doubles containing the nCats mean/median values. 
 * \param a the shape parameter of the gamma distribution. 
 * \param b the scale parameter of the gamma distribution. 
 * \param nCats the number of discrete categories. 
 * \param When true, the median of the category is calculated. Otherwise
          the mean is used. 
 * \return Does not return a value.
 * \throws Does not throw an error.
 */
void MbRandom::discretizeGamma(MbVector<double> &catRate, double a, double b, int nCats, bool median) {

	double factor = a / b * nCats;

	if (median) 
		{
		/* the median value for each category is used to represent all of the values
		   in that category */
		double interval = 1.0 / (2.0 * nCats);
		for (int i=0; i<nCats; i++) 
			catRate[i] = chiSquareQuantile((i * 2.0 + 1.0) * interval, 2.0 * a) / (2.0 * b);
		double t = 0.0;
		for (int i=0; i<nCats; i++) 
			t += catRate[i];
		for (int i=0; i<nCats; i++)     
			catRate[i] *= factor / t;
		}
	else 
		{
		/* the mean value for each category is used to represent all of the values
		   in that category */
		/* calculate the points in the gamma distribution */
		for (int i=0; i<nCats-1; i++) 
			catRate[i] = chiSquareQuantile((i + 1.0) / nCats, 2.0 * a) / (2.0 * b);
		/* calculate the cumulative values */
		double lnGammaValue = lnGamma(a + 1.0);
		for (int i=0; i<nCats-1; i++) 
			catRate[i] = incompleteGamma(catRate[i] * b, a + 1.0, lnGammaValue);
		catRate[nCats-1] = 1.0;
		/* calculate the relative values and rescale */
		for (int i=nCats-1; i>0; i--)
			{
			catRate[i] -= catRate[i-1];
			catRate[i] *= factor;
			}
		catRate[0] *= factor;
		}
	
}

#pragma mark Log Normal Distribution

/*!
 * This function returns the quantile of a normal probability 
 * distribution.
 *
 * \brief Natural quantile.
 * \param mu is the mean parameter of the normal. 
 * \param sigma is the variance parameter of the normal. 
 * \param p is the probability up to the quantile. 
 * \return Returns the quantile.
 * \throws Does not throw an error.
 */
double MbRandom::logNormalQuantile(double mu, double sigma, double p) {

  return exp( normalQuantile(mu, sigma, p) );
	
}

/*!
 * This function calculates the cumulative probability 
 * for a normally-distributed random variable.
 *
 * \brief Normal cumulative probability.
 * \param mu is the mean parameter of the normal. 
 * \param sigma is the variance parameter of the normal. 
 * \param x is the normal random variable. 
 * \return Returns the cumulative probability.
 * \see Adams, A. G. 1969. Areas under the normal curve. Cojputer J. 12:197-198.
 * \throws Does not throw an error.
 */
double MbRandom::normalCdf(double mu, double sigma, double x) {

	double cdf;
	double q;
	double z = (x - mu) / sigma;

	/* |X| <= 1.28 */
	if ( fabs(z) <= 1.28 )
		{
		double a1 = 0.398942280444;
		double a2 = 0.399903438504;
		double a3 = 5.75885480458;
		double a4 = 29.8213557808;
		double a5 = 2.62433121679;
		double a6 = 48.6959930692;
		double a7 = 5.92885724438;
		double y = 0.5 * z * z;
		q = 0.5 - fabs(z) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5 + a6 / ( y + a7 ) ) ) );
		}
	else if ( fabs(z) <= 12.7 )
		{
		double b0 = 0.398942280385;
		double b1 = 3.8052E-08;
		double b2 = 1.00000615302;
		double b3 = 3.98064794E-04;
		double b4 = 1.98615381364;
		double b5 = 0.151679116635;
		double b6 = 5.29330324926;
		double b7 = 4.8385912808;
		double b8 = 15.1508972451;
		double b9 = 0.742380924027;
		double b10 = 30.789933034;
		double b11 = 3.99019417011;
		double y = 0.5 * z * z;
		q = exp(-y) * b0 / (fabs(z) - b1 + b2 / (fabs(z) + b3 + b4 / (fabs(z) - b5 + b6 / (fabs(z) + b7 - b8 / (fabs(z) + b9 + b10 / (fabs(z) + b11))))));
		}
	else
		{
		q = 0.0;
		}
	if ( z < 0.0 )
		{
		/* negative x */
		cdf = q;
		}
	else
		{
		/* positive x */
		cdf = 1.0 - q;
		}
	return cdf;

}

/*!
 * This function returns the quantile of a normal probability 
 * distribution.
 *
 * \brief Natural quantile.
 * \param mu is the mean parameter of the normal. 
 * \param sigma is the variance parameter of the normal. 
 * \param p is the probability up to the quantile. 
 * \return Returns the quantile.
 * \throws Does not throw an error.
 */
double MbRandom::normalQuantile(double mu, double sigma, double p) {
	
	double z = pointNormal(p);
	double x = z * sigma + mu;
	return x;
	
}

#pragma mark Uniform Distribution

/*!
 * This function generates a uniformly-distributed random variable on the interval [0,1).
 * It is a version of the Marsaglia Multi-Carry.
 *
 * Taken from:
 *   Mathlib : A C Library of Special Functions
 *   Copyright (C) 2000, 2003  The R Development Core Team
 *
 * This random generator has a period of 2^60, which ensures it has the maximum
 * period of 2^32 for unsigned ints (32 bit ints).
 *
 * \brief Uniform[0,1) random variable.
 * \return Returns a uniformly-distributed random variable on the interval [0,1).
 * \throws Does not throw an error.
 * \see http://stat.fsu.edu/~geo/diehard.html
 */
double MbRandom::uniformRv(void) {

	// Returns a pseudo-random number between 0 and 1.
	double u = 0.0; // TAH: if I1 == 395220642 && I2 == 1045444223, this returns 0, so I put this check in...
	while(u == 0.0){
		I1 = 36969 * (I1 & 0177777) + (I1 >> 16);
		I2 = 18000 * (I2 & 0177777) + (I2 >> 16);
		u = ((I1 << 16)^(I2 & 0177777)) * 2.328306437080797e-10; 	/*!< in [0,1) */
	}
	return u;
}

/*!
 * This function calculates the cumulative probability  
 * for a uniform(0,1) random variable.
 *
 * \brief Uniform(0,1) cumulative probability.
 * \param x is the uniform random variable. 
 * \return Returns the cumulative probability.
 * \throws Does not throw an error.
 */
double MbRandom::uniformCdf(double x) {

	if ( x < 0.0 )
		return 0.0;
	else if ( x > 1.0 )
		return 1.0;
	else
		return x;
	
}

/*!
 * This function generates a discrete and uniformly-distributed random 
 * variable on the interval (a,b).
 *
 * \brief Discrete uniform(a,b) random variable.
 * \param a is the lower bound on the uniform. 
 * \param b is the upper bound on the uniform. 
 * \return Returns a discrete and uniformly-distributed random variable 
 *         on the interval (a,b).
 * \throws Does not throw an error.
 */
int MbRandom::discreteUniformRv(int a, int b) {

	return (int)((b-a+1) * uniformRv()) + a;

}

#pragma mark Beta Distribution

/*!
 * This function generates a Beta-distributed random variable.
 *
 * \brief Beta random variable.
 * \param a parameter of the Beta. 
 * \param b parameter of the Beta. 
 * \return Returns the random variable.
 * \throws Does not throw an error.
 */
double MbRandom::betaRv(double a, double b) {

	double z0 = rndGamma( a );
	double z1 = rndGamma( b );
	double sum = z0 + z1;
	double x = z0 / sum;

	/*double mu = ( a - 1.0 ) / ( a + b - 2.0 );
	double stdev = 0.5 / sqrt( a + b - 2.0 );
	double x;
	for (;;)
		{
		double y = normalRv();
		x = mu + stdev * y;
		if ( x < 0.0 || 1.0 < x )
			continue;
		double u = uniformRv();
		double test = (a - 1.0) * log(x / (a - 1.0)) + (b - 1.0) * log((1.0 - x) / (b - 1.0)) + (a + b - 2.0) * log(a + b - 2.0) + 0.5 * y * y;
		if ( log (u) <= test )
			break;
		}*/
	return x;

}

/*!
 * This function returns the probability density for a 
 * Beta-distributed random variable.
 *
 * \brief Beta probability density.
 * \param a parameter of the Beta. 
 * \param b parameter of the Beta. 
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
double MbRandom::betaPdf(double a, double b, double x) {

	double pdf;
	if ( x < 0.0 || 1.0 < x )
		pdf = 0.0;
	else
		pdf = pow(x, (a - 1.0)) * pow((1.0 - x), (b - 1.0)) / beta(a, b);
	return pdf;

}

/*!
 * This function returns the natural log of the probability density 
 * for a Beta-distributed random variable.
 *
 * \brief Beta log probability density.
 * \param a parameter of the Beta. 
 * \param b parameter of the Beta. 
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
double MbRandom::lnBetaPdf(double a, double b, double x) {

	return ( (lnGamma(a + b) - lnGamma(a) - lnGamma(b)) + (a - 1.0) * log(x) + (b - 1.0) * log(1.0 - x) );

}

/*!
 * This function returns the cumulative probability for a 
 * Beta-distributed random variable.
 *
 * \brief Beta cumulative probability.
 * \param a parameter of the Beta. 
 * \param b parameter of the Beta. 
 * \return Returns the cumulative probability.
 * \throws Does not throw an error.
 */
double MbRandom::betaCdf(double a, double b, double x) {

	double cdf;
	if ( x <= 0.0 )
		cdf = 0.0;
	else if ( x <= 1.0 )
		cdf = incompleteBeta(a, b, x);
	else
		cdf = 1.0;
	return cdf;

}

/*!
 * This function returns the quantile for a 
 * Beta-distributed random variable.
 *
 * \brief Beta quantile.
 * \param a parameter of the Beta. 
 * \param b parameter of the Beta. 
 * \param p is the probability up to the quantile. 
 * \return Returns the quantile.
 * \throws Does not throw an error.
 */
double MbRandom::betaQuantile(double a, double b, double p) {

#	define MAXK 20

	double bcoeff;
	double error = 0.0001;
	double errapp = 0.01;
	int j;

	/* estimate the solution */
	double x = a / ( a + b );

	double xOld = 0.0;
	int loopCnt = 2;
	double d[MAXK * (MAXK-1)];
	while ( errapp <= fabs ( ( x - xOld ) / x ) && loopCnt != 0 )
		{
		xOld = x;
		loopCnt--;
		/* cdfX = PROB { BETA(A,B) <= X }
		   q = ( cdf - cdfX ) / pdfX */
		double cdfX = betaCdf(a, b, x);
		double pdfX = betaPdf(a, b, x);
		double q = (p - cdfX) / pdfX;
		/* D(N,K) = C(N,K) * Q**(N+K-1) / (N-1)! */
		double t = 1.0 - x;
		double s1 = q * ( b - 1.0 ) / t;
		double s2 = q * ( 1.0 - a ) / x;
		d[2-1+0*MAXK] = s1 + s2;
		double tail = d[2-1+0*MAXK] * q / 2.0;
		x = x + q + tail;

		int k = 3;
		while ( error < fabs ( tail / x ) && k <= MAXK )
			{
			/* find D(2,K-2) */
			s1 = q * ((double)(k) - 2.0) * s1 / t;
			s2 = q * (2.0 - (double)(k)) * s2 / x;
			d[2-1+(k-2)*MAXK] = s1 + s2;
			/* find D(3,K-3), D(4,K-4), D(5,K-5), ... , D(K-1,1) */
			for (int i=3; i<=k-1; i++)
				{
				double sum2 = d[2-1+0*MAXK] * d[i-2+(k-i)*MAXK];
				bcoeff = 1.0;
				for ( j = 1; j <= k - i; j++ )
					{
					bcoeff = ( bcoeff * ( double ) ( k - i - j + 1 ) ) / ( double ) ( j );
					sum2 = sum2 + bcoeff * d[2-1+j*MAXK] * d[i-2+(k-i-j)*MAXK];
					}
				d[i-1+(k-i)*MAXK] = sum2 + d[i-2+(k-i+1)*MAXK] / (double)(i - 1);
				}
			/* compute D(K,0) and use it to expand the series */
			d[k-1+0*MAXK] = d[2-1+0*MAXK] * d[k-2+0*MAXK] + d[k-2+1*MAXK] / (double)(k - 1);
			tail = d[k-1+0*MAXK] * q / (double)(k);
			x += tail;
			/* check for divergence */
			if ( x <= 0.0 || 1.0 <= x )
				{
				cout << "Error in betaQuantile: The series has diverged" << endl;
				x = -1.0;
				return x;
				}
			k++;
			}
		}
	return x;
#	undef MAXK

}

#pragma mark Dirichlet Distribution

/*!
 * This function generates a Dirichlet-distributed random variable.
 *
 * \brief Dirichlet random variable.
 * \param *a is a pointer to a vector of doubles containing the parameters of the Dirichlet. 
 * \param n is an integer with the number of Dirichlet prameters. 
 * \param *z is a pointer to a vector of doubles containing the Dirichlet random variable. 
 * \return Does not return a value (the random variable is initialized in the parameter z).
 * \throws Does not throw an error.
 */
void MbRandom::dirichletRv(const MbVector<double> &a, MbVector<double> &z) {

	int n = a.size();
	double sum = 0.0;
	for(int i=0; i<n; i++)
		{
		/* z[i] = rndGamma(a[i]) / 1.0; */
		z[i] = rndGamma(a[i]);
		sum += z[i];
		}
	for(int i=0; i<n; i++)
		z[i] /= sum;

}

/*!
 * This function calculates the probability density 
 * for a Dirichlet-distributed random variable.
 *
 * \brief Dirichlet probability density.
 * \param *a is a pointer to a vector of doubles containing the Dirichlet parameters. 
 * \param *z is a pointer to a vector of doubles containing the random variables. 
 * \param n is the number of Dirichlet parameters/random variables.
 * \return Returns the probability density.
 * \throws Throws an MbException::ERROR.
 */
double MbRandom::dirichletPdf(const MbVector<double> &a, const MbVector<double> &z) {
	
	int n = a.size();
	double zSum = 0.0;
	for (int i=0; i<n; i++)
		zSum += z[i];

	double tol = 0.0001;
	if ( tol < fabs( zSum - 1.0 ) )
		{
		cout << "Fatal error in dirichletPdf" << endl;
		exit(1);
		//ui->error("Fatal error in dirichletPdf");
		//throw(MbException(MbException::ERROR));
		}

	double aSum = 0.0;
	for (int i=0; i<n; i++)
		aSum += a[i];

	double aProd = 1.0;
	for (int i=0; i<n; i++)
		aProd *= gamma(a[i]);

	double pdf = gamma(aSum) / aProd;

	for (int i=0; i<n; i++)
		pdf = pdf * pow( z[i], a[i] - 1.0 );

	return pdf;

}

/*!
 * This function calculates the natural log of the probability density 
 * for a Dirichlet-distributed random variable.
 *
 * \brief Natural log of Dirichlet probability density.
 * \param *a is a pointer to a vector of doubles containing the Dirichlet parameters. 
 * \param *z is a pointer to a vector of doubles containing the random variables. 
 * \param n is the number of Dirichlet parameters/random variables.
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
double MbRandom::lnDirichletPdf(const MbVector<double> &a, const MbVector<double> &z) {
	int n = a.size(); //!< we assume that a and z have the same size
		
	double alpha0 = 0.0;
	for (int i=0; i<n; i++)
		alpha0 += a[i];
	double lnP = lnGamma(alpha0);
	for (int i=0; i<n; i++)
		lnP -= lnGamma(a[i]);
	for (int i=0; i<n; i++)
		lnP += (a[i] - 1.0) * log(z[i]);	
	return lnP;

}

#pragma mark Poisson Distribution

/*!
 * This function generates a Poisson-distributed random 
 * variable with parameter lambda.
 *
 * \brief Poisson(lambda) random variable.
 * \param lambda the rate parameter of the Poisson. 
 * \return This function returns a Poisson-distributed integer.
 * \throws Does not throw an error.
 */
int MbRandom::poissonRv(double lambda) {

	if (lambda < 17.0)
		{
		if (lambda < 1.0e-6)
			{
			if (lambda == 0.0) 
				return 0;
			if (lambda < 0.0)
				{
				/* there should be an error here */
				cout << "Parameter negative in poisson function" << endl;
				}

			/* For extremely small lambda we calculate the probabilities of x = 1
			   and x = 2 (ignoring higher x). The reason for using this 
			   method is to prevent numerical inaccuracies in other methods. */
			return poissonLow(lambda);
			}
		else 
			{
			/* use the inversion method */
			return poissonInver(lambda);
			}
		}
	else 
		{
		if (lambda > 2.0e9) 
			{
			/* there should be an error here */
			cout << "Parameter too big in poisson function" <<endl;
			}
		/* use the ratio-of-uniforms method */
		return poissonRatioUniforms(lambda);
		}
		
}

/*!
 * This function calculates the cumulative probability for a
 * Poisson distribution. 
 *
 * \brief Poisson cumulative probability.
 * \param lambda is the rate parameter of the Poisson. 
 * \param x is the value of the random variable. 
 * \return Returns the cumulative probability. 
 * \throws Does not throw an error.
 */
double MbRandom::poissonCdf(double lambda, int x) {

	if ( x < 0 )
		return 0.0;
	double next = exp(-lambda);
	double cdf = next;
	for (int i=1; i<=x; i++)
		{
		double last = next;
		next = last * lambda / (double)i;
		cdf += next;
		}
	return cdf;

}

/*!
 * This function calculates the beta function.
 *
 * \brief Beta function.
 * \param a is an argument. 
 * \param b is an argument. 
 * \return Returns the value of the beta function. 
 * \throws Does not throw an error.
 */
double MbRandom::beta(double a, double b) {

  return ( exp(lnGamma(a) + lnGamma(b) - lnGamma(a + b)) );
  
}

/*!
 * This function calculates the gamma function for real x.
 *
 * \brief Gamma function.
 * \param x is the argument. 
 * \return Returns the value of the gamma function. 
 * \throws Does not throw an error.
 */
double MbRandom::gamma(double x) {

	double c[7] = { -1.910444077728E-03, 
	                8.4171387781295E-04, 
	                -5.952379913043012E-04, 
	                7.93650793500350248E-04, 
	                -2.777777777777681622553E-03, 
	                8.333333333333333331554247E-02, 
	                5.7083835261E-03 };
	double fact;
	int i;
	int n;
	double p[8] = { -1.71618513886549492533811, 
	                2.47656508055759199108314E+01, 
	                -3.79804256470945635097577E+02, 
	                6.29331155312818442661052E+02, 
	                8.66966202790413211295064E+02, 
	                -3.14512729688483675254357E+04, 
	                -3.61444134186911729807069E+04, 
	                6.64561438202405440627855E+04 };
	bool parity;
	double q[8] = { -3.08402300119738975254353E+01, 
	                3.15350626979604161529144E+02,
	                -1.01515636749021914166146E+03,
	                -3.10777167157231109440444E+03, 
	                2.25381184209801510330112E+04, 
	                4.75584627752788110767815E+03, 
	                -1.34659959864969306392456E+05, 
	                -1.15132259675553483497211E+05 };
	double sqrtpi = 0.9189385332046727417803297;
	double sum2;
	double value;
	double xbig = 35.040;
	double xden;
	double xminin = 1.18E-38;
	double xnum;
	double y;
	double y1;
	double ysq;
	double z;

	parity = false;
	fact = 1.0;
	n = 0;
	y = x;

	if ( y <= 0.0 )
		{
		/* argument negative */
		y = -x;
		y1 = ( double ) ( ( int ) ( y ) );
		value = y - y1;

		if ( value != 0.0 )
			{
			if ( y1 != ( double ) ( ( int ) ( y1 * 0.5 ) ) * 2.0 )
				parity = true;
			fact = -PI / sin(PI * value);
			y = y + 1.0;
			}
		else
			{
			//value = d_huge ( );
			value = HUGE_VAL;
			return value;
			}
		}
	if ( y < mbEpsilon() )
		{
		/* argument < EPS */
		if ( xminin <= y )
			{
			value = 1.0 / y;
			}
		else
			{
			//value = d_huge ( );
			value = HUGE_VAL;
			return value;
			}
		}
	else if ( y < 12.0 )
		{
		y1 = y;
		/* 0.0 < argument < 1.0 */
		if ( y < 1.0 )
			{
			z = y;
			y = y + 1.0;
			}
		/* 1.0 < argument < 12.0, reduce argument if necessary */
		else
			{
			n = int ( y ) - 1;
			y = y - ( double ) ( n );
			z = y - 1.0;
			}
		/* evaluate approximation for 1.0 < argument < 2.0 */
		xnum = 0.0;
		xden = 1.0;
		for ( i = 0; i < 8; i++ )
			{
			xnum = ( xnum + p[i] ) * z;
			xden = xden * z + q[i];
			}

		value = xnum / xden + 1.0;
		/* adjust result for case  0.0 < argument < 1.0 */
		if ( y1 < y )
			{
			value = value / y1;
			}
		/* adjust result for case  2.0 < argument < 12.0 */
		else if ( y < y1 )
			{
			for ( i = 1; i <= n; i++ )
				{
				value = value * y;
				y = y + 1.0;
				}
			}
		}
	else
		{
		/* evaluate for 12 <= argument */
		if ( y <= xbig )
			{
			ysq = y * y;
			sum2 = c[6];
			for ( i = 0; i < 6; i++ )
				{
				sum2 = sum2 / ysq + c[i];
				}
			sum2 = sum2 / y - y + sqrtpi;
			sum2 = sum2 + ( y - 0.5 ) * log ( y );
			value = exp ( sum2 );
			}
		else
			{
			//value = d_huge ( );
			value = HUGE_VAL;
			return value;
			}

		}
	/* final adjustments and return */
	if ( parity )
		{
		value = -value;
		}
	if ( fact != 1.0 )
		{
		value = fact / value;
		}

	return value;

}

#pragma mark Helper Functions

/*!
 * This function returns the incomplete beta function, which is
 *
 * BI(a,b,x) = Integral(0 <= T <= X) T**(A-1) (1-T)**(B-1) dt /
 *             Integral(0 <= T <= 1) T**(A-1) (1-T)**(B-1) dt
 *
 * \brief Incomplete beta function.
 * \param a is a beta parameter. 
 * \param b is a beta parameter. 
 * \param x is the upper limit of integration. 
 * \return Returns the incomplete beta function. 
 * \throws Does not throw an error.
 * \see Majumder & Bhattacharjee. 1973. Algorithm AS63. Applied 
 *      Statistics, 22.
 */
double MbRandom::incompleteBeta(double a, double b, double x) {

	double tol = 1.0E-07;

	double value;
	if ( x <= 0.0 )
		{
		value = 0.0;
		return value;
		}
	else if ( 1.0 <= x )
		{
		value = 1.0;
		return value;
		}

	/* change tail if necessary and determine S */
	double psq = a + b;

	double xx, cx, pp, qq;
	bool indx;
	if ( a < (a + b) * x )
		{
		xx = 1.0 - x;
		cx = x;
		pp = b;
		qq = a;
		indx = true;
		}
	else
		{
		xx = x;
		cx = 1.0 - x;
		pp = a;
		qq = b;
		indx = false;
		}

	double term = 1.0;
	int i = 1;
	value = 1.0;
	int ns = (int)(qq + cx * (a + b));

	/* use Soper's reduction formulas */
	double rx = xx / cx;

	double temp = qq - (double)i;
	if ( ns == 0 )
		rx = xx;

	int it = 0;
	int it_max = 1000;
	for (;;)
		{
		it++;
		if ( it_max < it )
			{
			cout << "Error in incompleteBeta: Maximum number of iterations exceeded!" << endl;
			exit(1);
			//ui->error("Error in incompleteBeta: Maximum number of iterations exceeded!");
			//throw(MbException(MbException::ERROR));
			}
		term = term * temp * rx / ( pp + ( double ) ( i ) );
		value = value + term;
		temp = fabs(term);
		if ( temp <= tol && temp <= tol * value )
			break;
		i++;
		ns--;
		if ( 0 <= ns )
			{
			temp = qq - (double)i;
			if ( ns == 0 )
				rx = xx;
			}
		else
			{
			temp = psq;
			psq = psq + 1.0;
			}
		}
		
	/* finish calculation */
	value = value * exp(pp * log(xx) + (qq - 1.0) * log(cx)) / (beta(a, b) * pp);
	if ( indx )
		value = 1.0 - value;
	return value;
	
}

/*!
 * This function returns the incomplete gamma ratio I(x,alpha) where x is
 * the upper limit of the integration and alpha is the shape parameter.
 *
 * \brief Incomplete gamma function.
 * \param alpha is the shape parameter of the gamma. 
 * \param x is the upper limit of integration. 
 * \return Returns -1 if in error and the incomplete gamma ratio otherwise. 
 * \throws Does not throw an error.
 * \see Bhattacharjee, G. P. 1970. The incomplete gamma integral. Applied 
 *      Statistics, 19:285-287.
 */
double MbRandom::incompleteGamma (double x, double alpha, double LnGamma_alpha) {

	int 			i;
	double 		p = alpha, g = LnGamma_alpha,
					accurate = 1e-8, overflow = 1e30,
					factor, gin = 0.0, rn = 0.0, a = 0.0, b = 0.0, an = 0.0, 
					dif = 0.0, term = 0.0, pn[6];

	if (x == 0.0) 
		return (0.0);
	if (x < 0 || p <= 0) 
		return (-1.0);

	factor = exp(p*log(x)-x-g);   
	if (x>1 && x>=p) 
		goto l30;
	gin = 1.0;  
	term = 1.0;  
	rn = p;
	l20:
		rn++;
		term *= x/rn;   
		gin += term;
		if (term > accurate) 
			goto l20;
		gin *= factor/p;
		goto l50;
	l30:
		a = 1.0-p;   
		b = a+x+1.0;  
		term = 0.0;
		pn[0] = 1.0;  
		pn[1] = x;  
		pn[2] = x+1;  
		pn[3] = x*b;
		gin = pn[2]/pn[3];
	l32:
		a++;  
		b += 2.0;  
		term++;   
		an = a*term;
		for (i=0; i<2; i++) 
			pn[i+4] = b*pn[i+2]-an*pn[i];
		if (pn[5] == 0) 
			goto l35;
		rn = pn[4]/pn[5];   
		dif = fabs(gin-rn);
		if (dif>accurate) 
			goto l34;
		if (dif<=accurate*rn) 
			goto l42;
	l34:
		gin = rn;
	l35:
		for (i=0; i<4; i++) 
			pn[i] = pn[i+2];
		if (fabs(pn[4]) < overflow) 
			goto l32;
		for (i=0; i<4; i++) 
			pn[i] /= overflow;
		goto l32;
	l42:
		gin = 1.0-factor*gin;
	l50:
		return (gin);

}

/*!
 * This function sets the two seeds for the random number generator, using the current time.
 *
 *
 * \brief Initializes random number seeds.
 * \return This function does not return anything. 
 * \throws Does not throw an error.
 */
void MbRandom::setSeed(void) {

	seedType x = (seedType)( time( 0 ) );
	I1 = x & 0xFFFF;
	I2 = x >> 16;
	
}

/*!
 * This function sets the two seeds for the random number generator.
 * If only one starting value is given we set the two seeds by 
 * using the two least signficant bytes as one seed
 * and the two most significant bytes shifted to the right as the second seed.
 *
 *
 * \brief Initializes random number seeds.
 * \return This function does not return anything. 
 * \throws Does not throw an error.
 */
void MbRandom::setSeed(seedType seed1, seedType seed2) { 
	if (seed1 == 0)
		setSeed();
    else if(seed2 == 0) 
		{
        I1 = seed1&0xFFFF;
        I2 = seed1>>16;
		}
	else 
		{
        I1 = seed1;
        I2 = seed2;
		}
	
}

/*!
 * This function gets the two seeds from the random number generator.
 *
 * \brief Return the random number seeds.
 * \param i1 [in/out] the first seed
 * \param i2 [in/out] the second seed
 * \return This function does not return anything. 
 * \throws Does not throw an error.
 */
void MbRandom::getSeed(seedType& i1, seedType& i2) {

	i1 = I1;
	i2 = I2;
	
}
 
/*!
 * This function calculates the log of the gamma function, which is equal to:
 * Gamma(alp) = {integral from 0 to infinity} t^{alp-1} e^-t dt
 * The result is accurate to 10 decimal places. Stirling's formula is used
 * for the central polynomial part of the procedure.
 *
 * \brief Natural log of the gamma function.
 * \param alp is the parameter of the gamma function. 
 * \return Returns the log of the gamma function. 
 * \throws Does not throw an error.
 * \see Pike, M. C. and I. D. Hill. 1966. Algorithm 291: Logarithm of the gamma
 *      function. Communications of the Association for Computing Machinery, 9:684.
 */
double MbRandom::lnGamma(double a) {

	double x = a;
	double f = 0.0;
	double z;
	if (x < 7) 
		{
		f = 1.0;  
		z = x - 1.0;
		while (++z < 7.0)  
			f *= z;
		x = z;   
		f = -log(f);
		}
	z = 1.0 / (x*x);
	return  (f + (x-0.5)*log(x) - x + 0.918938533204673 + 
			(((-0.000595238095238*z+0.000793650793651)*z-0.002777777777778)*z +
			0.083333333333333)/x);  

}

/*!
 * This function returns the round off unit for floating point arithmetic.
 * The returned value is a number, r, which is a power of 2 with the property
 * that, to the precision of the computer's arithmetic, 1 < 1 + r, but
 * 1 = 1 + r / 2. This function comes from John Burkardt.
 *
 * \brief Round off unity for floating point arithmetic.
 * \return Returns round off unity for floating point arithmetic. 
 * \throws Does not throw an error.
 */
double MbRandom::mbEpsilon(void) {

	double r = 1.0;
	while ( 1.0 < (double)(1.0 + r)  )
		r = r / 2.0;
	return 2.0 * r;
}

/*!
 * This function generates a normal(0,1) random variable.
 *
 * \brief Standard normal random variable.
 * \return Returns a standard normal random variable. 
 * \throws Does not throw an error.
 */
double MbRandom::normalRv(void) {

	double fac, rsq, v1, v2;
	
	if ( availableNormalRv == false )
		{
		do
			{
			v1 = 2.0 * uniformRv() - 1.0;
			v2 = 2.0 * uniformRv() - 1.0;
			rsq = v1 * v1 + v2 * v2;
			} while ( rsq >= 1.0 || rsq == 0.0 );
		fac = sqrt(-2.0 * log(rsq)/rsq);
		extraNormalRv = v1 * fac;
		availableNormalRv = true;
		return v2 * fac;
		}
	else
		{
		availableNormalRv = false;
		return extraNormalRv;
		}

}

/*!
 * This function quantiles of a standard normal distribution.
 *
 * \brief Quantile of a standard normal distribution.
 * \param prob is the probability up to the quantile. 
 * \return Returns quantile value. 
 * \throws Does not throw an error.
 * \see Odeh, R. E. and J. O. Evans. 1974. The percentage points of the normal
 *      distribution. Applied Statistics, 22:96-97.
 * \see Wichura, M. J.  1988. Algorithm AS 241: The percentage points of the
 *      normal distribution. 37:477-484.
 * \see Beasley, JD & S. G. Springer. 1977. Algorithm AS 111: The percentage
 *      points of the normal distribution. 26:118-121.
 */
double MbRandom::pointNormal(double prob) {

	double a0 = -0.322232431088;
	double a1 = -1.0;
	double a2 = -0.342242088547;
	double a3 = -0.0204231210245;
 	double a4 = -0.453642210148e-4;
 	double b0 = 0.0993484626060;
 	double b1 = 0.588581570495;
 	double b2 = 0.531103462366; 
 	double b3 = 0.103537752850; 
 	double b4 = 0.0038560700634;
 	double p = prob;
	double p1 = ( p < 0.5 ? p : 1.0 - p);
	if (p1 < 1e-20) 
	   return (-9999.0);
	double y = sqrt( log(1.0/(p1*p1)) );   
	double z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0);
	return ( p < 0.5 ? -z : z );

}

/*!
 * This function generates a Poisson-distributed random variable for
 * small values of lambda. The method is a simple calculation of the 
 * probabilities of x = 1 and x = 2. Larger values are ignored.
 *
 * \brief Poisson random variables for small lambda.
 * \param lambda is the rate parameter of the Poisson. 
 * \return Returns a Poisson-distributed random variable. 
 * \throws Does not throw an error.
 */
int MbRandom::poissonLow(double lambda) {

	double d, r;
	d = sqrt(lambda);
	if (uniformRv() >= d) 
		return 0;
	r = uniformRv() * d;
	if (r > lambda * (1.0 - lambda)) 
		return 0;
	if (r > 0.5 * lambda * lambda * (1.0 - lambda)) 
		return 1;
	return 2;

}

/*!
 * This function generates a Poisson-distributed random variable using
 * inversion by the chop down method. 
 *
 * \brief Poisson random variables using the chop down method.
 * \param lambda is the rate parameter of the Poisson. 
 * \return Returns a Poisson-distributed random variable. 
 * \throws Does not throw an error.
 */
int MbRandom::poissonInver(double lambda) {

	const int bound = 130;
	static double p_L_last = -1.0;
	static double p_f0;
	double r;
	double f;
	int x;

	if (lambda != p_L_last) 
		{
		p_L_last = lambda;
		p_f0 = exp(-lambda);
		} 

	while (1) 
		{  
		r = uniformRv();  
		x = 0;  
		f = p_f0;
		do 
			{
			r -= f;
			if (r <= 0.0) 
				return x;
			x++;
			f *= lambda;
			r *= x;
			} while (x <= bound);
		}

}

/*!
 * This function generates a Poisson-distributed random variable using
 * the ratio-of-uniforms rejectin method. 
 *
 * \brief Poisson random variables using the ratio-of-uniforms method.
 * \param lambda is the rate parameter of the Poisson. 
 * \return Returns a Poisson-distributed random variable. 
 * \throws Does not throw an error.
 * \see Stadlober, E. 1990. The ratio of uniforms approach for generating
 *      discrete random variates. Journal of Computational and Applied 
 *      Mathematics 31:181-189.
 */
int MbRandom::poissonRatioUniforms(double lambda) {

	static double p_L_last = -1.0;  /* previous L */
	static double p_a;              /* hat center */
	static double p_h;              /* hat width */
	static double p_g;              /* ln(L) */
	static double p_q;              /* value at mode */
	static int p_bound;             /* upper bound */
	int mode;                       /* mode */
	double u;                       /* uniform random */
	double lf;                      /* ln(f(x)) */
	double x;                       /* real sample */
	int k;                          /* integer sample */

	if (p_L_last != lambda) 
		{
		p_L_last = lambda;
		p_a = lambda + 0.5;
		mode = (int)lambda;
		p_g  = log(lambda);
		p_q = mode * p_g - lnFactorial(mode);
		p_h = sqrt(2.943035529371538573 * (lambda + 0.5)) + 0.8989161620588987408;
		p_bound = (int)(p_a + 6.0 * p_h);
		}

	while(1) 
		{
		u = uniformRv();
		if (u == 0.0) 
			continue;
		x = p_a + p_h * (uniformRv() - 0.5) / u;
		if (x < 0 || x >= p_bound) 
			continue;
		k = (int)(x);
		lf = k * p_g - lnFactorial(k) - p_q;
		if (lf >= u * (4.0 - u) - 3.0) 
			break;
		if (u * (u - lf) > 1.0) 
			continue;
		if (2.0 * log(u) <= lf) 
			break;
		}
	return(k);

}

/*
 * This function is used when generating gamma-distributed random variables.
 *
 * \brief Subfunction for gamma random variables.
 * \param s is the shape parameter of the gamma. 
 * \return Returns a gamma-distributed random variable. 
 * \throws Does not throw an error.
 */
double MbRandom::rndGamma(double s) {

	double r=0.0;
	if (s <= 0.0)      
		cout << "Gamma parameter less than zero" << endl;
	else if (s < 1.0)  
		r = rndGamma1(s);
	else if (s > 1.0)  
		r = rndGamma2(s);
	else           
		r =- log(uniformRv());
	return (r);
   
}

/*!
 * This function is used when generating gamma-distributed random variables.
 *
 * \brief Subfunction for gamma random variables.
 * \param s is the shape parameter of the gamma. 
 * \return Returns a gamma-distributed random variable. 
 * \throws Does not throw an error.
 */
double MbRandom::rndGamma1(double s) {

	double			r, x=0.0, small=1e-37, w;
	static double   a, p, uf, ss=10.0, d;
	
	if (s != ss) 
		{
		a  = 1.0 - s;
		p  = a / (a + s * exp(-a));
		uf = p * pow(small / a, s);
		d  = a * log(a);
		ss = s;
		}
	for (;;) 
		{
		r = uniformRv();
		if (r > p)        
			x = a - log((1.0 - r) / (1.0 - p)), w = a * log(x) - d;
		else if (r>uf)  
			x = a * pow(r / p, 1.0 / s), w = x;
		else            
			return (0.0);
		r = uniformRv();
		if (1.0 - r <= w && r > 0.0)
		if (r*(w + 1.0) >= 1.0 || -log(r) <= w)  
			continue;
		break;
		}
		
	return (x);
   
}

/*!
 * This function is used when generating gamma-distributed random variables.
 *
 * \brief Subfunction for gamma random variables.
 * \param s is the shape parameter of the gamma. 
 * \return Returns a gamma-distributed random variable. 
 * \throws Does not throw an error.
 */
double MbRandom::rndGamma2(double s) {

	double			r, d, f, g, x;
	static double	b, h, ss=0.0;
	
	if (s != ss) 
		{
		b  = s - 1.0;
		h  = sqrt(3.0 * s - 0.75);
		ss = s;
		}
	for (;;) 
		{
		r = uniformRv();
		g = r - r * r;
		f = (r - 0.5) * h / sqrt(g);
		x = b + f;
		if (x <= 0.0) 
			continue;
		r = uniformRv();
		d = 64.0 * r * r * g * g * g;
		if (d * x < x - 2.0 * f * f || log(d) < 2.0 * (b * log(x / b) - f))  
			break;
		}
		
	return (x);
   
}

/* log factorial ln(n!) */
/*!
 * This function calculates the natural log of the factorial of n.
 * The first time this function is called, it constructs a table 
 * of the natural log of the factorial up to 1023 (inclusive). Every
 * call of this function after this for values of n < 1024, then,
 * is simply a look-up. For values of n >= 1024, the log of the factorial
 * is calculated using the Stirling approximation.
 *
 * \brief Natural log of the factorial.
 * \param n is the number for the factorial (n!). 
 * \return Returns the natural log of the factorial of n. 
 * \throws Does not throw an error.
 */
double MbRandom::lnFactorial(int n) {

	static const double    
	C0 =  0.918938533204672722,
	C1 =  1.0/12.0, 
	C3 = -1.0/360.0;

	if (n < 1024) 
		{
		if (n <= 1) 
			{
			return 0.0;
			}
		if ( initializedFacTable == false ) 
			{
			/* make table of ln(n!) */
			double sum = facTable[0] = 0.0;
			for (int i=1; i<1024; i++) 
				{
				sum += log((double)i);
				facTable[i] = sum;
				}
			initializedFacTable = true;
			}
		return facTable[n];
		}

	/* not found in table. use Stirling approximation */
	double  n1, r;
	n1 = n;  r  = 1.0 / n1;
	return (n1 + 0.5) * log(n1) - n1 + C0 + r*(C1 + r*r*C3);
	
}



/*!
 * \return Returns an index from [0, nCats) where the probability of returning
 *		category i is prob[i]
 * Note that all rounding error is "given" to the last category (nCats - 1).
 */
unsigned MbRandom::categoricalRv(const double * prob, const unsigned nCats) {
	assert(nCats > 0);
	assert(prob);
	double u = this->uniformRv();
	for (unsigned j = 0; j < nCats; j++)
		{
		u -= prob[j];
		if (u < 0.0)
			return j;
		}
	assert(u < 1.0e-6); // we should only get here because of rounding error
	return nCats - 1;
}








