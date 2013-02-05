/*! 
 * \file
 * In this file we declare a class for random variables. It contains one
 * class, MbRandom, which when instantiated gives a stream of uniform
 * random variables. It also can transform the uniform random variable
 * to generate random variables with different distributions (such as
 * exponential, gamma, normal, etc.).
 *
 * \brief This file contains a definition for a class for random variables.
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
 * $Id: MbRandom.h,v 1.3 2006/09/11 17:29:51 paulvdm Exp $
 */

#ifndef MB_RANDOM_H
#define MB_RANDOM_H

#include <cmath>

#ifndef PI
#	define PI 3.141592653589793
#endif

/*! 
 * The Marsaglia Multi-Carry random number generator used in MrBayes 4 assumes a 32-bit data type 
 * for seeds. An unsigned int is 32 bits for both 32-bit and 64-bits processors. 
*/
typedef unsigned int seedType;

/*! 
 * MbRandom is a class that works with random variables. On creating an instance
 * of this class, a seed for a uniform random number is initialized. One can then
 * obtain random variables with different probability distributions, by transforming
 * a uniform(0,1) random variable in different ways. The class can either be
 * instantiated with a long integer, in which case the seed is set to the long
 * integers value (i.e., the class takes on a user-supplied seed) or the class
 * can be instantiated with no value passed in, in which case the seed is initialized
 * using the system time. The class also has functions for calculating the probability
 * or probability density for a random variable or for calculating the quantiles
 * of a probability distribution.
 *
 * \brief MbRandom is a class for generating random variables. 
*/
template <class T> class MbVector;
class MbRandom {

	public:
		                     MbRandom(void);                                                               /*!< constructor: initializes the seed with current time */
		                     MbRandom(seedType x);                                                         /*!< constructor: initializes the seed with supplied value */
					  void   getSeed(seedType &seed1, seedType &seed2);                                    /*!< retreives the seeds */
					  void   setSeed(void);                                                                /*!< initializes the seeds using the current time */
		              void   setSeed(seedType seed1, seedType seed2);                                      /*!< initializes the seeds */
					double   chiSquareRv(double v);                                       /* chi square */ /*!< Chi-square random variable */
                    double   chiSquarePdf(double v, double x);                                             /*!< the chi-square probability density */
                    double   lnChiSquarePdf(double v, double x);                                           /*!< natural log of the chi-square probability density */
                    double   chiSquareCdf(double v, double x);                                             /*!< the chi-square cumulative probability */
		            double   chiSquareQuantile(double prob, double v);                                     /*!< quantile of a chi square distribution */
		     inline double   exponentialRv(double lambda);                               /* exponential */ /*!< exponential random variable */
		     inline double   exponentialPdf(double lambda, double x);                                      /*!< Exponential probability density */
		     inline double   lnExponentialPdf(double lambda, double x);                                    /*!< natural log of Exponential probability density */
		     inline double   exponentialCdf(double lambda, double x);                                      /*!< Exponential cumulative probability */
		     inline double   exponentialQuantile(double lambda, double p);                                 /*!< quantile of an exponential distribution */
		            double   gammaRv(double a, double b);                                      /* gamma */ /*!< gamma random variable */
		            double   gammaPdf(double a, double b, double x);                                       /*!< Gamma probability density */
		            double   lnGammaPdf(double a, double b, double x);                                     /*!< natural log of Gamma probability density */
		            double   gammaCdf(double a, double b, double x);                                       /*!< Gamma cumulative probability */
		     inline double   gammaQuantile(double a, double b, double p);                                  /*!< quantile of gamma distribution */
 	 	     inline double   logNormalRv(double mu, double sigma);                        /* log normal */ /*!< log normal random variable */
		     inline double   logNormalPdf(double mu, double sigma, double x);                              /*!< log normal probability density */
		     inline double   lnLogNormalPdf(double mu, double sigma, double x);                            /*!< natural log of log normal probability density */
		     inline double   logNormalCdf(double mu, double sigma, double x);                              /*!< log normal cumulative probability */
		            double   logNormalQuantile(double mu, double sigma, double p);                         /*!< quantile of log normal distribution */
		     inline double   normalRv(double mu, double sigma);                               /* normal */ /*!< normal(mu,sigma) random variable */
		     inline double   normalPdf(double mu, double sigma, double x);                                 /*!< Normal probability density */
		     inline double   lnNormalPdf(double mu, double sigma, double x);                               /*!< natural log of Normal probability density */
		            double   normalCdf(double mu, double sigma, double x);                                 /*!< Normal cumulative probability */
		            double   normalQuantile(double mu, double sigma, double p);                            /*!< quantile of normal distribution */
		            double   uniformRv(void);                                           /* uniform(0,1) */ /*!< uniform(0,1) random variable */
		     inline double   uniformPdf(void);                                                             /*!< Uniform(0,1) probability density */
			 inline double   lnUniformPdf(void);                                                           /*!< natural log of Uniform(0,1) probability density */
		            double   uniformCdf(double x);                                                         /*!< Uniform(0,1) cumulative probability */
		     inline double   uniformQuantile(double p);                                                    /*!< quantile of a uniform(0,1) distribuiton */
		     inline double   uniformRv(double a, double b);                             /* uniform(a,b) */ /*!< uniform(a,b) random variable */
		     inline double   uniformPdf(double a, double b);                                               /*!< Uniform(a,b) probability density */
		     inline double   lnUniformPdf(double a, double b);                                             /*!< natural log of Uniform(a,b) probability density */
		            double   uniformCdf(double a, double b, double x);                                     /*!< Uniform(a,b) cumulative probability */
		     inline double   uniformQuantile(double a, double b, double p);                                /*!< quantile of a uniform(a,b) distribuiton */
		            double   betaRv(double a, double b);                                        /* beta */ /*!< Beta random variable */
		            double   betaPdf(double a, double b, double x);                                        /*!< Beta probability density */
		            double   lnBetaPdf(double a, double b, double x);                                      /*!< natural log of the Beta probability density */
		            double   betaCdf(double a, double b, double x);                                        /*!< Beta cumulative probability */
		            double   betaQuantile(double a, double b, double p);                                   /*!< quantile of the Beta distribution */
					  void   dirichletRv(const MbVector<double> &a, MbVector<double> &z);  /* dirichlet */ /*!< Dirichlet random variable */
		            double   dirichletPdf(const MbVector<double> &a, const MbVector<double> &z);           /*!< Dirichlet probability density */
		            double   lnDirichletPdf(const MbVector<double> &a, const MbVector<double> &z);         /*!< natural log of Dirichlet probability density */
		               int   discreteUniformRv(int a, int b);                       /* discrete uniform */ /*!< discrete uniform random variable */
		     inline double   discreteUniformProb(int a, int b);                                            /*!< discrete uniform probability */
		     inline double   lnDiscreteUniformProb(int a, int b);                                          /*!< natural log of discrete uniform probability */
		               int   poissonRv(double lambda);                                       /* poisson */ /*!< Poisson random variable */
		     inline double   poissonProb(double lambda, int x);                                            /*!< Poisson probability */
		     inline double   lnPoissonProb(double lambda, int x);                                          /*!< natural log of Poisson probability */
		            double   poissonCdf(double lambda, int x);                                             /*!< Poisson cumulative probability */
		            double   poissonQuantile(double lambda, double p);                                     /*!< quantile of a Poisson(lambda) distribution */
		              void   discretizeGamma(MbVector<double> &catRate, double a, double b, int nCats, bool median); /*!< calculates the average/median values for a discretized gamma distribution */
		            double   lnGamma(double a);                                                            /*!< log of the gamma function */
				  unsigned   categoricalRv(const double * prob, const unsigned nCats);
	private:
	                        /* private functions */
	                double   beta(double a, double b);                                                     /*!< calculates the beta function */
	                double   gamma(double x);                                                              /*!< calculates the gamma function */
	                double   incompleteBeta(double a, double b, double x);                                 /*!< calculates the incomplete beta function */
		            double   incompleteGamma (double x, double alpha, double LnGamma_alpha);               /*!< calculates the incomplete gamma ratio */
		            double   lnFactorial(int n);                                                           /*!< log of factorial [ln(n!)] */
		            double   mbEpsilon(void);                                                              /*!< round off unit for floating arithmetic */
		            double   normalRv(void);                                                               /*!< standard normal(0,1) random variable */
		            double   pointNormal(double prob);                                                     /*!< quantile of standard normal distribution */
		               int   poissonLow(double lambda);                                                    /*!< function used when calculating Poisson random variables */
		               int   poissonInver(double lambda);                                                  /*!< function used when calculating Poisson random variables */
		               int   poissonRatioUniforms(double lambda);                                          /*!< function used when calculating Poisson random variables */
		            double   rndGamma(double s);                                                           /*!< function used when calculating gamma random variable */
		            double   rndGamma1(double s);                                                          /*!< function used when calculating gamma random variable */
		            double   rndGamma2(double s);                                                          /*!< function used when calculating gamma random variable */
				   
		                    /* private data */
			      seedType   I1,I2;                                                                         /*!< seed values for the random number generator */
		              bool   initializedFacTable;                                                           /*!< a boolean which is false if the log factorial table has not been initialized */
		            double   facTable[1024];                                                                /*!< a table containing the log of the factorial up to 1024 */
		              bool   availableNormalRv;                                                             /*!< a boolean which is true if there is a normal random variable available */
		            double   extraNormalRv;                                                                 /*!< a normally-distributed random variable which */
		
};


/*!
 * This function generates an exponentially-distributed random variable.
 *
 * \brief Exponential random variable.
 * \param lambda is the rate parameter of the exponential. 
 * \return Returns an exponentially-distributed random variable.
 * \throws Does not throw an error.
 */
inline double MbRandom::exponentialRv(double lambda) {

	return -(1.0 / lambda) * std::log( uniformRv() );

}

/*!
 * This function calculates the probability density 
 * for an exponentially-distributed random variable.
 *
 * \brief Exponential probability density.
 * \param lambda is the rate parameter of the exponential. 
 * \param x is the exponential random variable. 
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
inline double MbRandom::exponentialPdf(double lambda, double x) {

	return lambda * exp(-lambda * x);
	
}


/*!
 * This function calculates the cumulative probability  
 * for an exponentially-distributed random variable.
 *
 * \brief Exponential cumulative probability.
 * \param lambda is the rate parameter of the exponential. 
 * \param x is the exponential random variable. 
 * \return Returns the cumulative probability.
 * \throws Does not throw an error.
 */
inline double MbRandom::exponentialCdf(double lambda, double x) {

	return 1.0 - exp(-lambda * x);
	
}


/*!
 * This function calculates the natural log of the probability density 
 * for an exponentially-distributed random variable.
 *
 * \brief Natural log of exponential probability density.
 * \param lambda is the rate parameter of the exponential. 
 * \param x is the exponential random variable. 
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
inline double MbRandom::lnExponentialPdf(double lambda, double x) {

	return (std::log(lambda) - lambda * x);
	
}

/*!
 * This function returns the quantile of a exponential probability 
 * distribution.
 *
 * \brief Exponential quantile.
 * \param lambda is the rate parameter. 
 * \param p is the probability up to the quantile. 
 * \return Returns the quantile.
 * \throws Does not throw an error.
 */
inline double MbRandom::exponentialQuantile(double lambda, double p) {

	return -(1.0 / lambda) * std::log(1.0 - p);
	
}


/*!
 * This function returns the quantile of a gamma probability 
 * distribution.
 *
 * \brief Gamma quantile.
 * \param a is the shape parameter. 
 * \param b is the scale parameter. 
 * \param p is the probability up to the quantile. 
 * \return Returns the quantile.
 * \throws Does not throw an error.
 */
inline double MbRandom::gammaQuantile(double a, double b, double p) {

	return chiSquareQuantile(p, 2.0 * a) / (2.0*b);

}

/*!
 * This function generates a log normally distributed random variable.
 *
 * \brief Log normal random variable.
 * \param mu is the mean parameter of the log normal. 
 * \param sigma is the variance parameter of the log normal. 
 * \return Returns a log normally distributed random variable.
 * \throws Does not throw an error.
 */
inline double MbRandom::logNormalRv(double mu, double sigma) {

	return exp( normalRv(mu, sigma) );
	
}


/*!
 * This function calculates the probability density 
 * for a log normally distributed random variable.
 *
 * \brief Log normal probability density.
 * \param mu is the mean parameter of the log normal. 
 * \param sigma is the variance parameter of the log normal. 
 * \param x is the log normal random variable. 
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
inline double MbRandom::logNormalPdf(double mu, double sigma, double x) {

	double pdf;
	if ( x <= 0.0 )
		{
		pdf = 0.0;
		}
	else
		{
		double y = ( std::log(x) - mu ) / sigma;
		pdf = exp( -0.5 * y * y ) / ( sigma * x * sqrt(2.0 * PI) );
		}
	return pdf;
	
}

/*!
 * This function calculates the cumulative probability 
 * for a log normally distributed random variable.
 *
 * \brief Log normal cumulative probability.
 * \param mu is the mean parameter of the log normal. 
 * \param sigma is the variance parameter of the log normal. 
 * \param x is the log normal random variable. 
 * \return Returns the cumulative probability.
 * \throws Does not throw an error.
 */
inline double MbRandom::logNormalCdf(double mu, double sigma, double x) {

	double cdf;
	if ( x <= 0.0 )
		{
		cdf = 0.0;
		}
	else
		{
		double logX = std::log(x);
		cdf = normalCdf( mu, sigma, logX );
		}
	return cdf;
	
}



/*!
 * This function calculates the natural log of the probability density 
 * for a log normally distributed random variable.
 *
 * \brief Natural log of the log normal probability density.
 * \param mu is the mean parameter of the log normal. 
 * \param sigma is the variance parameter of the log normal. 
 * \param x is the log normal random variable. 
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
inline double MbRandom::lnLogNormalPdf(double mu, double sigma, double x) {

	return ( - 0.5 * ( (std::log(x) - mu) / sigma ) * ( (std::log(x) - mu) / sigma ) ) - std::log( sigma * x * sqrt( 2.0 * PI ) );
	
}

/*!
 * This function generates a normally-distributed random variable.
 *
 * \brief Normal random variable.
 * \param mu is the mean of the normal. 
 * \param sigma is the variance of the normal. 
 * \return Returns a normally-distributed random variable.
 * \throws Does not throw an error.
 */
inline double MbRandom::normalRv(double mu, double sigma) {

	return ( mu + sigma * normalRv() );
	
}

/*!
 * This function calculates the probability density 
 * for a normally-distributed random variable.
 *
 * \brief Normal probability density.
 * \param mu is the mean parameter of the normal. 
 * \param sigma is the variance parameter of the normal. 
 * \param x is the normal random variable. 
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
inline double MbRandom::normalPdf(double mu, double sigma, double x) {

	double y = ( x - mu ) / sigma;
	return exp( -0.5 * y * y )  / ( sigma * sqrt ( 2.0 * PI ) );
	
}

/*!
 * This function calculates the natural log of the probability density 
 * for a normally-distributed random variable.
 *
 * \brief Natural log of normal probability density.
 * \param mu is the mean parameter of the normal. 
 * \param sigma is the variance parameter of the normal. 
 * \param x is the normal random variable. 
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
inline double MbRandom::lnNormalPdf(double mu, double sigma, double x) {

	return -0.5 * std::log(2.0 * PI * sigma) - 0.5 * (x - mu) * (x - mu) / (sigma * sigma);
	
}

/*!
 * This function calculates the probability density 
 * for a uniform(0,1) random variable.
 *
 * \brief Uniform(0,1) probability density.
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
inline double MbRandom::uniformPdf(void) {

	return 1.0;
	
}

/*!
 * This function calculates the natural log of the probability density 
 * for a uniform(0,1) random variable.
 *
 * \brief Natural log of uniform(0,1) probability density.
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
inline double MbRandom::lnUniformPdf(void) {

	return 0.0;
	
}

/*!
 * This function returns the quantile of a uniform(0,1) probability 
 * distribution.
 *
 * \brief Uniform(0,1) quantile.
 * \param p is the probability up to the quantile. 
 * \return Returns the quantile.
 * \throws Does not throw an error.
 */
inline double MbRandom::uniformQuantile(double p) {

	return p;

}

/*!
 * This function generates a uniformly-distributed random variable on the interval (a,b).
 *
 * \brief Uniform(a,b) random variable.
 * \param a is the lower bound on the uniform. 
 * \param b is the upper bound on the uniform. 
 * \return Returns a uniformly-distributed random variable on the interval (a,b).
 * \throws Does not throw an error.
 */
inline double MbRandom::uniformRv(double a, double b) {

	return ( a + uniformRv() * (b - a) );
	
}



/*!
 * This function calculates the probability density 
 * for a uniform(a,b) random variable.
 *
 * \brief Uniform(a,b) probability density.
 * \param a is the lower bound on the uniform. 
 * \param b is the upper bound on the uniform. 
 * \return Returns the probability density.
 * \throws Does not throw an error.
 */
inline double MbRandom::uniformPdf(double a, double b) {

	return 1.0 / (b-a);
	
}



/*!
 * This function calculates the cumulative probability 
 * for a uniform(a,b) random variable.
 *
 * \brief Uniform(a,b) cumulative probability.
 * \param a is the lower bound on the uniform. 
 * \param b is the upper bound on the uniform. 
 * \param x is the uniform random variable. 
 * \return Returns the cumulative probability.
 * \throws Does not throw an error.
 */
inline double MbRandom::uniformCdf(double a, double b, double x) {

	if ( x < a )
		return 0.0;
	else if ( x > b )
		return 1.0;
	else
		return (x-a) / (b-a);
	
}



/*!
 * This function calculates the natural log of the probability density 
 * for a uniform(a,b) random variable.
 *
 * \brief Natural log of uniform(a,b) probability density.
 * \param a is the lower bound on the uniform. 
 * \param b is the upper bound on the uniform. 
 * \return Returns the natural log of the probability density.
 * \throws Does not throw an error.
 */
inline double MbRandom::lnUniformPdf(double a, double b) {

	return ( -std::log(b-a) );
	
}



/*!
 * This function returns the quantile of a uniform(a,b) probability 
 * distribution.
 *
 * \brief Uniform(a,b) quantile.
 * \param a is the lower bound on the uniform. 
 * \param b is the upper bound on the uniform. 
 * \param p is the probability up to the quantile. 
 * \return Returns the quantile.
 * \throws Does not throw an error.
 */
inline double MbRandom::uniformQuantile(double a, double b, double p) {

	return a + (b - a) * p;

}

/*!
 * This function calculates the natural log of the probability for a
 * discrete uniform distribution. 
 *
 * \brief Natural log of discrete uniform probability.
 * \param a is the lower bound on the uniform. 
 * \param b is the upper bound on the uniform. 
 * \return Returns the natural log of the probability. 
 * \throws Does not throw an error.
 */
inline double MbRandom::discreteUniformProb(int a, int b) {

	return 1.0 / (b - a + 1);

}



/*!
 * This function calculates the natural log of the probability for a
 * discrete uniform distribution. 
 *
 * \brief Natural log of discrete uniform probability.
 * \param a is the lower bound on the uniform. 
 * \param b is the upper bound on the uniform. 
 * \return Returns the natural log of the probability. 
 * \throws Does not throw an error.
 */
inline double MbRandom::lnDiscreteUniformProb(int a, int b) {

	return std::log( 1.0 / (b - a + 1) );

}

/*!
 * This function calculates the natural log of the probability for a
 * Poisson distribution. 
 *
 * \brief Natural log of Poisson probability.
 * \param lambda is the rate parameter of the Poisson. 
 * \param x is the value of the random variable. 
 * \return Returns the natural log of the probability. 
 * \throws Does not throw an error.
 */
inline double MbRandom::poissonProb(double lambda, int x) {

	return exp(x * std::log(lambda) - lambda - lnFactorial(x));
	
}

/*!
 * This function calculates the natural log of the probability for a
 * Poisson distribution. 
 *
 * \brief Natural log of Poisson probability.
 * \param lambda is the rate parameter of the Poisson. 
 * \param x is the value of the random variable. 
 * \return Returns the natural log of the probability. 
 * \throws Does not throw an error.
 */
inline double MbRandom::lnPoissonProb(double lambda, int x) {

	return ( x * std::log(lambda) - lambda - lnFactorial(x) );
	
}



/*!
 * This function returns the quantile of a Poisson probability 
 * distribution.
 *
 * \brief Poisson(lambda) quantile.
 * \param lambda is the rate parameter of the Poisson. 
 * \param p is the probability up to the quantile. 
 * \return Returns the quantile.
 * \throws Does not throw an error.
 */
inline double MbRandom::poissonQuantile(double lambda, double p) {

	/* Starting with x = 0, find the first value for which
	   CDF(X-1) <= CDF <= CDF(X). */
	double sum = 0.0;
	int xmax = 100;
	for (int i=0; i<=xmax; i++)
		{
		double sumOld = sum;
		double newVal;
		if ( i == 0 )
			{
			newVal = exp(-lambda);
			sum = newVal;
			}
		else
			{
			double last = newVal;
			newVal = last * lambda / ( double ) ( i );
			sum += newVal;
			}
		if ( sumOld <= p && p <= sum )
			return i;
		}
	//cout << "Poisson quantile warning" << endl;
	return xmax;

}

#endif
