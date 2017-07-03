/***************************************************************************
 *   Copyright (C) 2009 by BUI Quang Minh   *
 *   minh.bui@univie.ac.at   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef RATEGAMMA_H
#define RATEGAMMA_H

#include "rateheterogeneity.h"

const int GAMMA_CUT_MEDIAN = 1; // 2 discrete Gamma approximations (mean or median) of Yang 1994
const int GAMMA_CUT_MEAN   = 2;

class PhyloTree;
/**
Discrete gamma distributed site-rate model from Yang 1994

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class RateGamma: virtual public RateHeterogeneity
{

	friend class RateGammaInvar;

public:
	/**
		constructor
		@param ncat number of rate categories
		@param shape Gamma shape parameter
		@param tree associated phylogenetic tree
	*/
    RateGamma(int ncat, double shape, bool median, PhyloTree *tree);

	/**
		destructor
	*/
    virtual ~RateGamma();

    /**
        start structure for checkpointing
    */
    virtual void startCheckpoint();

    /**
        save object into the checkpoint
    */
    virtual void saveCheckpoint();

    /**
        restore object from the checkpoint
    */
    virtual void restoreCheckpoint();

	/**
		@return true if this is a Gamma model (default: false)
	*/	
    virtual int isGammaRate() { 
        if (cut_median) return GAMMA_CUT_MEDIAN; 
        return GAMMA_CUT_MEAN;
    }

	virtual double getGammaShape() { return gamma_shape; }

	virtual void setGammaShape(double gs);

	/**
	 * @return model name with parameters in form of e.g. GTR{a,b,c,d,e,f}
	 */
	virtual string getNameParams();

	/**
		@return TRUE to use median rate for discrete categories, FALSE to use mean rate instead
        OBSOLETE, see isGammaRate()
	*/
//	bool isCutMedian() { return cut_median; }

	/**
		@return the number of rate categories
	*/
	virtual int getNRate() { return ncategory; }

	/**
		get the number of rate categories for site-specific category model
		@return the number of rate categories
	*/
	virtual int getNDiscreteRate() { return ncategory; }

	/**
		@param category category ID from 0 to #category-1
		@return the rate of the specified category
	*/
	virtual double getRate(int category) { return rates[category]; }

	/**
		set the rate of a specified category.
		@param category category ID from 0 to #category-1
		@param value the rate of the specified category
	*/
	virtual void setRate(int category, double value) { rates[category] = value; }

	/**
		get the proportion of sites under a specified category.
		@param category category ID from 0 to #category-1
		@return the proportion of the specified category
	*/
	virtual double getProp(int category) { return 1.0/ncategory; }

	/**
	 * 	return pointer to the rate array
	 */
	virtual double* getRates() { return rates; }

	/** discrete Gamma according to Yang 1994 (JME 39:306-314) and using median cutting point
		It takes 'ncategory' and 'gamma_shape' variables as input. On output, it write to 'rates' variable.
	*/
	void computeRates();

	/** discrete Gamma according to Yang 1994 (JME 39:306-314) and using mean of the portion of gamma distribution
		It takes 'ncategory' and 'gamma_shape' variables as input. On output, it write to 'rates' variable.
	*/
	void computeRatesMean ();

	/**
		Compute site-specific rates. Override this for Gamma model
		@param pattern_rates (OUT) pattern rates. Resizing if necesary
        @return total number of categories
	*/
	virtual int computePatternRates(DoubleVector &pattern_rates, IntVector &pattern_cat);

	/**
	 * setup the bounds for joint optimization with BFGS
	 */
	virtual void setBounds(double *lower_bound, double *upper_bound, bool *bound_check);

	/**
		the target function which needs to be optimized
		@param x the input vector x
		@return the function value at x
	*/
	virtual double targetFunk(double x[]);

	/**
		optimize parameters. Default is to optimize gamma shape
		@return the best likelihood
	*/
	virtual double optimizeParameters(double gradient_epsilon);

    /**
     *  Same as above but add parameters to control gamma bounds
     */
	virtual double optimizeParameters(double gradient_epsilon, double min_gamma, double max_gamma);

	/**
		override function from Optimization class, used by the minimizeOneDimen() to optimize
		gamma shape parameter
	*/
	virtual double computeFunction(double shape);


	/**
		return the number of dimensions
	*/
	virtual int getNDim() { return !fix_gamma_shape; }

	/**
		write information
		@param out output stream
	*/
	virtual void writeInfo(ostream &out);

	/**
		write parameters, used with modeltest
		@param out output stream
	*/
	virtual void writeParameters(ostream &out);

	virtual bool isFixGammaShape() const {
		return fix_gamma_shape;
	}

	virtual void setFixGammaShape(bool fixGammaShape) {
		fix_gamma_shape = fixGammaShape;
	}

    /**
        set number of rate categories
        @param ncat #categories
    */
	virtual void setNCategory(int ncat);

protected:

	/**
		this function is served for the multi-dimension optimization. It should pack the model parameters
		into a vector that is index from 1 (NOTE: not from 0)
		@param variables (OUT) vector of variables, indexed from 1
	*/
	virtual void setVariables(double *variables);

	/**
		this function is served for the multi-dimension optimization. It should assign the model parameters
		from a vector of variables that is index from 1 (NOTE: not from 0)
		@param variables vector of variables, indexed from 1
		@return TRUE if parameters are changed, FALSE otherwise (2015-10-20)
	*/
	virtual bool getVariables(double *variables);

	/**
		number of rate categories
	*/
	int ncategory;

	/**
		rates, containing ncategory elements
	*/
	double *rates;


	/**
		the gamma shape parameter 'alpha'
	*/
	double gamma_shape;

	/**
		TRUE to fix the gamma shape parameter
	*/
	bool fix_gamma_shape;

	/**
		TRUE to use median rate for discrete categories, FALSE to use mean rate instead
	*/
	bool cut_median;

public:

	//Normally, beta is assigned equal to alpha
	//double cmpPerPointGamma (const double prob, const double shape);

	/***********************************************************
	NUMERICAL SUBROUTINES
	THE FOLLOWING CODE COMES FROM tools.c in Yang's PAML package
	***********************************************************/
	/** returns ln(gamma(alpha)) for alpha>0, accurate to 10 decimal places.
	   Stirling's formula is used for the central polynomial part of the procedure.
	   Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
	   Communications of the Association for Computing Machinery, 9:684

	*/
	static double cmpLnGamma (double alpha);

	/** returns the incomplete gamma ratio I(x,alpha) where x is the upper
	   limit of the integration and alpha is the shape parameter.
	   returns (-1) if in error
	   (1) series expansion     if (alpha>x || x<=1)
	   (2) continued fraction   otherwise
	   RATNEST FORTRAN by
	   Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
	   19: 285-287 (AS32)
	*/
	static double cmpIncompleteGamma (double x, double alpha, double ln_gamma_alpha);

	/** functions concerning the CDF and percentage points of the gamma and
	   Chi2 distribution
	   returns z so that Prob{x<z}=prob where x ~ N(0,1) and (1e-12)<prob<1-(1e-12)
	   returns (-9999) if in error
	   Odeh RE & Evans JO (1974) The percentage points of the normal distribution.
	   Applied Statistics 22: 96-97 (AS70)

	   Newer methods:
	     Wichura MJ (1988) Algorithm AS 241: the percentage points of the
	       normal distribution.  37: 477-484.
	     Beasley JD & Springer SG  (1977).  Algorithm AS 111: the percentage
	       points of the normal distribution.  26: 118-121.
	*/
	static double cmpPointNormal (double prob);


	/** returns z so that Prob{x<z}=prob where x is Chi2 distributed with df=v
	   returns -1 if in error.   0.000002<prob<0.999998
	   RATNEST FORTRAN by
	       Best DJ & Roberts DE (1975) The percentage points of the
	       Chi2 distribution.  Applied Statistics 24: 385-388.  (AS91)
	   Converted into C by Ziheng Yang, Oct. 1993.
	*/
	static double cmpPointChi2 (double prob, double v);

	/* THE END OF THE CODES COMMING FROM tools.c in Yang's PAML package */

};




#endif
