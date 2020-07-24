/***************************************************************************
 *   Copyright (C) 2009-2015 by                                            *
 *   BUI Quang Minh <minh.bui@univie.ac.at>                                *
 *   Lam-Tung Nguyen <nltung@gmail.com>                                    *
 *                                                                         *
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
#include "tree/phylotree.h"
#include "rategamma.h"
#include <cmath>



RateGamma::RateGamma(int ncat, double shape, bool median, PhyloTree *tree) : RateHeterogeneity()
{
	ncategory = ncat;
	phylo_tree = tree;
	cut_median = median;
	//gamma_shape = MAX_GAMMA_SHAPE-1.0;
	gamma_shape = max(tree->params->min_gamma_shape, fabs(shape));
	fix_gamma_shape = false;
	rates = NULL;
	if (shape > 0.0) {
		// true unless -optfromgiven cmd line option
		fix_gamma_shape = !(Params::getInstance().optimize_from_given_params);
	} else if (shape == 0.0) {
		gamma_shape = max(tree->params->min_gamma_shape*5.0, random_double());
		cout << "Randomize initial gamma shape (alpha): " << gamma_shape << endl;
	}
	setNCategory(ncat);
}

void RateGamma::startCheckpoint() {
    checkpoint->startStruct("RateGamma");
}

void RateGamma::saveCheckpoint() {
    startCheckpoint();
    CKP_SAVE(gamma_shape);
//    CKP_SAVE(fix_gamma_shape);
//    CKP_SAVE(cut_median);
//    CKP_SAVE(ncategory);
    endCheckpoint();
    RateHeterogeneity::saveCheckpoint();
}

void RateGamma::restoreCheckpoint() {
    RateHeterogeneity::restoreCheckpoint();
    startCheckpoint();
    CKP_RESTORE(gamma_shape);
//    CKP_RESTORE(fix_gamma_shape);
//    CKP_RESTORE(cut_median);
//    CKP_RESTORE(ncategory);
    endCheckpoint();
    // necessary compute rates after restoring gamma_shape
	computeRates();
}

void RateGamma::setNCategory(int ncat) {
	ncategory = ncat;
	delete [] rates;
	rates = new double[ncategory];
 	for (int cat = 0; cat < ncategory; cat++) {
		rates[cat] = 1.0;
	}
	name = "+G" + convertIntToString(ncategory);
	full_name = "Gamma with " + convertIntToString(ncategory) + " categories";
	computeRates();
}


string RateGamma::getNameParams() {
	ostringstream str;
	str << "+G" << ncategory << '{' << gamma_shape << '}';
	return str.str();
}

RateGamma::~RateGamma()
{
	delete [] rates;
	rates = nullptr;
}

void RateGamma::computeRates() {
	int cat; /* category id */
	double sum_rates = 0.0;
	if (ncategory == 1) {
		rates[0] = 1.0;
		return;
	}
    double curScale = 0.0;
	for (cat = 0; cat < ncategory; cat++) {
		curScale += rates[cat];
	}
	if (!cut_median) {
		computeRatesMean();
	} else {
		for (cat = 0; cat < ncategory; cat ++) {
			double prob = ( 2.0 * cat + 1 ) / (2.0 * ncategory);
			double perPoint_ = cmpPointChi2 (prob, 2.0 * gamma_shape) / (2.0 * gamma_shape);
			perPoint_ = perPoint_ < 0.0 ? -perPoint_ : perPoint_;
			rates[ cat ] = perPoint_;
		}

		//rescale in order to make mean equal to 1.0


		for (cat = 0; cat < ncategory; cat ++)
			sum_rates += rates[ cat];

		for (cat = 0; cat < ncategory; cat ++)
			rates[ cat ] = rates[ cat ] * ncategory / sum_rates;
	}

	/* BQM 2015-02-25: Testing if RAxML forgot this rate rescaling step */
	if (phylo_tree && phylo_tree->params && phylo_tree->params->no_rescale_gamma_invar)
		return;

    double newScale = 0.0;
    for (cat = 0; cat < ncategory; cat++)
        newScale += rates[cat];

    if (newScale != curScale) {
        for (cat = 0; cat < ncategory; cat++)
            rates[cat] *= curScale/newScale;    
    }

	/* if invariable sites are present */
//	if (Params::getInstance().optimize_alg_gammai != "EM") {
//		double p_inv = getPInvar();
//		for (cat = 0; cat < ncategory; cat++)
//			rates[cat] = rates[cat]/(1.0 - p_inv);
//	}

	/* check for very small rates */
//	for (cat = 0; cat < ncategory; cat ++)
//		if (rates[cat] < MIN_GAMMA_RATE)
//			rates[cat] = MIN_GAMMA_RATE;
}

/*double RateGamma::cmpPerPointGamma (const double prob, const double shape) {
}*/

void RateGamma::computeRatesMean () {
	int i;
	double lnga1=cmpLnGamma(gamma_shape+1);
	double *freqK = new double[ncategory];
	for (i=0; i<ncategory-1; i++) /* cutting points, Eq. 9 */
		freqK[i]=cmpPointChi2((i+1.0)/ncategory, 2.0 * gamma_shape) / (2.0 * gamma_shape);
	for (i=0; i<ncategory-1; i++) /* Eq. 10 */
		freqK[i]=cmpIncompleteGamma(freqK[i]*gamma_shape, gamma_shape+1, lnga1);

	rates[0] = freqK[0]*ncategory;
	rates[ncategory-1] = (1-freqK[ncategory-2])*ncategory;
	for (i=1; i<ncategory-1; i++)  rates[i] = (freqK[i]-freqK[i-1])*ncategory;
	delete [] freqK;
}

void RateGamma::setGammaShape(double gs) {
	gamma_shape = gs;
    computeRates();
}

double RateGamma::computeFunction(double shape) {
	if (gamma_shape != shape) {
		gamma_shape = shape;
		computeRates();
		phylo_tree->clearAllPartialLH();
	}
	return -phylo_tree->computeLikelihood();
}

double RateGamma::targetFunk(double x[]) {
	getVariables(x);
	phylo_tree->clearAllPartialLH();
	return -phylo_tree->computeLikelihood();
}

void RateGamma::setBounds(double *lower_bound, double *upper_bound, bool *bound_check) {
	if (getNDim() == 0) return;
	lower_bound[1] = phylo_tree->params->min_gamma_shape;
	upper_bound[1] = MAX_GAMMA_SHAPE;
	bound_check[1] = false;
}

void RateGamma::setVariables(double *variables) {
	if (getNDim() == 0) return;
	variables[1] = gamma_shape;
}

bool RateGamma::getVariables(double *variables) {
	if (getNDim() == 0) return false;
    bool changed = (gamma_shape != variables[1]);
	gamma_shape = variables[1];
    if (changed)
        computeRates();
    return changed;
}

double RateGamma::optimizeParameters(double gradient_epsilon, double min_gamma, double max_gamma) {
	if (fix_gamma_shape)
		return phylo_tree->computeLikelihood();
	if (verbose_mode >= VB_MAX)
		cout << "Optimizing gamma shape..." << endl;
	double negative_lh;
	double current_shape = gamma_shape;
	double ferror, optx;
	optx = minimizeOneDimen(min_gamma, current_shape, max_gamma, max(gradient_epsilon, TOL_GAMMA_SHAPE), &negative_lh, &ferror);
//	if (gamma_shape != optx) {
//		gamma_shape = optx;
//		computeRates();
//		phylo_tree->clearAllPartialLH();
//	}
//	return phylo_tree->computeLikelihood();
	return -computeFunction(optx);
}

double RateGamma::optimizeParameters(double gradient_epsilon) {
	if (fix_gamma_shape)
		return phylo_tree->computeLikelihood();
	if (verbose_mode >= VB_MAX)
		cout << "Optimizing gamma shape..." << endl;
	double negative_lh;
	double current_shape = gamma_shape;
	double ferror, optx;
	optx = minimizeOneDimen(phylo_tree->params->min_gamma_shape, current_shape, MAX_GAMMA_SHAPE, max(gradient_epsilon, TOL_GAMMA_SHAPE), &negative_lh, &ferror);
//	gamma_shape = optx;
//	computeRates();
//	phylo_tree->clearAllPartialLH();
//	return -negative_lh;
	return -computeFunction(optx);
}

void RateGamma::writeInfo(ostream &out) {
	out << "Gamma shape alpha: " << gamma_shape << endl;
	//out << " (" << (cut_median ? "median" : "mean") << " rate per category)" << endl;
	//out << "Number of categories: " << ncategory << endl;
}

void RateGamma::writeParameters(ostream &out) {
	out << "\t" << gamma_shape;
}

int RateGamma::computePatternRates(DoubleVector &pattern_rates, IntVector &pattern_cat) {
	//cout << "Computing Gamma site rates by empirical Bayes..." << endl;

	phylo_tree->computePatternLhCat(WSL_RATECAT);

	int npattern = phylo_tree->aln->getNPattern();
	pattern_rates.resize(npattern);
	pattern_cat.resize(npattern);

	double *lh_cat = phylo_tree->_pattern_lh_cat;
	for (int i = 0; i < npattern; i++) {
		double sum_rate = 0.0, sum_lh = 0.0;
		int best = 0;
		for (int c = 0; c < ncategory; c++) {
			sum_rate += rates[c] * lh_cat[c];
			sum_lh += lh_cat[c];
			if (lh_cat[c] > lh_cat[best] || (lh_cat[c] == lh_cat[best] && random_double()<0.5))  // break tie at random
                best = c;
		}
		pattern_rates[i] = sum_rate / sum_lh;
		pattern_cat[i] = best;
		lh_cat += ncategory;
	}
    return ncategory;

//	pattern_rates.clear();
//	pattern_rates.insert(pattern_rates.begin(), ptn_rates, ptn_rates + npattern);
//	pattern_cat.resize(npattern, 0);
//	for (int i = 0; i < npattern; i++)
//		for (int j = 1; j < ncategory; j++)
//			if (fabs(rates[j] - ptn_rates[i]) < fabs(rates[pattern_cat[i]] - ptn_rates[i]))
//				pattern_cat[i] = j;
//	delete [] ptn_rates;
}


/*NUMERICAL SUBROUTINES
**************************************************************************************

**************************************************************************************
**************************************************************************************
**************************************************************************************
**************************************************************************************/

/* THE FOLLOWING CODE COMES FROM tools.c in Yang's PAML package */

//----------------------------------------------------------------------------------------
double RateGamma::cmpLnGamma (double alpha) {
	/* returns ln(gamma(alpha)) for alpha>0, accurate to 10 decimal places.
	   Stirling's formula is used for the central polynomial part of the procedure.
	   Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
	   Communications of the Association for Computing Machinery, 9:684
	*/
	double x=alpha, f=0, z;

	if (x<7) {
		f=1;  z=x-1;
		while (++z<7)  f*=z;
		x=z;   f=-log(f);
	}
	z = 1/(x*x);
	return  f + (x-0.5)*log(x) - x + .918938533204673
	        + (((-.000595238095238*z+.000793650793651)*z-.002777777777778)*z
	           +.083333333333333)/x;
} //end of function cmpLnGamma

//----------------------------------------------------------------------------------------
double RateGamma::cmpIncompleteGamma (double x, double alpha, double ln_gamma_alpha) {
	/* returns the incomplete gamma ratio I(x,alpha) where x is the upper
		   limit of the integration and alpha is the shape parameter.
	   returns (-1) if in error
	   (1) series expansion     if (alpha>x || x<=1)
	   (2) continued fraction   otherwise

	   RATNEST FORTRAN by
	   Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
	   19: 285-287 (AS32)
	*/

	int i;
	double p=alpha, g=ln_gamma_alpha;
	double accurate=1e-8, overflow=1e30;
	double factor, gin=0, rn=0, a=0,b=0,an=0,dif=0, term=0, pn[6];

	if (x==0) return (0);
	if (x<0 || p<=0) return (-1);

	factor=exp(p*log(x)-x-g);
	if (x>1 && x>=p) goto l30;
	/* (1) series expansion */
	gin=1;  term=1;  rn=p;
l20:
	rn++;
	term*=x/rn;   gin+=term;

	if (term > accurate) goto l20;

	gin*=factor/p;
	goto l50;
l30:

	/* (2) continued fraction */
	a=1-p;   b=a+x+1;  term=0;
	pn[0]=1;  pn[1]=x;  pn[2]=x+1;  pn[3]=x*b;
	gin=pn[2]/pn[3];
l32:
	a++;  b+=2;  term++;   an=a*term;
	for (i=0; i<2; i++) pn[i+4]=b*pn[i+2]-an*pn[i];
	if (pn[5] == 0) goto l35;
	rn=pn[4]/pn[5];   dif=fabs(gin-rn);
	if (dif>accurate) goto l34;
	if (dif<=accurate*rn) goto l42;
l34:
	gin=rn;
l35:
	for (i=0; i<4; i++) pn[i]=pn[i+2];
	if (fabs(pn[4]) < overflow) goto l32;
	for (i=0; i<4; i++) pn[i]/=overflow;
	goto l32;
l42:
	gin=1-factor*gin;

l50:
	return (gin);
} //end of function cmpIncompleteGamma


//----------------------------------------------------------------------------------------
/* functions concerning the CDF and percentage points of the gamma and
   Chi2 distribution
*/
double RateGamma::cmpPointNormal (double prob) {
	/* returns z so that Prob{x<z}=prob where x ~ N(0,1) and (1e-12)<prob<1-(1e-12)
	   returns (-9999) if in error
	   Odeh RE & Evans JO (1974) The percentage points of the normal distribution.
	   Applied Statistics 22: 96-97 (AS70)

	   Newer methods:
	     Wichura MJ (1988) Algorithm AS 241: the percentage points of the
	       normal distribution.  37: 477-484.
	     Beasley JD & Springer SG  (1977).  Algorithm AS 111: the percentage
	       points of the normal distribution.  26: 118-121.

	*/
	double a0=-.322232431088, a1=-1, a2=-.342242088547, a3=-.0204231210245;
	double a4=-.453642210148e-4, b0=.0993484626060, b1=.588581570495;
	double b2=.531103462366, b3=.103537752850, b4=.0038560700634;
	double y, z=0, p=prob, p1;

	p1 = (p<0.5 ? p : 1-p);

	if (p1<1e-20) return (-9999);

	y = sqrt (log(1/(p1*p1)));
	z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0);
	return (p<0.5 ? -z : z);
} //end of function cmpPointNormal



//----------------------------------------------------------------------------------------

double RateGamma::cmpPointChi2 (double prob, double v) {
	/* returns z so that Prob{x<z}=prob where x is Chi2 distributed with df=v
	   returns -1 if in error.   0.000002<prob<0.999998
	   RATNEST FORTRAN by
	       Best DJ & Roberts DE (1975) The percentage points of the
	       Chi2 distribution.  Applied Statistics 24: 385-388.  (AS91)
	   Converted into C by Ziheng Yang, Oct. 1993.
	*/
	double e=.5e-6, aa=.6931471805, p=prob, g;
	double xx, c, ch, a=0,q=0,p1=0,p2=0,t=0,x=0,b=0,s1,s2,s3,s4,s5,s6;

	if (p<.000002 || p>.999998 || v<=0) return (-1);

	g = cmpLnGamma (v/2);
	xx=v/2;   c=xx-1;
	if (v >= -1.24*log(p)) goto l1;

	ch=pow((p*xx*exp(g+xx*aa)), 1/xx);
	if (ch-e<0) return (ch);
	goto l4;
l1:
	if (v>.32) goto l3;
	ch=0.4;   a=log(1-p);
l2:
	q=ch;  p1=1+ch*(4.67+ch);  p2=ch*(6.73+ch*(6.66+ch));
	t=-0.5+(4.67+2*ch)/p1 - (6.73+ch*(13.32+3*ch))/p2;
	ch-=(1-exp(a+g+.5*ch+c*aa)*p2/p1)/t;
	if (fabs(q/ch-1)-.01 <= 0) goto l4;
	else                       goto l2;

l3:
	x=cmpPointNormal (p);
	p1=0.222222/v;   ch=v*pow((x*sqrt(p1)+1-p1), 3.0);
	if (ch>2.2*v+6)  ch=-2*(log(1-p)-c*log(.5*ch)+g);
l4:

	do {
		q=ch;   p1=.5*ch;
		if ((t=cmpIncompleteGamma (p1, xx, g))<0) {
			return (-1);
		}
		p2=p-t;
		t=p2*exp(xx*aa+g+p1-c*log(ch));
		b=t/ch;  a=0.5*t-b*c;

		s1=(210+a*(140+a*(105+a*(84+a*(70+60*a))))) / 420;
		s2=(420+a*(735+a*(966+a*(1141+1278*a))))/2520;
		s3=(210+a*(462+a*(707+932*a)))/2520;
		s4=(252+a*(672+1182*a)+c*(294+a*(889+1740*a)))/5040;
		s5=(84+264*a+c*(175+606*a))/2520;
		s6=(120+c*(346+127*c))/5040;
		ch+=t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
	} while (fabs(q/ch-1) > e);

	return (ch);
} //end of function cmpPointChi2


/* THE END OF THE CODES COMMING FROM tools.c in Yang's PAML package */
