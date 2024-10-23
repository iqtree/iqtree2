//
// C++ Implementation: optimization
//
// Description: 
//
//
// Author: BUI Quang Minh, Steffen Klaere, Arndt von Haeseler <minh.bui@univie.ac.at>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "optimization.h"
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "lbfgsb/lbfgsb_new.h"
#include "tools.h"


using namespace std;

const double ERROR_X = 1.0e-4;

double ran1(long *idum);
double *new_vector(long nl, long nh);
void free_vector(double *v, long nl, long nh);
double **new_matrix(long nrl, long nrh, long ncl, long nch);
void free_matrix(double **m, long nrl, long nrh, long ncl, long nch);
void fixBound(double x[], double lower[], double upper[], int n);

#define NR_END 1
#define FREE_ARG char*

#define GET_PSUM \
					for (n=1;n<=ndim;n++) {\
					for (sum=0.0,m=1;m<=mpts;m++) sum += p[m][n];\
					psum[n]=sum;}


/*
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran1(long *idum) {
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
*/

long idum = 123456;
double tt;

void nrerror(const char *error_text)
/* Numerical Recipes standard error handler */
{
	cerr << "NUMERICAL ERROR: " << error_text << endl;
	//exit(1);
	abort();
}

double *new_vector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

double **new_matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}


void free_vector(double *v, long nl, long nh)
/* free a double vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_matrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void fixBound(double x[], double lower[], double upper[], int n) {
	for (int i = 1; i <= n; i++) {
		if (x[i] < lower[i])
			x[i] = lower[i];
		else if (x[i] > upper[i])
			x[i] = upper[i];
	}
}

/**********************************************
	Optimization routines
**********************************************/
Optimization::Optimization()
{
}


Optimization::~Optimization()
{
}

/*****************************************************
	One dimensional optimization with Brent method
*****************************************************/
	
#define ITMAX 100
#define CGOLD 0.3819660
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#define TRY_AITKEN 0
/* Brents method in one dimension */
double Optimization::brent_opt (double ax, double bx, double cx, double tol,
                          double *foptx, double *f2optx, double fax, double fbx, double fcx) {
	int iter;
	double a,b,d=0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	double xw,wv,vx;
	double e=0.0;

	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=bx;
	fx=fbx;
	if (fax < fcx) {
		w=ax;
		fw=fax;
		v=cx;
		fv=fcx;
	} else {
		w=cx;
		fw=fcx;
		v=ax;
		fv=fax;
	}

#if TRY_AITKEN
    double delta = 0;     //change in U
    double lastDelta = 0; //previous change in U
    double lastU = 0;     //previous value of U
#endif
	for (iter=1;iter<=ITMAX;iter++) {
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			*foptx = fx;
			xw = x-w;
			wv = w-v;
			vx = v-x;
			*f2optx = 2.0*(fv*xw + fx*wv + fw*vx)/
			          (v*v*xw + x*x*wv + w*w*vx);
			return x;
		}

		if (fabs(e) > tol1) {
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0)
				p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			else {
				d=p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			}
		} else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}

		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu=computeFunction(u);
        
#if TRY_AITKEN
        //In a Nutshell, Aitken's sequence accelerator constructs
        //an auxiliary sequence, using
        //"the main sequence, minus the square of the delta divided
        // by the delta of the delta".
        
        double aitken = 0;
        bool   haveAitken = false;
        if (iter>=2) {
            lastDelta = delta;
            delta     = u - lastU;
            //For now. Later, use a calculated tolerance
            if (iter>=3 && 0.00001 < fabs(delta-lastDelta) && 0.00001 < fabs(delta) ) {
                haveAitken = true;
                aitken = u - ( delta * delta ) / (delta - lastDelta);
            }
        }
        cout.precision(6);
        if (!haveAitken) {
            cout << "f(" << u << ")=" << fu << endl;
        } else {
            cout << "f(" << u << ")=" << fu << ", A=" << aitken << ", f(A)=" << computeFunction(aitken) << endl;
        }
        lastU = u;
#endif

		if (fu <= fx) {
			if (u >= x)
				a=x;
			else
				b=x;

			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		} else {
			if (u < x)
				a=u;
			else
				b=u;
			if (fu <= fw || w == x) {
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			} else
				if (fu <= fv || v == x || v == w) {
					v=u;
					fv=fu;
				}
		}
	}

	*foptx = fx;
	xw = x-w;
	wv = w-v;
	vx = v-x;
	*f2optx = 2.0*(fv*xw + fx*wv + fw*vx)/(v*v*xw + x*x*wv + w*w*vx);

	return x;
}

#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT
#undef SIGN
#undef GOLD
#undef GLIMIT
#undef TINY


double Optimization::minimizeOneDimen(double xmin, double xguess, double xmax, double tolerance, double *fx, double *f2x) {
	double eps, optx, ax, bx, cx, fa, fb, fc;
	//int    converged;	/* not converged error flag */
		
	/* first attempt to bracketize minimum */
	if (xguess < xmin) xguess = xmin;
	if (xguess > xmax) xguess = xmax;
	eps = xguess*tolerance*50.0;
	ax = xguess - eps;
	if (ax < xmin) ax = xmin;
	bx = xguess;
	cx = xguess + eps;
	if (cx > xmax) cx = xmax;
	
	/* check if this works */
    // compute fb first to save some computation, if any
	fb = computeFunction(bx);
	fa = (ax<bx) ? computeFunction(ax) : fb;
	fc = (bx<cx) ? computeFunction(cx) : fb;

	/* if it works use these borders else be conservative */
	if ((fa < fb) || (fc < fb)) {
		if (ax != xmin) fa = computeFunction(xmin);
		if (cx != xmax) fc = computeFunction(xmax);
		ax = xmin;
		cx = xmax;
	}
	/*
	const int MAX_ROUND = 10;
	for (i = 0; ((fa < fb-tolerance) || (fc < fb-tolerance)) && (i<MAX_ROUND); i++) {
		// brent method require that fb is smaller than both fa and fc
		// find some random values until fb achieve this
			bx = (((double)rand()) / RAND_MAX)*(cx-ax) + ax;
			fb = computeFunction(bx);
	}*/

/*
	if ((fa < fb) || (fc < fb)) {
		if (fa < fc) { bx = ax; fb = fa; } else { bx = cx; fb = fc; }
		//cout << "WARNING: Initial value for Brent method is set at bound " << bx << endl;
	}*/
	//	optx = brent_opt(xmin, xguess, xmax, tolerance, fx, f2x, fa, fb, fc);
	//} else
	optx = brent_opt(ax, bx, cx, tolerance, fx, f2x, fa, fb, fc);
    if (*fx > fb) // if worse, return initial value
    {
        *fx = computeFunction(bx);
        return bx;
    }

	return optx; /* return optimal x */
}

double Optimization::minimizeOneDimenSafeMode(double xmin, double xguess, double xmax, double tolerance, double *f)
{
	double ferror;
	double optx = minimizeOneDimen(xmin, xguess, xmax, tolerance, f, &ferror);
	double fnew;
	// check value at the boundary
	if ((optx < xmax) && (fnew = computeFunction(xmax)) <= *f+tolerance) {
		//if (verbose_mode >= VB_MAX)
			//cout << "Note from Newton safe mode: " << optx << " (" << f << ") -> " << xmax << " ("<< fnew << ")" << endl;
		optx = xmax;
		*f = fnew;
	}
	if ((optx > xmin) && (fnew = computeFunction(xmin)) <= *f+tolerance) {
		//if (verbose_mode >= VB_MAX)
			//cout << "Note from Newton safe mode: " << optx << " -> " << xmin << endl;
		optx = xmin;
		*f = fnew;
	}
	return optx;
}

/*****************************************************
	One dimensional optimization with Newton Raphson 
	only applicable if 1st and 2nd derivatives are easy to compute
*****************************************************/


double Optimization::minimizeNewtonSafeMode(double xmin, double xguess, double xmax, double tolerance, double &f)
{
	double optx = minimizeNewton(xmin, xguess, xmax, tolerance, f);
	double fnew;
	// check value at the boundary
	if ((optx < xmax) && (fnew = computeFunction(xmax)) <= f+tolerance) {
		//if (verbose_mode >= VB_MAX)
			//cout << "Note from Newton safe mode: " << optx << " (" << f << ") -> " << xmax << " ("<< fnew << ")" << endl;
		optx = xmax;
		f = fnew;
	}
	if ((optx > xmin) && (fnew = computeFunction(xmin)) <= f+tolerance) {
		//if (verbose_mode >= VB_MAX)
			//cout << "Note from Newton safe mode: " << optx << " -> " << xmin << endl;
		optx = xmin;
		f = fnew;
	}
	return optx;
}

double Optimization::minimizeNewton(double x1, double xguess, double x2, double xacc, double &d2l, int maxNRStep)
{
	int j;
	double df,dx,dxold,f;
	double temp,xh,xl,rts, rts_old, xinit;

	rts = xguess;
	if (rts < x1) rts = x1;
	if (rts > x2) rts = x2;
	xinit = xguess;
//	finit = fold = fm = computeFuncDerv(rts,f,df);
    computeFuncDerv(rts,f,df);
	d2l = df;
	if (!isfinite(f) || !isfinite(df)) {
		nrerror("Wrong computeFuncDerv");
	}
	if (df >= 0.0 && fabs(f) < xacc) return rts;
	if (f < 0.0) {
		xl = rts;
		xh = x2;
	} else {
		xh = rts;
		xl = x1;	
	}

	dx=dxold=fabs(xh-xl);
	for (j=1;j<=maxNRStep;j++) {
		rts_old = rts;
		if (
			(df <= 0.0) // function is concave
//			|| (fm > fold + xacc) // increasing
			|| (((rts-xh)*df-f)*((rts-xl)*df-f) >= 0.0) // out of bound
			//|| (fabs(2.0*f) > fabs(dxold*df))  // converge too slow
			) {
			dxold=dx;
			dx=0.5*(xh-xl);
			rts=xl+dx;
            d2l = df;
			if (xl == rts) return rts;
		} else {
			dxold=dx;
			dx=f/df;
			temp=rts;
			rts -= dx;
			d2l = df;
			if (temp == rts) return rts;
		}
		if (fabs(dx) < xacc || (j == maxNRStep)) {
//			if (fm > finit) {
//				// happen in rare cases that it is worse than starting point: revert init value
//				fm = computeFunction(xinit);
//				return xinit;
//			}
			return rts_old;
//			return rts;
		}
//		fold = fm;
//		fm = computeFuncDerv(rts,f,df);
        computeFuncDerv(rts,f,df);
		if (!isfinite(f) || !isfinite(df)) nrerror("Wrong computeFuncDerv");
		if (df > 0.0 && fabs(f) < xacc) {
			d2l = df;
//			if (fm > finit) {
//				// happen in rare cases that it is worse than starting point: revert init value
//				fm = computeFunction(xinit);
//				return xinit;
//			}
			return rts;
		}
		if (f < 0.0)
			xl=rts;
		else if (f > 0.0)
			xh=rts;
	}
	nrerror("Maximum number of iterations exceeded in minimizeNewton");
	d2l = 0.0;
	return 0.0;
}

double Optimization::minimizeNewton(double x1, double xguess, double x2, double xacc, int maxNRStep)
{
	double var;
	double optx = minimizeNewton(x1, xguess, x2, xacc, var, maxNRStep);
	return optx;
}

/*****************************************************
	Multi dimensional optimization with Newton Raphson
	only applicable if 1st and 2nd derivatives are easy to compute
*****************************************************/

int matinv (double x[], int n, int m, double space[]);

double Optimization::minimizeNewtonMulti(double *x1, double *xguess, double *x2, double xacc, int N, int maxNRStep)
{
	int i, step;
    int N_2 = N*N;
	double df[N_2], dx[N], dxold[N], f[N+1], space[N];
	double temp, xh[N], xl[N], rts[N], rts_old[N];

    for (i = 0; i < N; i++) {
        rts[i] = xguess[i];
        if (rts[i] < x1[i]) rts[i] = x1[i];
        if (rts[i] > x2[i]) rts[i] = x2[i];
    }
    computeFuncDervMulti(rts, f, df);
    double score = f[N];
//    double score = targetFunk(rts-1);
//    cout << "Newton step 0: " << score << endl;

    bool stop = true;

    for (i = 0; i < N; i++) {
        if (!(fabs(f[i]) < xacc))
            stop = false;
        xl[i] = x1[i];
        xh[i] = x2[i];
//        if (f[i] < 0.0) {
//            xl[i] = rts[i];
//            xh[i] = x2[i];
//        } else {
//            xh[i] = rts[i];
//            xl[i] = x1[i];
//        }
        dx[i] = dxold[i] = fabs(xh[i]-xl[i]);
    }

    if (stop) {
        memcpy(xguess, rts, sizeof(double)*N);
        return score;
    }

	for (step = 1; step <= maxNRStep; step++) {
		memcpy(rts_old, rts, sizeof(double)*N);

        // first compute inverse of hessian matrix
        matinv(df, N, N, space);
        // now take the product
        for (i = 0; i < N; i++) {
            space[i] = 0.0;
            for (int k = 0; k < N; k++) {
                space[i] += df[i*N+k] * f[k];
            }
        }

         // identify alpha to be not out of bound
        double alpha;
        for (alpha = 1.0; alpha > xacc; alpha *= 0.5) {
            bool ok_bound = true;
            for (i = 0; i < N; i++)
                if (rts[i]-alpha*space[i] < xl[i] || rts[i]-alpha*space[i] > xh[i]) {
                    ok_bound = false;
                    break;
                }
            if (ok_bound) break;
        }

        do {
            stop = true;
            for (i = 0; i < N; i++) {
                dxold[i] = dx[i];
                dx[i] = alpha*space[i];
                temp=rts[i];
                rts[i] -= dx[i];
                if (fabs(dx[i]) >= xacc && (step != maxNRStep)) {
                    stop = false;
                } else {
                    rts[i] = rts_old[i];
                }
            }
            double new_score = targetFunk(rts-1);
            if (new_score <= score) break;
            // score increased, reduce step size by half
            memcpy(rts, rts_old, sizeof(double)*N);
            alpha *= 0.5;
        } while (!stop);
        
        if (stop) {
            break;
        }
        computeFuncDervMulti(rts, f, df);
        score = f[N];

//        double score = targetFunk(rts-1);
//        cout << "Newton step " << step << ": " << score << endl;
//        for (i = 0; i < N; i++) cout << "  " << rts[i];
//        cout << endl;

        stop = true;
        for (i = 0; i < N; i++) {
            if (!(fabs(f[i]) < xacc)) {
                stop = false;
            }
//            if (f[i] < 0.0)
//                xl[i] = rts[i];
//            else
//                xh[i] = rts[i];
        }
        if (stop)
            break;
	}

    // copy returned x-value
    memcpy(xguess, rts, sizeof(double)*N);

    if (stop)
        return score;

    cout << "Maximum number of iterations exceeded in minimizeNewtonMulti" << endl;
    return score;
}


/*****************************************************
	Multi dimensional optimization with BFGS method
*****************************************************/

#define ALF 1.0e-4
#define TOLX 1.0e-7
//static double maxarg1,maxarg2;
//#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
//        (maxarg1) : (maxarg2))

void Optimization::lnsrch(int n, double xold[], double fold, double g[], double p[], double x[],
                   double *f, double stpmax, int *check, double lower[], double upper[]) {
	int i;
	double a,alam,alam2=0,alamin,b,disc,f2=0,fold2=0,rhs1,rhs2,slope,sum,temp,
	test,tmplam;

	*check=0;
	for (sum=0.0,i=1;i<=n;i++) sum += p[i]*p[i];
	sum=sqrt(sum);
	if (sum > stpmax)
		for (i=1;i<=n;i++) p[i] *= stpmax/sum;
	for (slope=0.0,i=1;i<=n;i++)
		slope += g[i]*p[i];
	test=0.0;
	for (i=1;i<=n;i++) {
		temp=fabs(p[i])/max(fabs(xold[i]),1.0);
		if (temp > test) test=temp;
	}
	alamin=TOLX/test;
	alam=1.0;
	bool first_time = true;
	for (;;) {
//        for (; alam >= alamin;) {
//            bool out_of_bound = false;
//            for (i=1;i<=n;i++) {
//                x[i]=xold[i]+alam*p[i];
//                if (x[i] < lower[i] || x[i] > upper[i]) {
//                    out_of_bound = true;
//                    break;
//                }
//            }
//            if (out_of_bound)
//                alam *= 0.5;
//            else
//                break;
//        }

        for (i=1;i<=n;i++) x[i]=xold[i]+alam*p[i];

		fixBound(x, lower, upper, n);
		//checkRange(x);
		*f=targetFunk(x);
		if (alam < alamin) {
			for (i=1;i<=n;i++) x[i]=xold[i];
			*check=1;
			return;
		} else if (*f <= fold+ALF*alam*slope) return;
		else {
			if (first_time)
				tmplam = -slope/(2.0*(*f-fold-slope));
			else {
				rhs1 = *f-fold-alam*slope;
				rhs2=f2-fold2-alam2*slope;
				a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
				b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
				if (a == 0.0) tmplam = -slope/(2.0*b);
				else {
					disc=b*b-3.0*a*slope;
					if (disc<0.0) //nrerror("Roundoff problem in lnsrch.");
						tmplam = 0.5 * alam;
					else if (b <= 0.0) tmplam=(-b+sqrt(disc))/(3.0*a);
					else tmplam = -slope/(b+sqrt(disc));
				}
				if (tmplam>0.5*alam)
					tmplam=0.5*alam;
			}
		}
		alam2=alam;
		f2 = *f;
		fold2=fold;
		alam=max(tmplam,0.1*alam);
		first_time = false;
	}
}
#undef ALF
#undef TOLX


const int MAX_ITER = 3;
extern double random_double(int *rstream);

bool Optimization::restartParameters(double guess[], int ndim, double lower[], double upper[], bool bound_check[], int iteration) {
    int i;
    bool restart = false;
    if (iteration <= MAX_ITER) {
        for (i = 1; i <= ndim; i++) {
	    if (bound_check[i]) {
	        if (fabs(guess[i]-lower[i]) < 1e-4 || fabs(guess[i]-upper[i]) < 1e-4) {
		     restart = true;
	             break;
		} // if (fabs...
	    } // if bound_check
	} // for
    } // if iteration
    if (restart) {
      // TODO: Dominik observed that this warning is printed even during model
      // testing. A quieter solution might be preferred.
	cout << "Restart estimation at the boundary... " << std::endl;
        for (i = 1; i <= ndim; i++) {
	    guess[i] = random_double() * (upper[i] - lower[i])/3 + lower[i];
	}
    }	
    return (restart);
}

double Optimization::minimizeMultiDimen(double guess[], int ndim, double lower[], double upper[], bool bound_check[], double gtol, double *hessian) {
    int i, iter;
    double fret, minf = 1e+12;
    double *minx = new double [ndim+1];
    int count = 0;
    bool restart;
    do {
        dfpmin(guess, ndim, lower, upper, gtol, &iter, &fret, hessian);
        if (fret < minf) {
            minf = fret;
            for (i = 1; i <= ndim; i++) {
                minx[i] = guess[i];
            }
        }
        count++;
        // restart the search if at the boundary
        // it's likely to end at a local optimum at the boundary
        restart = restartParameters(guess,ndim,lower,upper,bound_check,count);
    } while (restart);
    if (count > 1) {
        for (i = 1; i <= ndim; i++) {
            guess[i] = minx[i];
        }
        fret = minf;
    }
    delete [] minx;
    
    return fret;
}


#define ITMAX 200
#define SQR(a) ((a)*(a))
#define EPS 3.0e-8
#define TOLX (4*EPS)
#define STPMX 100.0

#define FREEALL free_vector(xi,1,n);free_vector(pnew,1,n); \
free_matrix(hessin,1,n,1,n);free_vector(hdg,1,n);free_vector(g,1,n); \
free_vector(dg,1,n);



void Optimization::dfpmin(double p[], int n, double lower[], double upper[]
                          , double gtol, int *iter, double *fret, double *hessian) {
    int check,i,its,j;
    double den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test;
    double *dg,*g,*hdg,**hessin,*pnew,*xi;
    
    dg=new_vector(1,n);
    g=new_vector(1,n);
    hdg=new_vector(1,n);
    hessin=new_matrix(1,n,1,n);
    pnew=new_vector(1,n);
    xi=new_vector(1,n);
    fp = derivativeFunk(p,g);
    for (i=1;i<=n;++i)
    {
        xi[i] = -g[i];
        sum += p[i]*p[i];
    }
    if (hessian) {
        int pos = 0;
        for (i=1; i<=n; i++) {
            for (j=1; j<=n; j++, ++pos) {
                hessin[i][j] = hessian[pos];
            }
        }
    } else {
        for (i=1;i<=n;i++) {
            for (j=1;j<=n;j++) {
                hessin[i][j]=0.0;
            }
            hessin[i][i]=1.0;
        }
    }
	//checkBound(p, xi, lower, upper, n);
	//checkDirection(p, xi);

	stpmax=STPMX*max(sqrt(sum),(double)n);
	for (its=1;its<=ITMAX;its++) {
		*iter=its;
		lnsrch(n,p,fp,g,xi,pnew,fret,stpmax,&check, lower, upper);
		fp = *fret;
		for (i=1;i<=n;i++) {
			xi[i]=pnew[i]-p[i];
			p[i]=pnew[i];
		}
		test=0.0;
		for (i=1;i<=n;i++) {
			temp=fabs(xi[i])/max(fabs(p[i]),1.0);
			if (temp > test) test=temp;
		}
		if (test < TOLX) {
			FREEALL
			return;
		}
		for (i=1;i<=n;i++) dg[i]=g[i];
		derivativeFunk(p,g);
		test=0.0;
		den=max(fabs(*fret),1.0); // fix bug found by Tung, as also suggested by NR author
		for (i=1;i<=n;i++) {
			temp=fabs(g[i])*max(fabs(p[i]),1.0)/den;
			if (temp > test) test=temp;
		}
		if (test < gtol) {
			FREEALL
			return;
		}
        for (i=1;i<=n;i++) {
            dg[i]=g[i]-dg[i];
        }
        for (i=1;i<=n;i++) {
            hdg[i]=0.0;
            for (j=1;j<=n;j++) {
                hdg[i] += hessin[i][j]*dg[j];
            }
        }
		fac=fae=sumdg=sumxi=0.0;
		for (i=1;i<=n;i++) {
			fac += dg[i]*xi[i];
			fae += dg[i]*hdg[i];
			sumdg += SQR(dg[i]);
			sumxi += SQR(xi[i]);
		}
		if (fac*fac > EPS*sumdg*sumxi)
		{
			fac=1.0/fac;
			fad=1.0/fae;
            for (i=1;i<=n;i++) {
                dg[i]=fac*xi[i]-fad*hdg[i];
            }
			for (i=1;i<=n;i++) {
				for (j=1;j<=n;j++) {
					hessin[i][j] += fac*xi[i]*xi[j]
					                -fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
				}
			}
		}
		for (i=1;i<=n;i++) {
			xi[i]=0.0;
			for (j=1;j<=n;j++) xi[i] -= hessin[i][j]*g[j];
		}
		//checkBound(p, xi, lower, upper, n);
		//checkDirection(p, xi);
		//if (*iter > 200) cout << "iteration=" << *iter << endl;
	}
	// BQM: disable this message!
	//nrerror("too many iterations in dfpmin");
	FREEALL
}
#undef ITMAX
#undef SQR
#undef EPS
#undef TOLX
#undef STPMX
#undef FREEALL
#undef FMAX


/**
	the approximated derivative function
	@param x the input vector x
	@param dfx the derivative at x
	@return the function value at x
*/
double Optimization::derivativeFunk(double x[], double dfx[]) {
	/*
	if (!checkRange(x))
		return INFINITIVE;
	*/
	int ndim = getNDim();
	double *h = new double[ndim+1];
    double temp;
    int dim;
	double fx = targetFunk(x);
	for (dim = 1; dim <= ndim; dim++ ){
		temp = x[dim];
		h[dim] = ERROR_X * fabs(temp);
		if (h[dim] == 0.0) h[dim] = ERROR_X;
		x[dim] = temp + h[dim];
		h[dim] = x[dim] - temp;
		dfx[dim] = (targetFunk(x));
		x[dim] = temp;
	}
	for (dim = 1; dim <= ndim; dim++ )
        dfx[dim] = (dfx[dim] - fx) / h[dim];
    delete [] h;
	return fx;
}


/*#define NRANSI
#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

double Optimization::brent(double ax, double bx, double cx, double tol,
	double *xmin)
{
	int iter;
	double a,b,d=0.0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	double e=0.0;

	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=computeFunction(x);
	for (iter=1;iter<=ITMAX;iter++) {
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			*xmin=x;
			return fx;
		}
		if (fabs(e) > tol1) {
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			else {
				d=p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			}
		} else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu=computeFunction(u);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			} else if (fu <= fv || v == x || v == w) {
				v=u;
				fv=fu;
			}
		}
	}
	nrerror("Too many iterations in brent");
	*xmin=x;
	return fx;
}

#undef SIGN
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT
#undef NRANSI*/

/*#define JMAX 20

double Optimization::minimizeNewton(double xmin, double xguess, double xmax, double tolerance, double &f)
{
	return rtsafe(xmin, xguess, xmax, tolerance, f);
	//double fe;
	//return minimizeOneDimen(xmin, rtn, xmax, tolerance, &f, &fe);

	int j;
	double df,ddf,dx,rtn,rtnold, fstart=0, fnew;

	rtn=xguess;
	if (rtn < xmin) rtn = xmin;
	if (rtn > xmax) rtn = xmax;
	

	for (j=1;j<=JMAX;j++) {
		f = computeFuncDerv(rtn,df,ddf);
		if (!isfinite(f)) 
			return 0;
		if (j == 1) fstart = f;
		if (ddf == 0.0) break;
		dx=(df/fabs(ddf));
		if (fabs(dx) <= tolerance) break;

		rtnold = rtn; rtn = rtn-dx;
		if (rtn < xmin) rtn = xmin;
		if (rtn > xmax) rtn = xmax;
		dx = rtnold-rtn;

		while (fabs(dx) > tolerance && (fnew = computeFunction(rtn)) > f + tolerance) {
			dx /= 2;
			rtn = rtnold - dx;
		}
		if (fabs(dx) <= tolerance) { rtn = rtnold; break; }
	}
	//if (j > JMAX)
		//nrerror("Maximum number of iterations exceeded in Newton-Raphson");
	if (f <= fstart && j <= JMAX && (j > 1 || xguess > xmin+tolerance)) 
		return rtn;
	// Newton does not work, turn to other method
	double fe;
	return minimizeOneDimen(xmin, xguess, xmax, tolerance, &f, &fe);
}*/

/*
double Optimization::minimizeNewton(double xmin, double xguess, double xmax, double tolerance, double &f)
{
	int j;
	double df,ddf,dx,dxold,rtn,temp, fold, fstart;
	double xmin_orig = xmin, xmax_orig = xmax;

	rtn=xguess;
	dx=dxold=(xmax-xmin);
	fstart = fold = f = computeFuncDerv(rtn,df,ddf);

	for (j=1;j<=JMAX;j++) {
		if (ddf <= 0.0) break;
		if ((((rtn-xmax)*ddf-df) * ((rtn-xmin)*ddf-df) > 0) || // run out of range
			(fabs(2.0*df) > fabs(dxold*ddf)) // dx not decreasing fast enough
			) // f even increase
		{
			dxold = dx;
			dx = 0.5*(xmax-xmin);
			rtn = xmin+dx;
			if (xmin == rtn) break;
		} else {
			dxold=dx;
			if (ddf == 0.0)
				nrerror("2nd derivative is zero");
			dx=df/ddf;
			temp=rtn;
			//if (f > fold) dx /= 2.0;
			//if (ddf < 0) dx = -dx;
			rtn -= dx;
			//if (rtn < xmin) rtn = xmin;
			//if (rtn > xmax) rtn = xmax;
			dx = temp - rtn;
			if (temp == rtn) break;
		}
		if (fabs(dx) < tolerance) break;
		fold = f;
		f = computeFuncDerv(rtn,df,ddf);
		//if (f > fold) break; // Does not decrease function, escape
		if (df < 0.0) 
			xmin = rtn;
		else
			xmax = rtn;
	}
	if (j > JMAX)
	nrerror("Maximum number of iterations exceeded in Newton-Raphson");
	if (f <= fstart) return rtn;
	// Newton does not work (find a max instead of min), turn to other method
	double fe;
	return minimizeOneDimen(xmin_orig, xguess, xmax_orig, tolerance, &f, &fe);
	//return 0.0;
}
#undef JMAX*/


double Optimization::L_BFGS_B(int n, double* x, double* l, double* u, double pgtol, int maxit) {
	int i;
	double Fmin;
	int fail;
	int fncount;
	int grcount;
	char msg[100];

	int m = 10;          // number of BFGS updates retained in the "L-BFGS-B" method. It defaults to 5.

	int *nbd;           // 0: unbounded; 1: lower bounded; 2: both lower & upper; 3: upper bounded
	nbd = new int[n];
	for (i=0; i<n; i++)
		nbd[i] = 2;

	double factr = 1e+7; // control the convergence of the "L-BFGS-B" method.
	// Convergence occurs when the reduction in the object is within this factor
	// of the machine tolerance.
	// Default is 1e7, that is a tolerance of about 1e-8

//	double pgtol = 0;   // helps control the convergence of the "L-BFGS-B" method.
//    pgtol = 0.0;
	// It is a tolerance on the projected gradient in the current search direction.
	// Default is zero, when the check is suppressed

	int trace = 0;      // non-negative integer.
	if (verbose_mode >= VB_MAX) {
		trace = 1;
	}
	// If positive, tracing information on the progress of the optimization is produced.
	// Higher values may produce more tracing information.

	int nREPORT = 10;   // The frequency of reports for the "L-BFGS-B" methods if "trace" is positive.
	// Defaults to every 10 iterations.

/*#ifdef USE_OLD_PARAM
	lbfgsb(n, m, x, l, u, nbd, &Fmin, fn, gr1, &fail, ex,
			factr, pgtol, &fncount, &grcount, maxit, msg, trace, nREPORT);
#else*/

	lbfgsb(n, m, x, l, u, nbd, &Fmin, &fail,
			factr, pgtol, &fncount, &grcount, maxit, msg, trace, nREPORT);
//#endif

    if (fail == 51 || fail == 52) {
        cout << msg << endl;
    }

	delete[] nbd;
    
    return Fmin;
}

void Optimization::lbfgsb(int n, int m, double *x, double *l, double *u, int *nbd,
		double *Fmin, int *fail,
		double factr, double pgtol,
		int *fncount, int *grcount, int maxit, char *msg,
		int trace, int nREPORT)
{
	char task[60];
	double f, *g, dsave[29], *wa;
	int tr = -1, iter = 0, *iwa, isave[44], lsave[4];

	/* shut up gcc -Wall in 4.6.x */

	for(int i = 0; i < 4; i++) lsave[i] = 0;

	if(n == 0) { /* not handled in setulb */
		*fncount = 1;
		*grcount = 0;
		*Fmin = optimFunc(n, u);
		strcpy(msg, "NOTHING TO DO");
		*fail = 0;
		return;
	}
	if (nREPORT <= 0) {
		cerr << "REPORT must be > 0 (method = \"L-BFGS-B\")" << endl;
		exit(1);
	}
	switch(trace) {
	case 2: tr = 0; break;
	case 3: tr = nREPORT; break;
	case 4: tr = 99; break;
	case 5: tr = 100; break;
	case 6: tr = 101; break;
	default: tr = -1; break;
	}

	*fail = 0;
	g = (double*) malloc (n * sizeof(double));
	/* this needs to be zeroed for snd in mainlb to be zeroed */
	wa = (double *) malloc((2*m*n+4*n+11*m*m+8*m) * sizeof(double));
	iwa = (int *) malloc(3*n * sizeof(int));
	strcpy(task, "START");
	while(1) {
		/* Main workhorse setulb() from ../appl/lbfgsb.c : */
		setulb(n, m, x, l, u, nbd, &f, g, factr, &pgtol, wa, iwa, task,
				tr, lsave, isave, dsave);
		/*    Rprintf("in lbfgsb - %s\n", task);*/
		if (strncmp(task, "FG", 2) == 0) {
			f = optimGradient(n, x, g);
			if (!isfinite(f)) {
				cerr << "L-BFGS-B needs finite values of 'fn'" << endl;
				exit(1);
			}
			
		} else if (strncmp(task, "NEW_X", 5) == 0) {
			iter++;
			if(trace == 1 && (iter % nREPORT == 0)) {
				cout << "iter " << iter << " value " << f << endl;
			}
			if (iter > maxit) {
				*fail = 1;
				break;
			}
		} else if (strncmp(task, "WARN", 4) == 0) {
			*fail = 51;
			break;
		} else if (strncmp(task, "CONV", 4) == 0) {
			break;
		} else if (strncmp(task, "ERROR", 5) == 0) {
			*fail = 52;
			break;
		} else { /* some other condition that is not supposed to happen */
			*fail = 52;
			break;
		}
	}
	*Fmin = f;
	*fncount = *grcount = isave[33];
	if (trace) {
		cout << "final value " << *Fmin << endl;
		if (iter < maxit && *fail == 0)
			cout << "converged" << endl;
		else
			cout << "stopped after " << iter << " iterations\n";
	}
	strcpy(msg, task);
	free(g);
	free(wa);
	free(iwa);
}

double Optimization::optimFunc(int nvar, double *vars) { 
    return targetFunk(vars-1);
}


double Optimization::optimGradient(int nvar, double *x, double *dfx) {
    return derivativeFunk(x-1, dfx-1);
//    const double ERRORX = 1e-5;
//	double fx = optimFunc(nvar, x);
//	double h, temp;
//	for (int dim = 0; dim <= nvar; dim++ ){
//		temp = x[dim];
//		h = ERRORX * fabs(temp);
//		if (h == 0.0) h = ERRORX;
//		x[dim] = temp + h;
//		h = x[dim] - temp;
//		dfx[dim] = (optimFunc(nvar, x) - fx) / h;
//		x[dim] = temp;
//	}
//	return fx;
}
