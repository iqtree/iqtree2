/*
 *
 * lbfgsb_new.h
 * HAL_HAS
 *
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2014, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * All rights reserved. CSIRO is willing to grant you a license to HAL-HAS on the terms of the GNU General Public
 * License version 3 as published by the Free Software Foundation (http://www.gnu.org/licenses/gpl.html), except
 * where otherwise indicated for third party material.
 * The following additional terms apply under clause 7 of that license:
 * EXCEPT AS EXPRESSLY STATED IN THIS AGREEMENT AND TO THE FULL EXTENT PERMITTED BY APPLICABLE LAW, THE SOFTWARE
 * IS PROVIDED "AS-IS". CSIRO MAKES NO REPRESENTATIONS, WARRANTIES OR CONDITIONS OF ANY KIND, EXPRESS OR IMPLIED,
 * INCLUDING BUT NOT LIMITED TO ANY REPRESENTATIONS, WARRANTIES OR CONDITIONS REGARDING THE CONTENTS OR ACCURACY
 * OF THE SOFTWARE, OR OF TITLE, MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NON-INFRINGEMENT, THE ABSENCE
 * OF LATENT OR OTHER DEFECTS, OR THE PRESENCE OR ABSENCE OF ERRORS, WHETHER OR NOT DISCOVERABLE.
 * TO THE FULL EXTENT PERMITTED BY APPLICABLE LAW, IN NO EVENT SHALL CSIRO BE LIABLE ON ANY LEGAL THEORY (INCLUDING,
 * WITHOUT LIMITATION, IN AN ACTION FOR BREACH OF CONTRACT, NEGLIGENCE OR OTHERWISE) FOR ANY CLAIM, LOSS, DAMAGES
 * OR OTHER LIABILITY HOWSOEVER INCURRED.  WITHOUT LIMITING THE SCOPE OF THE PREVIOUS SENTENCE THE EXCLUSION OF
 * LIABILITY SHALL INCLUDE: LOSS OF PRODUCTION OR OPERATION TIME, LOSS, DAMAGE OR CORRUPTION OF DATA OR RECORDS;
 * OR LOSS OF ANTICIPATED SAVINGS, OPPORTUNITY, REVENUE, PROFIT OR GOODWILL, OR OTHER ECONOMIC LOSS; OR ANY SPECIAL,
 * INCIDENTAL, INDIRECT, CONSEQUENTIAL, PUNITIVE OR EXEMPLARY DAMAGES, ARISING OUT OF OR IN CONNECTION WITH THIS
 * AGREEMENT, ACCESS OF THE SOFTWARE OR ANY OTHER DEALINGS WITH THE SOFTWARE, EVEN IF CSIRO HAS BEEN ADVISED OF
 * THE POSSIBILITY OF SUCH CLAIM, LOSS, DAMAGES OR OTHER LIABILITY.
 * APPLICABLE LEGISLATION SUCH AS THE AUSTRALIAN CONSUMER LAW MAY APPLY REPRESENTATIONS, WARRANTIES, OR CONDITIONS,
 * OR IMPOSES OBLIGATIONS OR LIABILITY ON CSIRO THAT CANNOT BE EXCLUDED, RESTRICTED OR MODIFIED TO THE FULL EXTENT
 * SET OUT IN THE EXPRESS TERMS OF THIS CLAUSE ABOVE "CONSUMER GUARANTEES".  TO THE EXTENT THAT SUCH CONSUMER
 * GUARANTEES CONTINUE TO APPLY, THEN TO THE FULL EXTENT PERMITTED BY THE APPLICABLE LEGISLATION, THE LIABILITY
 * OF CSIRO UNDER THE RELEVANT CONSUMER GUARANTEE IS LIMITED (WHERE PERMITTED AT CSIRO’S OPTION) TO ONE OF FOLLOWING
 * REMEDIES OR SUBSTANTIALLY EQUIVALENT REMEDIES:
 * (a)               THE REPLACEMENT OF THE SOFTWARE, THE SUPPLY OF EQUIVALENT SOFTWARE, OR SUPPLYING RELEVANT
 *                   SERVICES AGAIN;
 * (b)               THE REPAIR OF THE SOFTWARE;
 * (c)               THE PAYMENT OF THE COST OF REPLACING THE SOFTWARE, OF ACQUIRING EQUIVALENT SOFTWARE, HAVING THE
 *                   RELEVANT SERVICES SUPPLIED AGAIN, OR HAVING THE SOFTWARE REPAIRED.
 * IN THIS CLAUSE, CSIRO INCLUDES ANY THIRD PARTY AUTHOR OR OWNER OF ANY PART OF THE SOFTWARE OR MATERIAL DISTRIBUTED
 * WITH IT.  CSIRO MAY ENFORCE ANY RIGHTS ON BEHALF OF THE RELEVANT THIRD PARTY.
 * Third Party Components
 * The following third party components are distributed with the Software.  You agree to comply with the license
 * terms for these components as part of accessing the Software.  Other third party software may also be identified
 * in separate files distributed with the Software.
 * ___________________________________________________________________
 * 
 * R : A Computer Language for Statistical Data Analysis version 3.0.1 (http://cran.r-project.org/src/base/R-3/R-3.0.1.tar.gz)
 * Copyright (C) 2000-2004 The R Core Team
 * This software is licensed under GNU GPL
 * 
 * JACOBI_EIGENVALUE.C (http://people.sc.fsu.edu/~jburkardt/c_src/jacobi_eigenvalue/jacobi_eigenvalue.c)
 * Copyright (C) 2003-2013 John Burkardt
 * This software is licensed under GNU LGPL (http://www.gnu.org/licenses/lgpl.html)
 * ___________________________________________________________________
 */


#ifndef __RAL_RAS__lbfgsb_new__
#define __RAL_RAS__lbfgsb_new__

#include <cstring>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <float.h>
//#include "gradient.h"

using namespace std;

// Function to access the L-BFGS-B function
// 1. int n : The number of the variables
// 2. double* x : initial values of the variables
// 3. double* l : lower bounds of the variables
// 4. int maxit : max # of iterations
// 5. void* ex  : the wrapped variables for objective function
// After the function is invoked, the values of x will be updated
void lbfgsb_R(int n, double* x, double* l, int maxit, void* ex);

// Function to access the L-BFGS-B function
// 1. int n : The number of the variables
// 2. double* x : initial values of the variables
// 3. double* l : lower bounds of the variables
// 4. double* u : upper bounds of the variables
// 4. int maxit : max # of iterations
// 5. void* ex  : the wrapped variables for objective function
// After the function is invoked, the values of x will be updated
void lbfgsb_R2(int n, double* x, double* l, double* u, int maxit, void* ex);



// ========================================================= //
// FUNCTIONS converted from R v 3.0.1
// ========================================================= //

typedef double optimfn(int, double *, void *);
typedef void optimgr(int, double *, double *, void *);

void lbfgsb(int n, int m, double *x, double *l, double *u, int *nbd,
		double *Fmin, optimfn fminfn, optimgr fmingr, int *fail,
		void *ex, double factr, double pgtol,
		int *fncount, int *grcount, int maxit, char *msg,
		int trace, int nREPORT);


void setulb(int n, int m, double *x, double *l, double *u, int *nbd,
		double *f, double *g, double factr, double *pgtol,
		double *wa, int * iwa, char *task, int iprint,
		int *lsave, int *isave, double *dsave);

void mainlb(int n, int m, double *x,
		double *l, double *u, int *nbd, double *f, double *g,
		double factr, double *pgtol, double *ws, double * wy,
		double *sy, double *ss, double *wt, double *wn,
		double *snd, double *z, double *r, double *d,
		double *t, double *wa, int *indx, int *iwhere,
		int *indx2, char *task, int iprint,
		char *csave, int *lsave, int *isave, double *dsave);

void errclb(int n, int m, double factr, double *l, double *u,
		int *nbd, char *task, int *info, int *k);

void prn3lb(int n, double *x, double *f, char *task, int iprint,
		int info, int iter, int nfgv, int nintol, int nskip,
		int nact, double sbgnrm, int nint,
		char *word, int iback, double stp, double xstep,
		int k);

void prn1lb(int n, int m, double *l, double *u, double *x,
		int iprint, double epsmch);

void active(int n, double *l, double *u,
		int *nbd, double *x, int *iwhere, int iprint,
		int *prjctd, int *cnstnd, int *boxed);

void projgr(int n, double *l, double *u,
		int *nbd, double *x, double *g, double *sbgnrm);

void timer(double * ttime);

void cauchy(int n, double *x, double *l, double *u, int *nbd,
		double *g, int *iorder, int * iwhere, double *t,
		double *d, double *xcp, int m,
		double *wy, double *ws, double *sy, double *wt,
		double *theta, int *col, int *head, double *p,
		double *c, double *wbp, double *v, int *nint,
		int iprint, double *sbgnrm, int *info, double * epsmch);

void freev(int n, int *nfree, int *indx,
		int *nenter, int *ileave, int *indx2, int *iwhere,
		int *wrk, int *updatd, int *cnstnd, int iprint,
		int *iter);

void formk(int n, int *nsub, int *ind, int * nenter, int *ileave,
		int *indx2, int *iupdat, int * updatd, double *wn,
		double *wn1, int m, double *ws, double *wy, double *sy,
		double *theta, int *col, int *head, int *info);

void cmprlb(int n, int m, double *x,
		double *g, double *ws, double *wy, double *sy,
		double *wt, double *z, double *r, double *wa,
		int *indx, double *theta, int *col, int *head,
		int *nfree, int *cnstnd, int *info);

void subsm(int n, int m, int *nsub, int *ind,
		double *l, double *u, int *nbd, double *x,
		double *d, double *ws, double *wy, double *theta,
		int *col, int *head, int *iword, double *wv,
		double *wn, int iprint, int *info);

void lnsrlb(int n, double *l, double *u,
		int *nbd, double *x, double *f, double *fold,
		double *gd, double *gdold, double *g, double *d,
		double *r, double *t, double *z, double *stp,
		double *dnorm, double *dtd, double *xstep,
		double *stpmx, int *iter, int *ifun, int *iback, int *nfgv,
		int *info, char *task, int *boxed, int *cnstnd,
		char *csave, int *isave, double *dsave);

void matupd(int n, int m, double *ws,
		double *wy, double *sy, double *ss, double *d,
		double *r, int *itail, int *iupdat, int *col,
		int *head, double *theta, double *rr, double *dr,
		double *stp, double *dtd);

void prn2lb(int n, double *x, double *f, double *g, int iprint,
		int iter, int nfgv, int nact, double sbgnrm,
		int nint, char *word, int iword, int iback,
		double stp, double xstep);

void pvector(char *title, double *x, int n);

void formt(int m, double *wt, double *sy, double *ss,
		int *col, double *theta, int *info);

void bmv(int m, double *sy, double *wt,
		int *col, double *v, double *p, int *info);

void hpsolb(int n, double *t, int *iorder, int iheap);

void dcsrch(double *f, double *g, double *stp,
		/*Chgd: the next five are no longer pointers:*/
		double ftol, double gtol, double xtol,
		double stpmin, double stpmax,
		char *task, int *isave, double *dsave);

void dcstep(double *stx, double *fx, double *dx,
		double *sty, double *fy, double *dy, double *stp,
		double *fp, double *dp, int *brackt, double *stpmin,
		double *stpmax);

// ========================================================= //
// Other fortan functions
// ========================================================= //

#define DY(I) dy[(I)-1]
#define DX(I) dx[(I)-1]

int dcopy(int *n, double *dx, int *incx,
		double *dy, int *incy);

int dscal(int *n, double *da, double *dx,
		int *incx);

double ddot(int *n, double *dx, int *incx, double *dy,
		int *incy);

int daxpy(int *n, double *da, double *dx,
		int *incx, double *dy, int *incy);

int dpofa(double *a, int *lda, int *n, int *info);

int dtrsl(double *t, int *ldt, int *n,
		double *b, int *job, int *info);

#endif /* defined(__RAL_RAS__lbfgsb_new__) */
