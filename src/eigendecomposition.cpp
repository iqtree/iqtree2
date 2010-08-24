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
#include "eigendecomposition.h"
#include "optimization.h"
#include <math.h>
#include <string.h>
#include <iostream>

const double ZERO = 0.000001;
using namespace std;

EigenDecomposition::EigenDecomposition()
{
	total_num_subst = 1.0;
}

void EigenDecomposition::eigensystem_sym(double **rate_params, double *state_freq, 
	double *eval, double **evec, double **inv_evec, int num_state) 
{
	double forg[num_state], new_forg[num_state], forg_sqrt[num_state], off_diag[num_state], eval_new[num_state];
	double *a[num_state];
	double *b[num_state];
	int i, j, k, error, new_num, inew, jnew;
	double zero;

	for (i=0; i < num_state; i++)
		a[i] = new double[num_state];
	for (i=0; i < num_state; i++)
		b[i] = new double[num_state];

	/* get relative transition matrix and frequencies */
	memcpy(forg, state_freq, num_state * sizeof(double));
	for (i = 0; i < num_state; i++)
		memcpy(a[i], rate_params[i], num_state * sizeof(double)); 

	//rtfdata(a, forg, num_state); 
	//    write (a, forg);

	computeRateMatrix(a, forg, num_state); /* make 1 PAM rate matrix */

	/* copy a to b */
	for (i = 0; i < num_state; i++)
		for (j = 0; j < num_state; j++)
			b[i][j] = a[i][j];

	eliminateZero(b, forg, num_state, a, new_forg, new_num);

	transformRateMatrix(a, new_forg, forg_sqrt, new_num); 

	// make this matrix tridiagonal
	tred2(a, new_num, eval_new, off_diag);
	// compute eigenvalues and eigenvectors
	tqli(eval_new, off_diag, new_num, a);

	// now get back eigen
	//for (i = 0,inew = 0; i < num_state; i++)
	for (i = num_state-1,inew = new_num-1; i >= 0; i--)
		eval[i] = (forg[i] > ZERO) ? eval_new[inew--] : 0;
		//eval[i] = (forg[i] > ZERO) ? eval_new[inew++] : 0;

	// calculate the actual eigenvectors of Q and its inverse matrix
	//for (i = 0, inew = 0; i < num_state; i++)
	for (i = num_state-1,inew = new_num-1; i >= 0; i--)
		if (forg[i] > ZERO) {
// 			for (j = 0, jnew = 0; j < num_state; j++) 
			for (j = num_state-1, jnew = new_num-1; j >= 0; j--) 
				if (forg[j] > ZERO) {
					evec[i][j] = a[inew][jnew] / forg_sqrt[inew];
					inv_evec[i][j] = a[jnew][inew] * forg_sqrt[jnew];
					//jnew++;
					jnew--;
				} else {
					evec[i][j] = (i == j);
					inv_evec[i][j] = (i == j);
				}
// 			inew++;
 			inew--;
		} else 
		for (j=0; j < num_state; j++) {
			evec[i][j] = (i==j);
			inv_evec[i][j] = (i==j);
		}



	/* check eigenvalue equation */
	error = 0;
	for (j = 0; j < num_state; j++) {
		for (i = 0, zero = 0.0; i < num_state; i++) {
			for (k = 0; k < num_state; k++) zero += b[i][k] * evec[k][j];
			zero -= eval[j] * evec[i][j];
			if (fabs(zero) > 1.0e-5) {
				error = 1;
				break;
			}
		}
	}
	if (error) {
		std::cout << "\nWARNING: Eigensystem doesn't satisfy eigenvalue equation!\n";
		std::cout << "Rate matrix R: " << endl;
		for (i = 0; i < num_state; i++) {
			for (j = 0; j < num_state; j++) std::cout << rate_params[i][j] << " ";
			cout << endl;
		}
		std::cout << "State frequencies: " << endl;
		for (i = 0; i < num_state; i++) std::cout << state_freq[i] << " ";
		std::cout << endl;
	}

	for (i=num_state-1; i>= 0; i--)
		delete [] b[i];

	for (i=num_state-1; i>= 0; i--)
		delete [] a[i];

} // eigensystem_new

EigenDecomposition::~EigenDecomposition()
{
}

/* make rate matrix with 0.01 expected substitutions per unit time */
void EigenDecomposition::computeRateMatrix(double **a, double *stateFrqArr_, int num_state) {
	
/*
	if (myrate.isNsSyHeterogenous())
		return;
*/
	int i, j;
	double delta, temp, sum;
	double m[num_state];

	
	for (i = 0; i < num_state; i++) {
		for (j = 0; j < num_state; j++) {
			a[i][j] = stateFrqArr_[j]*a[i][j];
		}
	}

	for (i = 0, sum = 0.0; i < num_state; i++) {
		for (j = 0, temp = 0.0; j < num_state; j++)
			temp += a[i][j];
		m[i] = temp; /* row sum */
		sum += temp*stateFrqArr_[i]; /* exp. rate */
	}

	delta = total_num_subst / sum; /* 0.01 subst. per unit time */

	for (i = 0; i < num_state; i++) {
		for (j = 0; j < num_state; j++) {
			if (i != j)
				a[i][j] = delta * a[i][j];
			else
				a[i][j] = delta * (-m[i]);
		}
	}
} /* onepamratematrix */

void EigenDecomposition::eliminateZero(double **mat, double *forg, int num, 
	double **new_mat, double *new_forg, int &new_num) {
	int i, j, inew, jnew;
	new_num = 0;
	for (i = 0; i < num; i++)
		if (forg[i] > ZERO) 
			new_forg[new_num++] = forg[i];
	if (new_num == num) return;
	//writeDouble(forg, num);
	//writeMat(mat, num);
	for (i = 0, inew = 0; i < num; i++)
		if (forg[i] > ZERO) {
			for (j = 0, jnew = 0; j < num; j++) 
				if (forg[j] > ZERO) {
					new_mat[inew][jnew] = mat[i][j];
					jnew++;
				}
			inew++;
		}
	//writeMat(new_mat, new_num);
	//writeDouble(new_forg, new_num);
}

void EigenDecomposition::transformRateMatrix(double **a, double *stateFrq, double *stateFrq_sqrt, int num_state) {
	int i, j;

	for (i = 0; i < num_state; i++)
		stateFrq_sqrt[i] = sqrt(stateFrq[i]);
	for (i = 0; i < num_state; i++)
		for (j = 0; j < i; j++) {
			a[i][j] *= stateFrq_sqrt[i] / stateFrq_sqrt[j];
			a[j][i] = a[i][j];
		}
}


void EigenDecomposition::tred2(double **a, int n, double *d, double *e)
{
	int l,k,j,i;
	double scale,hh,h,g,f;

	for (i=n-1;i>0;i--) {
		l=i-1;
		h=scale=0.0;
		if (l > 0) {
			for (k=0;k<=l;k++)
				scale += fabs(a[i][k]);
			if (scale == 0.0)
				e[i]=a[i][l];
			else {
				for (k=0;k<=l;k++) {
					a[i][k] /= scale;
					h += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i]=scale*g;
				h -= f*g;
				a[i][l]=f-g;
				f=0.0;
				for (j=0;j<=l;j++) {
					a[j][i]=a[i][j]/h;
					g=0.0;
					for (k=0;k<=j;k++)
						g += a[j][k]*a[i][k];
					for (k=j+1;k<=l;k++)
						g += a[k][j]*a[i][k];
					e[j]=g/h;
					f += e[j]*a[i][j];
				}
				hh=f/(h+h);
				for (j=0;j<=l;j++) {
					f=a[i][j];
					e[j]=g=e[j]-hh*f;
					for (k=0;k<=j;k++)
						a[j][k] -= (f*e[k]+g*a[i][k]);
				}
			}
		} else
			e[i]=a[i][l];
		d[i]=h;
	}
	d[0]=0.0;
	e[0]=0.0;
	/* Contents of this loop can be omitted if eigenvectors not
			wanted except for statement d[i]=a[i][i]; */
	for (i=0;i<n;i++) {
		l=i;
		if (d[i] != 0.0) {
			for (j=0;j<l;j++) {
				g=0.0;
				for (k=0;k<l;k++)
					g += a[i][k]*a[k][j];
				for (k=0;k<l;k++)
					a[k][j] -= g*a[k][i];
			}
		}
		d[i]=a[i][i];
		a[i][i]=1.0;
		for (j=0;j<l;j++) a[j][i]=a[i][j]=0.0;
	}
}

/**
	return a^2
*/
inline double sqr(double a) {
	return (a == 0.0) ? 0.0 : a*a;
}

double pythag(double a, double b)
{
	double absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+sqr(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+sqr(absa/absb)));
}

#define NRANSI
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

void EigenDecomposition::tqli(double *d, double *e, int n, double **z)
{
	int m,l,iter,i,k;
	double s,r,p,g,f,dd,c,b;

	for (i=1;i<n;i++) e[i-1]=e[i];
	e[n-1]=0.0;
	for (l=0;l<n;l++) {
		iter=0;
		do {
			for (m=l;m<n-1;m++) {
				dd=fabs(d[m])+fabs(d[m+1]);
				if ((double)(fabs(e[m])+dd) == dd) break;
			}
			if (m != l) {
				if (iter++ == 30) 
					nrerror("Too many iterations in tqli");
				g=(d[l+1]-d[l])/(2.0*e[l]);
				r=pythag(g,1.0);
				g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
				s=c=1.0;
				p=0.0;
				for (i=m-1;i>=l;i--) {
					f=s*e[i];
					b=c*e[i];
					e[i+1]=(r=pythag(f,g));
					if (r == 0.0) {
						d[i+1] -= p;
						e[m]=0.0;
						break;
					}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=(d[i]-g)*s+2.0*c*b;
					d[i+1]=g+(p=s*r);
					g=c*r-b;
					for (k=0;k<n;k++) {
						f=z[k][i+1];
						z[k][i+1]=s*z[k][i]+c*f;
						z[k][i]=c*z[k][i]-s*f;
					}
				}
				if (r == 0.0 && i >= l) continue;
				d[l] -= p;
				e[l]=g;
				e[m]=0.0;
			}
		} while (m != l);
	}
}
#undef SIGN
#undef NRANSI
