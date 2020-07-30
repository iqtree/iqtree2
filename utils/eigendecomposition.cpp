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
#include <cmath>
#include <string.h>
#include <iostream>
#include <stdlib.h>
#include "tools.h"

using namespace std;

EigenDecomposition::EigenDecomposition()
{
	total_num_subst = 1.0;
	normalize_matrix = true;
    ignore_state_freq = false;
}

void EigenDecomposition::eigensystem(
	double **rate_params, double *state_freq, 
	double *eval, double **evec, double **inv_evec, int num_state) 
{
	double *forg = new double[num_state];
	double *evali = new double[num_state];
	double *new_forg = new double[num_state];
	double *eval_new = new double[num_state];
	double **a = (double**)new double[num_state];
	double **b = (double**)new double[num_state];
	double **evec_new = (double**)new double[num_state];
	int *ordr = new int[num_state + 1];
	int i, j, k, error, new_num, inew, jnew;
	double zero;


	for (i=0; i < num_state; i++)
		a[i] = new double[num_state];
	for (i=0; i < num_state; i++)
		b[i] = new double[num_state];
	for (i=0; i < num_state; i++)
		evec_new[i] = new double[num_state];

	/* get relative transition matrix and frequencies */
	memcpy(forg, state_freq, num_state * sizeof(double));
    // BQM 2015-09-07: normalize state frequencies to 1
    double sum = 0.0;
	for (i = 0; i < num_state; i++) {
		sum += forg[i];
	}
    sum = 1.0/sum;
	for (i = 0; i < num_state; i++) {
		forg[i] *= sum;
	}
	for (i = 0; i < num_state; i++) {
		memcpy(a[i], rate_params[i], num_state * sizeof(double));
	}
	//rtfdata(a, forg, num_state); 
	//    write (a, forg);

	computeRateMatrix(a, forg, num_state); /* make 1 PAM rate matrix */

	/* copy a to b */
	for (i = 0; i < num_state; i++) {
		for (j = 0; j < num_state; j++) {
			b[i][j] = a[i][j];
		}
	}

	eliminateZero(b, forg, num_state, a, new_forg, new_num);

	elmhes(a, ordr, new_num); /* compute eigenvalues and eigenvectors */
	//    writeInt (ordr);

	eltran(a, evec_new, ordr, new_num);

	//  writeMat (evec);

	hqr2(new_num, 1, new_num, a, evec_new, eval_new, evali);


	// now get back eigen
	//for (i = 0,inew = 0; i < num_state; i++)
	for (i = num_state-1,inew = new_num-1; i >= 0; i--)
		eval[i] = (forg[i] > ZERO_FREQ) ? eval_new[inew--] : 0;

	// calculate the actual eigenvectors of Q and its inverse matrix
	//for (i = 0, inew = 0; i < num_state; i++)
	for (i = num_state-1,inew = new_num-1; i >= 0; i--)
		if (forg[i] > ZERO_FREQ) {
			for (j = num_state-1, jnew = new_num-1; j >= 0; j--)
				if (forg[j] > ZERO_FREQ) {
					evec[i][j] = evec_new[inew][jnew];
					jnew--;
				} else {
					evec[i][j] = (i == j);
				}
 			inew--;
		} else 
		for (j=0; j < num_state; j++) {
			evec[i][j] = (i==j);
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
		cout << "\nWARNING: Eigensystem doesn't satisfy eigenvalue equation!\n";
		cout << "Rate matrix R: " << endl;
		for (i = 0; i < num_state; i++) {
			for (j = 0; j < num_state; j++) cout << rate_params[i][j] << " ";
			cout << endl;
		}
		cout << "State frequencies: " << endl;
		for (i = 0; i < num_state; i++) cout << state_freq[i] << " ";
		cout << endl;
	}

	for (i=num_state-1; i>= 0; i--)
		delete [] evec_new[i];
	for (i=num_state-1; i>= 0; i--)
		delete [] b[i];
	for (i=num_state-1; i>= 0; i--)
		delete [] a[i];
	delete [] ordr;
	delete [] evec_new;
	delete [] b;
	delete [] a;
	delete [] eval_new;
	delete [] new_forg;
	delete [] evali;
	delete [] forg;
	
	luinverse(evec, inv_evec, num_state); /* inverse eigenvectors are in Ievc */
	//checkevector(evec, inv_evec, num_state); /* check whether inversion was OK */

} /* eigensystem */


void EigenDecomposition::eigensystem_sym(double **rate_params, double *state_freq, 
	double *eval, double *evec, double *inv_evec, int num_state)
{
	double *forg = new double[num_state];
	double *new_forg = new double[num_state];
	double *forg_sqrt = new double[num_state];
	double *off_diag = new double[num_state];
	double *eval_new = new double[num_state];
	double **a = (double**)new double[num_state];
	double **b = (double**)new double[num_state];
	int i, j, k, new_num, inew, jnew;
	double error = 0.0;
	double zero;

	for (i = 0; i < num_state; i++) {
		a[i] = new double[num_state];
	}
	for (i = 0; i < num_state; i++) {
		b[i] = new double[num_state];
	}
	/* get relative transition matrix and frequencies */
	memcpy(forg, state_freq, num_state * sizeof(double));
    
    // BQM 2015-09-07: normalize state frequencies to 1
    double sum = 0.0;
	for (i = 0; i < num_state; i++) {
		sum += forg[i];
	}
    sum = 1.0/sum;
	for (i = 0; i < num_state; i++) {
		forg[i] *= sum;
	}
    
	for (i = 0; i < num_state; i++) {
		memcpy(a[i], rate_params[i], num_state * sizeof(double));
	}
	//rtfdata(a, forg, num_state); 
	//    write (a, forg);

	computeRateMatrix(a, forg, num_state); /* make 1 PAM rate matrix */

	/* copy a to b */
	for (i = 0; i < num_state; i++) {
		for (j = 0; j < num_state; j++) {
			b[i][j] = a[i][j];
		}
	}

	eliminateZero(b, forg, num_state, a, new_forg, new_num);

	symmetrizeRateMatrix(a, new_forg, forg_sqrt, new_num); 

	// make this matrix tridiagonal
	tred2(a, new_num, eval_new, off_diag);
	// compute eigenvalues and eigenvectors
	tqli(eval_new, off_diag, new_num, a);
    
    // make sure that all eval are non-positive
    
    for (i = 0; i < new_num; i++)
        ASSERT(eval_new[i] <= 0.01);

	// now get back eigen
	//for (i = 0,inew = 0; i < num_state; i++)
	for (i = num_state-1,inew = new_num-1; i >= 0; i--)
		eval[i] = (forg[i] > ZERO_FREQ) ? eval_new[inew--] : 0;

	ASSERT(inew == -1);
	// calculate the actual eigenvectors of Q and its inverse matrix
	//for (i = 0, inew = 0; i < num_state; i++)
	for (i = num_state-1,inew = new_num-1; i >= 0; i--)
		if (forg[i] > ZERO_FREQ) {
			for (j = num_state-1, jnew = new_num-1; j >= 0; j--)
				if (forg[j] > ZERO_FREQ) {
					evec[i*num_state+j] = a[inew][jnew] / forg_sqrt[inew];
					inv_evec[i*num_state+j] = a[jnew][inew] * forg_sqrt[jnew];
					jnew--;
				} else {
					evec[i*num_state+j] = (i == j);
					inv_evec[i*num_state+j] = (i == j);
				}
 			inew--;
		} else 
		for (j=0; j < num_state; j++) {
			evec[i*num_state+j] = (i==j);
			inv_evec[i*num_state+j] = (i==j);
		}

    // Only print eigenvalues and eigenvectors if state space is manageable.
	if ((verbose_mode >= VB_MAX) && num_state < 30) {
		cout << "eigenvalues:";
		for (i = 0; i < num_state; i++)
			cout << "\t" << eval[i];
		cout << endl;
		cout << "eigenvectors:" << endl;
		for (i = 0; i < num_state; i++) {
			for (j = 0; j < num_state; j++) {
				cout << "\t" << evec[i*num_state+j];
			}
			cout << endl;
		}
		cout << "inv_eigenvectors:" << endl;
		for (i = 0; i < num_state; i++) {
			for (j = 0; j < num_state; j++) {
				cout << "\t" << inv_evec[i*num_state+j];
			}
			cout << endl;
		}
		cout << endl;
	}

	/* check eigenvalue equation */
	error = 0.0;
	for (j = 0; j < num_state; j++) {
		for (i = 0; i < num_state; i++) {
			for (k = 0, zero = 0.0; k < num_state; k++) zero += b[i][k] * evec[k*num_state+j];
            zero -= eval[j] * evec[i*num_state+j];
            error = max(error, fabs(zero));
		}
	}
	if (error >= 0.1) {
		cout.precision(5);
        cout.unsetf(ios::fixed);
		cout << "\nWARNING: Eigensystem doesn't satisfy eigenvalue equation! (gap=" << error << ")" << endl;
        
		cout << " State frequencies (might be un-normalized): " << new_num << " states freq > " << ZERO_FREQ << endl;
        double sum = 0.0;
        cout.precision(7);
		for (i = 0; i < num_state; i++) {
            cout << state_freq[i] << " ";
            sum += state_freq[i];
        }
		cout << endl;
        cout << "sum = " << sum << endl;
		ASSERT(0);
	}

	for (i=num_state-1; i>= 0; i--)
		delete [] b[i];

	for (i=num_state-1; i>= 0; i--)
		delete [] a[i];

	delete [] b;
	delete [] a;
	delete [] eval_new;
	delete [] off_diag;
	delete [] forg_sqrt;
	delete [] new_forg;
	delete [] forg;
	
} // eigensystem_new


void EigenDecomposition::eigensystem_nonrev(
	double *rate_matrix, double *state_freq, double *eval, double *eval_imag,
	double *evec, double *inv_evec, int num_state)
{
	double *forg = new double[num_state];
	double *evali = new double[num_state];
	double *new_forg = new double[num_state];
	double *eval_new = new double[num_state];
	double **a = (double**)new double[num_state];
	double **b = (double**)new double[num_state];
	double **evec_new = (double**)new double[num_state];
	double **inv_evec_new = (double**)new double[num_state];
	int *ordr = new int[num_state + 1];
	int i, j, error, new_num, inew, jnew;
//	double zero;


	for (i=0; i < num_state; i++)
		a[i] = new double[num_state];
	for (i=0; i < num_state; i++)
		b[i] = new double[num_state];
	for (i=0; i < num_state; i++)
		evec_new[i] = new double[num_state];
	for (i=0; i < num_state; i++)
		inv_evec_new[i] = new double[num_state];

	/* get relative transition matrix and frequencies */
	memcpy(forg, state_freq, num_state * sizeof(double));
    // BQM 2015-09-07: normalize state frequencies to 1
    double sum = 0.0;
	for (i = 0; i < num_state; i++) {
		sum += forg[i];
	}
    sum = 1.0/sum;
	for (i = 0; i < num_state; i++) {
		forg[i] *= sum;
	}
	for (i = 0; i < num_state; i++) {
		memcpy(a[i], &rate_matrix[i * num_state], num_state * sizeof(double));
	}
	//rtfdata(a, forg, num_state);
	//    write (a, forg);

	computeRateMatrix(a, forg, num_state); /* make 1 PAM rate matrix */

	/* copy a to b */
	for (i = 0; i < num_state; i++)
		for (j = 0; j < num_state; j++)
			b[i][j] = a[i][j];

	eliminateZero(b, forg, num_state, a, new_forg, new_num);

	elmhes(a, ordr, new_num); /* compute eigenvalues and eigenvectors */
	//    writeInt (ordr);

	eltran(a, evec_new, ordr, new_num);

	//  writeMat (evec);

	hqr2(new_num, 1, new_num, a, evec_new, eval_new, evali);
    
    // check that complex eigenvalues are conjugated
	for (i = 0; i < new_num; i++) {
		if (evali[i] != 0.0) {
			ASSERT(i < new_num - 1 && evali[i + 1] != 0.0);
			i++;
		}
	}
	luinverse(evec_new, inv_evec_new, num_state); /* inverse eigenvectors are in Ievc */

	// now get back eigen
	//for (i = 0,inew = 0; i < num_state; i++)
	for (i = num_state-1,inew = new_num-1; i >= 0; i--)
        if (forg[i] > ZERO_FREQ) {
            eval[i] = eval_new[inew];
            eval_imag[i] = evali[inew];
            inew--;
        } else {
            eval[i] = 0.0;
            eval_imag[i] = 0.0;
        }
//		eval[i] = (forg[i] > ZERO) ? eval_new[inew--] : 0;
		//eval[i] = (forg[i] > ZERO) ? eval_new[inew++] : 0;

	// calculate the actual eigenvectors of Q and its inverse matrix
	//for (i = 0, inew = 0; i < num_state; i++)
	for (i = num_state-1,inew = new_num-1; i >= 0; i--)
		if (forg[i] > ZERO_FREQ) {
// 			for (j = 0, jnew = 0; j < num_state; j++)
			for (j = num_state-1, jnew = new_num-1; j >= 0; j--)
				if (forg[j] > ZERO_FREQ) {
					evec[i*num_state+j] = evec_new[inew][jnew];
					inv_evec[i*num_state+j] = inv_evec_new[inew][jnew];
					//jnew++;
					jnew--;
				} else {
					evec[i*num_state+j] = (i == j);
					inv_evec[i*num_state+j] = (i == j);
				}
// 			inew++;
 			inew--;
		} else
		for (j=0; j < num_state; j++) {
			evec[i*num_state+j] = (i==j);
			inv_evec[i*num_state+j] = (i==j);
		}

	/* check eigenvalue equation */
	error = 0;
//	for (j = 0; j < num_state; j++) {
//		for (i = 0, zero = 0.0; i < num_state; i++) {
//			for (k = 0; k < num_state; k++) zero += b[i][k] * evec[k][j];
//			zero -= eval[j] * evec[i][j];
//			if (fabs(zero) > 1.0e-5) {
//				error = 1;
//				break;
//			}
//		}
//	}
	if (error) {
		cout << "\nWARNING: Eigensystem doesn't satisfy eigenvalue equation!\n";
		cout << "Rate matrix R: " << endl;
		for (i = 0; i < num_state; i++) {
			for (j = 0; j < num_state; j++) cout << rate_matrix[i*num_state+j] << " ";
			cout << endl;
		}
		cout << "State frequencies: " << endl;
		for (i = 0; i < num_state; i++) cout << state_freq[i] << " ";
		cout << endl;
	}

	for (i=num_state-1; i>= 0; i--)
		delete [] inv_evec_new[i];
	for (i=num_state-1; i>= 0; i--)
		delete [] evec_new[i];
	for (i=num_state-1; i>= 0; i--)
		delete [] b[i];
	for (i=num_state-1; i>= 0; i--)
		delete [] a[i];
	delete [] ordr;
	delete [] inv_evec_new;
	delete [] evec_new;
	delete [] b;
	delete [] a;
	delete [] eval_new;
	delete [] new_forg;
	delete [] evali;
	delete [] forg;

	//checkevector(evec, inv_evec, num_state); /* check whether inversion was OK */

} /* eigensystem */


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
	double *m = new double[num_state];

	if (!ignore_state_freq)
	for (i = 0; i < num_state; i++) {
		for (j = 0; j < num_state; j++) {
			a[i][j] = stateFrqArr_[j]*a[i][j];
		}
	}

	for (i = 0, sum = 0.0; i < num_state; i++) {
		for (j = 0, temp = 0.0; j < num_state; j++)
            if (j != i)
			temp += a[i][j];
		m[i] = temp; /* row sum */
		sum += temp*stateFrqArr_[i]; /* exp. rate */
	}

	if (normalize_matrix) {
		delta = total_num_subst / sum; /* 0.01 subst. per unit time */

		for (i = 0; i < num_state; i++) {
			for (j = 0; j < num_state; j++) {
				if (i != j)
					a[i][j] = delta * a[i][j];
				else
					a[i][j] = delta * (-m[i]);
			}
		}
	} else {
		for (i = 0; i < num_state; i++)
			a[i][i] = -m[i];
	}
	delete [] m;
} /* onepamratematrix */

void EigenDecomposition::eliminateZero(double **mat, double *forg, int num, 
	double **new_mat, double *new_forg, int &new_num) {
	int i, j, inew, jnew;
	new_num = 0;
	for (i = 0; i < num; i++) {
		if (forg[i] > ZERO_FREQ)
			new_forg[new_num++] = forg[i];
    }
	if (new_num == num) return;
	//writeDouble(forg, num);
	//writeMat(mat, num);
	for (i = 0, inew = 0; i < num; i++)
		if (forg[i] > ZERO_FREQ) {
			for (j = 0, jnew = 0; j < num; j++) 
				if (forg[j] > ZERO_FREQ) {
					new_mat[inew][jnew] = mat[i][j];
					jnew++;
				}
			inew++;
		} else {
            for (j = 0; j < num; j++)
                mat[i][j] = 0.0;
            for (j = 0; j < num; j++)
                mat[j][i] = 0.0;
        }
	if (verbose_mode >= VB_MED) {
		cout << "new_num_states = " << new_num << endl;
        for (i = 0; i < new_num; i++) {
            cout << new_mat[i][i] << " ";
        }
        cout << endl;
    }
	//writeMat(new_mat, new_num);
	//writeDouble(new_forg, new_num);
}

void EigenDecomposition::symmetrizeRateMatrix(double **a, double *stateFrq, double *stateFrq_sqrt, int num_state) {
	int i, j;

	for (i = 0; i < num_state; i++)
		stateFrq_sqrt[i] = sqrt(stateFrq[i]);
	for (i = 0; i < num_state; i++) {
        double tmp = 1.0/stateFrq_sqrt[i];
		for (j = 0; j < i; j++) {
            a[j][i] *= stateFrq_sqrt[j]*tmp;
            a[i][j] = a[j][i];
            
//			a[i][j] *= stateFrq_sqrt[i] / stateFrq_sqrt[j];
//            a[j][i] = a[i][j];

//            a[j][i] *= stateFrq_sqrt[j] / stateFrq_sqrt[i];
//			if (fabs(a[j][i] - a[i][j]) > 1e-5) {
//                cout << a[i][j] << "  " << a[j][i];
//                assert(0);
//            }
		}
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
                if (iter++ == 100) {
					outWarning("Too many iterations in tqli");
                    break;
                }
 
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

void EigenDecomposition::elmhes(double **a, int *ordr, int n) {
	int m, j, i;
	double y, x;


	for (i = 0; i < n; i++)
		ordr[i] = 0;
	for (m = 2; m < n; m++) {
		x = 0.0;
		i = m;
		for (j = m; j <= n; j++) {
			if (fabs(a[j - 1][m - 2]) > fabs(x)) {
				x = a[j - 1][m - 2];
				i = j;
			}
		}
		ordr[m - 1] = i;      /* vector */

		if (i != m) {
			for (j = m - 2; j < n; j++) {
				y = a[i - 1][j];
				a[i - 1][j] = a[m - 1][j];
				a[m - 1][j] = y;
			}
			for (j = 0; j < n; j++) {
				y = a[j][i - 1];
				a[j][i - 1] = a[j][m - 1];
				a[j][m - 1] = y;
			}
		}
		if (x != 0.0) {
			for (i = m; i < n; i++) {
				y = a[i][m - 2];
				if (y != 0.0) {
					y /= x;
					a[i][m - 2] = y;
					for (j = m - 1; j < n; j++)
						a[i][j] -= y * a[m - 1][j];
					for (j = 0; j < n; j++)
						a[j][m - 1] += y * a[j][i];
				}
			}
		}
	}
} /* elmhes */


void EigenDecomposition::eltran(double **a, double **zz, int *ordr, int n) {
	int i, j, m;


	for (i = 0; i < n; i++) {
		for (j = i + 1; j < n; j++) {
			zz[i][j] = 0.0;
			zz[j][i] = 0.0;
		}
		zz[i][i] = 1.0;
	}
	if (n <= 2)
		return;
	for (m = n - 1; m >= 2; m--) {
		for (i = m; i < n; i++)
			zz[i][m - 1] = a[i][m - 2];
		i = ordr[m - 1];
		if (i != m) {
			for (j = m - 1; j < n; j++) {
				zz[m - 1][j] = zz[i - 1][j];
				zz[i - 1][j] = 0.0;
			}
			zz[i - 1][m - 1] = 1.0;
		}
	}
} /* eltran */


void EigenDecomposition::mcdiv(double ar, double ai, double br, double bi,
                  double *cr, double *ci) {
	double s, ars, ais, brs, bis;


	s = fabs(br) + fabs(bi);
	ars = ar / s;
	ais = ai / s;
	brs = br / s;
	bis = bi / s;
	s = brs * brs + bis * bis;
	*cr = (ars * brs + ais * bis) / s;
	*ci = (ais * brs - ars * bis) / s;
} /* mcdiv */


void EigenDecomposition::hqr2(int n, int low, int hgh, double **h,
                 double **zz, double *wr, double *wi) {
	int i, j, k, l=0, m, en, na, itn, its;
	double p=0, q=0, r=0, s=0, t, w, x=0, y, ra, sa, vi, vr, z=0, norm, tst1, tst2;
	int notlas; /* boolean */


	norm = 0.0;
	k = 1;
	/* store isolated roots and compute matrix norm */
	for (i = 0; i < n; i++) {
		for (j = k - 1; j < n; j++)
			norm += fabs(h[i][j]);
		k = i + 1;
		if (i + 1 < low || i + 1 > hgh) {
			wr[i] = h[i][i];
			wi[i] = 0.0;
		}
	}
	en = hgh;
	t = 0.0;
	itn = n * 30;
	while (en >= low) {    /* search for next eigenvalues */
		its = 0;
		na = en - 1;
		while (en >= 1) {
			/* look for single small sub-diagonal element */
			for (l = en; l > low; l--) {
				s = fabs(h[l - 2][l - 2]) + fabs(h[l - 1][l - 1]);

				if (s == 0.0)
					s = norm;
				tst1 = s;
				tst2 = tst1 + fabs(h[l - 1][l - 2]);
				if (tst2 == tst1)
					goto L100;
			}
			l = low;
		L100:
			x = h[en - 1][en - 1];    /* form shift */
			if (l == en || l == na)
				break;
			if (itn == 0) {
				/* all eigenvalues have not converged */
				cout << "\n\n\nHALT: PLEASE REPORT ERROR B TO DEVELOPERS\n\n\n";
				exit(1);
			}
			y = h[na - 1][na - 1];
			w = h[en - 1][na - 1] * h[na - 1][en - 1];
			/* form exceptional shift */
			if (its == 10 || its == 20) {
				t += x;
				for (i = low - 1; i < en; i++)
					h[i][i] -= x;
				s = fabs(h[en - 1][na - 1]) + fabs(h[na - 1][en - 3]);
				x = 0.75 * s;
				y = x;
				w = -0.4375 * s * s;
			}
			its++;
			itn--;
			/* look for two consecutive small sub-diagonal elements */
			for (m = en - 2; m >= l; m--) {
				z = h[m - 1][m - 1];
				r = x - z;
				s = y - z;
				p = (r * s - w) / h[m][m - 1] + h[m - 1][m];
				q = h[m][m] - z - r - s;
				r = h[m + 1][m];
				s = fabs(p) + fabs(q) + fabs(r);
				p /= s;
				q /= s;
				r /= s;
				if (m == l)
					break;
				tst1 = fabs(p) *
				       (fabs(h[m - 2][m - 2]) + fabs(z) + fabs(h[m][m]));
				tst2 = tst1 + fabs(h[m - 1][m - 2]) * (fabs(q) + fabs(r));
				if (tst2 == tst1)
					break;
			}
			for (i = m + 2; i <= en; i++) {
				h[i - 1][i - 3] = 0.0;
				if (i != m + 2)
					h[i - 1][i - 4] = 0.0;
			}
			for (k = m; k <= na; k++) {
				notlas = (k != na);
				if (k != m) {
					p = h[k - 1][k - 2];
					q = h[k][k - 2];
					r = 0.0;
					if (notlas)
						r = h[k + 1][k - 2];
					x = fabs(p) + fabs(q) + fabs(r);
					if (x != 0.0) {
						p /= x;
						q /= x;
						r /= x;
					}
				}
				if (x != 0.0) {
					if (p < 0.0) /* sign */
						s = - sqrt(p * p + q * q + r * r);
					else
						s = sqrt(p * p + q * q + r * r);
					if (k != m)
						h[k - 1][k - 2] = -s * x;
					else {
						if (l != m)
							h[k - 1][k - 2] = -h[k - 1][k - 2];
					}
					p += s;
					x = p / s;
					y = q / s;
					z = r / s;
					q /= p;
					r /= p;
					if (!notlas) {
						for (j = k - 1; j < n; j++) {    /* row modification */
							p = h[k - 1][j] + q * h[k][j];
							h[k - 1][j] -= p * x;
							h[k][j] -= p * y;
						}
						j = (en < (k + 3)) ? en : (k + 3); /* min */
						for (i = 0; i < j; i++) {    /* column modification */
							p = x * h[i][k - 1] + y * h[i][k];
							h[i][k - 1] -= p;
							h[i][k] -= p * q;
						}
						/* accumulate transformations */
						for (i = low - 1; i < hgh; i++) {
							p = x * zz[i][k - 1] + y * zz[i][k];
							zz[i][k - 1] -= p;
							zz[i][k] -= p * q;
						}
					} else {
						for (j = k - 1; j < n; j++) {    /* row modification */
							p = h[k - 1][j] + q * h[k][j] + r * h[k + 1][j];
							h[k - 1][j] -= p * x;
							h[k][j] -= p * y;
							h[k + 1][j] -= p * z;
						}
						j = (en < (k + 3)) ? en : (k + 3); /* min */
						for (i = 0; i < j; i++) {    /* column modification */
							p = x * h[i][k - 1] + y * h[i][k] + z * h[i][k + 1];
							h[i][k - 1] -= p;
							h[i][k] -= p * q;
							h[i][k + 1] -= p * r;
						}
						/* accumulate transformations */
						for (i = low - 1; i < hgh; i++) {
							p = x * zz[i][k - 1] + y * zz[i][k] +
							    z * zz[i][k + 1];
							zz[i][k - 1] -= p;
							zz[i][k] -= p * q;
							zz[i][k + 1] -= p * r;
						}
					}
				}
			}           /* for k */
		}               /* while infinite loop */
		if (l == en) {           /* one root found */
			h[en - 1][en - 1] = x + t;
			wr[en - 1] = h[en - 1][en - 1];
			wi[en - 1] = 0.0;
			en = na;
			continue;
		}
		y = h[na - 1][na - 1];
		w = h[en - 1][na - 1] * h[na - 1][en - 1];
		p = (y - x) / 2.0;
		q = p * p + w;
		z = sqrt(fabs(q));
		h[en - 1][en - 1] = x + t;
		x = h[en - 1][en - 1];
		h[na - 1][na - 1] = y + t;
		if (q >= 0.0) {           /* real pair */
			if (p < 0.0) /* sign */
				z = p - fabs(z);
			else
				z = p + fabs(z);
			wr[na - 1] = x + z;
			wr[en - 1] = wr[na - 1];
			if (z != 0.0)
				wr[en - 1] = x - w / z;
			wi[na - 1] = 0.0;
			wi[en - 1] = 0.0;
			x = h[en - 1][na - 1];
			s = fabs(x) + fabs(z);
			p = x / s;
			q = z / s;
			r = sqrt(p * p + q * q);
			p /= r;
			q /= r;
			for (j = na - 1; j < n; j++) {    /* row modification */
				z = h[na - 1][j];
				h[na - 1][j] = q * z + p * h[en - 1][j];
				h[en - 1][j] = q * h[en - 1][j] - p * z;
			}
			for (i = 0; i < en; i++) {    /* column modification */
				z = h[i][na - 1];
				h[i][na - 1] = q * z + p * h[i][en - 1];
				h[i][en - 1] = q * h[i][en - 1] - p * z;
			}
			/* accumulate transformations */
			for (i = low - 1; i < hgh; i++) {
				z = zz[i][na - 1];
				zz[i][na - 1] = q * z + p * zz[i][en - 1];
				zz[i][en - 1] = q * zz[i][en - 1] - p * z;
			}
		} else {           /* complex pair */
			wr[na - 1] = x + p;
			wr[en - 1] = x + p;
			wi[na - 1] = z;
			wi[en - 1] = -z;
		}
		en -= 2;
	}                   /* while en >= low */
	/* backsubstitute to find vectors of upper triangular form */
	if (norm != 0.0) {
		for (en = n; en >= 1; en--) {
			p = wr[en - 1];
			q = wi[en - 1];
			na = en - 1;
			if (q == 0.0) {/* real vector */
				m = en;
				h[en - 1][en - 1] = 1.0;
				if (na != 0) {
					for (i = en - 2; i >= 0; i--) {
						w = h[i][i] - p;
						r = 0.0;
						for (j = m - 1; j < en; j++)
							r += h[i][j] * h[j][en - 1];
						if (wi[i] < 0.0) {
							z = w;
							s = r;
						} else {
							m = i + 1;
							if (wi[i] == 0.0) {
								t = w;
								if (t == 0.0) {
									tst1 = norm;
									t = tst1;
									do {
										t = 0.01 * t;
										tst2 = norm + t;
									} while (tst2 > tst1);
								}
								h[i][en - 1] = -(r / t);
							} else {    /* solve real equations */
								x = h[i][i + 1];
								y = h[i + 1][i];
								q = (wr[i] - p) * (wr[i] - p) + wi[i] * wi[i];
								t = (x * s - z * r) / q;
								h[i][en - 1] = t;
								if (fabs(x) > fabs(z))
									h[i + 1][en - 1] = (-r - w * t) / x;
								else
									h[i + 1][en - 1] = (-s - y * t) / z;
							}
							/* overflow control */
							t = fabs(h[i][en - 1]);
							if (t != 0.0) {
								tst1 = t;
								tst2 = tst1 + 1.0 / tst1;
								if (tst2 <= tst1) {
									for (j = i; j < en; j++)
										h[j][en - 1] /= t;
								}
							}
						}
					}
				}
			} else if (q > 0.0) {
				m = na;
				if (fabs(h[en - 1][na - 1]) > fabs(h[na - 1][en - 1])) {
					h[na - 1][na - 1] = q / h[en - 1][na - 1];
					h[na - 1][en - 1] = (p - h[en - 1][en - 1]) /
					                    h[en - 1][na - 1];
				} else
					mcdiv(0.0, -h[na - 1][en - 1], h[na - 1][na - 1] - p, q,
					      &h[na - 1][na - 1], &h[na - 1][en - 1]);
				h[en - 1][na - 1] = 0.0;
				h[en - 1][en - 1] = 1.0;
				if (en != 2) {
					for (i = en - 3; i >= 0; i--) {
						w = h[i][i] - p;
						ra = 0.0;
						sa = 0.0;
						for (j = m - 1; j < en; j++) {
							ra += h[i][j] * h[j][na - 1];
							sa += h[i][j] * h[j][en - 1];
						}
						if (wi[i] < 0.0) {
							z = w;
							r = ra;
							s = sa;
						} else {
							m = i + 1;
							if (wi[i] == 0.0)
								mcdiv(-ra, -sa, w, q, &h[i][na - 1],
								      &h[i][en - 1]);
							else {    /* solve complex equations */
								x = h[i][i + 1];
								y = h[i + 1][i];
								vr = (wr[i] - p) * (wr[i] - p);
								vr = vr + wi[i] * wi[i] - q * q;
								vi = (wr[i] - p) * 2.0 * q;
								if (vr == 0.0 && vi == 0.0) {
									tst1 = norm * (fabs(w) + fabs(q) + fabs(x) +
									               fabs(y) + fabs(z));
									vr = tst1;
									do {
										vr = 0.01 * vr;
										tst2 = tst1 + vr;
									} while (tst2 > tst1);
								}
								mcdiv(x * r - z * ra + q * sa,
								      x * s - z * sa - q * ra, vr, vi,
								      &h[i][na - 1], &h[i][en - 1]);
								if (fabs(x) > fabs(z) + fabs(q)) {
									h[i + 1]
									[na - 1] = (q * h[i][en - 1] -
									            w * h[i][na - 1] - ra) / x;
									h[i + 1][en - 1] = (-sa - w * h[i][en - 1] -
									                    q * h[i][na - 1]) / x;
								} else
									mcdiv(-r - y * h[i][na - 1],
									      -s - y * h[i][en - 1], z, q,
									      &h[i + 1][na - 1], &h[i + 1][en - 1]);
							}
							/* overflow control */
							t = (fabs(h[i][na - 1]) > fabs(h[i][en - 1])) ?
							    fabs(h[i][na - 1]) : fabs(h[i][en - 1]);
							if (t != 0.0) {
								tst1 = t;
								tst2 = tst1 + 1.0 / tst1;
								if (tst2 <= tst1) {
									for (j = i; j < en; j++) {
										h[j][na - 1] /= t;
										h[j][en - 1] /= t;
									}
								}
							}
						}
					}
				}
			}
		}
		/* end back substitution. vectors of isolated roots */
		for (i = 0; i < n; i++) {
			if (i + 1 < low || i + 1 > hgh) {
				for (j = i; j < n; j++)
					zz[i][j] = h[i][j];
			}
		}
		/* multiply by transformation matrix to give vectors of
		 * original full matrix. */
		for (j = n - 1; j >= low - 1; j--) {
			m = ((j + 1) < hgh) ? (j + 1) : hgh; /* min */
			for (i = low - 1; i < hgh; i++) {
				z = 0.0;
				for (k = low - 1; k < m; k++)
					z += zz[i][k] * h[k][j];
				zz[i][j] = z;
			}
		}
	}
	return;
} /* hqr2 */

void EigenDecomposition::luinverse(double **inmat, double **imtrx, int size) {
	double eps = 1.0e-20; /* ! */
	int i, j, k, l, maxi=0, idx, ix, jx;
	double sum, tmp, maxb, aw;
	int *index = new int[size];
	double *wk;
	double **omtrx = (double**) new double[size];

	for (i = 0; i < size; i++)
		omtrx[i] = new double[size];

	/* copy inmat to omtrx */
	for (i = 0; i < size; i++)
		for (j = 0; j < size; j++)
			omtrx[i][j] = inmat[i][j];

	wk = (double *) calloc((size_t)size, sizeof(double));
	aw = 1.0;
	for (i = 0; i < size; i++) {
		maxb = 0.0;
		for (j = 0; j < size; j++) {
			if (fabs(omtrx[i][j]) > maxb)
				maxb = fabs(omtrx[i][j]);
		}
		if (maxb == 0.0) {
			/* Singular matrix */
			ASSERT(0 && "\n\n\nHALT: PLEASE REPORT ERROR C TO DEVELOPERS\n\n\n");
		}
		wk[i] = 1.0 / maxb;
	}
	for (j = 0; j < size; j++) {
		for (i = 0; i < j; i++) {
			sum = omtrx[i][j];
			for (k = 0; k < i; k++)
				sum -= omtrx[i][k] * omtrx[k][j];
			omtrx[i][j] = sum;
		}
		maxb = 0.0;
		for (i = j; i < size; i++) {
			sum = omtrx[i][j];
			for (k = 0; k < j; k++)
				sum -= omtrx[i][k] * omtrx[k][j];
			omtrx[i][j] = sum;
			tmp = wk[i] * fabs(sum);
			if (tmp >= maxb) {
				maxb = tmp;
				maxi = i;
			}
		}
		if (j != maxi) {
			for (k = 0; k < size; k++) {
				tmp = omtrx[maxi][k];
				omtrx[maxi][k] = omtrx[j][k];
				omtrx[j][k] = tmp;
			}
			aw = -aw;
			wk[maxi] = wk[j];
		}
		index[j] = maxi;
		if (omtrx[j][j] == 0.0)
			omtrx[j][j] = eps;
		if (j != size - 1) {
			tmp = 1.0 / omtrx[j][j];
			for (i = j + 1; i < size; i++)
				omtrx[i][j] *= tmp;
		}
	}
	for (jx = 0; jx < size; jx++) {
		for (ix = 0; ix < size; ix++)
			wk[ix] = 0.0;
		wk[jx] = 1.0;
		l = -1;
		for (i = 0; i < size; i++) {
			idx = index[i];
			sum = wk[idx];
			wk[idx] = wk[i];
			if (l != -1) {
				for (j = l; j < i; j++)
					sum -= omtrx[i][j] * wk[j];
			} else if (sum != 0.0)
				l = i;
			wk[i] = sum;
		}
		for (i = size - 1; i >= 0; i--) {
			sum = wk[i];
			for (j = i + 1; j < size; j++)
				sum -= omtrx[i][j] * wk[j];
			wk[i] = sum / omtrx[i][i];
		}
		for (ix = 0; ix < size; ix++)
			imtrx[ix][jx] = wk[ix];
	}
	free((char *)wk);
	wk = NULL;
	for (i = size-1; i >= 0; i--)
		delete [] omtrx[i];
	delete [] omtrx;
	delete [] index;
} /* luinverse */

void EigenDecomposition::checkevector(double *evec, double *ivec, int nn) {
	int i, j, ia, ib, ic, error;
	double **matx = (double**) new double [nn];
	double sum;

	for (i = 0; i < nn; i++)
		matx[i] = new double[nn];

	/* multiply matrix of eigenvectors and its inverse */
	for (ia = 0; ia < nn; ia++) {
		for (ic = 0; ic < nn; ic++) {
			sum = 0.0;
			for (ib = 0; ib < nn; ib++) sum += evec[ia*nn+ib] * ivec[ib*nn+ic];
			matx[ia][ic] = sum;
		}
	}
	/* check whether the unitary matrix is obtained */
	error = 0;
	for (i = 0; i < nn; i++) {
		for (j = 0; j < nn; j++) {
			if (i == j) {
				if (fabs(matx[i][j] - 1.0) > 1.0e-5)
					error = 1;
			} else {
				if (fabs(matx[i][j]) > 1.0e-5)
					error = 1;
			}
		}
	}
	if (error) {
		cout << "\nWARNING: Inversion of eigenvector matrix not perfect!\n";
	}

	for (i = nn-1; i >= 0; i--)
		delete [] matx[i];
	delete [] matx;
} /* checkevector */
