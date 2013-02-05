/*! 
 * \file
 * This file defines math utility functions in the namespace
 * MbMath. Access these functions by using MbMath::<function>.
 *  
 * \brief Definitions of math utility functions
 *
 * MrBayes version 4.0 beta
 *
 * (c) Copyright 2005.
 * \version 4.0 Beta
 * \date Last modified: $Date: 2006/09/01 21:51:54 $
 * \author John Huelsenbeck (1)
 * \author Bret Larget (2)
 * \author Paul van der Mark (3)
 * \author Fredrik Ronquist (3)
 * \author Donald Simon (4)
 * \author (authors listed in alphabetical order)
 * (1) Division of Biological Science, University of California, San Diego
 * (2) Departments of Botany and of Statistics, University of Wisconsin - Madison
 * (3) School of Computational Science, Florida State University
 * (4) Department of Mathematics/Computer Science, Duquesne University
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License (the file gpl.txt included with this
 * distribution or http: *www.gnu.org/licenses/gpl.txt) for more
 * details.
 *
 * $Id: MbMath.cpp,v 1.11 2006/09/01 21:51:54 ronquist Exp $
 */

#include <cmath>
#include <iostream>

#include "MbMath.h"
#include "MbMatrix.h"
#include "MbVector.h"
#include "MbRandom.h"

using namespace std;

/*!
 * Back-substitution of Gaussian elimination
 * 
 * \brief Back-substitution
 * \param u Matrix to back substitute
 * \param b Solution vector
 * \return Returns nothing
 */
void MbMath::backSubstitutionRow(MbMatrix<double> &u, MbVector<double> &b) {

	int n = u.dim1();
	b[n-1] /= u[n-1][n-1];
	for (int i=n-2; i>=0; i--) 
		{
		double dotProduct = 0.0;
		for (int j=i+1; j<n; j++)
			dotProduct += u[i][j] * b[j];
		b[i] = (b[i] - dotProduct) / u[i][i];
		}

}

/*!
 * This function computes the L and U decomposition of a matrix. Basically,
 * we find matrices lMat and uMat such that: lMat * uMat = aMat
 *
 * \brief Compute LU decomposition
 * \param aMat The matrix to LU decompose (destroyed)
 * \param lMat The L matrix
 * \param uMat The U matrix
 * \return Returns nothing
 */
void MbMath::computeLandU(MbMatrix<double> &aMat, MbMatrix<double> &lMat, MbMatrix<double> &uMat) {

	int n = aMat.dim1();
	for (int j=0; j<n; j++) 
		{
		for (int k=0; k<j; k++)
			for (int i=k+1; i<j; i++)
				aMat[i][j] = aMat[i][j] - aMat[i][k] * aMat[k][j];

		for (int k=0; k<j; k++)
			for (int i=j; i<n; i++)
				aMat[i][j] = aMat[i][j] - aMat[i][k] * aMat[k][j];

		for (int m=j+1; m<n; m++)
	  		aMat[m][j] /= aMat[j][j]; 
		}

	for (int row=0; row<n; row++)
		{
		for (int col=0; col<n; col++) 
			{
			if ( row <= col ) 
				{
				uMat[row][col] = aMat[row][col];
				lMat[row][col] = (row == col ? 1.0 : 0.0);
				}
			else 
				{
				lMat[row][col] = aMat[row][col];
				uMat[row][col] = 0.0;
				}
			}
		}

}

/*!
 * This function approximates the matrix exponential, f = e^a, using
 * the Pade method, which has the advantage of error control. The error
 * is controlled by setting qValue appropriately (using the function SetQValue).
 *
 * \brief Pade approximation of Matrix exponential
 * \param a [in] Input matrix
 * \param f [out] Output matrix, e^a
 * \return Returns nothing
 * \see
 * Golub, G. H., and C. F. Van Loan. 1996. Matrix Computations, Third Edition.
 *    The Johns Hopkins University Press, Baltimore, Maryland. [algorithm 11.3.1]
 * \todo See if ldexp is faster than regular matrix division by scalar
 */
int MbMath::expMatrixPade(MbMatrix<double> &A, MbMatrix<double> &F, int qValue) {

	int dim = A.dim1();
	if (dim != A.dim2())
		return (1);
	
	// create identity matrices
	MbMatrix<double> D(dim,dim,0.0);
	for (int i=0; i<dim; i++)
		D[i][i] = 1.0;
	MbMatrix<double> N(D.copy()), X(D.copy());

	// create uninitialized matrix
	MbMatrix<double> cX(dim, dim);
	
	// We assume that we have a rate matrix where rows sum to zero
	// Then the infinity-norm is twice the maximum absolute value
	// of the diagonal cells.
	double normA = 0.0;
	for (int i=0; i<dim; i++) {
		double x = fabs (A[i][i]);
		if (x > normA)
			normA = x;
	}
	normA *= 2.0;

	// Calculate 1 + floor (log2(normA))
	int y;	
	frexp(normA, &y);	// this will give us the floor(log2(normA)) part in y
	y++;

	// Get max(0,y)
	int j = 0;
	if (y > 0)
		j = y;

	// divide A by scalar 2^j
	A /= ldexp (1.0, j);
	
	double c = 1.0;
	for (int k=1; k<=qValue; k++) {
		c = c * (qValue - k + 1.0) / ((2.0 * qValue - k + 1.0) * k);

		/* X = AX */
		X = A * X;

		/* N = N + cX */
		cX = c * X;
		N = N + cX;

		/* D = D + (-1)^k*cX */
		if (k % 2 == 0)
			D = D + cX;
		else
			D = D - cX;
		}

	MbMath::gaussianElimination(D, N, F);

	for (int k=0; k<j; k++)
		F = F * F;
	
	for (int i=0; i<dim; i++)
		{
		for (j=0; j<dim; j++)
			{
			if (F[i][j] < 0.0)
				F[i][j] = fabs(F[i][j]);
			}
		}
	return (0);

}

/*!
 * This function returns the factorial of x, x!
 *
 * \brief Return x!
 * \param x The x value
 * \return The factorial x!
 */
double MbMath::factorial(int x) {
	
	double fac = 1.0;
	for (int i=1; i<=x; i++)
		fac *= i;
	return (fac);
		
}

/*!
 * Forward substitution of Gaussian elimination
 *
 * \brief Forward substitution
 * \param L [in/out] Matrix for forward substitution
 * \param b [in/out] Solution vector
 * \return Returns nothing
 */
void MbMath::forwardSubstitutionRow(MbMatrix<double> &L, MbVector<double> &b) {

	int n = L.dim1();
	b[0] = b[0] / L[0][0];
	for (int i=1; i<n; i++) 
		{
		double dotProduct = 0.0;
		for (int j=0; j<i; j++)
	      	dotProduct += L[i][j] * b[j];
		b[i] = (b[i] - dotProduct) / L[i][i];
		}
}
	
/*!
 * Gaussian elimination
 *
 * \brief Gaussian elimination
 * \param a ??
 * \param bMat ??
 * \param xMat ??
 * \return Returns nothing
 */
void MbMath::gaussianElimination (MbMatrix<double> &a, MbMatrix<double> &bMat, MbMatrix<double> &xMat) {

	int n = a.dim1();
	MbMatrix<double> lMat(n, n);
	MbMatrix<double> uMat(n, n);
	MbVector<double> bVec(n);

	computeLandU (a, lMat, uMat);

	for (int k=0; k<n; k++) 
		{
		for (int i=0; i<n; i++)
			bVec[i] = bMat[i][k];

		/* Answer of Ly = b (which is solving for y) is copied into b. */
		forwardSubstitutionRow (lMat, bVec);

		/* Answer of Ux = y (solving for x and the y was copied into b above) 
			is also copied into b. */
		backSubstitutionRow(uMat, bVec);
		for (int i=0; i<n; i++)
			xMat[i][k] = bVec[i];
		}
	
}

/*!
 * This function returns the hypotenuse of a
 * triangle with the legs being a and b
 *
 * \brief Return hypotenuse
 * \param a First leg
 * \param b Second leg
 * \return Hypotenuse
 */
double MbMath::hypotenuse(double a, double b) {
	
	double r;
	if ( fabs(a) > fabs(b) ) 
		{
		r = b / a;
		r = fabs(a) * sqrt(1+r*r);
		} 
	else if ( b != 0.0 ) 
		{
		r = a / b;
		r = fabs(b) * sqrt(1+r*r);
		} 
	else 
		{
		r = 0.0;
		}
	return r;
		
}

/*!
 * This function returns the natural logarithm
 * of the factorial of x, ln(x!)
 *
 * \brief Return ln(x!)
 * \param x The x value
 * \return The ln factorial, ln(x!)
 */
double MbMath::lnFactorial(int x) {
	
	double lnFac = 0.0;
	for (int i=1; i<=x; i++)
		lnFac += log( (double)(i) );
	return (lnFac);

}

/*!
 * Calculates the log of the gamma function. The Gamma function is equal
 * to:
 *
 *	  Gamma(alp) = {integral from 0 to infinity} t^{alp-1} e^-t dt
 *
 * The result is accurate to 10 decimal places. Stirling's formula is used
 * for the central polynomial part of the procedure.
 *
 * \brief Calculate ln of gamma function
 * \param alp Input value
 * \return lnGamma(alp)
 * \see
 * Pike, M. C. and I. D. Hill.  1966.  Algorithm 291: Logarithm of the gamma
 *    function.  Communications of the Association for Computing
 *    Machinery, 9:684.
 */
double MbMath::lnGamma(double alp) {

	double x = alp;
	double f = 0.0;
	double z;
	if ( x < 7 ) 
		{
		f = 1.0;  
		z = x-1.0;
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
 *
 * This function calculates the p and q values needed to control the error of the
 * Pade approximation for calculating the matrix exponential, P = e^{Q * v}. 
 * The error, e(p,q) is:
 *
 *    e(p,q) = 2^(3-(p+q)) * ((p!*q!) / (p+q)! * (p+q+1)!)
 *
 * Setting p = q will minimize the error for a given amount of work. This function 
 * assumes that p = q. The function takes in as a parameter the desired tolerance
 * for the accuracy of the matrix exponentiation, and returns qV = p = q, that
 * will achieve the tolerance.
 * 
 * \brief Calculate p=q needed to control error of Pade approximation
 * \param tolerance The desired tolerance
 * \return The int value giving the desired tolerance
 * \see
 * Golub, G. H., and C. F. Van Loan. 1996. Matrix Computations, Third Edition.
 *    The Johns Hopkins University Press, Baltimore, Maryland.
 */
int MbMath::findPadeQValue(const double tolerance) {

	// Here we want to calculate
	// double x = pow(2.0, 3.0 - (0 + 0)) * MbMath::factorial(0) * MbMath::factorial(0) / (MbMath::factorial(0+0) * MbMath::factorial(0+0+1));
	// that is, the expression below for qV = 0. However, we can simplify that to
	double x = 8.0;
	int qV = 0;
	while (x > tolerance) {
		qV++;
		x = pow(2.0, 3.0 - (qV + qV)) * MbMath::factorial(qV) * MbMath::factorial(qV) / (MbMath::factorial(qV+qV) * MbMath::factorial(qV+qV+1));
	}
	return (qV);

}

/*!
 * Transpose the matrix a. The matrix a should be m X n whereas the
 * matrix t should be n X m. If not, we return 1. On success, 0 is
 * returned.
 *
 * \brief Transpose a matrix
 * \param a [in] Matrix to transpose
 * \param t [out] Transposed matrix
 * \return 0 on success, 1 on failure
 */
int MbMath::transposeMatrix(const MbMatrix<double> &a, MbMatrix<double> &t) {
	
	int m = a.dim1();
	int n = a.dim2();
	
	if ( m != t.dim2() || n != t.dim1() )
		return (1);

	for (int i=0; i<m; i++)
		for (int j=0; j<n; j++)
			t[j][i] = a[i][j];
	return (0);

}


