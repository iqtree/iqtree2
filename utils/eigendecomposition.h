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
#ifndef EIGENDECOMPOSITION_H
#define EIGENDECOMPOSITION_H

//const double ZERO_FREQ = 0.000001;
const double ZERO_FREQ = 1e-10;


/**
Eigenvalues, eigenvectors decomposition

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class EigenDecomposition{
public:
    EigenDecomposition();

    ~EigenDecomposition();

	/**
		EigenSystem for symmetric matrix
		@param rate_params rate parameters (not the rate matrix)
		@param state_freq state frequencies
		@param eval (OUT) eigenvalues
		@param evec (OUT) eigenvectors
		@param inv_evec (OUT) inverse matrix of eigenvectors
		@param num_state (IN) number of states
	*/
	void eigensystem_sym(double **rate_params, double *state_freq, 
	double *eval, double *evec, double *inv_evec, int num_state);

	/**
		EigenSystem for general non-symmetric matrix
		@param rate_params rate parameters (not the rate matrix)
		@param state_freq state frequencies
		@param eval (OUT) eigenvalues
		@param evec (OUT) eigenvectors
		@param inv_evec (OUT) inverse matrix of eigenvectors
		@param num_state (IN) number of states
	*/
	void eigensystem(double **rate_params, double *state_freq, 
	double *eval, double **evec, double **inv_evec, int num_state);

	/**
		EigenSystem for general non-symmetric matrix without state frequencies
		@param rate_matrix rate matrix
		@param eval (OUT) real part of eigenvalues
		@param eval_imag (OUT) imaginary part of eigenvalues
		@param evec (OUT) eigenvectors
		@param inv_evec (OUT) inverse matrix of eigenvectors
		@param num_state (IN) number of states
	*/
    void eigensystem_nonrev(double *rate_matrix, double *state_freq, double *eval, double *eval_imag,
    		double *evec, double *inv_evec, int num_state);


	/** TRUE to normalize rate matrix to 1.0 subst per unit time */
	bool normalize_matrix;

	/**
		the total number of substitutions per unit time
	*/
	double total_num_subst;
	
    /** TRUE to ignore state_freq in computation, default: FALSE */
    bool ignore_state_freq;


protected:

	/**
		compute the rate matrix and then normalize it such that the total number of substitutions is 1.
		@param rate_matrix (IN/OUT) As input, it contains rate parameters. On output it is filled with rate matrix entries
		@param state_freq state frequencies
		@param num_state number of states
	*/
	virtual void computeRateMatrix(double **rate_matrix, double *state_freq, int num_state);

	/**
		Eliminate zero entries in the rate matrix. 
		Return the new non-zero matrix with possibly reduced dimension.
		@param mat input rate matrix
		@param forg state frequencies
		@param num number of states
		@param new_mat (OUT) the new rate matrix
		@param new_forg (OUT) new state frequencies
		@param new_num (OUT) new number of states
	*/
	void eliminateZero(double **mat, double *forg, int num, 
		double **new_mat, double *new_forg, int &new_num);

/*********************************************************
* aided function for symmetric matrix
*********************************************************/

	/**
		transform the rate matrix into symmetric form, used for subsequent eigen decomposition
		@param a (IN/OUT) rate matrix
		@param stateFrq state frequencies
		@param stateFrq_sqrt square root of state frequencies
		@param num_state number of states
	*/
	void symmetrizeRateMatrix(double **a, double *stateFrq, double *stateFrq_sqrt, int num_state);


	/**
		Householder transformation of symmetric matrix A
		to tridiagonal form 
		@param a the input matrix, must be symmetric. On output,
			a is replaced by the orthogonal matrix effecting the transformation
		@param  n the size of matrix a
		@param d [0..n-1] returned the diagonal elements of the tridiagonal matrix
		@param e [0..n-1] returned the off-diagonal elements with e[0]=0
	*/
	void tred2(double **a, int n, double *d, double *e);

	/**
		QL algorithm with implicit shifts to determine eigenvalues and
		eigenvectors of a real tridiagonal symmetric matrix.
		@param d [0..n-1] diagonal elements of the tridiagonal matrix. 
			On output d return the eigenvalues.
		@param e [0..n-1] off-diagonal elements of the tridiagonal matrix, e[0] arbitrary.
			On output e is destroyed.
		@param n matrix size
		@param z must be input as the matrix returned by tred2
			z[k] return the normalized eigenvector corresponding to d[k]
	*/
	void tqli(double *d, double *e, int n, double **z);

/*********************************************************
* aided function for non-symmetric matrix
*********************************************************/

	/**
		convert a non-symmetric matrix into Hessenberg form with zeros everywhere
		below the diagonal except for the first sub-diagonal row
		@param a (IN-OUT) the matrix
		@param ordr (OUT) the order of columns
		@param n (IN) size of matrix 
	*/
	void elmhes(double **a, int *ordr, int n);

	/*
		something here
	*/
	void eltran(double **a, double **zz, int *ordr, int n);

	/*
		something here
	*/
	void mcdiv(double ar, double ai, double br, double bi,
	           double *cr, double *ci);

	/**
		QR algorithm for non-symmetric matrix to calculate eigenvectors and eigenvalues
		of a Hessenberg matrix (should be preceded by elmhes function)
		@param n (IN) size of matrix 
	*/
	void hqr2(int n, int low, int hgh, double **h, double **zz, double *wr, double *wi);

	/**
		compute the inverse of a square matrix
		@param inmat (IN) the matrix
		@param imtrx (OUT) the inverse of the input matrix
		@param size the size of matrix
	*/
	void luinverse(double **inmat, double **imtrx, int size);

	void checkevector(double *evec, double *ivec, int nn);

};

#endif
