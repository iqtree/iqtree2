/*!
 * \file
 * This file contains the declaration of MbEigensystem
 *
 * \brief Declaration of MbEigensystem
 *
 * MrBayes version 4.0 beta
 *
 * (c) Copyright 2005.
 * \version 4.0 Beta
 * \date Last modified: $Date: 2006/06/04 16:25:45 $
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
 * $Id: MbEigensystem.h,v 1.1 2006/06/04 16:25:45 ronquist Exp $
 */

#ifndef MbEigensystem_H
#define MbEigensystem_H

#include "MbMatrix.h"
#include "MbVector.h"
#include <complex>

using namespace std;
typedef complex<double> mbcompl;

/*!
 * This class constructs and stores the eigenvalue decomposition of
 * a square real (i.e., non-complex) matrix. The inverse of the
 * eigenvectors are also calculated and stored.
 * 
 * \see
 * -# Peters, G., and J.H. Wilkinson. 1970. Eigenvectors of real
 *       and complex matrices by LR and QR triangularisations.
 *       Numer. Math. 16:184-204.
 * -# Martin, R.S., and J.H. Wilkinson. 1968. Similarity reduction
 *       of a general matrix to Hessenberg form. Numer. Math.
 *       12:349-368.
 * -# Parlett, B.N., and C. Reinsch. 1969. Balancing a matrix for
 *       calculation of eigenvalues and eigenvectors. Numer.
 *       Math. 13:292-304.
 * \brief Construct, store and calculate eigenvalue decomposition
 * \todo See if some functions can be replaced with or moved to MbMath:: functions (Someone)
 * \todo Better documentation of algorithms
 * \todo Check if complex functions work
 */

class MbEigensystem {

public:
                        MbEigensystem(const MbMatrix<double> &m);                                                    //!< construct the eigenvalue decomposition
                        ~MbEigensystem(void);                                                                        //!< destructor
               double   getDeterminant (void);                                                                       //!< return determinant
     MbMatrix<double>   &getEigenvectors(void) { return eigenvectors; }                                              //!< return the eigenvector matrix
     MbMatrix<double>   &getInverseEigenvectors(void) { return inverseEigenvectors; }                                //!< return the inverse eigenvector matrix
	 MbVector<double>   &getRealEigenvalues(void) { return realEigenvalues; }                                        //!< return the real parts of the eigenvalues
	 MbVector<double>   &getImagEigenvalues(void) { return imaginaryEigenvalues; }                                   //!< return the imaginary parts of the eigenvalues
    MbMatrix<mbcompl>   &getComplexEigenvectors(void) { return complexEigenvectors; }                                //!< return the eigenvector matrix
    MbMatrix<mbcompl>   &getComplexInverseEigenvectors() { return complexInverseEigenvectors; }                      //!< return the inverse eigenvector matrix
                 bool   getIsComplex(void) { return isComplex; }                                                     //!< returns 'true' if there are complex eigenvalues
                  int   update(const MbMatrix<double> &m);                                                           //!< update the eigensystem for matrix m

private:
                  int   n;                                                                                           //!< row and column dimension (square matrix)
     MbMatrix<double>   eigenvectors;                                                                                //!< matrix for internal storage of eigenvectors
     MbMatrix<double>   inverseEigenvectors;                                                                         //!< matrix for internal storage of the inverse eigenvectors
    MbMatrix<mbcompl>   complexEigenvectors;                                                                         //!< matrix for internal storage of complex eigenvectors
    MbMatrix<mbcompl>   complexInverseEigenvectors;                                                                  //!< matrix for internal storage of the inverse of the complex eigenvectors
     MbVector<double>   realEigenvalues;                                                                             //!< vector for internal storage of the eigenvalues (real part)
     MbVector<double>   imaginaryEigenvalues;                                                                        //!< vector for internal storage of the eigenvalues (imaginary part)
                 bool   isComplex;                                                                                   //!< flag whether there are complex eigenvalues
                 void   allocateComplexEigenvectors(void);                                                           //!< allocate space for complex eigenvectors
	 			 void   balance(MbMatrix<double> &A, MbVector<double> &scale, int *low, int *high);                  //!< balances a matrix
                 void   balback(int low, int high, MbVector<double> &scale, MbMatrix<double> &eivec);                //!< reverses the balancing
                 bool   checkForComplexEigenvalues(void);                                                            //!< returns 'true' if there are complex eigenvalues
                 void   complexLUBackSubstitution(MbMatrix<mbcompl> &a, int *indx, MbVector<mbcompl> &b);            //!< back-substitutes a complex LU-decomposed matrix
                  int   complexLUDecompose(MbMatrix<mbcompl> &a, double *vv, int *indx, double *pd);                 //!< calculates the LU-decomposition of a complex matrix
                 void   elmhes(int low, int high, MbMatrix<double> &a, MbVector<int> &perm);                         //!< reduces matrix to upper Hessenberg form
                 void   elmtrans(int low, int high, MbMatrix<double> &a, MbVector<int> &perm, MbMatrix<double> &h);  //!< copies the Hessenberg matrix
                  int   hqr2(int low, int high, MbMatrix<double> &h, MbVector<double> &wr, MbVector<double> &wi, MbMatrix<double> &eivec); //!< computes eigenvalues and eigenvectors
                 void   initializeComplexEigenvectors(void);                                                         //!< sets up the complex eigenvector matrix
                  int   invertMatrix(MbMatrix<double> &a, MbMatrix<double> &aInv);                                   //!< inverts a matrix
                  int   invertComplexMatrix(MbMatrix<mbcompl> &a, MbMatrix<mbcompl> &aInv);                          //!< inverts a complex matrix
                 void   luBackSubstitution (MbMatrix<double> &a, int *indx, double *b);                              //!< back-substitutes an LU-decomposed matrix
                  int   luDecompose(MbMatrix<double> &a, double *vv, int *indx, double *pd);                         //!< calculates the LU-decomposition of a matrix

};

#endif

