/*! 
 * \file
 * This file declares math utility functions in the namespace
 * MbMath. Access these functions by using MbMath::<function>.
 *  
 * \brief Declaration of math utility functions
 *
 * MrBayes version 4.0 beta
 *
 * (c) Copyright 2005.
 * \version 4.0 Beta
 * \date Last modified: $Date: 2006/09/11 21:05:15 $
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
 * $Id: MbMath.h,v 1.10 2006/09/11 21:05:15 ronquist Exp $
 */

#ifndef MbMath_H
#define MbMath_H

#include "MbMatrix.h"
#include "MbVector.h"
#include <cmath>

using namespace std;


/*!
 * This namespace contains math utility functions. Call these functions
 * by using MbMath::<function>. An alternative is to declare 'using
 * namespace MbMath;', after which functions can be accessed without the
 * MbMath:: prefix. Some functions return 0 on success and 1 on failure;
 * it is up to the calling function to handle the error.
 *
 * \brief Math utility functions declared in namespace MbMath
 * \todo Remove divideMatrixByPowerOfTwo and divideMatrixByPower of two
 * and use overloaded operators in MbMatrix instead (Someone)
 */
namespace MbMath {

      void   backSubstitutionRow(MbMatrix<double> &u, MbVector<double> &b);                               //!< back substitution of row
      void   computeLandU(MbMatrix<double> &aMat, MbMatrix<double> &lMat, MbMatrix<double> &uMat);        //!< LU decomposition
       int   expMatrixPade(MbMatrix<double> &a, MbMatrix<double> &f, int q);                              //!< exponentiate matrix using Pade approximation
    double   factorial(int x);                                                                            //!< return x! (x factorial)
       int   findPadeQValue(const double tolerance);                                                      //!< set p and q of the Pade method to achieve desired tolerance
      void   forwardSubstitutionRow(MbMatrix<double> &L, MbVector<double> &b);                            //!< forward substitution of row
      void   gaussianElimination (MbMatrix<double> &a, MbMatrix<double> &bMat, MbMatrix<double> &xMat);   //!< gaussian elimination
    double   hypotenuse(double a, double b);                                                              //!< return hypotenuse of triangle with legs a and b
    double   lnFactorial(int x);                                                                          //!< return ln(x!)
	double   lnGamma(double alp);                                                                         //!< return lnGamma(alp)
       int   transposeMatrix(const MbMatrix<double> &a, MbMatrix<double> &t);                             //!< transpose matrix

}

#endif

