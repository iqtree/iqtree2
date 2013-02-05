/*! 
 * \file
 * This file declares MbTransitionMatrix, which is used to hold information about
 * general continuous-time Markov models describing how characters change
 * on a phylogenetic tree. In particular, MbTransitionMatrix stores the transition
 * rate matrix and calculates transition probabilities and stationary frequencies
 * from that rate matrix.
 *  
 * \brief Declaration of MbTransitionMatrix
 *
 * MrBayes version 4.0 beta
 *
 * (c) Copyright 2005.
 * \version 4.0 Beta
 * \date Last modified: $Date: 2006/08/25 21:00:28 $
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
 * $Id: MbTransitionMatrix.h,v 1.9 2006/08/25 21:00:28 ronquist Exp $
 */

#ifndef MbTransitionMatrix_H
#define MbTransitionMatrix_H

#include <complex>

#include "MbEigensystem.h"
#include "MbVector.h"
#include "MbMatrix.h"

using namespace std;


/*!
 * This class sets up a continuous-time Markov model describing how discrete characters (DNA, RNA, amino acid,
 * morphological etc.) change on a phylogenetic tree. A continuous-time Markov chain is defined by a matrix describing
 * the infinitessimal rates of change from one state to another. This rate matrix, also known as the Q matrix or the
 * instantaneous rate matrix, allows one to calculate the probability of observing a change from state i to state
 * j over a branch length (evolutionary time) of v. The probabilities of change can be expressed in a matrix known
 * as the P matrix or the transition probability matrix. The instantaneous rate matrix also allows one to determine
 * the long-term behavior of the chain, that is, the probability of finding the process in state i after infinite
 * time. This is also known as the stationary frequency of i.
 *
 * The class accommodates both reversible and irreversible models. An irreversible model is created by a constructor
 * taking a rate set parameter and a data type parameter (defaulting to DNA). A reversible model is created by a
 * separate constructor, which also takes a state frequency parameter. Both constructors also take a bool indicating
 * whether the Eigensystem or the Pade approximation will be used to calculate the Q matrix from the P matrix.
 *
 * The rates in an reversible model are assumed to be in the following order (illustrated for the DNA data type):
 *
 *          to
 *           A      C      G      T
 * from A  -----  r0*f1  r1*f2  r2*f3
 *      C  r0*f0  -----  r3*f2  r4*f3
 *      G  r1*f0  r3*f1  -----  r5*f3
 *      T  r2*f0  r4*f1  r5*f2  -----
 *
 * That is, the rates are in the order of the upper diagonal: A<->C, A<->G, A<->T, C<->G, C<->T, G<->T.
 *
 * For an irreversible model, the rates are assumed to be in the order (for a DNA character):
 *
 *          to
 *          A    C    G    T
 * from A  ---  r0   r1   r2
 *      C  r3   ---  r4   r5
 *      G  r6   r7   ---  r8
 *      T  r9   r10  r11  ---
 *
 * That is, the off-diagonal rates are given one row at a time: A->C, A->G, A->T, C->A, C->G,
 * C->T, G->A, G->C, G->T, T->A, T->C, T->G.
 *
 */

class MbTransitionMatrix {

	public:
                      MbTransitionMatrix(const MbVector<double> &rate, bool useEigen = true);   //!< constructor for irreversible model
                      MbTransitionMatrix(const MbVector<double> &rate, const MbVector<double> &pi, bool useEigen = true);   //!< constructor for reversible model
	                  ~MbTransitionMatrix ();                                                   //!< destructor

			  double  getPadeTolerance(void) { return padeTolerance; }                     //!< set tolerance of Pade approximation
	MbMatrix<double>  getQ(void) { return Q.copy(); }                                      //!< return copy of Q matrix
	MbVector<double>  getStationaryFreqs(void) {return pi.copy(); }                        //!< get stationary frequencies (pi)
	            bool  getUseEigens(void) { return useEigens; }                             //!< is eigensystem used?
	            bool  isReversible(double tolerance = 1E-6);                               //!< is Q time-reversible?
	            void  restoreQ(void);                                                      //!< restore Q matrix and eigensystem
				void  setPadeTolerance(const double tol);                                  //!< set tolerance of Pade approximation
				void  setUseEigens (const bool flag=true);                                 //!< use eigensystem (true) or Pade approx (false)
    MbMatrix<double>  &tiProbs(const double v, MbMatrix<double> &P);                       //!< calculate transition probabilities (P) for length v
	             int  updateQ(const MbVector<double> &rate);                               //!< update Q matrix (and eigensystem if used)
	             int  updateQ(const MbVector<double> &rate, const MbVector<double> &pi);   //!< update Q matrix (and eigensystem if used)

	private:
	          double  *c_ijk;                                                              //!< vector of precalculated product of eigenvectors and their inverse
	 complex<double>  *cc_ijk;                                                             //!< vector of precalculated product of eigenvectors and their inverse
	 complex<double>  *ceigenvalue;                                                        //!< vector holding complex eigenvalues
	 complex<double>  *ceigValExp;                                                         //!< working space for calculating exp(-lambda*v) from complex eigensystem
	   MbEigensystem  *eigens;                                                             //!< pointer to eigensystem object
	          double  *eigenvalue;                                                         //!< vector holding eigenvalues
              double  *eigValExp;                                                          //!< working space for calculating exp(-lambda*v)
                bool  isComplex;                                                           //!< does Q have complex eigensystem?
                bool  isOldComplex;                                                        //!< did last Q (oldQ) have complex eigensystem?
                bool  isRev;                                                               //!< is Q reversible?
	             int  numStates;                                                           //!< number of states
	          double  *oldC_ijk;                                                           //!< old precalculated product of eigenvectors and their inverse
	          double  *oldEigenvalue;                                                      //!< old eigenvalues
	 complex<double>  *oldCC_ijk;                                                          //!< old precalculated product of complex eigenvectors and their inverse
	 complex<double>  *oldCEigenvalue;                                                     //!< old complex eigenvalues
	MbMatrix<double>  Q;                                                                   //!< the Q (rate) matrix
	MbMatrix<double>  oldQ;                                                                //!< old Q (rate) matrix
    MbVector<double>  pi;                                                                  //!< stationary frequencies (for irrev matrix)
	             int  padeQValue;                                                          //!< integer value used to control error in Pade approximation
	          double  padeTolerance;                                                       //!< tolerance for Pade approximation of matrix exponential
                bool  useEigens;                                                           //!< use eigensystem (true) or Pade approximation (false)
	            void  allocateComplexEigens(void);                                         //!< allocate space for complex eigensystem calculations
	            void  allocateEigens(void);                                                //!< allocate space for eigensystem calculations
	            void  calcCijk(void);                                                      //!< precalculations for matrix exponentiation using eigensystem
	            void  calcComplexCijk(void);                                               //!< precalculations for matrix exponentiation using complex eigensystem
	            void  calcStationaryFreq(void);                                            //!< calculate the stationary probabilites
	            void  freeComplexEigens(void);                                             //!< free space for complex eigensystem calculations
	            void  freeEigens(void);                                                    //!< free space for eigensystem calculations
	            void  initializeEigenVariables(void);                                      //!< initialize local variables for eigensystem calculations
	            void  rescaleQ(void);                                                      //!< rescale Q matrix (using pi)
	            void  tiProbsComplexEigens(const double v, MbMatrix<double> &P);           //!< calculates transition probabilities using complex eigensystem
	            void  tiProbsEigens(const double v, MbMatrix<double> &P);                  //!< calculates transition probabilities using eigensystem
	            void  tiProbsPade(const double v, MbMatrix<double> &P);                    //!< calculates transition probabilities using Pade approximation

};

#endif
