/*! 
 * \file
 * This file contains the functions in the MbTransitionMatrix class implementing
 * general continuous-time Markov models describing how discrete characters
 * change on a phylogenetic tree. 
 *  
 * \brief Function definitions for MbTransitionMatrix
 *
 * MrBayes version 4.0 beta
 *
 * (c) Copyright 2005.
 * \version 4.0 Beta
 * \date Last modified: $Date: 2006/08/25 21:00:27 $
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
 * $Id: MbTransitionMatrix.cpp,v 1.11 2006/08/25 21:00:27 ronquist Exp $
 */

#include <complex>

#include "MbTransitionMatrix.h"
#include "MbMath.h"

using namespace std;

/*!
 * This constructor builds a non-reversible rate matrix
 * for discrete characters.
 * 
 * \brief Creates non-reversible rate matrix
 * \param rate The irrev rates 
 * \param useEigen Use eigensystem (true) or Pade approximation (false)
 * 
 */
MbTransitionMatrix::MbTransitionMatrix(const MbVector<double> &rate, bool useEigen)
    : c_ijk(0), cc_ijk(0), ceigenvalue(0), ceigValExp(0), eigens(0), eigenvalue(0), eigValExp(0), 
	isComplex(false), isOldComplex(false), isRev(false), numStates(0), oldC_ijk(0), oldEigenvalue(0), 
	oldCC_ijk(0), oldCEigenvalue(0), useEigens(useEigen) {

	// Find number of rates
	int nSt = (int) (floor(sqrt((double)rate.dim()))) + 1;
	
	// If the number of rates is incorrect, return an empty transition matrix object
	if ( rate.dim() != nSt*(nSt-1) )
		return;
	
	// initialize
	numStates = nSt;

	// Allocate space for the stationary frequencies
	pi = MbVector<double>(numStates, 0.0);

	// Allocate space for variables that hold Q
	Q = MbMatrix<double>(numStates, numStates, 0.0);
	oldQ = MbMatrix<double>(numStates, numStates, 0.0);

	// Allocate space for eigensystem calculations
	if (useEigens) {
		allocateEigens();
	}

	setPadeTolerance(1E-6);
	
	// Update the Q matrix and the precalculated values
	// Do it twice to make sure restore always works
	updateQ(rate);
	updateQ(rate);

}

/*!
 * This constructor builds a time-reversible transition matrix to describe
 * evolutionary substitutions along a phylogenetic tree for discrete characters.
 * Such a model is often used to describe substitutions for DNA or amino acid data
 * but it works equally well for other characters with discrete states.
 * 
 * \brief Creates time-reversible transition matrix
 * \param rate Vector of rates
 * \param pi Vector of stationary state frequencies
 * \param useEigen Use eigensystem (true) or Pade approximation (false)
 * 
 */
MbTransitionMatrix::MbTransitionMatrix(const MbVector<double> &rate, const MbVector<double> &pi, bool useEigen)
    : c_ijk(0), cc_ijk(0), ceigenvalue(0), ceigValExp(0), eigens(0), eigenvalue(0), eigValExp(0), 
	isComplex(false), isOldComplex(false), isRev(true), numStates(0), oldC_ijk(0), oldEigenvalue(0), 
	oldCC_ijk(0), oldCEigenvalue(0), useEigens(useEigen) {

	// Check for consistency
	int nSt = pi.dim();
	if ( 2*rate.dim() != nSt*(nSt-1) )
		return;	// return empty transition matrix

	// Initialize
	numStates = nSt;

	// Store the stationary frequencies in case someone asks for them
	this->pi = pi.copy();

	// Allocate space for variables that hold Q
	Q = MbMatrix<double>(numStates, numStates, 0.0);
	oldQ = MbMatrix<double>(numStates, numStates, 0.0);

	// Allocate space for eigensystem calculations
	if (useEigens)
		allocateEigens();
	
	setPadeTolerance(1E-6);
	
	// Update the Q matrix and the precalculated values
	// Do it twice to ensure that restore always works
	updateQ(rate, pi);
	updateQ(rate, pi);
	
}

/*!
 * Destructor. Deallocates memory used for Q matrix
 * and eigensystem.
 *
 * \brief Destructor
 *
 */
MbTransitionMatrix::~MbTransitionMatrix(void) {

	if (useEigens) {
		freeEigens();
		if (isComplex || isOldComplex)
			freeComplexEigens();
	}

}

/*!
 * This function allocates space for variables that are needed in
 * calculation of transition probabilities using a complex eigensystem
 * of the Q matrix.
 *
 * \brief Allocate space for complex eigensystem
 */
void MbTransitionMatrix::allocateComplexEigens(void) {

	cc_ijk         = new complex<double>[numStates*numStates*numStates];
    oldCC_ijk      = new complex<double>[numStates*numStates*numStates];
	ceigValExp     = new complex<double>[numStates];
	ceigenvalue    = new complex<double>[numStates];
	oldCEigenvalue = new complex<double>[numStates];

}

/*!
 * This function allocates space for variables that are needed in
 * calculation of transition probabilities using the eigensystem
 * of the Q matrix. It also initializes the eigensystem calculator
 * with the current Q matrix so that we can immediately retrieve
 * eigenvalues.
 *
 * \brief Allocate space for eigensystem
 */
void MbTransitionMatrix::allocateEigens(void) {

	eigens        = new MbEigensystem(Q);
	c_ijk         = new double[numStates*numStates*numStates];
    oldC_ijk      = new double[numStates*numStates*numStates];
	eigValExp     = new double[numStates];
	eigenvalue    = new double[numStates];
	oldEigenvalue = new double[numStates];

}

/*!
 * This function precalculates the product of the eigenvectors and their
 * inverse for faster calculation of transition probabilities. The output
 * is a vector of precalculated values (c_ijk). The input is the eigenvector
 * matrix and the inverse of the eigenvector matrix. This function also
 * fetches the eigenvalues from the eigensystem and stores them in an array
 * of doubles in this class.
 *
 * \brief Calculate c_ijk and get eigenvalues from eigensystem
 */
void MbTransitionMatrix::calcCijk(void) {

	// keep a copy of the eigenvalues
	double *p = (double *)(eigens->getRealEigenvalues());
	double *q = eigenvalue;
	memcpy(q, p, numStates*sizeof(double));

	// calculate c_ijk
	MbMatrix<double> ev  = eigens->getEigenvectors();
	MbMatrix<double> iev = eigens->getInverseEigenvectors();
	double *pc = c_ijk;
	for (int i=0; i<numStates; i++)
		for (int j=0; j<numStates; j++)
			for (int k=0; k<numStates; k++)
			 	*pc++ = ev[i][k] * iev[k][j];

}

/*!
 * This function precalculates the product of the eigenvectors and their
 * inverse for faster calculation of transition probabilities when we have
 * at least one complex eigenvalue. The output is a vector of precalculated
 * complex values (cc_ijk). The input is the complex eigenvector matrix
 * and the inverse of the complex eigenvector matrix. This function also
 * fetches the real and imaginary eigenvalues from the eigensystem and stores
 * them in two arrays of double values (eigenvalues and ieigenvalues) in this
 * class.
 *
 * \brief Calculate cc_ijk and get eigenvalues from complex eigensystem
 */
void MbTransitionMatrix::calcComplexCijk(void) {

	// keep a copy of the complex eigenvalues
	double *p = (double *)(eigens->getRealEigenvalues());
	double *q = (double *)(eigens->getImagEigenvalues());
	for (int i=0; i<numStates; i++)
		ceigenvalue[i] = complex<double>(*p++, *q++);

	// calculate cc_ijk
	MbMatrix<complex<double> > cev = eigens->getComplexEigenvectors();
	MbMatrix<complex<double> > ciev = eigens->getComplexInverseEigenvectors();
	complex<double> *pc = cc_ijk;
	for (int i=0; i<numStates; i++)
		for (int j=0; j<numStates; j++)
			for (int k=0; k<numStates; k++)
			 	*(pc++) = cev[i][k] * ciev[k][j];

}

/*!
 * This function calculates the stationary frequencies of the rate matrix. The
 * rate matrix, Q, is the infinitesimal generator of the Markov chain. It is an
 * n X n matrix whose off-diagonal elements are q_ij >= 0 and whose diagonal elements
 * are specified such that each row sums to zero. The rate matrix is finite (has
 * a fixed number of states) and we assume that the input matrix is irreducible, as
 * is the usual case for substitution models. Because Q is irreducible and finite,
 * it has a stationary distribution, pi, which is a row vector of n probabilities.
 * The stationary probabilities can be calculated by solving the homogeneous system
 * of equations, pi*Q = 0, where 0 is a vector of zeros.
 *
 * We do the following to calculate the stationary frequencies. 
 *
 * 1. We perform an LU decomposition of the transpose of the matrix Q.
 *
 *    Q' = LU
 *
 * 2. Now we set Ux = z (x will eventually hold the stationary probabilities). 
 *    Because L is nonsingular, we have z = 0. We proceed to back substitute on
 *    Ux = z = 0. When u_nn = 0, we can put in any solution for x. Here, we put
 *    in x_n = 1. We then solve the other values of x through back substitution.
 *
 * 3. The solution obtained in 2 is not a probability vector. We normalize the
 *    vector such that the sum of the elements is 1.
 *
 * Note that the only time we need to use this function is when we don't
 * know the stationary frequencies of the rate matrix beforehand. For most
 * substitution models used in molecular evolution, the stationary frequencies
 * are built into the rate matrix itself. These models are time-reversible.
 * This function is useful for the non-reversible models.
 *
 * \brief Calculate stationary frequencies of non-reversible model.
 * 
 * \see
 * Stewart, W. J. 1999. Numerical methods for computing stationary distributions of 
 *    finite irreducible Markov chains. In "Advances in Computational
 *    Probability", W. Grassmann, ed. Kluwer Academic Publishers. 
 *
 */
void MbTransitionMatrix::calcStationaryFreq(void) {

	// transpose the rate matrix (qMatrix) and put into QT
	MbMatrix<double> QT(numStates, numStates);
	MbMath::transposeMatrix(Q, QT);

	// compute the LU decomposition of the transposed rate matrix
	MbMatrix<double> L(numStates, numStates);
	MbMatrix<double> U(numStates, numStates);
	MbMath::computeLandU(QT, L, U);
	
	// back substitute into z = 0 to find un-normalized stationary frequencies
	// start with x_n = 1.0
	pi[numStates-1] = 1.0;
	for (int i=numStates-2; i>=0; i--)
		{
		double dotProduct = 0.0;
		for (int j=i+1; j<numStates; j++)
			dotProduct += U[i][j] * pi[j];
		pi[i] = (0.0 - dotProduct) / U[i][i];
		}
		
	// normalize the solution vector
	double sum = 0.0;
	for (int i=0; i<numStates; i++)
		sum += pi[i];
	for (int i=0; i<numStates; i++)
		pi[i] /= sum;

}

/*!
 * This function frees up the space allocated for variables needed in
 * the calculation of transition probabilities using a complex eigensystem
 * of the Q matrix.
 *
 * \brief Free space for complex eigensystem
 */
void MbTransitionMatrix::freeComplexEigens(void) {

	delete [] cc_ijk;
    delete [] oldCC_ijk;
	delete [] ceigValExp;
	delete [] ceigenvalue;
	delete [] oldCEigenvalue;

}

/*!
 * This function frees up the space allocated for variables needed in
 * the calculation of transition probabilities using the eigensystem
 * of the Q matrix. It also deletes the eigensystem itself.
 *
 * \brief Free space for eigensystem
 */
void MbTransitionMatrix::freeEigens(void) {

	delete eigens;
	delete [] c_ijk;
	delete [] oldC_ijk;
	delete [] eigValExp;
	delete [] eigenvalue;
	delete [] oldEigenvalue;

}

/*!
 * Initialize local variables for eigensystem calculations. Perform precalcuations
 * using the eigenvector and inverse eigenvector matrix, taking into account that
 * these may be, or may have been, complex.
 *
 * \brief Initialize local eigensystem variables
 */
void MbTransitionMatrix::initializeEigenVariables(void) {
	
	if (useEigens) {
		bool wasSecondLastComplex = isComplex;
		
		// Make sure the eigensystem is up to date
		eigens->update(Q);

		// Check if any eigenvalues are complex
		isComplex = eigens->getIsComplex();

		// Precalculate the product of the eigenvectors and their inverse
		if (isComplex == false) {
			if (wasSecondLastComplex && !isOldComplex)
				freeComplexEigens();
			calcCijk();
		}
		else {
			if (!wasSecondLastComplex && !isOldComplex)
				allocateComplexEigens();
			calcComplexCijk();
		}
	}
}

/*!
 * This function checks that the rate matrix, a, is time reversible. It takes as
 * input the rate matrix, a, and the stationary frequencies of the process, f. 
 * It checks that f[i] * q[i][j] = f[j] * q[j][i] for all i != j. It does this
 * by accumulating the difference | f[i] * q[i][j] - f[j] * q[j][i] | for all
 * off-diagonal comparisons. If this difference is less than tolerance,
 * it reports that the rate matrix is time-reversible. If the flag isRev
 * is set to true, then we do not need to check because then we have determined
 * previously that the rate matrix is reversible.
 *
 * \brief Check reversibility
 * \param tolerance The tolerated deviation
 */ 
bool MbTransitionMatrix::isReversible(double tolerance) {

	if ( tolerance <= 0.0 || isRev )
		return (isRev);
	
	double diff = 0.0;
	for (int i=0; i<numStates; i++)
		for (int j=i+1; j<numStates; j++)
			diff += fabs(pi[i] * Q[i][j] - pi[j] * Q[j][i]);
	if (diff < tolerance)
		isRev = true;
	
	return isRev;

}

/*!
 * Rescale the rate matrix so that the mean rate at stationarity
 * is 1.0. Requires that stationary frequencies have been
 * calculated first.
 *
 * \brief Rescale Q given pi vector
 */
void MbTransitionMatrix::rescaleQ (void) {

	// Calculate the scaler, a factor by which all elements of Q
	// are multiplied such that the mean rate of substitution is 1
	double scaler = 0.0;
	for (int i=0; i<numStates; i++)
		scaler += pi[i] * Q[i][i];
			
	// Rescale rate matrix
	scaler = -1.0 / scaler;
	for (int i=0; i<numStates; i++)
		for (int j=0; j<numStates; j++)
			Q[i][j] *= scaler;

}
			
/*!
 * This function restores the Q matrix and, if necessary, the 
 * eigensystem and the c_ijk precalculated values.
 *
 * \brief Restore Q matrix and precalculated values
 *
 */
void MbTransitionMatrix::restoreQ(void) {

	// Quick swapping relies on the shallow inline assignment
	// in MbMatrix and MbVector classes
	MbMatrix<double> tempM = Q;
	Q = oldQ;
	oldQ = tempM;
	
	if (useEigens) {
		bool tempB = isComplex;
		isComplex = isOldComplex;
		isOldComplex = tempB;

		double *tempD = eigenvalue;
		eigenvalue = oldEigenvalue;
		oldEigenvalue = tempD;

		tempD = c_ijk;
		c_ijk = oldC_ijk;
		oldC_ijk = tempD;
		
		if ( isComplex || isOldComplex ) {
			complex<double> *tempC = ceigenvalue;
			ceigenvalue = oldCEigenvalue;
			oldCEigenvalue = tempC;

			tempC = cc_ijk;
			cc_ijk = oldCC_ijk;
			oldCC_ijk = tempC;
		}
	}

}

/*!
 * This functions sets the tolerance for the Pade approximation
 * and finds the associated q value (number of cycles) for the 
 * Pade algorithm for matrix exponentiation. The tolerance value
 * roughly gives the maximum percentage error in the exponentiated
 * matrix. For more details, see the algorithm citation for
 * expMatrixPade in MbMath.cpp.
 *
 * \brief Find q value for Pade algorithm
 * \param tol Tolerance (epsilon(p,q))
 */
void  MbTransitionMatrix::setPadeTolerance(const double tol) {
	
	if (tol > 0.0) {
		padeTolerance = tol;
		padeQValue = MbMath::findPadeQValue(tol);
	}
	else
		padeQValue = 4;

}

/*!
 * This functions sets the calculation of transition probabilities to use
 * either the eigensystem of the rate matrix or the Pade approximation.
 * Allocation or freeing of memory space is also taken care of here.
 *
 * \brief Set method for calculation of transition probabilities
 * \param flag Use eigensystem (true) or Pade approximation (false)
 */
void MbTransitionMatrix::setUseEigens (const bool flag) {

	if (useEigens && flag)
		return;
	else if (useEigens == false && flag == false)
		return;
	
	if (flag) {
		useEigens = true;
		// precalculate twice to make sure that calls to
		// restoreQ() will always give reasonable results
		// first swap back to old matrix
		restoreQ();
		// allocate and construct eigensystem for this Q
		allocateEigens();
		// initialize local eigensystem variables
		initializeEigenVariables();
		// now swap back to the current rate matrix and
		// repeat the procedure
		restoreQ();
		initializeEigenVariables();
	}
	else {
		// we can now free the the eigensystem variables
		if (isComplex || isOldComplex)
			freeComplexEigens();
		freeEigens();
		useEigens = isComplex = isOldComplex = false;
	}
	useEigens = flag;

}

/*!
 * This function returns transition probabilities in the MbMatrix P
 *
 * \brief Calculate transition probabilities
 * \param v [in] Branch length (times any rate multiplier)
 * \param P [out] Matrix of transition probabilities
 * \return Returns reference to P
 */
MbMatrix<double> &MbTransitionMatrix::tiProbs(const double v, MbMatrix<double> &P) {
	
	if (!useEigens)
		tiProbsPade(v, P);
	else if (!isComplex)
		tiProbsEigens(v, P);
	else
		tiProbsComplexEigens(v, P);

	return P;
}

/*!
 * This function calculates transition probabilities using
 * complex eigenvalues and eigenvectors.
 *
 * \brief Calculate P matrix using complex eigensystem
 * \param v Branch length (times any rate multiplier)
 * \param P Matrix of transition probabilities
 * \todo Speed up by storing ceigValExp and cc_ijk in doubles
 */
void MbTransitionMatrix::tiProbsComplexEigens(const double v, MbMatrix<double> &P) {
	
	for (int s=0; s<numStates; s++)
		ceigValExp[s] = exp(ceigenvalue[s] * v);

	const complex<double> *ptr = cc_ijk;
	for (int i=0; i<numStates; i++) {
		for (int j=0; j<numStates; j++) {
			complex<double> sum = complex<double>(0.0, 0.0);
			for(int s=0; s<numStates; s++)
				sum += (*ptr++) * ceigValExp[s];
			P[i][j] = (sum.real() < 0.0) ? 0.0 : sum.real();
		}
	}

}

/*!
 * This function calculates transition probabilities using
 * eigenvalues and eigenvectors.
 *
 * \brief Calculate P matrix using eigensystem
 * \param v Branch length (times any rate multiplier)
 * \param P Matrix of transition probabilities
 */
void MbTransitionMatrix::tiProbsEigens(const double v, MbMatrix<double> &P) {
	
	for (int s=0; s<numStates; s++)
		eigValExp[s] = exp(eigenvalue[s] * v);

	double *ptr = c_ijk;
	for (int i=0; i<numStates; i++) {
		for (int j=0; j<numStates; j++) {
			double sum = 0.0;
			for(int s=0; s<numStates; s++)
				sum += (*ptr++) * eigValExp[s];
			P[i][j] = (sum < 0.0) ? 0.0 : sum;
		}
	}

}

/*!
 * This function calculates the transition probabilities using the Pade
 * approximation to the matrix exponential.
 *
 * \brief Calculate P matrix using Pade approximation
 *
 * \param v Branch length (times any rate multiplier)
 * \param *tiProbs Pointer to result array of transition probabilities
 *
 */
void MbTransitionMatrix::tiProbsPade(const double v, MbMatrix<double> &P) {

	MbMatrix<double> A = Q * v;
	MbMath::expMatrixPade(A, P, padeQValue);

}

/*!
 * This function initializes the Q matrix for a non-reversible
 * model, rescales it, and then performs the relevant precalculations
 * for the computation of transition probabilities. Note that the
 * Q matrix resulting from the input rateset can, in principle, be
 * time-reversible but we really do not gain anything by checking
 * whether it is. If you are interested in finding out whether
 * the matrix is time-reversible, call isReversible()
 *
 * \brief Update Q matrix and precalculations (time-irreversible)
 * \param rate Irreversible rates
 * \return 0 for success, 1 for failure (wrong dimension on rate)
 */
int MbTransitionMatrix::updateQ(const MbVector<double> &rate) {
    
	if ( rate.dim() != numStates*(numStates-1) )
		return (1);

	// Mark Q here as tentatively non-reversible
	// If someone calls us using isReversible(), we
	// will check to make sure
	isRev = false;

	// Swap values so that we can restore
	restoreQ();
	
	// Initialize the Q matrix
	int index = 0;
	for (int i=0; i<numStates; i++)
		Q[i][i] = 0.0;
	for (int i=0; i<numStates; i++) {
		for (int j=0; j<numStates; j++) {
			if (i != j)
				Q[i][i] -= (Q[i][j] = rate[index++]);
		}
	}
	calcStationaryFreq();
	
	// Rescale the Q matrix
	rescaleQ();

	// Initialize local variables for eigensystem calculations
	if (useEigens)
		initializeEigenVariables();

	return (0);

}

/*!
 * This function initializes the Q matrix for a time-reversible
 * model, rescales it, and then performs the relevant precalculations
 * for the computation of transition probabilities.
 *
 * \brief Update Q matrix and precalculations (time-reversible)
 * \param rate Reversible rates
 * \param pi Stationary state frequencies
 * \return 0 for success, 1 for failure (wrong dimensions on rate or pi)
 */
int MbTransitionMatrix::updateQ(const MbVector<double> &rate, const MbVector<double> &pi) {
    
	if ( (2*rate.dim() != numStates*(numStates-1)) || (pi.dim() != numStates) )
		return (1);

	// Swap values so that we can restore
	restoreQ ();

	// Store copy of the stationary frequencies if somebody wants them
	this->pi.inject(pi);

	// Initialize the Q matrix
	int index = 0;
	for (int i=0; i<numStates; i++)
		Q[i][i] = 0.0;
	for (int i=0; i<numStates; i++) {
		double freqi = pi[i];
		for (int j=i+1; j<numStates; j++) {
			double freqj = pi[j];
			double r = rate[index++];
			Q[i][i] -= (Q[i][j] = r * freqj);
			Q[j][j] -= (Q[j][i] = r * freqi);
		}
	}
	
	// Rescale the Q matrix
	rescaleQ();

	// Initialize local variables for eigensystem calculations
	if (useEigens)
		initializeEigenVariables();

	return (0);

}
