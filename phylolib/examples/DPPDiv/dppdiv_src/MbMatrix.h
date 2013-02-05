/*! 
 * \file
 * This file declares the class MbMatrix, a templated class for
 * matrices and common operations on matrices (2D arrays). It
 * also contains definitions for MbMatrix and related functions
 * operating on templated matrices, such as printing and reading.
 *  
 * \brief Declaration and definitions for MbMatrix
 *
 * MrBayes version 4.0 beta
 *
 * (c) Copyright 2005.
 * \version 4.0 Beta
 * \date Last modified: $Date: 2006/08/25 14:52:04 $
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
 * $Id: MbMatrix.h,v 1.9 2006/08/25 14:52:04 paulvdm Exp $
 */

#ifndef MbMatrix_H 
#define MbMatrix_H 
 
#include <iostream> 
#include <iomanip> 
#include <cstdlib> 
#include <cstring>
 

/*! 
 * MrBayes templated matrix type. We used the Template Numerical 
 * Toolkit (TNT) code as a model for this class. TNT is similar to the 
 * LAPACK matrix type. The TNT code comes with the following 
 * disclaimer: 
 *  
 * "Template Numerical Toolkit (TNT) 
 * 
 *    Mathematical and Computational Sciences Division 
 *    National Institute of Technology, 
 *    Gaithersburg, MD USA 
 * 
 * This software was developed at the National Institute of Standards and 
 * Technology (NIST) by employees of the Federal Government in the course 
 * of their official duties. Pursuant to title 17 Section 105 of the 
 * United States Code, this software is not subject to copyright protection 
 * and is in the public domain. NIST assumes no responsibility whatsoever for 
 * its use by other parties, and makes no guarantees, expressed or implied, 
 * about its quality, reliability, or any other characteristic." 
 * 
 * Storage corresponds to C (row-major) ordering. 
 * Elements are accessed via A[i][j] notation.  
 * 
 * Array assignment is by reference (i.e. shallow assignment). 
 * That is, B=A implies that the A and B point to the 
 * same matrix, so modifications to the elements of A 
 * will be reflected in B. If an independent copy 
 * is required, then B = A.copy() can be used.  Note 
 * that this facilitates returning matrices from functions 
 * without relying on compiler optimizations to eliminate 
 * extensive data copying. 
 * 
 * The indexing and layout of this matrix object makes 
 * it compatible with C and C++ algorithms that utilize 
 * the familiar C[i][j] notation.  This includes numerous 
 * textbooks, such as Numerical Recipes, and various 
 * public domain codes. 
 * 
 * This class employs its own garbage collection via 
 * the use of reference counts.  That is, whenever 
 * an internal array storage no longer has any references 
 * to it, it is destroyed. 
 * 
 * Note that the multiplication operator is overloaded to
 * do matrix multiplication rather than element-wise
 * multiplication.
 *
 * \brief Templated matrices and matrix operations 
 */ 
template <class T> 
class MbMatrix { 
 
public: 
                            MbMatrix(void);                                   //!< null constructor (0 X 0 matrix) 
                            MbMatrix(int m, int n);                           //!< creates a m X n matrix without initialization 
                            MbMatrix(int m, int n,  T *a);                    //!< creates a m X n matrix as a view of a one-dimensional array 
                            MbMatrix(int m, int n, const T &a);               //!< creates a m X n matrix initializing all elements to a constant 
                   inline   MbMatrix(const MbMatrix &A);                      //!< creates a m X n matrix with elements shared by another matrix, A 
                            ~MbMatrix(void);                                  //!< destructor 
 
							operator T**() { return &(v[0]); }                //!< type cast operator
							operator const T**() { return &(v[0]); }          //!< type cast operator for const
          inline MbMatrix   &operator=(const T &a);                           //!< assignment operator (all elements have the value a) 
		         MbMatrix   &operator=(const MbMatrix &A) { return ref(A); }  //!< assignment operator (shallow copy, elements share data) 
                     bool   operator==(const MbMatrix &A) const;              //!< equality operator 
                        T   *operator[](int i) { return v[i]; }               //!< indexing operator (allows reference of data using [][] notation) 
                  const T   *operator[](int i) const { return v[i]; }         //!< indexing operator (allows reference of data using [][] notation)
		  inline MbMatrix   &ref(const MbMatrix &A);                          //!< creates a reference to another matrix (shallow copy) 
                 MbMatrix   copy(void) const;                                 //!< creates a copy of another matrix (deep copy, with separate data elements) 
                 MbMatrix   &inject(const MbMatrix &A);                       //!< copy the values of elements from one matrix to another 
                      int   dim1(void) const { return m; }                    //!< number of rows 
                      int   dim2(void) const { return n; }                    //!< number of columns 
					  int   getRefCount(void) const { return *refCount; }     //!< get the number of matrices that share the same data 
 
	private: 
                        T   **v; 
                      int   m; 
                      int   n; 
                      int   *refCount; 
                     bool   cArray; 
 
                     void   destroy(void); 
 
}; 
 

template <class T> std::ostream& operator<<(std::ostream &s, const MbMatrix<T> &A);  //!< operator << 
template <class T> std::istream& operator>>(std::istream &s, MbMatrix<T> &A);        //!< operator >> 

template <class T> MbMatrix<T> operator+(const MbMatrix<T> &A, const MbMatrix<T> &B);     //!< operator + 
template <class T> MbMatrix<T> operator-(const MbMatrix<T> &A, const MbMatrix<T> &B);     //!< operator - 
template <class T> MbMatrix<T> operator*(const MbMatrix<T> &A, const MbMatrix<T> &B);     //!< operator * (matrix multiplication) 
template <class T> MbMatrix<T> &operator+=(const MbMatrix<T> &A, const MbMatrix<T> &B);   //!< operator += 
template <class T> MbMatrix<T> &operator-=(const MbMatrix<T> &A, const MbMatrix<T> &B);   //!< operator -= 
template <class T> MbMatrix<T> &operator*=(const MbMatrix<T> &A, const MbMatrix<T> &B);   //!< operator *= (matrix multiplication)
template <class T> MbMatrix<T> operator+(const T &a, const MbMatrix<T> &B);               //!< operator + for scalar + matrix 
template <class T> MbMatrix<T> operator-(const T &a, const MbMatrix<T> &B);               //!< operator - for scalar - matrix 
template <class T> MbMatrix<T> operator*(const T &a, const MbMatrix<T> &B);               //!< operator * for scalar * matrix 
template <class T> MbMatrix<T> operator/(const T &a, const MbMatrix<T> &B);               //!< operator / for scalar / matrix 
template <class T> MbMatrix<T> operator+(const MbMatrix<T> &A, const T &b);               //!< operator + for matrix + scalar 
template <class T> MbMatrix<T> operator-(const MbMatrix<T> &A, const T &b);               //!< operator - for matrix - scalar 
template <class T> MbMatrix<T> operator*(const MbMatrix<T> &A, const T &b);               //!< operator * for matrix * scalar 
template <class T> MbMatrix<T> operator/(const MbMatrix<T> &A, const T &b);               //!< operator / for matrix / scalar 
template <class T> MbMatrix<T> &operator+=(const MbMatrix<T> &A, const T &b);             //!< operator += for scalar 
template <class T> MbMatrix<T> &operator-=(const MbMatrix<T> &A, const T &b);             //!< operator -= for scalar 
template <class T> MbMatrix<T> &operator*=(const MbMatrix<T> &A, const T &b);             //!< operator *= for scalar 
template <class T> MbMatrix<T> &operator/=(const MbMatrix<T> &A, const T &b);             //!< operator /= for scalar 

template <class T> bool        operator!=(const MbMatrix<T> &A, const MbMatrix<T> &B);   //!< inequality (operator !=) 

// Definitions of inlined member functions

/*!
 * Copy constructor, which creates a shallow copy of the
 * MbMatrix argument. Matrix data are not copied but shared.
 * Thus, in MbMatrix B(A), subsequent changes to A will be
 * reflected by changes in B. For an independent copy, use
 * MbMatrix B(A.copy()), or B = A.copy(), instead. Note
 * the use of garbage collection in this class, through the
 * reference counter refCount.
 *
 * \brief Shallow copy constructor
 * \param A Matrix to copy
 */
template <class T>
inline MbMatrix<T>::MbMatrix(const MbMatrix &A)
    : v(A.v), m(A.m), n(A.n), refCount(A.refCount), cArray(A.cArray) {

	(*refCount)++;
}

/*!
 * Assign all elements of the matrix the value of
 * the constant scalar a.
 *
 * \brief Assign scalar to all elements
 * \param a Scalar used in assignment
 * \return Assigned matrix
 */
template <class T>
inline MbMatrix<T> &MbMatrix<T>::operator=(const T &a) {

	T *p = &(v[0][0]);
	T *end = p + m*n;
	for (; p<end; p++)
		*p = a;
	return *this;

}

/*!
 * Create a reference (shallow assignment) to another existing matrix.
 * In B.ref(A), B and A share the same data and subsequent changes
 * to the matrix elements of one will be reflected in the other. Note that
 * the reference counter is always allocated, even for null matrices,
 * so we need not test whether refCount is NULL.
 *
 * This is what operator= calls, and B=A and B.ref(A) are equivalent
 * operations.
 *
 * \brief Make this reference to A
 * \param A Matrix to take reference of
 * \return Matrix to which the reference was assigned
 */
template <class T>
inline MbMatrix<T> &MbMatrix<T>::ref(const MbMatrix &A) {

	if (this != &A)
		{
		(*refCount)--;
		if (*refCount < 1)
			destroy();
		m = A.m;
		n = A.n;
		v = A.v;
		refCount = A.refCount;
		cArray = A.cArray;
        (*refCount)++;
		}
	return *this;
	
}


// Defintions of member functions that are not inlined

/*!
 * Null constructor. Creates a (0 X O) ('null') matrix.
 * Note that the reference count will be set to 1 for this
 * null matrix. This is to simplify the rest of the
 * code at the cost of allocating and deleting an int
 * everytime a null matrix is needed.
 *
 * \brief Null constructor
 */
template <class T>
MbMatrix<T>::MbMatrix(void)
    : v(0), m(0), n(0), refCount(0), cArray(false) {

	refCount = new int;
	*refCount = 1;
}

/*!
 * Create a new (m X n) matrix, without initializing matrix 
 * elements. If m or n are not positive, a null matrix is created. Note
 * that the reference count will be set to 1 regardless of whether
 * we create a null matrix. This is to simplify the rest of the
 * code at the cost of allocating and deleting an int every time
 * a null matrix is needed.
 *
 * This version avoids the O(m*n) initialization overhead.
 *
 * \brief Constructor of uninitialized matrix
 * \param m The first (row) dimension of the matrix
 * \param n The second (column) dimension of the matrix
 */
template <class T>
MbMatrix<T>::MbMatrix(int m, int n)
    : v(0), m(0), n(0), refCount(0), cArray(false) {

	if (m > 0 && n > 0) {
		// allocate and initialize pointers
		T *p = new T[m * n];
		v = new T*[m];
		for (int i=0; i<m; i++) {
			v[i] = p;
			p += n;
		}
		this->m = m;
		this->n = n;
	}    
	refCount = new int;
	*refCount = 1;
}

/*!
 * Create a new (m X n) matrix as a view of an existing 'two-dimensional'
 * C array stored in C order, i.e. right-most dimension varying fastest,  
 * often referred to as "row-major" ordering. The storage for this pre-existing
 * array will never be destroyed by the MbMatrix class because of the cArray
 * flag, which is set to true for matrix space allocated in this way. However,
 * garbage collection will make sure that row pointers are deallocated when
 * appropriate.
 *
 * \param m The first (row) dimension of the matrix
 * \param n The second (column) dimension of the matrix
 * \param a Pointer to C array used as data storage
 */
template <class T>
MbMatrix<T>::MbMatrix(int m, int n, T *a)
    : v(0), m(0), n(0) , refCount(0), cArray(false) {

	if (m > 0 && n > 0) {
		// initialize pointers
		T *p = a;
		v = new T*[m];
		for (int i=0; i<m; i++) {
			v[i] = p;
			p += n;
		}
		this->m = m;
		this->n = n;
	}
	refCount = new int;
	*refCount = 1;
	cArray = true;
}

/*!
 * Create a new (m X n) matrix, initializing matrix elements
 * to the constant value specified by the third argument. Most
 * often used to create a matrix of zeros, as in MbMatrix A(m, n, 0.0).
 * 
 * \brief Constructor of initialized matrix.
 * \param m The first (row) dimension of the matrix
 * \param n The second (column) dimension of the matrix
 * \param a Value for initialization.
 */
template <class T>
MbMatrix<T>::MbMatrix(int m, int n, const T &a)
    : v(0), m(0), n(0), refCount(0), cArray(false) {

	if (m > 0 && n > 0) {
		// allocate and initialize pointers
		T *p = new T[m * n];
		v = new T*[m];
		for (int i=0; i<m; i++) {
			v[i] = p;
			p += n;
		}
		this->m = m;
		this->n = n;
		// set values
		T *end = &(v[0][0]) + m*n;
		for (p=&(v[0][0]); p<end; p++)
			*p = a;
	}
	refCount = new int;
	*refCount = 1;
}

/*!
 * Destructor. Note that refCount is decreased and only if
 * refCount reaches 0 do we delete allocated memory. This is
 * garbage collection as implemented in JAVA and other languages.
 * Note that null matrices also have a reference count allocated,
 * so that we can always access the value of refCount.
 *
 * \brief Destructor with garbage collection
 */
template <class T>
MbMatrix<T>::~MbMatrix(void) {

	(*refCount)--;
	if (*refCount < 1)
		destroy();

}

/*!
 * Equality operator. The dimensions of the two matrices
 * are first compared. If they are not the same, then false
 * is returned. Second, all elements are compared. If
 * they are the same, true is returned, otherwise false
 * is returned. Note that this operator is not useful
 * for float and double matrices, but it is handy for int
 * and bool matrices, as well as for matrices of other types
 * that have a sensible operator!= defined.
 *
 * \brief Equality operator
 * \param A Matrix to compare (*this) to
 * \return True if (*this)==A, false otherwise.
 */
template <class T>
bool MbMatrix<T>::operator==(const MbMatrix &A) const {

	if (m != A.m || n != A.n)
		return false;
	for (int i=0; i<m; i++) {
		for (int j=0; j<n; j++)
			if (v[i][j] != A.v[i][j])
				return false;
	}
	return true;

}

/*!
 * Create a new version of an existing matrix.  Used in B = A.copy()
 * or in the construction of a new matrix that does not share
 * data with the copied matrix, e.g. in MbMatrix B(A.copy()).
 *
 * \brief Create independent copy
 * \return Copy of this.
 */
template <class T>
MbMatrix<T> MbMatrix<T>::copy(void) const {

	MbMatrix A(m, n);
	memcpy (&(A.v[0][0]), &(v[0][0]), m*n*sizeof(T));
	return A;

}

/*
 * Copy the elements from one matrix to another, in place.
 * That is, if you call B.inject(A), both A and B must conform
 * (i.e. have the same row and column dimensions).
 *
 * This differs from B = A.copy() in that references to B
 * before this assignment are also affected.  That is, if
 * we have 
 *
 * MbMatrix A(n);
 * MbMatrix C(n);
 * MbMatrix B(C);        (elements of B and C are shared) 
 *
 * then B.inject(A) affects both B and C, while B=A.copy() creates
 * a new matrix B which shares no data with C or A.
 *
 * A is the matrix from which elements will be copied.
 * The function returns an instance of the modifed matrix. That is, in 
 * B.inject(A), it returns B.  If A and B are not conformant, no 
 * modifications to B are made.
 *
 * \brief Inject elements of A into (*this)
 * \param A Matrix with elements to inject
 * \return Injected matrix
 */
template <class T>
MbMatrix<T> &MbMatrix<T>::inject(const MbMatrix &A) {

	if (A.m == m && A.n == n)
		memcpy(&(v[0][0]), &(A.v[0][0]), m*n*sizeof(T));
	return *this;
	
}

/*!
 * This is a garbage collector, which is
 * called only when there is no more element
 * referencing this matrix.
 *
 * \brief Garbage collection
 */
template <class T>
void MbMatrix<T>::destroy(void) {

	if (v != 0) {
		if (cArray == false) 
			{
			delete [] (v[0]);
			}
		delete [] (v);
	}
	if (refCount != 0)
		{
		delete refCount;
		}
	cArray = false;
}


// Definitions of related templated functions on matrices

/*!
 * Printing of a matrix to an ostream object.
 * We use the format:
 * [<m>,<n>]
 * ((v_11,v_12,v_13,...,v_1n),
 * (v_i1,v_i2,v_i3,...,v_in),
 * (v_m1,v_m2,v_m3,...,v_mn)) 
 *
 * \brief operator<<
 * \param A Matrix to output
 * \param s ostream to output to
 * \return ostream object (for additional printing)
 */
template <class T>
std::ostream& operator<<(std::ostream &s, const MbMatrix<T> &A) {

	int M = A.dim1();
	int N = A.dim2();
	s << "[" << M << "," << N << "]\n";
	s << "(";
	for (int i=0; i<M; i++) {
		s << "(";
		for (int j=0; j<N; j++) {
			s << A[i][j];
			if (j != N-1)
				s << ",";
		}
		if (i != M-1)
			s << "),\n";
		else
			s << ")";
	}
	s << ")";
	return s;
	
}

/*!
 * Reading of a matrix from an istream object.
 * We expect the format:
 * [<m>,<n>]
 * ((v_11,v_12,v_13,...,v_1n),
 * (v_i1,v_i2,v_i3,...,v_in),
 * (v_m1,v_m2,v_m3,...,v_mn)) 
 * On failure, a null matrix is returned.
 *
 * \brief operator>>
 * \param A Matrix to receive input
 * \param s istream to read from
 * \return istream object (for additional reading)
 */
template <class T>
std::istream& operator>>(std::istream &s, MbMatrix<T> &A) {

	A = MbMatrix<T>();	// make sure we return null matrix on failure
	int M, N;
	char c;
	s >> c;
	if (c != '[')
		return s;
	s >> M;
	s >> c;
	if (c != ',')
		return s;
	s >> N;
	MbMatrix<T> B(M,N);
	s >> c;
	if (c != ']')
		return s;
	s.ignore();  // ignore the newline
	s >> c;
	if (c != '(')
		return s;
	for (int i=0; i<M; i++) {
		s >> c;
		if (c != '(')
			return s;
		for (int j=0; j<N; j++) {
			s >> B[i][j];
			if (j < N-1) {
				s >> c;
				if (c != ',') return s;
			}
		}
		s >> c;
		if (c != ')')
			return s;
		if (i != M-1)
			s.ignore();  // ignore newline
	}
	s >> c;
	if (c != ')')
		return s;
	A = B;
	return s;

}

/*!
 * This function performs addition of a scalar to
 * each element of a matrix and returns the
 * resulting matrix.
 *
 * \brief operator+ (scalar)
 * \param A Matrix
 * \param b Scalar
 * \return A + b
 */
template <class T>
MbMatrix<T> operator+(const MbMatrix<T> &A, const T &b) {

	MbMatrix<T> B(A.copy());
	for (int i=0; i<B.dim1(); i++)
		for (int j=0; j<B.dim2(); j++)
			B[i][j] = A[i][j] + b;
	return B;

}

/*!
 * This function performs subtraction of a scalar from
 * each element of a matrix and returns the
 * resulting matrix.
 *
 * \brief operator- (scalar)
 * \param A Matrix
 * \param b Scalar
 * \return A - b
 */
template <class T>
MbMatrix<T> operator-(const MbMatrix<T> &A, const T &b) {

	MbMatrix<T> B(A.copy());
	for (int i=0; i<B.dim1(); i++)
		for (int j=0; j<B.dim2(); j++)
			B[i][j] = A[i][j] - b;
	return B;

}

/*!
 * This function performs multiplication of a scalar to
 * each element of a matrix and returns the
 * resulting matrix.
 *
 * \brief operator* (scalar)
 * \param A Matrix
 * \param b Scalar
 * \return A * b
 */
template <class T>
MbMatrix<T> operator*(const MbMatrix<T> &A, const T &b) {

	MbMatrix<T> B(A.copy());
	for (int i=0; i<B.dim1(); i++)
		for (int j=0; j<B.dim2(); j++)
			B[i][j] = A[i][j] * b;
	return B;

}

/*!
 * This function performs division with a scalar of
 * each element of a matrix and returns the
 * resulting matrix.
 *
 * \brief operator/ (scalar)
 * \param A Matrix
 * \param b Scalar
 * \return A / b
 */
template <class T>
MbMatrix<T> operator/(const MbMatrix<T> &A, const T &b) {

	MbMatrix<T> B(A.copy());
	for (int i=0; i<B.dim1(); i++)
		for (int j=0; j<B.dim2(); j++)
			B[i][j] = A[i][j] / b;
	return B;

}

/*!
 * This function performs addition of a scalar to
 * each element of a matrix and returns the
 * resulting matrix.
 *
 * \brief operator+ (scalar first)
 * \param a Scalar
 * \param B Matrix
 * \return a + B
 */
template <class T>
MbMatrix<T> operator+(const T &a, const MbMatrix<T> &B) {

	MbMatrix<T> A(B.copy());
	for (int i=0; i<A.dim1(); i++)
		for (int j=0; j<A.dim2(); j++)
			A[i][j] = a + B[i][j];
	return A;

}

/*!
 * This function subtracts each element of a
 * a matrix from a scalar and returns the
 * resulting matrix.
 *
 * \brief operator- (scalar first)
 * \param a Scalar
 * \param B Matrix
 * \return a - B
 */
template <class T>
MbMatrix<T> operator-(const T &a, const MbMatrix<T> &B) {

	MbMatrix<T> A(B.copy());
	for (int i=0; i<A.dim1(); i++)
		for (int j=0; j<A.dim2(); j++)
			A[i][j] = a - B[i][j];
	return A;

}

/*!
 * This function performs multiplication of a scalar to
 * each element of a matrix and returns the
 * resulting matrix.
 *
 * \brief operator* (scalar first)
 * \param a Scalar
 * \param B Matrix
 * \return a * B
 */
template <class T>
MbMatrix<T> operator*(const T &a, const MbMatrix<T> &B) {

	MbMatrix<T> A(B.copy());
	for (int i=0; i<A.dim1(); i++)
		for (int j=0; j<A.dim2(); j++)
			A[i][j] = a * B[i][j];
	return A;

}

/*!
 * This function performs division of a scalar by
 * each element of a matrix and returns the
 * resulting matrix.
 *
 * \brief operator/ (scalar first)
 * \param a Scalar
 * \param B Matrix
 * \return a / B
 */
template <class T>
MbMatrix<T> operator/(const T &a, const MbMatrix<T> &B) {

	MbMatrix<T> A(B.copy());
	for (int i=0; i<A.dim1(); i++)
		for (int j=0; j<A.dim2(); j++)
			A[i][j] = a / B[i][j];
	return A;

}

/*!
 * This function performs addition of a scalar to
 * each element of a matrix in place and returns the
 * resulting matrix.
 *
 * \brief operator+= (scalar)
 * \param A Matrix
 * \param b Scalar
 * \return A += b
 */
template <class T>
MbMatrix<T> &operator+=(MbMatrix<T> &A, const T &b) {

	for (int i=0; i<A.dim1(); i++)
		for (int j=0; j<A.dim2(); j++)
			A[i][j] += b;
	return A;

}

/*!
 * This function performs subtraction of a scalar from
 * each element of a matrix in place and returns the
 * resulting matrix.
 *
 * \brief operator-= (scalar)
 * \param A Matrix
 * \param b Scalar
 * \return A -= b
 */
template <class T>
MbMatrix<T> &operator-=(MbMatrix<T> &A, const T &b) {

	for (int i=0; i<A.dim1(); i++)
		for (int j=0; j<A.dim2(); j++)
			A[i][j] -= b;
	return A;

}

/*!
 * This function performs multiplication of a scalar to
 * each element of a matrix in place and returns the
 * resulting matrix.
 *
 * \brief operator*= (scalar)
 * \param A Matrix
 * \param b Scalar
 * \return A *= b
 */
template <class T>
MbMatrix<T> &operator*=(MbMatrix<T> &A, const T &b) {

	for (int i=0; i<A.dim1(); i++)
		for (int j=0; j<A.dim2(); j++)
			A[i][j] *= b;
	return A;

}

/*!
 * This function performs division with a scalar of
 * each element of a matrix in place and returns the
 * resulting matrix.
 *
 * \brief operator/= (scalar)
 * \param A Matrix
 * \param b Scalar
 * \return A /= b
 */
template <class T>
MbMatrix<T> &operator/=(MbMatrix<T> &A, const T &b) {

	for (int i=0; i<A.dim1(); i++)
		for (int j=0; j<A.dim2(); j++)
			A[i][j] /= b;
	return A;

}

/*!
 * This function performs elementwise addition of two
 * matrices and returns the resulting matrix. If the
 * matrices are not conformant, a null matrix is returned.
 *
 * \brief operator+
 * \param A Matrix 1
 * \param B Matrix 2
 * \return A + B, null matrix on failure
 */
template <class T>
MbMatrix<T> operator+(const MbMatrix<T> &A, const MbMatrix<T> &B) {

	int m = A.dim1();
	int n = A.dim2();
	if (B.dim1() != m ||  B.dim2() != n)
		return MbMatrix<T>();
	else {
		MbMatrix<T> C(m,n);
		for (int i=0; i<m; i++) {
			for (int j=0; j<n; j++)
				C[i][j] = A[i][j] + B[i][j];
		}
		return C;
	}

}

/*!
 * This function performs elementwise subtraction of two
 * matrices and returns the resulting matrix. If the
 * matrices are not conformant, a null matrix is returned.
 *
 * \brief operator-
 * \param A Matrix 1
 * \param B Matrix 2
 * \return A - B, null matrix on failure
 */
template <class T>
MbMatrix<T> operator-(const MbMatrix<T> &A, const MbMatrix<T> &B) {

	int m = A.dim1();
	int n = A.dim2();
	if (B.dim1() != m ||  B.dim2() != n)
		return MbMatrix<T>();
	else {
		MbMatrix<T> C(m,n);
		for (int i=0; i<m; i++) {
			for (int j=0; j<n; j++)
				C[i][j] = A[i][j] - B[i][j];
		}
		return C;
	}

}

/*!
 * Compute C = A*B, where C[i][j] is the dot-product of 
 * row i of A and column j of B. Note that this operator
 * does not perform elementwise multiplication. If the 
 * matrices do not have the right dimensions for matrix
 * multiplication (that is, if the number of columns of A
 * is different from the number of rows of B), the function
 * returns a null matrix.
 *
 * \brief Matrix multiplication
 * \param A An (m X n) matrix
 * \param B An (n X k) matrix
 * \return A * B, an (m X k) matrix, or null matrix on failure
 */
template <class T>
MbMatrix<T> operator*(const MbMatrix<T> &A, const MbMatrix<T> &B) {

	if ( A.dim2() != B.dim1() )
		return MbMatrix<T>();
	int M = A.dim1();
	int N = A.dim2();
	int K = B.dim2();
	MbMatrix<T> C(M,K);
	for (int i=0; i<M; i++) {
		for (int j=0; j<K; j++) {
			T sum = 0;
			for (int k=0; k<N; k++)
				sum += A[i][k] * B [k][j];
			C[i][j] = sum;
		}
	}
	return C;

}

/*!
 * This function performs elementwise addition on two
 * matrices and puts the result in the first matrix.
 * If the two matrices are nonconformant, the first
 * matrix is left intact.
 *
 * \brief operator+=
 * \param A Matrix 1
 * \param B Matrix 2
 * \return A += B, A unmodified on failure
 */
template <class T>
MbMatrix<T>&  operator+=(MbMatrix<T> &A, const MbMatrix<T> &B) {

	int m = A.dim1();
	int n = A.dim2();
	if (B.dim1() == m && B.dim2() == n) {
		for (int i=0; i<m; i++) {
			for (int j=0; j<n; j++)
				A[i][j] += B[i][j];
		}
	}
	return A;

}

/*!
 * This function performs elementwise subtraction on two
 * matrices and puts the result in the first matrix.
 * If the two matrices are nonconformant, the first
 * matrix is left intact.
 *
 * \brief operator-=
 * \param A Matrix 1
 * \param B Matrix 2
 * \return A -= B, A unmodified on failure
 */
template <class T>
MbMatrix<T>&  operator-=(MbMatrix<T> &A, const MbMatrix<T> &B) {

	int m = A.dim1();
	int n = A.dim2();
	if (B.dim1() == m && B.dim2() == n) {
		for (int i=0; i<m; i++) {
			for (int j=0; j<n; j++)
				A[i][j] -= B[i][j];
		}
	}
	return A;

}

/*!
 * Compute C = A*B, where C[i][j] is the dot-product of 
 * row i of A and column j of B. Then assign the result to
 * A. Note that this operator does not perform elementwise
 * multiplication. If the matrices are not both square of the
 * same dimension, then the operation is not possible to
 * perform and we return an unomidified A.
 *
 * \brief Matrix multiplication with assignment (operator *=)
 * \param A An (n X n) matrix
 * \param B An (n X n) matrix
 * \return A = A * B, an (n X n) matrix, or unmodified A on failure
 */
template <class T>
MbMatrix<T> &operator*=(MbMatrix<T> &A, const MbMatrix<T> &B) {

	if (A.dim1()==A.dim2() && B.dim1()==B.dim2() && A.dim1()==B.dim1()) {
		int N = A.dim1();
		MbMatrix<T> C(N,N);
		for (int i=0; i<N; i++) {
			for (int j=0; j<N; j++) {
				T sum = 0;
				for (int k=0; k<N; k++)
					sum += A[i][k] * B [k][j];
				C[i][j] = sum;
			}
		}
		A = C;
	}
	return A;

}

/*!
 * Inequality operator. It calls operator== and returns
 * the reverse of the bool result. Note that this operator
 * is not useful for float and double matrices, but it is handy
 * for int and bool matrices, as well as for matrices of other types
 * that have a sensible operator!= defined.
 *
 * \brief Inequality operator
 * \param A Matrix 1
 * \param B Matrix 2
 * \return True if A != B, false otherwise
 */
template <class T>
bool operator!=(const MbMatrix<T> &A, const MbMatrix<T> &B) {

	if (A == B)
		return false;
	else
		return true;

}

#endif
