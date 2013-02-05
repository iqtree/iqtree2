/*! 
 * \file
 * This file declares MbVector, a templated class for vectors
 * (1D arrays). It also contains definitions for MbVector
 * and related functions operating on templated vectors, such
 * as printing and reading.
 *  
 * \brief Declaration and definitions for MbVector
 *
 * MrBayes version 4.0 beta
 *
 * (c) Copyright 2005.
 * \version 4.0 Beta
 * \date Last modified: $Date: 2006/09/26 21:47:23 $
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
 * $Id: MbVector.h,v 1.12 2006/09/26 21:47:23 paulvdm Exp $
 */

#ifndef MbVector_H
#define MbVector_H

#include <iostream>
#include <iomanip>
#include <cstdlib>

/*!
 * MrBayes templated vector type. We used the Template Numerical Toolkit
 * (TNT) code as a model for this class. TNT vectors are similar to 
 * the LAPACk vector type. The TNT code comes with the following
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
 * This header file provides a one-dimensional, numerical array 
 * which looks like a conventional C array. Elements are 
 * accessed via the familiar A[i] notation. 
 *
 * Array assignment is by reference (i.e. shallow assignment).
 * That is, B=A implies that the A and B point to the
 * same array, so modifications to the elements of A
 * will be reflected in B. If an independent copy
 * is required, then B = A.copy() can be used.  Note
 * that this facilitates returning arrays from functions
 * without relying on compiler optimizations to eliminate
 * extensive data copying.
 *
 * The indexing and layout of this array object makes
 * it compatible with C and C++ algorithms that utilize
 * the familiar C[i] notation.  This includes numerous
 * textbooks, such as Numercial Recipes, and various
 * public domain codes.
 *
 * This class employs its own garbage collection via
 * the use of reference counts.  That is, whenever
 * an internal array storage no longer has any references
 * to it, it is destroyed.
 *
 * \brief Vector class
 */
template <class T>
class MbVector {

public:
                            MbVector(void);                                  //!< null constructor (zero-length vector)
                 explicit   MbVector(int x);                                 //!< creates vector of length n without initializing values
                            MbVector(int x,  T *a);                          //!< creates vector of length n as a view of vector a
                            MbVector(int x, const T &a);                     //!< creates vector of length n initializing all elements with value a
                   inline   MbVector(const MbVector &A);                     //!< creates vector of length n sharing data with vector A
                            ~MbVector(void);                                 //!< destructor

                            operator T*() { return &(v[0]); }                //!< type cast to T pointer
				            operator const T*() { return &(v[0]); }          //!< type cast to T pointer for const
          inline MbVector   &operator=(const T &a);                          //!< assignment operator (all elements have the value a)
		         MbVector   &operator=(const MbVector &A) { return ref(A); } //!< assignment operator (shallow copy, elements share data)
                     bool   operator==(const MbVector &A) const;             //!< equality operator
                 inline T   &operator[](int i) { return v[i]; }              //!< indexing operator (allows reference of data using [] notation)
           inline const T   &operator[](int i) const { return v[i]; };       //!< indexing operator (allows reference of data using [] notation) (const)
          inline MbVector   &ref(const MbVector &A);                         //!< creates a reference to another array (shallow copy)
                 MbVector   copy(void) const;                                //!< creates a copy of another array (deep copy, with separate data elements)
                 MbVector   &inject(const MbVector & A);                     //!< copy the values of elements from one array to another
               inline int   dim1(void) const { return n; }                   //!< get the dimensions of the vector (number of elements of vector)
               inline int   dim(void) const { return n; }                    //!< get the dimensions of the vector (number of elements of vector)
			   inline int   getRefCount(void) const { return *refCount; }    //!< get the number of vectors that share the same data
			   inline int   size() const { return n; }                       //!< get the number of elements of the vector

                        T   *v;                                            //!< pointer to values
                      int   n;                                             //!< number of elements
                      int   *refCount;                                     //!< number of references to the vector
private:


                     void   destroy(void);                                 //!< garbage collector

};

template <class T> std::ostream& operator<<(std::ostream &s, const MbVector<T> &A);        //!< operator <<
template <class T> std::istream& operator>>(std::istream &s, MbVector<T> &A);              //!< operator >>

template <class T> MbVector<T>   operator+(const MbVector<T> &A, const MbVector<T> &B);    //!< operator +
template <class T> MbVector<T>   operator-(const MbVector<T> &A, const MbVector<T> &B);    //!< operator -
template <class T> MbVector<T>   operator*(const MbVector<T> &A, const MbVector<T> &B);    //!< operator *
template <class T> MbVector<T>   operator/(const MbVector<T> &A, const MbVector<T> &B);    //!< operator /
template <class T> MbVector<T>   &operator+=(const MbVector<T> &A, const MbVector<T> &B);  //!< operator +=
template <class T> MbVector<T>   &operator-=(const MbVector<T> &A, const MbVector<T> &B);  //!< operator -=
template <class T> MbVector<T>   &operator*=(const MbVector<T> &A, const MbVector<T> &B);  //!< operator *=
template <class T> MbVector<T>   &operator/=(const MbVector<T> &A, const MbVector<T> &B);  //!< operator /=

template <class T> bool          operator!=(const MbVector<T> &A, const MbVector<T> &B);  //!< operator /=

template <class T> MbVector<T>   operator*(const MbVector<T> &A, const T &b); //!< operator *

// Defintions of inlined member functions

/*!
 * Copy constructor, which creates a shallow copy of the
 * MbVector argument. Vector data are not copied but shared.
 * Thus, in MbVector B(A), subsequent changes to A will be
 * reflected by changes in B. For an independent copy, use
 * MbVector B(A.copy()), or B = A.copy(), instead. Note
 * the use of garbage collection in this class, through the
 * reference counter refCount.
 *
 * \brief Shallow copy constructor
 * \param A Vector to copy
 */
template <class T>
inline MbVector<T>::MbVector(const MbVector &A)
    : v(A.v), n(A.n), refCount(A.refCount) {

	(*refCount)++;
	
}

/*!
 * Assign all elements of the vector the value of
 * the constant scalar a.
 *
 * \brief Assign scalar to all elements
 * \param a Scalar used in assignment
 * \return Assigned vector
 */
template <class T>
inline MbVector<T> &MbVector<T>::operator=(const T &a) {

	T *end = v + n;
	for (T *p=v; p<end; p++)
		*p = a;
	return *this;

}

/*!
 * Create a reference (shallow assignment) to another existing vector.
 * In B.ref(A), B and A share the same data and subsequent changes
 * to the vector elements of one will be reflected in the other. Note that
 * the reference counter is always allocated, even for null vectors,
 * so we need not test whether refCount is NULL.
 *
 * This is what operator= calls, and B=A and B.ref(A) are equivalent
 * operations.
 *
 * \brief Make this reference to A
 * \param A Vector to take reference of
 * \return Vector for which reference was set to A
 */
template <class T>
inline MbVector<T> &MbVector<T>::ref(const MbVector &A) {

	if (this != &A)
		{
		(*refCount)--;
		if ( *refCount < 1)
			destroy();
		n = A.n;
		v = A.v;
		refCount = A.refCount;
		(*refCount)++;
		}
	return *this;
	
}


// Defintions of member functions that are not inlined

/*!
 * Null constructor. Creates a 0-length (NULL) vector.
 * Note that the reference count will be set to 1 for this
 * null vector. This is to simplify the rest of the
 * code at the cost of allocating and deleting an int
 * everytime a null vector is needed.
 *
 * \brief Null constructor
 */
template <class T>
MbVector<T>::MbVector(void)
    : v(0), n(0), refCount(0) {

	refCount = new int;
	*refCount = 1;
	
}

/*!
 * Create a new vector of length n, without initializing vector 
 * elements. If x is not positive, a null vector is created. Note
 * that the reference count will be set to 1 regardless of whether
 * we create a null vector. This is to simplify the rest of the
 * code at the cost of allocating and deleting an int every time
 * a null vector is needed.
 *
 * This version avoids the O(n) initialization overhead.
 *
 * \brief Constructor of uninitialized vector
 * \param x Dimension (length) of new vector
 */
template <class T>
MbVector<T>::MbVector(int x)
    : v(0), n(0), refCount(0) {

	if (x > 0) {
		v = new T[x];
		n = x;
	}    
	refCount = new int;
	*refCount = 1;

}

/*!
 * Create a new x-length vector as a view of an existing one-dimensional
 * C array. The storage for this pre-existing array will never be destroyed
 * by the MbVector class since the reference count is set to 2. When the
 * field is lost by MbVector, the reference count will still be 1 and the
 * field is not garbage collected.
 *
 * \param x The dimension (length) of the new vector
 * \param a Pointer to C array used as data storage
 */
template <class T>
MbVector<T>::MbVector(int x, T *a)
    : v(0), n(0) , refCount(0) {

	if (x > 0) {
		v = a;
		n = x;
	}
	refCount = new int;
	*refCount = 2;

}

/*!
 * Constructor, which creates a vector with x elements.
 * The elements will be initialized to the constant
 * specified by the second argument. Most often used to
 * create a vector of zeros, as in MbVector A(n, 0.0).
 * 
 * \brief Constructor of initialized vector.
 * \param x Number of elements.
 * \param a Value for initialization.
 */
template <class T>
MbVector<T>::MbVector(int x, const T &a)
    : v(0), n(0), refCount(0) {

	if (x > 0) {
		v = new T[x];
		n = x;
		T *end = v+n;
		for (T *p = v; p<end; p++)
			*p = a;
	}
	refCount = new int;
	*refCount = 1;

}

/*!
 * Destructor. Note that refCount is decreased and only if
 * refCount reaches 0 do we delete allocated memory. This is
 * garbage collection as implemented in JAVA and other languages.
 * Note that null vectors also have a reference count allocated,
 * so that we can always access the value of refCount.
 *
 * \brief Destructor with garbage collection
 */
template <class T>
MbVector<T>::~MbVector(void) {

	(*refCount)--;
	if (*refCount < 1)
		destroy();

}

/*!
 * Equality operator. The dimensions of the two vectors
 * are first compared. If it is not the same, then false
 * is returned. Second, all elements are compared. If
 * they are the same, true is returned, otherwise false
 * is returned. Note that this operator is not useful
 * for float and double vectors, but it is handy for int
 * and bool vectors, as well as for vectors of other types
 * that have a sensible operator!= defined.
 *
 * \brief Equality operator.
 * \param A Vector to compare this to.
 * \return True if this==A, false otherwise.
 */
template <class T>
bool MbVector<T>::operator==(const MbVector &A) const {

	if (n != A.n)
		return false;
	for (int i=0; i<n; i++)
		if (v[i] != A.v[i])
			return false;
	return true;

}

/*!
 * Create a new version of an existing vector.  Used in B = A.copy()
 * or in the construction of a new vector that does not share
 * data with the copied vector, e.g. in MbVector B(A.copy()).
 *
 * \brief Create independent copy
 * \return Copy of this.
 */
template <class T>
MbVector<T> MbVector<T>::copy(void) const {

	MbVector A(n);
	memcpy (A.v, v, n*sizeof(T));
	return A;

}

/*
 * Copy the elements from one vector to another, in place.
 * That is, if you call B.inject(A), both A and B must conform
 * (i.e. have the same dimension).
 *
 * This differs from B = A.copy() in that references to B
 * before this assignment are also affected.  That is, if
 * we have 
 *
 * MbVector A(n);
 * MbVector C(n);
 * MbVector B(C);        (elements of B and C are shared) 
 *
 * then B.inject(A) affects both B and C, while B=A.copy() creates
 * a new vector B which shares no data with C or A.
 *
 * A is the vector from which elements will be copied.
 * The function returns an instance of the modifed vector. That is, in 
 * B.inject(A), it returns B.  If A and B are not conformant, no 
 * modifications to B are made.
 *
 * \brief Inject elements of A into this
 * \param A Vector with elements to inject
 * \return Injected vector
 */
template <class T>
MbVector<T> &MbVector<T>::inject(const MbVector &A) {

	if (A.n == n)
		memcpy(v, A.v, n*sizeof(T));
	return *this;
	
}

/*!
 * This is a garbage collector, which is
 * called only when there is no more element
 * referencing this vector.
 *
 * \brief Garbage collection
 */
template <class T>
void MbVector<T>::destroy(void) {

	if (v != 0)
		delete [] (v);
	if (refCount != 0)
		delete refCount;

}


// Definitions of related templated functions on vectors

/*!
 * Printing of a vector to an ostream object.
 * We use the format '[<dim>] (v_1,v_2,v_3,...,v_n)'.
 *
 * \brief operator<<
 * \param A Vector to output
 * \param s ostream to output to
 * \return ostream object (for additional printing)
 */
template <class T>
std::ostream& operator<<(std::ostream &s, const MbVector<T> &A) {

	int N = A.dim();
	s << "(";
	for (int i=0; i<A.dim(); i++) {
		s << A[i];
		if (i != N-1)
			s << ",";
	}
	s << ")";
	return s;

}

/*!
 * Reading of a vector from an istream object.
 * We expect the format:
 * [<dim>] (v_1,v_2,v_3,...,v_n)
 * On failure, a null vector is returned.
 * 
 * \brief operator>>
 * \param A Vector to receive input
 * \param s istream to read from
 * \return istream object (for additional reading)
 */
template <class T>
std::istream& operator>>(std::istream &s, MbVector<T> &A) {

	A = MbVector<T>();	// make sure we return null vector on failure
	int N;
	char c;
	s >> c;
	if (c != '[')
		return s;
	s >> N;
	MbVector<T> B(N);
	s >> c;
	if (c != ']')
		return s;
	s.ignore();  // ignore the space
	s >> c;
	if (c != '(')
		return s;
	for (int i=0; i<N; i++) {
		s >> B[i];
		if (i < N-1) {
			s >> c;
			if (c != ',') return s;
		}
	}
	s >> c;
	if (c != ')')
		return s;
	A = B;
	return s;

}

/*!
 * This function performs elementwise addition on two
 * vectors and returns the resulting vector. If the
 * vectors are not conformant, a null vector is returned.
 *
 * \brief operator+
 * \param A Vector 1
 * \param B Vector 2
 * \return A + B, null vector on failure
 */
template <class T>
MbVector<T> operator+(const MbVector<T> &A, const MbVector<T> &B) {

	int n = A.dim1();
	if (B.dim1() != n )
		return MbVector<T>();
	else {
		MbVector<T> C(n);
		for (int i=0; i<n; i++)
			C[i] = A[i] + B[i];
		return C;
	}

}

/*!
 * This function performs elementwise subtraction on two
 * vectors and returns the resulting vector. If the
 * vectors are not conformant, a null vector is returned.
 *
 * \brief operator-
 * \param A Vector 1
 * \param B Vector 2
 * \return A - B, null vector on failure
 */
template <class T>
MbVector<T> operator-(const MbVector<T> &A, const MbVector<T> &B) {

	int n = A.dim1();
	if (B.dim1() != n )
		return MbVector<T>();
	else {
		MbVector<T> C(n);
		for (int i=0; i<n; i++)
			C[i] = A[i] - B[i];
		return C;
	}

}

/*!
 * This function performs elementwise multiplication on two
 * vectors and returns the resulting vector. If the
 * vectors are not conformant, a null vector is returned.
 *
 * \brief operator*
 * \param A Vector 1
 * \param B Vector 2
 * \return A * B, null vector on failure
 */
template <class T>
MbVector<T> operator*(const MbVector<T> &A, const MbVector<T> &B) {

	int n = A.dim1();
	if (B.dim1() != n )
		return MbVector<T>();
	else {
		MbVector<T> C(n);
		for (int i=0; i<n; i++)
			C[i] = A[i] * B[i];
		return C;
	}

}

/*!
 * This function performs elementwise multiplication on a
 * vector and a scalar value and returns the resulting vector. 
 *
 * \brief operator*
 * \param A Vector 1
 * \param B Scalar value
 * \return A * B, null vector on failure
 */
template <class T>
MbVector<T> operator*(const MbVector<T> &A, const T &b) {

    int n = A.dim1();
    MbVector<T> C(n);
    for(int i=0; i<n; i++)
        C[i] = A[i] * b;
    return C;
}
 
/*!
 * This function performs elementwise division on two
 * vectors and returns the resulting vector. If the
 * vectors are not conformant, a null vector is returned.
 *
 * \brief operator/
 * \param A Vector 1
 * \param B Vector 2
 * \return A / B, null vector on failure.
 */
template <class T>
MbVector<T> operator/(const MbVector<T> &A, const MbVector<T> &B) {

	int n = A.dim1();
	if (B.dim1() != n )
		return MbVector<T>();
	else {
		MbVector<T> C(n);
		for (int i=0; i<n; i++)
			C[i] = A[i] / B[i];
		return C;
	}

}

/*!
 * This function performs elementwise addition on two
 * vectors and puts the result in the first vector.
 * If the two vectors are nonconformant, the first
 * vector is left intact.
 *
 * \brief operator+=
 * \param A Vector 1
 * \param B Vector 2
 * \return A += B, A unmodified on failure
 */
template <class T>
MbVector<T>&  operator+=(MbVector<T> &A, const MbVector<T> &B) {

	int n = A.dim1();
	if (B.dim1() == n) {
		for (int i=0; i<n; i++)
			A[i] += B[i];
	}
	return A;

}

/*!
 * This function performs elementwise subtraction on two
 * vectors and puts the result in the first vector.
 * If the two vectors are nonconformant, the first
 * vector is left intact.
 *
 * \brief operator-=
 * \param A Vector 1
 * \param B Vector 2
 * \return A -= B, A unmodified on failure
 */
template <class T>
MbVector<T>&  operator-=(MbVector<T> &A, const MbVector<T> &B) {

	int n = A.dim1();
	if (B.dim1() == n) {
		for (int i=0; i<n; i++)
			A[i] -= B[i];
	}
	return A;

}

/*!
 * This function performs elementwise multiplication on two
 * vectors and puts the result in the first vector.
 * If the two vectors are nonconformant, the first
 * vector is left intact.
 *
 * \brief operator*=
 * \param A Vector 1
 * \param B Vector 2
 * \return A *= B, A unmodified on failure
 */
template <class T>
MbVector<T>&  operator*=(MbVector<T> &A, const MbVector<T> &B) {

	int n = A.dim1();
	if (B.dim1() == n) {
		for (int i=0; i<n; i++)
			A[i] *= B[i];
	}
	return A;

}

/*!
 * This function performs elementwise division on two
 * vectors and puts the result in the first vector.
 * If the two vectors are nonconformant, the first
 * vector is left intact.
 *
 * \brief operator/=
 * \param A Vector 1
 * \param B Vector 2
 * \return A /= B, A unmodified on failure
 */
template <class T>
MbVector<T>&  operator/=(MbVector<T> &A, const MbVector<T> &B) {

	int n = A.dim1();
	if (B.dim1() == n) {
		for (int i=0; i<n; i++)
			A[i] /= B[i];
	}
	return A;

}

/*!
 * Inequality operator. It calls operator== and returns
 * the reverse of the bool result. Note that this operator
 * is not useful for float and double vectors, but it is handy
 * for int and bool vectors, as well as for vectors of other types
 * that have a sensible operator!= defined.
 *
 * \brief Inequality operator
 * \param A Vector 1
 * \param B Vector 2
 * \return True if A != B, false otherwise
 */
template <class T>
bool operator!=(const MbVector<T> &A, const MbVector<T> &B) {

	if (A == B)
		return false;
	else
		return true;

}

#endif
