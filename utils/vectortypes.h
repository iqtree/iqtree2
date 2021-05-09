//
//  vectortypes.h
//  iqtree
//
//  Created by James Barbetti on 5/2/21.
//

#ifndef vectortypes_h
#define vectortypes_h
#include <vector>
#include <algorithm>

/**
        vector of double number
 */
typedef std::vector<double> DoubleVector;

/**
        vector of int
 */
typedef std::vector<int> IntList;


/**
        vector of int
 */
typedef std::vector<int> IntVector;

/**
        vector of bool
 */


template <class T, class S> class CastingVector: public S {
    //
    //A subclass of vector (or other container) class S,
    //where the entries in S are treated as though
    //they are of type T (via implicit casting).
    //
public:
    typedef S super;
    typedef typename S::size_type size_type;
    typedef CastingVector<T,S> this_type;
    
    CastingVector(): super() {}
    explicit CastingVector(size_type initialSize): super(initialSize) {}
    CastingVector(size_type initialSize, const T initialValue)
        : super(initialSize, initialValue) {}
    CastingVector(const CastingVector& rhs): super(rhs) {}
    
    class const_iterator: public super::const_iterator {
        public:
            const_iterator( typename super::const_iterator i) : super::const_iterator(i) {};
            inline T operator*() { return super::const_iterator::operator*(); }
    };

    class iterator : public super::iterator {
        public:
            iterator( typename super::iterator i) : super::iterator(i) {};
            inline T operator*() { return super::iterator::operator*(); }
    };
    
    const_iterator begin() const { return const_iterator ( super::begin() ); }
    iterator       begin()       { return iterator( super::begin() ); }
    const_iterator end()   const { return const_iterator ( super::end() ); }
    iterator       end()         { return iterator( super::end() ); }
    inline T operator[] (typename super::size_type i) const {
        return super::operator[] (i);
    }
    
    class reference {
        private:
            S& to_vector;
            size_type at_index;
        public:
            reference(S& vector, size_type index):
                to_vector(vector), at_index(index) {}
            operator T() { return to_vector[at_index]; }
            reference& operator= (const T new_value) {
                to_vector[at_index] = new_value;
                return *this;
            }
            reference& operator= (const reference& elsewhere) {
                to_vector[at_index] = elsewhere.to_vector[elsewhere.at_index];
                return *this;
            }
    };
    inline reference operator[] (typename super::size_type i) {
        return reference(*(dynamic_cast<S*>(this)), i);
    }
    this_type& reverseAll() {
        std::reverse ( super::begin(), super::end());
        return *this;
    }
};

typedef CastingVector<bool, std::vector<char>> BoolVector;

/**
        vector of char
 */
typedef std::vector<char> CharVector;

/**
        vector of string
 */
typedef std::vector<std::string> StrVector;


/**
        matrix of double number
 */
#define mmatrix(T) std::vector< std::vector<T> >

/**
        matrix of double
 */

typedef mmatrix(double) DoubleMatrix;


#endif /* vectortypes_h */
