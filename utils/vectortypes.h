//
//  vectortypes.h
//  iqtree
//
//  Created by James Barbetti on 5/2/21.
//

#ifndef vectortypes_h
#define vectortypes_h
#include <vector>
#include <algorithm> //for std::transform
#include <string>    //for std::string

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
        vector of intptr_t
 */
typedef std::vector<intptr_t> Int64Vector;

/**
 * @brief  A subclass of vector (or other container) class S, where the
 *         elements in (or members of) S are treated as though they are 
 *         of type T.
 * @tparam T - the type we want to pretend is stored in the container
 * @tparam S - the container superclass, that contains members of a type
 *             (which should be interconvertible with T, via implicit
 *              casting).
 * @note   S must implement an S::size_type constructor.
 * @note   S must have S::const_iterator and S::iterator derived types
 * @note   S must implement const and non-const implementations 
 *           of begin() and end() that return S::const_iterator 
 *           and S::iterator.
 * @note   S must implement const and non-const operator[].
 * @note   All other member functions of S are exposed just as they are,
 *         for CastingVector<T,S> inherits publicly from S.
 */
template <class T, class S> class CastingVector: public S {
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
            explicit const_iterator( typename super::const_iterator i) 
                : super::const_iterator(i) {};
            inline T operator*() { return super::const_iterator::operator*(); }
    };

    class iterator : public super::iterator {
        public:
            explicit iterator( typename super::iterator i) 
                : super::iterator(i) {};
            inline T operator*() { return T(super::iterator::operator*()); }
    };
    
    const_iterator begin() const { return const_iterator ( super::begin() ); }
    iterator       begin()       { return iterator( super::begin() ); }
    const_iterator end()   const { return const_iterator ( super::end() ); }
    iterator       end()         { return iterator( super::end() ); }
    inline T operator[] (typename super::size_type i) const {
        return T(super::at(i));
    }
    
    class reference {
        private:
            S& to_vector;
            size_type at_index;
        public:
            reference(S& vector, size_type index):
                to_vector(vector), at_index(index) {}
            operator T() { return T(to_vector[at_index]); }
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

/**
 * @brief In single-threaded code, this is equivalent to std::vector<bool>
 *        (but uses about 8 times as much memory!).
 *        In multi-threaded code, this is NOT equivalent.  Because the
 *        elements are char rather than individual bits in bit-fields.
 *        So reads and writes of individual elements, in separate threads,
 *        are less likely to interfere with one another (but concurrent write
 *        access is still potentially problematic, if two or more threads 
 *        write to and/or read from the same element).
 * 
 */
typedef CastingVector<bool, std::vector<char> > BoolVector;

/**
        vector of char
 */
typedef std::vector<char> CharVector;

/**
        vector of std::string
 */
class StrVector: public std::vector<std::string> {
public:
    typedef std::vector<std::string> super;
    StrVector() = default;
    StrVector(const StrVector&) = default;
    explicit StrVector(size_t size);
    ~StrVector() = default;

    /**
     * @brief  Append string literals, from an array of const char*.
     * @tparam COUNT - the number of string literals in the array.
     * @param  literals - the array of const char* string literals.
     * @param  force_upper_case - true if the strings are to be forced to
     *         upper case.
     * @return StrVector& - reference to *this.
     * @note   upper-case is done character-by-character via toupper()
     *         (which may well NOT be appropriate in languages other 
     *          than English) -James B.
     */
    template <int COUNT> StrVector& appendLiterals(const char* (&literals)[COUNT], 
                                                   bool force_upper_case = false) {
        for (int i=0; i<COUNT; ++i) {
            std::string literal = literals[i];
            if (force_upper_case) {
                std::transform(literal.begin(), literal.end(), literal.begin(),
                               []( char c){ return toupper(c); });
            }
            emplace_back(literal);
        }
        return *this;
    }

    /**
     * @brief  Set to the string literals, read from an array of const char*.
     * @tparam COUNT - the number of string literals in the array.
     * @param  literals - the array of const char* string literals.
     * @param  force_upper_case - true if the strings are to be forced to
     *         upper case. false (the default), if not.
     * @return StrVector& - reference to *this.
     * @note   implemented in terms of appendLiterals().
     */
    template <int COUNT> StrVector& setToLiterals(const char* (&literals)[COUNT], 
                                                  bool force_upper_case = false ) {
        clear();
        return appendLiterals(literals, force_upper_case);
    }

    /**
     * @brief  indicate whether a given string literal is currently
     *         contained by *this.
     * @param  find_me a (const char*) string literal 
     * @return true  - if it (find_me) is contained in *this.
     * @return false - if it is not.
     * @note   Comparison is done via std::string's operator==.
     */
    bool contains(const char* find_me)        const;

    /**
     * @brief  indicate whether a given string literal is currently
     *         contained by *this.
     * @param  find_me the string (a const std::string reference) to
     *         look for.
     * @return true  - if it (find_me) is contained in *this.
     * @return false - if it is not.
     * @note   Comparison is done via std::string's operator==.
     */
    bool contains(const std::string& find_me) const;

    /**
     * @brief  construct a string, by concating all the strings in *this,
     *         in order, separating (not terminating) strings from *this
     *         with a specified separator.
     * @param  separator - the (const char* string literal) separator to use
     * @return std::string - the constructed string
     * @note   if size() yields zero, returns an empty string.
     */
    std::string join(const char* separator)   const;

    /**
     * @brief  construct a string, by concating all the strings in *this,
     *         in order, separating (not terminating) strings from *this
     *         with a specified separator.
     * @param  separator - the (const char* string literal) separator to use
     * @return std::string - the constructed string
     * @note   if size() yields zero, returns an empty string.
     */
    std::string join(const std::string& separator) const;

    /**
     * @brief  sort the strings in *this into lexicographic order.
     * @note   sorting uses std::sort, and comparison is done using std::less.
     */

    void sort();
};

/**
        matrix of double number
 */
#define mmatrix(T) std::vector< std::vector<T> >

/**
        matrix of double
 */

typedef mmatrix(double) DoubleMatrix;

#endif /* vectortypes_h */
