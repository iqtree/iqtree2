/*

    Abstract class for 64-bit floating point numbers
    Based on vectorclass (VCL), mainly for templating with double

    @author: minh
    @date:   2016-09-24


*/

#ifndef VECTORF64_H
#define VECTORF64_H


//typedef int64_t Vec1db;
//typedef bool Vec1db;


/*****************************************************************************
*
*          Vec1db: Vector of 1 Booleans for use with Vec1d
*
*****************************************************************************/

class Vec1db {
public:
    bool xmm; // Double vector
    // Default constructor:
    Vec1db() {
    }
    // Constructor to broadcast scalar value:
    Vec1db(bool b) {
        xmm = b;
    }
    // Assignment operator to broadcast scalar value:
    Vec1db & operator = (bool b) {
        *this = Vec1db(b);
        return *this;
    }
private: // Prevent constructing from int, etc.
    Vec1db(int b);
    Vec1db & operator = (int x);
public:
    // Member function to change a single element in vector
    // Note: This function is inefficient. Use load function if changing more than one element
    Vec1db const & insert(uint32_t index, bool value) {
        xmm = value;
        return *this;
    }
    // Member function extract a single element from vector
    bool extract(uint32_t index) const {
        return xmm;
    }
    // Extract a single element. Operator [] can only read an element, not write.
    bool operator [] (uint32_t index) const {
        return extract(index);
    }
    static int size() {
        return 1;
    }
};


/*****************************************************************************
*
*          Operators for Vec1db
*
*****************************************************************************/

// vector operator & : bitwise and
static inline Vec1db operator & (Vec1db const & a, Vec1db const & b) {
    return Vec1db(a.xmm && b.xmm);
}
static inline Vec1db operator && (Vec1db const & a, Vec1db const & b) {
    return Vec1db(a.xmm && b.xmm);
}

// vector operator &= : bitwise and
static inline Vec1db & operator &= (Vec1db & a, Vec1db const & b) {
    a = a & b;
    return a;
}




/*****************************************************************************
*
*          Vec1d: Vector of 1 double precision floating point values
*
*****************************************************************************/

class Vec1d {
public:
    double xmm; // double vector
    // Default constructor:
    Vec1d() {
    }
    // Constructor to broadcast the same value into all elements:
    Vec1d(double d) {
        xmm = d;
    }

    // Member function to load from array (unaligned)
    Vec1d & load(double const * p) {
        xmm = *p;
        return *this;
    }
    // Member function to load from array, aligned by 8
    Vec1d const & load_a(double const * p) {
        xmm = *p;
        return *this;
    }
    // Partial load. Load n elements and set the rest to 0
    Vec1d & load_partial(int n, double const * p) {
        switch (n) {
        case 1:
            xmm = *p; break;
        default:
            xmm = 0.0;
        }
        return *this;
    }
    // Member function to store into array (unaligned)
    void store(double * p) const {
        *p = xmm;
    }
    // Member function to store into array, aligned by 8
    void store_a(double * p) const {
        *p = xmm;
    }

    // cut off vector to n elements. The last 4-n elements are set to zero
    Vec1d & cutoff(int n) {
        if (n == 0)
            xmm = 0.0;
        return *this;
    }

    // Member function to change a single element in vector
    // Note: This function is inefficient. Use load function if changing more than one element
    Vec1d const & insert(uint32_t index, double value) {
        xmm = value;
        return *this;
    };
    // Member function extract a single element from vector
    double extract(uint32_t index) const {
        return xmm;
    }
    // Extract a single element. Use store function if extracting more than one element.
    // Operator [] can only read an element, not write.
    double operator [] (uint32_t index) const {
        return extract(index);
    }

    static int size() {
        return 1;
    }
};

/*****************************************************************************
*
*          Operators for Vec1d
*
*****************************************************************************/

// vector operator + : add element by element
static inline Vec1d operator + (Vec1d const & a, Vec1d const & b) {
    return Vec1d(a.xmm + b.xmm);
}

// vector operator + : add vector and scalar
static inline Vec1d operator + (Vec1d const & a, double b) {
    return a + Vec1d(b);
}
static inline Vec1d operator + (double a, Vec1d const & b) {
    return Vec1d(a) + b;
}

// vector operator += : add
static inline Vec1d & operator += (Vec1d & a, Vec1d const & b) {
    a = a + b;
    return a;
}

// postfix operator ++
static inline Vec1d operator ++ (Vec1d & a, int) {
    Vec1d a0 = a;
    a = a + 1.0;
    return a0;
}

// prefix operator ++
static inline Vec1d & operator ++ (Vec1d & a) {
    a = a + 1.0;
    return a;
}

// vector operator - : subtract element by element
static inline Vec1d operator - (Vec1d const & a, Vec1d const & b) {
    return Vec1d(a.xmm - b.xmm);
}

// vector operator - : subtract vector and scalar
static inline Vec1d operator - (Vec1d const & a, double b) {
    return a - Vec1d(b);
}
static inline Vec1d operator - (double a, Vec1d const & b) {
    return Vec1d(a) - b;
}

// vector operator - : unary minus
// Change sign bit, even for 0, INF and NAN
static inline Vec1d operator - (Vec1d const & a) {
    return Vec1d(-a.xmm);
}

// vector operator -= : subtract
static inline Vec1d & operator -= (Vec1d & a, Vec1d const & b) {
    a = a - b;
    return a;
}

// postfix operator --
static inline Vec1d operator -- (Vec1d & a, int) {
    Vec1d a0 = a;
    a = a - 1.0;
    return a0;
}

// prefix operator --
static inline Vec1d & operator -- (Vec1d & a) {
    a = a - 1.0;
    return a;
}

// vector operator * : multiply element by element
static inline Vec1d operator * (Vec1d const & a, Vec1d const & b) {
    return Vec1d(a.xmm * b.xmm);
}

// vector operator * : multiply vector and scalar
static inline Vec1d operator * (Vec1d const & a, double b) {
    return a * Vec1d(b);
}
static inline Vec1d operator * (double a, Vec1d const & b) {
    return Vec1d(a) * b;
}

// vector operator *= : multiply
static inline Vec1d & operator *= (Vec1d & a, Vec1d const & b) {
    a = a * b;
    return a;
}

// vector operator / : divide all elements by same integer
static inline Vec1d operator / (Vec1d const & a, Vec1d const & b) {
    return Vec1d(a.xmm/b.xmm);
}

// vector operator / : divide vector and scalar
static inline Vec1d operator / (Vec1d const & a, double b) {
    return a / Vec1d(b);
}
static inline Vec1d operator / (double a, Vec1d const & b) {
    return Vec1d(a) / b;
}

// vector operator /= : divide
static inline Vec1d & operator /= (Vec1d & a, Vec1d const & b) {
    a = a / b;
    return a;
}

// vector operator == : returns true for elements for which a == b
static inline Vec1db operator == (Vec1d const & a, Vec1d const & b) {
    return Vec1db(a.xmm == b.xmm);
}

// vector operator != : returns true for elements for which a != b
static inline Vec1db operator != (Vec1d const & a, Vec1d const & b) {
    return Vec1db(a.xmm != b.xmm);
}

// vector operator < : returns true for elements for which a < b
static inline Vec1db operator < (Vec1d const & a, Vec1d const & b) {
    return Vec1db(a.xmm < b.xmm);
}

// vector operator <= : returns true for elements for which a <= b
static inline Vec1db operator <= (Vec1d const & a, Vec1d const & b) {
    return Vec1db(a.xmm <= b.xmm);
}

// vector operator > : returns true for elements for which a > b
static inline Vec1db operator > (Vec1d const & a, Vec1d const & b) {
    return b < a;
}

// vector operator >= : returns true for elements for which a >= b
static inline Vec1db operator >= (Vec1d const & a, Vec1d const & b) {
    return b <= a;
}

// General arithmetic functions, etc.

// Horizontal add: Calculates the sum of all vector elements.
static inline double horizontal_add (Vec1d const & a) {
    return a.xmm;
}

// function max: a > b ? a : b
static inline Vec1d max(Vec1d const & a, Vec1d const & b) {
    return max(a.xmm,b.xmm);
}

// function min: a < b ? a : b
static inline Vec1d min(Vec1d const & a, Vec1d const & b) {
    return min(a.xmm,b.xmm);
}


// function abs: absolute value
// Removes sign bit, even for -0.0f, -INF and -NAN
static inline Vec1d abs(Vec1d const & a) {
    return Vec1d(fabs(a.xmm));
}

// function log: logarithm
// Removes sign bit, even for -0.0f, -INF and -NAN
static inline Vec1d log(Vec1d const & a) {
    return Vec1d(log(a.xmm));
}

// Fused multiply and add functions

// Multiply and add
static inline Vec1d mul_add(Vec1d const & a, Vec1d const & b, Vec1d const & c) {
    return a * b + c;
}

// Multiply and subtract
static inline Vec1d mul_sub(Vec1d const & a, Vec1d const & b, Vec1d const & c) {
    return a * b - c;
}

// Multiply and inverse subtract
static inline Vec1d nmul_add(Vec1d const & a, Vec1d const & b, Vec1d const & c) {
    return c - a * b;
}


/*****************************************************************************
*
*          Horizontal Boolean functions
*
*****************************************************************************/

// horizontal_and. Returns true if all bits are 1
static inline bool horizontal_and (Vec1db const & a) {
    return a.xmm;
}

// horizontal_or. Returns true if at least one bit is 1
static inline bool horizontal_or (Vec1db const & a) {
    return a.xmm;
}

// instances of exp_d template
static inline Vec1d exp(Vec1d const & x) {
    return Vec1d(exp(x.xmm));
}



#endif //VECTORF64_H
