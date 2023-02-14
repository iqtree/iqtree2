//
//  distancematrix.h - Template classes, Matrix and SquareMatrix,
//                     used in distance matrix tree construction
//                     algorithms.
//  Copyright James Barbetti (2020)
//
//  LICENSE:
//* This program is free software; you can redistribute it and/or modify
//* it under the terms of the GNU General Public License as published by
//* the Free Software Foundation; either version 2 of the License, or
//* (at your option) any later version.
//*
//* This program is distributed in the hope that it will be useful,
//* but WITHOUT ANY WARRANTY; without even the implied warranty of
//* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//* GNU General Public License for more details.
//*
//* You should have received a copy of the GNU General Public License
//* along with this program; if not, write to the
//* Free Software Foundation, Inc.,
//* 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//

#ifndef distancematrix_h
#define distancematrix_h

#include <iostream>    //for std::istream
#include <sstream>     //for std::stringstream
#include <vector>      //for std::vector
#include <math.h>      //for floor()
#if USE_GZSTREAM
#include <utils/gzstream.h>  //for igzstream and pigzstream
#else
#include <fstream>     //for std::ifstream
#endif
#include <utils/safe_io.h>   //for safeGetTrimmedLineAsStream
#include <utils/progress.h>  //for progress_display

#define MATRIX_ALIGNMENT 64
    //MUST be a power of 2 (else x & MATRIX_ALIGNMENT_MASK
    //would be no good and x % MATRIX_ALIGNMENT would be needed).
    //Assumed: sizeof(NJFloat) divides MATRIX_ALIGNMENT
    //Vectorized versions run faster (particularly on older
    //hardware), if rows are aligned.
#define MATRIX_ALIGNMENT_MASK (MATRIX_ALIGNMENT - 1)

template <class P> inline P* matrixAlign(P* p) {
    //If we've got an array that mighnt't be MATRIX_ALIGNMENT-byte
    //aligned, but we've got MATRIX_ALIGNMENT/sizeof(P) extra items
    //in it, we can point to the first item in the array that *is*
    //MATRIX_ALIGNMENT-byte aligned.  This function returns the
    //address of that item.
    uintptr_t address = reinterpret_cast<uintptr_t>(p);
    auto offset = address & MATRIX_ALIGNMENT_MASK;
    if (0<offset) {
        return p + (MATRIX_ALIGNMENT - offset)/sizeof(P);
    } else {
        return p;
    }
}

/**
 * @brief  A matrix, stored in row major order (but with
 *         rows not necessarily adjacent to one another
 *         in memory, if columns have been deleted).
 * @tparam T the type of the cells in the matrix
 * @note   I resorted to declaring the data, rows, and
 *         rowTotals members public, because of problems
 *         I had accessing them from BoundingMatrix. - James B.
 */
template <class T=double> class Matrix
{
protected:
    intptr_t row_count;    //the current number of rows.
    intptr_t column_count; //the current number of columns (need not
                           //be the same; but will be in SquareMatrix).
    intptr_t shrink_r;     //if row_count reaches *this*, pack the array

public:
    typedef T cell_type;
    T*     data;           //The single one-big-block memory allocation
                           //for all of the cells in the matrix
    T**    rows;           //An array (of size at least row_count)
                           //of pointers to rows of the matrix,
                           //as it is now.
    
    Matrix() : row_count(0), column_count(0), shrink_r(0)
             , data(nullptr), rows(nullptr) {
    }
    Matrix(intptr_t r, intptr_t c): Matrix() {
        setDimensions( r, c);
    }
    /**
     * @brief Release any allocated memory
     */
    void actual_clear() {
        delete [] data;
        delete [] rows;
        data         = nullptr;
        rows         = nullptr;
        row_count    = 0;
        column_count = 0;
    }
    virtual void clear() {
        actual_clear();
    }
    /**
     * @brief  get a const reference to the pointer to the const
     *         first element in the matrix, in the (r)th row.
     * @param  r  
     * @return const T*& 
     * @note   for performance rasons, there is no range checking.
     */
    const T*& getRow(size_t r) const {    
        return rows[r];
    }
    /**
     * @brief  get a reference to the pointer to the 
     *         first element in the matrix, in the (r)th row.
     * @param  r  
     * @return const T*& 
     * @note   for performance rasons, there is no range checking.
     */
    T*& getRow(size_t r) {
        return rows[r];
    }
    /**
     * @brief Append the values, found in the (c)th column of the matrix
     *        to a (supplied) vector of the cell type, T, of the matrix.
     * @param c    the column
     * @param dest a reference to a vector of T.
     */
    void appendColumnToVector(size_t c, std::vector<T>& dest) {
        dest.resize(row_count);
        T* dest_ptr = dest.data();
        //Not parallelized: called from within other for-loops
        //that are already parallelized.
        for (intptr_t r=0; r<row_count; ++r) {
            dest_ptr[r] = rows[r][c];
        }
    }
    /**
     * @brief provides read-only access to a matrix element, by row and column
     * @param r row number
     * @param c column number
     * @return  const T& - read-only reference to the element in row r
     *          and column c.
     * @note    for performance reasons, no range checking is performed.
     */
    const T& cell (size_t r, size_t c) const {
        return rows[r][c];
    }
    /**
     * @brief provides read/write access to a matrix element, by row and column
     * @param r row number
     * @param c column number
     * @return  T& reference to the element in row r and column c
     * @note    for performance reasons, no range checking is performed.
     */
    T& cell (size_t r, size_t c) {
        return rows[r][c];
    }
    /**
     * @brief Set the dimensions (number of rows, number of columns)
     *        of the matrix.
     * @param r the desired number of rows
     * @param c the desired number of columns
     * @note  All elements in the matrix will be zeroed. Existing 
     *        memory allocated for the matrix (if any) will be released
     *        before the matrix is resized.
     * @note  potentially throws an out-of-memory exception.
     *        if this happens, the matrix will be empty.
     * @note  It is assumed that sizeof(T) divides MATRIX_ALIGNMENT.     
     * @note  Allocations are aligned, as are rows, to MATRIX_ALIGNMENT
     *        byte boundaries (this means that there may be unused 
     *        memory between rows, if column_count*sizeof(T) is not 
     *        divisible by MATRIX_ALIGNMENT.
     */
    void setDimensions(size_t r, size_t c) {
        actual_clear();
        if (0==c || 0==r) {
            return;
        }
        try {
            size_t w     = widthNeededFor(c);
            row_count    = r;
            column_count = c;
            shrink_r     = (r+r)/3;
            if (shrink_r<100) {
                shrink_r=0;
            }
            data        = new T[ r * w + MATRIX_ALIGNMENT/sizeof(T)];
            try {
                rows    = new T*[r];
            }
            catch (...) {
                delete [] data;
                data = nullptr;
                throw;
            }
            T *rowStart = matrixAlign(data);
            for (intptr_t row=0; row<row_count; ++row) {
                rows[row]      = rowStart;
                rowStart    += w;
            }
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for (intptr_t row=0; row<row_count; ++row) {
                zeroRow(row);
            }
        }
        catch (...) {
            clear();
            throw;
        }
    }
    /**
     * @brief Set *this matrix to match another. Used in assignment operators
     *        and copy constructors.
     * @param rhs The matrix that *this is to match.
     * @note  any memory currently allocated for *this matrix will be 
     *        deallocated.  If an out-of-memory (or other) exception is
     *        thrown, *this will be an empty matrix.
     * @note  the copying of the elements from the source matrix 
     *        is parallelized across rows.
     */
    void assign(const Matrix<T> &rhs ) {
        if (this==&rhs) {
            return;
        }
        setDimensions(rhs.row_count, rhs.column_count);
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (intptr_t r=0; r<row_count; ++r) {
            T *             destRow      = rows[r];
            T const *       sourceRow    = rhs.rows[r];
            T const * const endSourceRow = sourceRow + column_count;
            for (; sourceRow<endSourceRow; ++destRow, ++sourceRow) {
                *destRow = *sourceRow;
            }
        }
    }
    /**
     * @brief  returns width, rounded up so that each row will
     *         have a starting address that is MATRIX_ALIGNMENT-byte 
     *         aligned.  
     * @param  width the width (to round up)
     * @return size_t the width (rounded up)
     * @note   it is assumed that , if sizeof(T)<MATRIX_ALIGNMENT,
     *         then sizeof(T) evenly divides MATRIX_ALIGNMENT.
     */
    size_t widthNeededFor(size_t width) {
        if (MATRIX_ALIGNMENT<=sizeof(T)) {
            return width;
        }
        size_t leftOver  = (width * sizeof(T)) & MATRIX_ALIGNMENT_MASK;
        if (leftOver==0) {
            return width;
        }
        return width + (MATRIX_ALIGNMENT-leftOver) / sizeof(T);
    }
    /**
     * @brief Zero all the elements in a single row of the matrix
     * @param r 
     * @note  Not prallelized (it's called from code that parallelizes
     *        across rows).
     */
    void zeroRow(size_t r) {
        T* rowStart = rows[r];
        T* rowStop  = rowStart + column_count;
        for (T* rowZap=rowStart; rowZap<rowStop; ++rowZap) {
            *rowZap = 0;
        }
    }
    /**
     * @brief remove a column (by copying the content of the
     *        last column, into the column, and then removing 
     *        the last column).
     * @param c the column to remove
     * @note  c is not range checked
     * @note  parallelized across rows.
     */
    virtual void removeColumn(intptr_t c) {
        if ( c + 1 < column_count) {
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for (intptr_t r=0; r<row_count; ++r) {
                T* rowData = rows[r];
                rowData[c] = rowData[column_count-1]; //U-R
            }
        }
        --column_count;
    }
    /**
     * @brief remove a row (by copying the content of the
     *        last row, into the row, and then removing 
     *        the last row).
     * @param row_to_remove the row to remove
     * @note  row_to_remove is NOT range checked
     * @note  parallelized across columns
     * @note  This could also have been implemented more simply, by setting
     *        rows[rowNum] = rows[n] (instead of copying the content of the
     *        last row).  Originally it was.  But in practice row copying 
     *        seemed to result in faster UPGMA and NJ implementations (!).
     *        Not what I eas expecting, if I'm honest. -James B.
     * @note  The packing code (that packs the array rows into a smaller
     *        part of the array, when the number of rows falls), should
     *        really be executed as part of removeColumn!  Because the
     *        point is to move rows that have "drifted apart" in memory,
     *        due to column deletions, closer together, to improve cache
     *        (or memory page) utilisation.
     */
    virtual void removeRow(size_t row_to_remove) {
        --row_count;
        T*       destRow   = rows[row_to_remove];
        const T* sourceRow = rows[row_count];
        rows[row_count] = nullptr;
        if (destRow!=sourceRow) {
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for (intptr_t c=0; c<column_count; ++c) {
                destRow[c] = sourceRow[c];
            }
        }

        if ( row_count == shrink_r && 0 < shrink_r) {
            //Move the data in the array closer to the front.
            //This also helps (but: only very slightly. 5%ish?).
            size_t w = widthNeededFor(column_count);
            destRow  = data;
            for (intptr_t r=1; r<row_count; ++r) {
                destRow += w;
                sourceRow = rows[r];
                #ifdef _OPENMP
                #pragma omp parallel for
                #endif
                for (intptr_t c=0; c<column_count; ++c) {
                    destRow[c] = sourceRow[c];
                }
                rows[r] = destRow;
            }
            shrink_r    = (row_count+row_count)/3;
            if (shrink_r<100) shrink_r=0;
        }
    }
    /**
     * @brief load distances, from a flat array of doubles 
     *        (in row-major order), into the matrix.
     * @param matrix the flat array (it assumed that the array contains
     *               a matrix of - at least - row_count - rows, each of
     *               exactly column_count columns, stored in row-major
     *               order).
     * @note  it is assumed that a const double, x can be cast to T
     *        via a C-Style cast: (T) x.
     * @note  it would be better if this didn't make use of C-style casts.
     *        The problem, though is that I use SquareMatrix (which is
     * .      a susclass) to represent I-matrices (matrices of ints
     *        that indentify clusters).
     */     
    virtual void loadDistancesFromFlatArray(const double* matrix) {
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (intptr_t row=0; row<row_count; ++row) {
            const double* sourceStart = matrix + row * column_count;
            const double* sourceStop  = sourceStart + column_count;
            T*            dest        = rows[row];
            for (const double* source=sourceStart;
                 source<sourceStop; ++source, ++dest ) {
                *dest = (T) *source;
            }
        }
    }
};

/**
 * @brief  A Matrix where row_count is always the same as column_count
 * @tparam T the cell type
 * @note   This class can be, and is, used for variance (V) and rectangular
 *         sorted distance (S) and index (I) matrices, not just square 
 *         distance (D) matrices. Lines that access the upper-right triangle
 *         of the matrix are tagged with U-R.
 */
template <class T=double> class SquareMatrix: public Matrix<T>
{
public:
    typedef Matrix<T> super;
public:
    using super::rows;
    using super::removeColumn;
    using super::removeRow;
    using super::row_count;
    using super::column_count;
    T*     rowTotals; //The U vector (e.g. in the NJ, UNJ, and UPGMA algorithms)
    
    /**
     * @brief  return the rank of the matrix
     * @return size_t the rank (the number of rows, and also of columns)
     */
    size_t getSize() {
        return row_count;
    }
    /**
     * @brief set the rank of the matrix (and zero its content)
     * @param rank 
     * @note  implemented in terms of super::setDimensions
     */
    virtual void setSize(intptr_t rank) {
        delete [] rowTotals;
        rowTotals = nullptr;
        super::setDimensions(rank,rank);
        rowTotals = new T[rank];
        for (int r=0; r<rank; ++r) {
            rowTotals[r] = (T)0.0;
        }
    }
    /**
     * @brief assign *this as a copy of another SquareMatrix
     * @param rhs the SquareMatrix to copy
     */
    void assign(const SquareMatrix& rhs) {
        if (this==&rhs) {
            return;
        }        
        delete [] rowTotals;
        rowTotals = nullptr;
        super::assign(rhs);
        try {
            rowTotals = new T[row_count];
        }
        catch (...) {
            super::clear();
            throw;
        }
        for (intptr_t r=0; r<row_count; ++r) {
            rowTotals[r] = rhs.rowTotals[r];
        }
    }
    SquareMatrix(): super(), rowTotals(nullptr) {
    }
    SquareMatrix(const SquareMatrix& rhs): rowTotals(nullptr) {
        assign(rhs);
    }
    virtual ~SquareMatrix() {
        actual_clear();
    }
    void actual_clear() {
        delete [] rowTotals;
        rowTotals = nullptr;
    }
    virtual void clear() override {
        super::clear();
        actual_clear();
    }
    SquareMatrix& operator=(const SquareMatrix& rhs) {
        assign(rhs);
        return *this;
    }
    /**
     * @brief Recalculate row totals
     * @note  Row total calculation is parallelized across rows.
     * @note  Lines that access the lower left triangle are tagged L-L.
     *        Lines that access the upper right traingle are tagged U-R.
     */
    virtual void calculateRowTotals() const {
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (intptr_t r=0; r<row_count; ++r) {
            T total = (T)0.0;
            const T* rowData = rows[r];
            for (intptr_t c=0; c<r; ++c) {
                total += rowData[c]; //L-L
            }
            for (intptr_t c=r+1; c<column_count; ++c) {
                total += rowData[c]; //U-R
            }
            rowTotals[r] = total;
        }
    }

    /**
     * @brief recalculate total for row, a, excluding column b (a<=b).
     * @note  At present, this member function is not called
     *        (it turns out to be faster to "roll-in" row total
     *        recalculations into existing column for-loops elsewhere).
     * @param a the row to recalculate for
     * @param b the column to exclude
     */
    void recalculateTotalForOneRow(size_t a, size_t b) {
        T replacementRowTotal = 0;
        auto rowA = rows[a];
        #ifdef _OPENMP
        #pragma omp parallel for reduction(+:replacementRowTotal)
        #endif
        for (size_t i=0; i<column_count; ++i) {
            replacementRowTotal += rowA[i];
        }
        if (a<=column_count) {
            replacementRowTotal -= rowA[a];
        }
        if (a!=b && b<=column_count) {
            replacementRowTotal -= rowA[b];
        }
        rowTotals[a] = replacementRowTotal;
    }
    /**
     * @brief Remove row (and matching column) from a square matrix, 
     *        by swapping the last row (and column) into its place.
     * @param rowNum the row (and column) to remove
     */
    void removeRowAndColumn(size_t rowNum)  {
        removeColumn(rowNum);
        removeRowOnly(rowNum);
    }
    void removeRowOnly(size_t rowNum) {
        //Remove row from a rectangular matrix.
        //Don't touch the columns in the row
        //(Used for the S and I matrices in BIONJ).
        removeRow(rowNum);
        rowTotals[rowNum] = rowTotals[row_count];
    }
    virtual void addCluster(const std::string &name) {
    }
};

template <class M> inline void copyLowerTriangleOfMatrixToUpper(M& matrix) {
    intptr_t n = matrix.getSize();
    for (intptr_t r = 0; r < n; ++r) {
        for (intptr_t c = r + 1; c < n; ++c) {
            matrix.cell(r, c) = matrix.cell(c, r);
        }
    }
}

/**
 * @brief  Load (part of) a row of a matrix, from a line of a stream
 * @tparam F the stream or fyle type (usually stringstream)
 * @tparam M the matrix type
 * @param line    the stream from which a line is to be read
 * @param matrix  the matrix that is being read
 * @param square  true if loading from a square matrix
 *                (AND it is assumed that the matrix will be symmetric
 *                 around its diagonal, with the upper-right triangle
 *                 the transpose of the lower-left triangle)
 * @param r       the row number
 * @param cStart  the starting column number 
 * @param cStop   the stopping column number
 * @return intptr_t the column number where the line or file ran out
 */
template <class F=std::stringstream, class M> intptr_t loadPartialRowFromLine
    (F& line,  M& matrix, bool square,
     intptr_t r, intptr_t cStart, intptr_t cStop) {
    intptr_t c      = cStart;
    for (; line.tellg() != -1 && c < cStop; ++c) {
        line >> matrix.cell(r, c);
        //Ensure matrix is symmetric (as it is read!)
        //Note: I'd rather throw if there's a disagreement!
        //But I want backwards compatibility with BIONJ2009.
        if (square && c < r && matrix.cell(r, c) != matrix.cell(c, r)) {
            //If loading from a square matrix file,
            //and in lower triangle, and the corresponding entry
            //in the upper triangle was different, average them both.
            typename M::cell_type v = (matrix.cell(r, c) + matrix.cell(c, r)) * (typename M::cell_type)0.5;
            matrix.cell(c, r) = v; //U-R
            matrix.cell(r, c) = v;
        }
    }
    return c;
}

/**
 * @brief  if it hasn't been determined already, determine - from a single matrix
 *         row that has just been read in, whether the matrix representation 
 *         in the file being read is square, upper triangle, or lower triangle.
 * 
 * @tparam M the matrix type
 * @tparam P the type used to report progress (expected to implement
 *           hide() and show() members).
 * @param r  the row number (of a row that has just been loaded)
 * @param c  the number of columns that were read from the row
 * @param cStop 
 * @param square reference: true if the in-file array is believed to be square 
 *               will be set to false if it is determined the in-file matrix
 *               data is assumed to be in a lower-triangle or an upper-triangle
 *               format.  Otherwise, left as is.
 * @param lower  reference: will be set to true if lower-triangle format detected.
 * @param upper  reference: will be set to false if upper-triangle format detected.
 * @param matrix   a reference to the matrix being loaded.
 * @param progress a status-reporting object that implements hide() and show(). 
 * @note  this really only does anything for the first row read.
 *        if no column were read, we assume lower-triangle format.
 *        if row_count-1 were read, we assume upper-triangle format.
 *        otherwise, we continue to assume square format.
 * @note  if upper-triangle format is detected, the items read into
 *        the first row of the matrix will be offset (to the left)
 *        by one place, because there won't have been data for the 
 *        matrix entry, in the first row, that belongs on the diagonal.
 *        So, they'll all have to be shuffled one place to the right.
 * @note  this function exists at all only because lizard complains
 *        about other (higher-level) functions if I don't move this
 *        logic out into a lower level function.  -James B.
 */
template <class M, class P> 
void recognizeFormat(intptr_t r, intptr_t c, intptr_t cStop,
                     bool& square, bool& lower, bool& upper, 
                     M& matrix, P& progress) {
    const char* format_comment = "";
    if (r!=0) {
        return;
    }
    if (!square) {
        return;
    }
    if (c == 0) {
        //Implied lower-triangle format
        square         = false;
        lower          = true;
        format_comment = "Input appears to be in lower-triangle format" ;
    }
    else if (c + 1 == cStop) {
        //Implied upper-triangle format
        square         = false;
        upper          = true;
        format_comment = "Input appears to be in upper-triangle format" ;
        //Shuffle everything read in this row, one place to the right,
        //to make room for the zero to be written onto the top-left cell
        //(which is on the diagonal).
        for (size_t shift = cStop - 1; 0 < shift; --shift) {
            matrix.cell(0, shift) = matrix.cell(0, shift - 1);
        }
        matrix.cell(0, 0) = 0;
    }
    else {
        return;
    }

    #if USE_PROGRESS_DISPLAY
        progress.hide();
    #endif

    std::cout << format_comment << std::endl;

    #if USE_PROGRESS_DISPLAY
        progress.show();
    #endif
}

template <class P> inline void incrementProgress(P& progress) {
    ++progress;
}

/**
 * @brief  Load sequence names (and a distance matrix) from an open file
 * @tparam F the file type
 * @tparam M the matrix type
 * @param  in             reference to the input file (of type F)
 * @param  reportProgress true, if progress is to be reported. false otherwise.
 * @param  matrix         reference to the matrix to load (of type M)
 * @note   sequence names may not contain white space characters
 *         (spaces, tabs, or linefeeds).
 * @note   throws strings (not pretty! should be throwing std::exception!)
 *         if it can't make sense of the file format (specifically, if it
 * .       sees more columns than it expects, in any row of the input file)
 */
template <class F=std::stringstream, class M> 
    void loadDistanceMatrixFromOpenFile
        (F& in, bool reportProgress, M& matrix) {
    intptr_t rank;
    bool lower  = false;
    bool upper  = false;
    bool square = true;
    std::stringstream first_line;
    safeGetTrimmedLineAsStream<F>(in, first_line);
    first_line >> rank;
    matrix.setSize(rank);
#if USE_PROGRESS_DISPLAY
    const char* taskDescription = reportProgress ? "Loading distance matrix" : "";
    progress_display progress(rank, taskDescription, "loaded", "row");
#else
    intptr_t progress=0;
#endif
    for (intptr_t r = 0; r < rank; ++r) {
        std::stringstream line;
        safeGetTrimmedLineAsStream(in, line);
        std::string name;
        line >> name;
        matrix.addCluster(name);
        if (upper) {
            //Copy column r from upper triangle (above the diagonal)
            //to the left part of row r (which is in the lower triangle).
            for (intptr_t c = 0; c < r; ++c) {
                matrix.cell(r, c) = matrix.cell(c, r);
            }
        }
        matrix.cell(r, r) = 0;
        intptr_t cStart = (upper) ? (r + 1) : 0;
        intptr_t cStop  = (lower) ? r : rank;
        intptr_t c      = loadPartialRowFromLine(line, matrix, square, 
                                                 r,    cStart, cStop);

        if (line.tellg() == -1 && c < cStop) {
            recognizeFormat(r,c,cStop,square,lower,upper,matrix,progress);
        }
        else if (line.tellg() != -1) {
            std::stringstream problem;
            problem 
                << "Expected to see columns"
                << " [" << (cStart + 1) << ".." << cStop << "]"
                << " in row " << r << " of the distance matrix,"
                << " but there were more distances.";
            throw problem.str();
        }
        incrementProgress(progress);
    }
    if (lower) {
        copyLowerTriangleOfMatrixToUpper(matrix);
    }
}

/**
 * @brief  Loads squarence names, and distance matrix from a file
 *         (when given the path to the file)
 * @tparam M the distance
 * @param  distanceMatrixFilePath the path of the file containing the matrix
 * @param  reportProgress         true if progress is to be reported. Else, false.
 * @param  matrix                 the matrix to load
 * @return true   if it worked.
 * @return false  if it failed (an error message explaining what the problem
 *                was will have been written to std::cerr).
 * @note   Will, if USE_GZSTREAM is defined, and non-zero, use an igzstream
 *         (gzip-format-capable file stream).  Will otherwise use a 
 *         std::ifstream (standard template library input filestream).
 * @note   Does not write a message to standard output (or standard error)
 *         if the matrix is not symmetric.
 */
template <class M> bool loadDistanceMatrixInto
    (const std::string& distanceMatrixFilePath, 
     bool reportProgress, M& matrix) {        
    #if USE_GZSTREAM
    igzstream     in;
    #else
    std::ifstream in;
    #endif        
    try {
        in.exceptions(std::ios::failbit | std::ios::badbit);
        in.open(distanceMatrixFilePath.c_str(), std::ios_base::in);
        loadDistanceMatrixFromOpenFile(in, reportProgress, matrix);
        in.close();
        return true;
    } catch (std::ios::failure &) {
        std::cerr << "Load matrix failed: IO error"
            << " reading file: " << distanceMatrixFilePath << std::endl;
        return false;
    } catch (const char *str) {
        std::cerr << "Load matrix failed: " << str << std::endl;
        return false;
    } catch (const std::string &str) {
        std::cerr << "Load matrix failed: " << str << std::endl;
        return false;
    }
    return true;
}

#endif /* distancematrix_h */
