//
//  distancematrix.h
//  iqtree
//
//  Created by James Barbetti on 12/8/20.
//

#ifndef distancematrix_h
#define distancematrix_h

#include "gzstream.h"                //for igzstream
#include <fstream>
#include <iostream>                  //for std::istream
#include <sstream>                   //for std::stringstream
#include <vector>                    //for std::vector
#include "io.h"                      //for safeGetTrimmedLineAsStream
#include "progress.h"                //for progress_display


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
    if (0<offset)
    {
        return p + (MATRIX_ALIGNMENT - offset)/sizeof(P);
    } else {
        return p;
    }
}

template <class T=double> class Matrix
{
protected:
    size_t row_count;
    size_t column_count;
    size_t shrink_r; //if row_count reaches *this*, pack the array
public:
    T*     data;
    T**    rows;
    
    Matrix() : row_count(0), column_count(0), shrink_r(0)
             , data(nullptr), rows(nullptr) {
    }
    virtual void clear() {
        delete [] data;
        delete [] rows;
        data         = nullptr;
        rows         = nullptr;
        row_count    = 0;
        column_count = 0;
    }
    const T*& getRow(size_t r) const {
        return rows[r];
    }
    T*& getRow(size_t r) {
        return rows[r];
    }
    void appendColumnToVector(size_t c, std::vector<T>& dest) {
        for (size_t r=0; r<row_count; ++r) {
            dest.emplace_back( rows[r][c] );
        }
    }
    const T& cell (size_t r, size_t c) const {
        return rows[r][c];
    }
    T& cell (size_t r, size_t c) {
        return rows[r][c];
    }
    void setDimensions(size_t r, size_t c) {
        clear();
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
            rows        = new T*[r];
            T *rowStart = matrixAlign(data);
            for (size_t row=0; row<row_count; ++row) {
                rows[row]      = rowStart;
                rowStart    += w;
            }
            #pragma omp parallel for
            for (size_t row=0; row<row_count; ++row) {
                zeroRow(row);
            }
        }
        catch (...) {
            clear();
            throw;
        }
    }
    void assign(const Matrix<T> &rhs ) {
        if (this==&rhs) {
            return;
        }
        setDimensions(rhs.row_count, rhs.column_count);
        #pragma omp parallel for
        for (size_t r=0; r<row_count; ++r) {
            T *             destRow      = rows[r];
            T const *       sourceRow    = rhs.rows[r];
            T const * const endSourceRow = sourceRow + column_count;
            for (; sourceRow<endSourceRow; ++destRow, ++sourceRow) {
                *destRow = *sourceRow;
            }
        }
    }
    
    size_t widthNeededFor(size_t width) {
        //
        //returns width, rounded up so that each row will
        //have a starting address that is MATRIX_ALIGNMENT-byte aligned.
        //
        if (MATRIX_ALIGNMENT<=sizeof(T))
        {
            return width;
        }
        size_t leftOver  = (width * sizeof(T)) & MATRIX_ALIGNMENT_MASK;
        if (leftOver==0) {
            return width;
        }
        return width + (MATRIX_ALIGNMENT-leftOver) / sizeof(T);
    }
    void zeroRow(size_t r) {
        T* rowStart = rows[r];
        T* rowStop  = rowStart + column_count;
        for (T* rowZap=rowStart; rowZap<rowStop; ++rowZap) {
            *rowZap = 0;
        }
    }
    virtual void removeColumn(size_t c) {
        if ( c + 1 < column_count) {
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for (size_t r=0; r<row_count; ++r) {
                T* rowData = rows[r];
                rowData[c] = rowData[column_count-1]; //U-R
            }
        }
        --column_count;
    }
    virtual void removeRow(size_t row_to_remove) {
        --row_count;
        //was rows[rowNum] = rows[n];... but let's copy
        //instead.  On average it seems (very slightly) faster.
        T*       destRow   = rows[row_to_remove];
        const T* sourceRow = rows[row_count];
        rows[row_count] = nullptr;
        if (destRow!=sourceRow) {
            #pragma omp parallel for
            for (size_t c=0; c<column_count; ++c) {
                destRow[c] = sourceRow[c];
            }
        }
        if ( row_count == shrink_r && 0 < shrink_r) {
            //Move the data in the array closer to the front.
            //This also helps (but: only very slightly. 5%ish?).
            size_t w = widthNeededFor(column_count);
            destRow  = data;
            for (size_t r=1; r<row_count; ++r) {
                destRow += w;
                sourceRow = rows[r];
                #pragma omp parallel for
                for (size_t c=0; c<column_count; ++c) {
                    destRow[c] = sourceRow[c];
                }
                rows[r] = destRow;
            }
            shrink_r    = (row_count+row_count)/3;
            if (shrink_r<100) shrink_r=0;
        }
    }
};

template <class T=double> class SquareMatrix: public Matrix<T>
{
    //Note 1: This is a separate class so that it can be
    //        used for square variance (V) and rectangular
    //        sorted distance (S) and index (I) matrices,
    //        not just square distance (D) matrices.
    //        Lines that access the upper-right triangle
    //        of the matrix are tagged with U-R.
    //Note 2: I resorted to declaring the data, rows, and
    //        rowTotals members public, because of problems
    //        I had accessing them from BoundingMatrix.
    //Note 3: Perhaps there should be separate SquareMatrix
    //        and RectangularMatrix classes?
public:
    typedef Matrix<T> super;
public:
    using super::rows;
    using super::removeColumn;
    using super::removeRow;
    using super::row_count;
    using super::column_count;
    T*     rowTotals; //The U vector
    
    size_t getSize() {
        return row_count;
    }
    virtual void setSize(size_t rank) {
        super::setDimensions(rank,rank);
        delete [] rowTotals;
        rowTotals = new T[rank];
        for (int r=0; r<rank; ++r) {
            rowTotals[r] = 0.0;
        }
    }
    void assign(const SquareMatrix& rhs) {
        if (this==&rhs) {
            return;
        }
        super::assign(rhs);
        delete [] rowTotals;
        rowTotals = new T[row_count];
        for (size_t r=0; r<row_count; ++r) {
            rowTotals[r] = rhs.rowTotals[r];
        }
    }
    SquareMatrix(): super(), rowTotals(nullptr) {
    }
    SquareMatrix(const SquareMatrix& rhs): rowTotals(nullptr) {
        assign(rhs);
    }
    virtual ~SquareMatrix() {
        clear();
    }
    virtual void clear() {
        super::clear();
        delete [] rowTotals;
        rowTotals = nullptr;
    }
    SquareMatrix& operator=(const SquareMatrix& rhs) {
        assign(rhs);
        return *this;
    }
    virtual void calculateRowTotals() const {
        //Note: Although this isn't currently in use,
        //it's been kept, in case it is needed
        //(after, say, every 200 iterations of
        //neighbour-joining) to deal with accumulated
        //rounding error.  It might be.
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (size_t r=0; r<row_count; ++r) {
            T total(0.0);
            const T* rowData = rows[r];
            for (size_t c=0; c<r; ++c) {
                total += rowData[c];
            }
            for (size_t c=r+1; c<column_count; ++c) {
                total += rowData[c]; //U-R
            }
            rowTotals[r] = total;
        }
    }
    void recalculateTotalForOneRow(size_t a, size_t b) {
        //recalculate total for row, a, excluding
        //column b (a<=b).
        //Note: At present, this member function is not called
        //(it turns out to be faster to "roll-in" row total
        // recalculations into existing column for-loops elsewhere).
        //
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

    void removeRowAndColumn(size_t rowNum)  {
        //Remove row (and matching column) from a
        //square matrix, by swapping the last row
        //(and column) into its place.
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
    
template <class M> bool loadDistanceMatrixInto
    (const std::string distanceMatrixFilePath, bool reportProgress, M& matrix) {
    size_t rank;
    igzstream in;
    bool lower  = false;
    bool upper  = false;
    bool square = true;
    try {
        in.exceptions(std::ios::failbit | std::ios::badbit);
        in.open(distanceMatrixFilePath.c_str(), std::ios_base::in);
        std::stringstream first_line;
        safeGetTrimmedLineAsStream(in, first_line);
        first_line >> rank;
        matrix.setSize(rank);
        const char* taskDescription = reportProgress ? "Loading distance matrix" : "";
        progress_display progress(rank, taskDescription, "loaded", "row");
        for (size_t r=0; r<matrix.getSize(); ++r) {
            std::stringstream line;
            safeGetTrimmedLineAsStream(in, line);
            std::string name;
            line >> name;
            matrix.addCluster(name);
            if (upper) {
                //Copy column r from upper triangle (above the diagonal)
                //to the left part of row r (which is in the lower triangle).
                for (size_t c=0; c<r; ++c) {
                    matrix.cell(r,c) = matrix.cell(c,r);
                }
                matrix.cell(r,r) = 0;
            }
            int cStart = (upper) ? (r+1) : 0;
            int cStop  = (lower) ? r     : matrix.getSize();
            size_t c = cStart;
            for (; line.tellg()!=-1 && c<cStop; ++c) {
                line >> matrix.cell(r,c);
                //Ensure matrix is symmetric (as it is read!)
                //Note: I'd rather throw if there's a disagreement!
                //But I want backwards compatibility with BIONJ2009.
                if (square && c<r && matrix.cell(r,c) != matrix.cell(c,r)) {
                    //If loading from a square matrix file,
                    //and in lower triangle, and the corresponding entry
                    //in the upper triangle was different, average them both.
                    auto v = ( matrix.cell(r,c) + matrix.cell(c,r) ) * 0.5;
                    matrix.cell(c,r) = v; //U-R
                    matrix.cell(r,c) = v;
                }
            }
            if (line.tellg()==-1 && c<cStop)
            {
                if (square && r==0 && c==0) {
                    //Implied lower-triangle format
                    square = false;
                    lower  = true;
                    progress.hide();
                    std::cout << "Input appears to be in lower-triangle format" << std::endl;
                    progress.show();
                }
                else if (square && r==0 && c+1==cStop) {
                    //Implied upper-triangle format
                    square = false;
                    upper  = true;
                    progress.hide();
                    std::cout << "Input appears to be in upper-triangle format" << std::endl;
                    progress.show();
                    for (size_t shift=cStop-1; 0<shift; --shift) {
                        matrix.cell(0,shift) = matrix.cell(0, shift-1);
                    }
                    matrix.cell(0,0)=0;
                }
            }
            else if (line.tellg()!=-1) {
                std::stringstream problem;
                problem << "Expected to see columns [" << (cStart+1) << ".." << cStop << "]"
                    << " in row " << r << " of the distance matrix, but there were more distances.";
                throw problem.str();
            }
            ++progress;
        }
        if (lower) {
            size_t n = matrix.getSize();
            for (size_t r=0; r<n; ++r) {
                for (size_t c=r+1; c<n; ++c) {
                    matrix.cell(r,c) = matrix.cell(c,r);
                }
            }
        }
        in.close();
    } catch (std::ios::failure &) {
        std::cerr << "Load matrix failed: IO error"
            << " reading file: " << distanceMatrixFilePath << std::endl;
        return false;
    } catch (const char *str) {
        std::cerr << "Load matrix failed: " << str << std::endl;
        return false;
    } catch (std::string &str) {
        std::cerr << "Load matrix failed: " << str << std::endl;
        return false;
    }
    //Note: The old code wrote a message to standard output,
    //      if the matrix was not symmetric.  This code doesn't.
    return true;
}

#endif /* distancematrix_h */
