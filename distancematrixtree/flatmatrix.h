//
// flatmatrix.h
// Defines FlatMatrix (a distance matrix of double precision
// distances, stored sequentially in row-major order).
// Copyright (C) 2020-22, James Barbetti.
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

#ifndef flatmatrix_h
#define flatmatrix_h
#include <stdlib.h> //for size_t
#include <vector>   //for std::vector
#include <string>   //for std::string
#include <utils/vectortypes.h> //for StrVector

class FlatMatrix {
protected:
    StrVector sequenceNames;   //The names of the sequences
    intptr_t  rowCount;        //The rank of the matrix
    double*   distanceMatrix;  //The data in the matrix, in row-major order
    bool      borrowed;        //if true, *this instance does not own the matrix
                               //if false, *this instance owns it, and must 
                               //delete[] it, in its destructor.
public:
    typedef double cell_type;
    /**
     * @brief Construct an empty FlatMatrix
     * @note  The rank will be zero.  No memory will be allocated.
     */
    FlatMatrix();
    /**
     * @brief Construct a new FlatMatrix
     * @param sequence_names the names of the sequences
     * @param distance_data  and n*n matrix, in row-major order, of the distances
     *                       between the sequences.
     */
    FlatMatrix(const StrVector& sequence_names,
               double* distance_data);
    virtual ~FlatMatrix();
    
    const StrVector&   getSequenceNames()           const;
    size_t             getMaxSeqNameLength()        const;
    const std::string& getSequenceName(intptr_t i)  const;
    const std::string& sequenceName(intptr_t i)     const;
    std::string&       sequenceName(intptr_t i);
    void               setSequenceName(intptr_t i, 
                                       const std::string& new_name);
    virtual void       setSize(intptr_t rows);
    intptr_t           getSize()                    const;
    /**
     * @brief  get a pointer to the start of the matrix data
     *         (a flat matrix, storing rowCount * rowCount entries in
     *         row-major order).
     * @return const double* - pointer to the start of the matrix data
     */
    const double*      getDistanceMatrix()          const;
    /**
     * @brief  read-only matrix-entry access
     * @param  r row number    (0<=r<rowCount) 
     * @param  c column number (0<=c<rowCount) 
     * @return double content of the entry in row r and column c
     */
    double             cell(intptr_t r, intptr_t c) const;
    /**
     * @brief  read-write matrix-entry access
     * @param  r row number    (0<=r<rowCount) 
     * @param  c column number (0<=c<rowCount) 
     * @return double& - a reference to the entry, in the matrix, in
     *         row r and column c.
     */
    double&            cell(intptr_t r, intptr_t c);
    /**
     * @brief  add a sequence (or a cluster) to the list of sequences
     * @param  clusterName 
     */
    void               addCluster(const std::string& clusterName);
    /**
     * @brief  write the sequence names, and the distance matrix,
     *         in one of the Phylip distance matrix formats, to the 
     *         file with the specified file path
     * @param  format    a string indicating the format "square", "upper", 
     *                   or "lower".
     * @param  precision the precision (the number of digits to display,
     *                   after the decimal point, for entries in the matrix)
     * @param  compression_level how much to compress the distance matrix
     * @param  report_progress   if true, report progress as the write happens.
     *                           if false, don't.
     * @param  file_name the file path, to write the file to
     * @return true  if the file write succeeds
     * @return false if the file write fails (an error will be written to std::cerr)
     */
    bool               writeToDistanceFile(const std::string& format,
                                           int precision,
                                           int compression_level,
                                           bool report_progress,
                                           const std::string& file_name) const;

    /**
     * @brief 
     * 
     * @tparam S the stream type 
     * @tparam P the progress display class
     * @param  format    the format ("square", "upper", or "lower")
     * @param  precision the precision (number of digits to use after decimal point)
     * @param  out       the output stream
     * @param  progress  pointer to the class instance used to report progress
     */
    template <class S /*stream type*/, class P /*progress display*/>
    void               writeDistancesToOpenFile(const std::string& format,
                                                int precision, S &out,
                                                P* progress) const;

    /**
    * @brief Write distances, from one row of the matrix, to a string.
    *        Each distance is prefixed by a single white space, ' '.
    * @param nseqs    The number of sequences
    * @param seq1     The sequence (the row number)
    * @param rowStart The column number within the row, to start with
    * @param rowStop  The column number after the last column to output
    * @param line     A stringstream to which a description is to be written
    */
    virtual void appendRowDistancesToLine(intptr_t nseqs,    intptr_t seq1, 
                                         intptr_t rowStart, intptr_t rowStop,
                                         std::stringstream& line) const; 
};

#endif /* flatmatrix_h */
