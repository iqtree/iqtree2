//
// flatmatrix.h
// Defines FlatMatrix (a distance matrix of double precision
// distances, stored sequentially in row-major order).
// Copyright (C) 2020, James Barbetti.
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

class FlatMatrix {
private:
    std::vector<std::string> sequenceNames;
    size_t                   rowCount;
    double*                  distanceMatrix;
    bool                     borrowed;
public:
    typedef double cell_type;
    FlatMatrix();
    FlatMatrix(const std::vector<std::string>& sequence_names,
               double* distance_data);
    virtual ~FlatMatrix();
    
    const std::vector<std::string>& getSequenceNames() const;
    size_t             getMaxSeqNameLength()    const;
    const std::string& sequenceName(size_t i)   const;
    std::string&       sequenceName(size_t i);
    void               setSize(size_t rows);
    size_t             getSize();
    const double*      getDistanceMatrix()      const;
    double             cell(size_t r, size_t c) const;
    double&            cell(size_t r, size_t c);
    void               addCluster(const std::string& clusterName);
    bool               writeToDistanceFile(const std::string& format,
                                           int precision,
                                           int compression_level,
                                           const std::string& file_name) const;
    template <class S>
    void          writeDistancesToOpenFile(const std::string& format,
                                           int precision, S &out) const;
};

#endif /* flatmatrix_h */
