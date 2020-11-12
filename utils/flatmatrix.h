//
// flatmatrix.h
// Defines FlatMatrix (a distance matrix of double precision
// distances, stored sequentially in row-major order).
// Copyright (C) 2020, James Barbetti.
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
    FlatMatrix();
    FlatMatrix(const std::vector<std::string>& sequence_names,
               double* distance_data);
    virtual ~FlatMatrix();
    
    const std::vector<std::string>& getSequenceNames() const;
    void          setSize(size_t rows);
    size_t        getSize();
    const double* getDistanceMatrix() const;
    double        cell(size_t r, size_t c) const;
    double&       cell(size_t r, size_t c);
    void          addCluster(const std::string& clusterName);
    size_t        getMaxSeqNameLength() const;
    bool          writeToDistanceFile(const std::string& format,
                                      int compression_level,
                                      const std::string& file_name) const;
    template <class S>
    void          writeDistancesToOpenFile(const std::string& format,
                                           S &out) const;
};

#endif /* flatmatrix_h */
