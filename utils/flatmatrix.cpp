//
// flatmatrix.cpp
// Distance matrix stored sequentially in row-major order
//
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
//
// This file, created by James Barbetti on 12-Nov-2020.
// (but: most of the code was moved here from decenttree.cpp,
//  where it had been since 24-Sep-2020).
// The only bits from 12-Nov-2020 are writeToDistanceFile,
// writeDistancesToOpenFile, and getMaxSeqNameLength.
//

#include "flatmatrix.h"
#include <math.h>       //for log10
#include <iostream>     //for std::fstream
#include <sstream>      //for std::stringstream
#if USE_GZSTREAM
#include "gzstream.h"   //for ogzstream
#else
#include <fstream>      //for std::ofstream
#endif

FlatMatrix::FlatMatrix(): rowCount(0), distanceMatrix(nullptr), borrowed(false) {
}

FlatMatrix::FlatMatrix(const StrVector& sequence_names,
                       double* distance_data)
    : sequenceNames(sequence_names), rowCount(sequence_names.size()),
      distanceMatrix(distance_data), borrowed(true){
}

FlatMatrix::~FlatMatrix() {
    if (!borrowed) {
        delete [] distanceMatrix;
    }
    distanceMatrix = nullptr;
}

const StrVector& FlatMatrix::getSequenceNames() const {
    return sequenceNames;
}

const std::string& FlatMatrix::sequenceName(intptr_t i) const {
    return sequenceNames[i];
}

std::string& FlatMatrix::sequenceName(intptr_t i) {
    return sequenceNames[i];
}

void FlatMatrix::setSequenceName(intptr_t i, 
                                 const std::string& new_name) {
    sequenceNames[i] = new_name;
}

void FlatMatrix::setSize(intptr_t rows) {
    if (!borrowed) {
        delete [] distanceMatrix;
    }
    borrowed = false;
    rowCount = rows;
    distanceMatrix = new double [ rowCount * rowCount ];
    memset(distanceMatrix, 0, rowCount*rowCount*sizeof(double));
}

intptr_t FlatMatrix::getSize() {
    return rowCount;
}

const double* FlatMatrix::getDistanceMatrix() const {
    return distanceMatrix;
}

double FlatMatrix::cell(intptr_t r, intptr_t c) const {
    return distanceMatrix[r * rowCount + c];
}

double& FlatMatrix::cell(intptr_t r, intptr_t c) {
    return distanceMatrix[r * rowCount + c];
}

void FlatMatrix::addCluster(const std::string& clusterName) {
    sequenceNames.emplace_back(clusterName);
}

bool FlatMatrix::writeToDistanceFile(const std::string& format,
                                     int precision,
                                     int compression_level,
                                     const std::string& file_name) const {
    try {
        if (format.find(".gz") == std::string::npos) {
            std::ofstream out;
            out.exceptions(std::ios::failbit | std::ios::badbit);
            out.open(file_name.c_str());
            writeDistancesToOpenFile(format, precision, out);
            out.close();
        } else {
            //Todo: Decide. Should we be insisting the file name ends with .gz too?
            #if USE_GZSTREAM
            ogzstream out;
            #else
            std::ofstream out;
            #endif
            out.exceptions(std::ios::failbit | std::ios::badbit);
            #if USE_GZSTREAM
            out.open(file_name.c_str(), std::ios::out, compression_level);
            #else
            out.open(file_name.c_str(), std::ios::out);
            #endif
            writeDistancesToOpenFile(format, precision, out);
            out.close();
        }
    } catch (std::ios::failure) {
        return false;
    }
    return true;
}

template <class S>
void FlatMatrix::writeDistancesToOpenFile(const std::string& format,
                                          int precision, S &out) const {
    intptr_t nseqs   = sequenceNames.size();
    size_t max_len = getMaxSeqNameLength();
    if (max_len < 10) {
        max_len = 10;
    }
    out << nseqs << std::endl;
    out.precision(precision);
    bool lower = (format.substr(0,5) == "lower");
    bool upper = (format.substr(0,5) == "upper");
    for (intptr_t seq1 = 0; seq1 < nseqs; ++seq1)  {
        std::stringstream line;
        line.width(max_len);
        line << std::fixed << std::left << sequenceNames[seq1];
        line.precision(precision);
        size_t rowStart = upper ? (seq1+1) : 0;
        size_t rowStop  = lower ? (seq1)   : nseqs;

        appendRowDistancesToLine(nseqs, seq1, rowStart, rowStop, line);
        line << "\n";
        out << line.str();
    }
    out.flush();
}

void FlatMatrix::appendRowDistancesToLine(intptr_t nseqs,    intptr_t seq1, 
                                          intptr_t rowStart, intptr_t rowStop,
                                          std::stringstream& line) const {
    intptr_t pos = seq1 * nseqs + rowStart;
    for (intptr_t seq2 = rowStart; seq2 < rowStop; ++seq2, ++pos) {
        if (distanceMatrix[pos] <= 0) {
            line << " 0";
        } else {
            line << " " << distanceMatrix[pos];
        }
    }
}

size_t FlatMatrix::getMaxSeqNameLength() const {
    size_t   len   = 0;
    intptr_t nseqs = sequenceNames.size();
    for (intptr_t i = 0; i < nseqs; ++i) {
        if (sequenceNames[i].length() > len) {
            len = sequenceNames[i].length();
        }
    }
    return len;
}
