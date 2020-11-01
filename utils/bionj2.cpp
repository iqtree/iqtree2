//
//  bionj2.cpp - Implementations of NJ and BIONJ algorithms
//               (that work in terms of .mldist inputs and
//                NEWICK outputs).
//
//  Copyright (C) 2020, James Barbetti.
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
// Short names used for matrices and vectors (all indices start at 0).
//   D distance matrix (input, read from an .mldist file)
//   V estimated variance matrix (used in BIONJ, but not in NJ)
//   S bottom left triangle of the distance matrix,
//     with each row sorted in ascending order (see [SMP2011])
//   I index matrix, indicating (for each row of S, which cluster
//     each distance was calculated for).
//   Q connection-cost matrix (doesn't actually exist, cells of it
//     are calculated).
//   U vector of row totals (each row, the sum of the corresponding
//     row of the D matrix).  BIONJ implementations also use vectors
//     (indexed by cluster) of cluster totals.
//
// Created by James Barbetti on 18/6/2020.
// Notes:
// 1. Before 12-Aug-2020, the Matrix and ClusterTree template
//    classes were declared in this file (rather than in
//    separate headers)
// 2. Before 31-Oct-2020,
//    (a) UPGMA_Matrix, VectorizedUPGMA_Matrix were here
//        (they are now in upgma.h)
//    (b) NJMatrix, UNJMatrix, BIONJMatrix, VectorizedMatrix,
//        VectorNJ, VectorBIONJ were here (now in nj.h)
//    (c) BoundingMatrix, RapidNJ, RapidBIONJ were here
//        (they are now in rapidnj.h).
//

#include "starttree.h"
#include "upgma.h"
#include "nj.h"
#include "rapidnj.h"

namespace StartTree
{

#define ADVERTISE(type, shortName, longName) \
    f.advertiseTreeBuilder( new Builder<type>( shortName, longName))

void addBioNJ2020TreeBuilders(Factory& f) {
    const char* defaultName = "RapidNJ";
    ADVERTISE(NJMatrix<NJFloat>,   "NJ",      "Neighbour Joining (Saitou, Nei [1987])");
    ADVERTISE(UNJMatrix<NJFloat>,  "UNJ",     "Unweighted Neighbour Joining (Gascel [1997])");
    ADVERTISE(RapidNJ,             "NJ-R",    "Rapid Neighbour Joining"
                                              " (Simonsen, Mailund, Pedersen [2011])");
    ADVERTISE(VectorNJ,            "NJ-V",    "Vectorized Neighbour Joining (Saitou, Nei [1987])");
    ADVERTISE(BIONJMatrix<NJFloat>,"BIONJ",   "BIONJ (Gascuel, Cong [2009])");
    ADVERTISE(RapidBIONJ,          "BIONJ-R", "Rapid BIONJ (Saitou, Nei [1987], Gascuel [2009],"
                                              " Simonson Mailund Pedersen [2011])");
    ADVERTISE(VectorBIONJ,         "BIONJ-V", "Vectorized BIONJ (Gascuel, Cong [2009])");
    ADVERTISE(UPGMA_Matrix<NJFloat>,"UPGMA",  "UPGMA (Sokal, Michener [1958])");
    ADVERTISE(VectorizedUPGMA_Matrix<NJFloat>,"UPGMA-V", "Vectorized UPGMA (Sokal, Michener [1958])");
    ADVERTISE(BoundingMatrix<double>,"NJ-R-D","Double precision Rapid Neighbour Joining");
    ADVERTISE(RapidNJ,           defaultName, "Rapid Neighbour Joining (Simonsen, Mailund, Pedersen [2011]) (default)");
    f.setNameOfDefaultTreeBuilder(defaultName);
}
}; //end of namespace
