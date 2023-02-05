//
//  bionj2.cpp - Advertisment of NJ and BIONJ algorithms
//               (that work in terms of .mldist inputs and
//                NEWICK outputs); this class defines functions
//               that make those algorithms available, via a 
//               StartTree::Registry.
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
// Created by James Barbetti on 18-Jun-2020.
// Notes:
// 1. Before 12-Aug-2020, the Matrix and ClusterTree template
//    classes were declared in this file (rather than in
//    separate headers) (see clustertree.h, distancematrix.h)
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
#include "auctionmatrix.h"
#include "fancyrapidnj.h"

namespace StartTree
{

#define ADVERTISE(type, shortName, longName) \
    f.advertiseTreeBuilder( new Builder<type>( shortName, longName))

/**
 * @brief Register tree-building implementations, in a registry of same,
 *        that are defined in the upgma.h, nj.h, rapidnj.h, auctionmatrix.h,
 *        and fancyrapidnj.h headers.
 * @param f - the "starting tree" distance-matrix phylogenetic inference 
 *            algorithm registry to which the implementations are to be
 *            added.
 * @note  The vectorized implementations are only included if the
 *        USE_VECTORCLASS_LIBRARY is defined and is set to a non-zero value.
 */
void addBioNJ2020TreeBuilders(Registry& f) {
    const char* defaultName = "RapidNJ";
    ADVERTISE(NJMatrix<NJFloat>,     "NJ",      "Neighbour Joining (Saitou, Nei [1987])");
    ADVERTISE(UNJMatrix<NJFloat>,    "UNJ",     "Unweighted Neighbour Joining (Gascel [1997])");
    ADVERTISE(RapidNJ,               "NJ-R",    "Rapid Neighbour Joining"
                                                " (Simonsen, Mailund, Pedersen [2011])");

    ADVERTISE(RapidNJ,           defaultName,   "Rapid Neighbour Joining"
                                                " (Simonsen, Mailund, Pedersen [2011])");

    ADVERTISE(FancyNJMatrix<NJFloat>,"ONJ-R",   "Rapid Neighbour Joining (a rival version)"
                                                " (Simonsen, Mailund, Pedersen [2011])");

    #if USE_VECTORCLASS_LIBRARY
    ADVERTISE(Vectorized_RapidNJ,    "NJ-R-V",  "Rapid Neighbour Joining (Vectorized)"
                                                " (Simonsen, Mailund, Pedersen [2011])");
    ADVERTISE(VectorizedFancyNJMatrix<NJFloat>, "ONJ-R-V",  "Rapid Neighbour Joining (a rival version)"
                                                " (Simonsen, Mailund, Pedersen [2011]) (Vectorized)");
    #endif

    #if USE_VECTORCLASS_LIBRARY
    ADVERTISE(VectorNJ,             "NJ-V",    "Vectorized Neighbour Joining (Saitou, Nei [1987])");
    #endif

    ADVERTISE(BIONJMatrix<NJFloat>, "BIONJ",   "BIONJ (Gascuel, Cong [2009])");
    ADVERTISE(RapidBIONJ,           "BIONJ-R", "Rapid BIONJ (Saitou, Nei [1987], Gascuel [2009],"
                                               " Simonson Mailund Pedersen [2011])");
    #if USE_VECTORCLASS_LIBRARY
    ADVERTISE(VectorBIONJ,          "BIONJ-V", "Vectorized BIONJ (Gascuel, Cong [2009])");
    #endif

    ADVERTISE(UPGMA_Matrix<NJFloat>,"UPGMA",   "UPGMA (Sokal, Michener [1958])");

    #if USE_VECTORCLASS_LIBRARY
    ADVERTISE(VectorizedUPGMA_Matrix<NJFloat>, "UPGMA-V", "Vectorized UPGMA (Sokal, Michener [1958])");
    #endif

    ADVERTISE(BoundingMatrix<double>,"NJ-R-D", "Double precision Rapid Neighbour Joining");
    f.setNameOfDefaultTreeBuilder(defaultName);
    ADVERTISE(DistanceAuctionMatrix, "AUCTION",    "Auction Joining");
}
}; //end of namespace
