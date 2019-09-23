/***************************************************************************
 *   Copyright (C) 2018 by BUI Quang Minh   *
 *   m.bui@anu.edu.au   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef SUPERALIGNMENTUNLINKED_H
#define SUPERALIGNMENTUNLINKED_H

#include "superalignment.h"

class SuperAlignmentUnlinked : public SuperAlignment
{
public:
    /** constructor initialize from a supertree */
    SuperAlignmentUnlinked(Params &params);
    
    /** constructor initialize empty alignment */
    SuperAlignmentUnlinked();
    
    /**
     initialize seq_names, taxon_index, buildPattern
     */
    virtual void init(StrVector *sequence_names = NULL);
    
    /**
     * build all patterns of super alignent from partitions and taxa_index
     * it is in form of a binary alignment, where 0 means absence and 1 means presence
     * of a gene in a sequence
     */
    virtual void buildPattern();

    /**
     determine if the pattern is constant. update the is_const variable.
     */
    virtual void computeConst(Pattern &pat);

    /* build seq_states containing set of states per sequence
     * @param add_unobs_const TRUE to add all unobserved constant states (for +ASC model)
     */
//    void buildSeqStates(bool add_unobs_const = false);

    /** TRUE if all taxon sets are separate */
    bool unlinked_taxa;
    
};

#endif

