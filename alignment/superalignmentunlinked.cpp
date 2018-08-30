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

#include "superalignmentunlinked.h"

/** constructor initialize from a supertree */
SuperAlignmentUnlinked::SuperAlignmentUnlinked(Params &params)
: SuperAlignment()
{
    readFromParams(params);
    init();
}

/** constructor initialize empty alignment */
SuperAlignmentUnlinked::SuperAlignmentUnlinked()
: SuperAlignment()
{
}

void SuperAlignmentUnlinked::init(StrVector *sequence_names) {
    // start original code
    
    max_num_states = 0;
    // first build taxa_index and partitions
    int part, seq, npart = partitions.size();

    ASSERT(!sequence_names);
    
    map<string, int> name2part;
    vector<Alignment*>::iterator it;
    size_t total_seqs = 0;
    for (it = partitions.begin(); it != partitions.end(); it++)
        total_seqs += (*it)->seq_names.size();
    taxa_index.resize(total_seqs, IntVector(npart, -1));
    
    part = 0;
    for (auto it = partitions.begin(); it != partitions.end(); it++, part++) {
        // Make sure that all partitions have different seq names
        for (auto sit = (*it)->seq_names.begin(); sit != (*it)->seq_names.end(); sit++) {
            if (name2part.find(*sit) != name2part.end())
                outError("Duplicate taxon name " + (*sit) + " in partitions " + partitions[name2part[*sit]]->name + " and " + (*it)->name);
            name2part[*sit] = it - partitions.begin();
        }
        int nseq = (*it)->getNSeq();
        for (seq = 0; seq < nseq; seq++) {
            taxa_index[seq_names.size()+seq][part] = seq;
        }
        seq_names.insert(seq_names.end(), (*it)->seq_names.begin(), (*it)->seq_names.end());
    }
    // now the patterns of sequence-genes presence/absence
    buildPattern();
}
