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
#include "utils/timeutil.h"

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
    unlinked_taxa = true;
}

void SuperAlignmentUnlinked::init(StrVector *sequence_names) {
    // start original code
    
    max_num_states = 0;
    // first build taxa_index and partitions
    map<string, int> name2part;
    unlinked_taxa = true;
    for (auto it = partitions.begin(); it != partitions.end(); it++) {
        // Make sure that all partitions have different seq names
        for (auto sit = (*it)->seq_names.begin(); sit != (*it)->seq_names.end(); sit++) {
            if (name2part.find(*sit) != name2part.end()) {
                unlinked_taxa = false;
                break;
            }
            name2part[*sit] = (it) - partitions.begin();
        }
    }

    if (!unlinked_taxa) {
        // if some taxon sets are overlapping
        SuperAlignment::init(sequence_names);
        cout << "Linked " << seq_names.size() << " total sequences" << endl;
        return;
    }
    
    for (auto it = partitions.begin(); it != partitions.end(); it++) {
        seq_names.insert(seq_names.end(), (*it)->seq_names.begin(), (*it)->seq_names.end());
    }

    cout << "Unlinked " << seq_names.size() << " total sequences" << endl;
    
    /*
    taxa_index.resize(total_seqs, IntVector(npart, -1));
    for (auto it = partitions.begin(), part = 0, seq = 0; it != partitions.end(); it++, part++) {
        int part_nseq = (*it)->getNSeq();
        for (int part_seq = 0; part_seq < part_nseq; part_seq++, seq++) {
            taxa_index[seq][part] = part_seq;
        }
    }
    ASSERT(seq == total_seqs);
    */
    // now the patterns of sequence-genes presence/absence
    buildPattern();
}

void SuperAlignmentUnlinked::buildPattern() {
    if (!unlinked_taxa) {
        SuperAlignment::buildPattern();
        return;
    }
    int part, npart = partitions.size();
    seq_type = SEQ_BINARY;
    num_states = 2; // binary type because the super alignment presents the presence/absence of taxa in the partitions
    STATE_UNKNOWN = 2;
    site_pattern.resize(npart, -1);
    clear();
    pattern_index.clear();
    /*
    VerboseMode save_mode = verbose_mode;
    verbose_mode = min(verbose_mode, VB_MIN); // to avoid printing gappy sites in addPattern
    size_t nseq = getNSeq();
    int start_seq = 0;
    resize(npart, Pattern(nseq));
    for (part = 0; part < npart; part++) {
        Pattern *pat = &at(part);
        for (int seq = 0; seq < partitions[part]->getNSeq(); seq++)
            pat->at(start_seq + seq) = 1;
        //addPattern(pat, part);
        computeConst(*pat);
        
        // NOT USED FOR TOPO_UNLINKED
        //pattern_index[*pat] = part;
        site_pattern[part] = part;
        start_seq += partitions[part]->getNSeq();
    }
    ASSERT(start_seq == nseq);
    verbose_mode = save_mode;
    */
    resize(1, Pattern(getNSeq(), npart));
    computeConst(at(0));
    for (part = 0; part < npart; part++) {
        site_pattern[part] = 0;
    }
    
    countConstSite();
//    buildSeqStates();
}

void SuperAlignmentUnlinked::computeConst(Pattern &pat) {
    if (!unlinked_taxa) {
        SuperAlignment::computeConst(pat);
        return;
    }
    bool is_const = (partitions.size() == 1);
    bool is_invariant = (partitions.size() == 1);
    bool is_informative = (partitions.size() > 1);
    pat.const_char = (is_const) ? 1 : (STATE_UNKNOWN+1);
    
    pat.num_chars = (is_const) ? 1 : 2; // number of states with >= 1 appearance

    pat.flag = 0;
    if (is_const) pat.flag |= PAT_CONST;
    if (is_invariant) pat.flag |= PAT_INVARIANT;
    if (is_informative) pat.flag |= PAT_INFORMATIVE;
}

/*
void SuperAlignmentUnlinked::buildSeqStates(bool add_unobs_const) {
    if (!unlinked_taxa) {
        SuperAlignment::buildSeqStates(add_unobs_const);
        return;
    }
    seq_states.clear();
    if (add_unobs_const) {
        seq_states.resize(getNSeq(), IntVector({0,1}));
    } else {
        if (partitions.size() == 1)
            seq_states.resize(getNSeq(), IntVector({1}));
        else
            seq_states.resize(getNSeq(), IntVector({0,1}));
    }
}
*/
