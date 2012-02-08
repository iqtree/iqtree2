/***************************************************************************
 *   Copyright (C) 2009 by BUI Quang Minh   *
 *   minh.bui@univie.ac.at   *
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
#include "phylotree.h"
#include "superalignment.h"
#include "phylosupertree.h"

SuperAlignment::SuperAlignment(PhyloSuperTree *super_tree)
 : Alignment()
{
	int site, seq, nsite = super_tree->size();
	PhyloSuperTree::iterator it;
	map<string,int> name_map;
	for (site = 0, it = super_tree->begin(); it != super_tree->end(); it++, site++) {
		partitions.push_back((*it)->aln);
		int nseq = (*it)->aln->getNSeq();
		for (seq = 0; seq < nseq; seq++) {
			int id = getSeqID((*it)->aln->getSeqName(seq));
			if (id < 0) {
				seq_names.push_back((*it)->aln->getSeqName(seq));
				id = seq_names.size()-1;
				IntVector vec(nsite, -1);
				vec[site] = seq;
				taxa_index.push_back(vec);
			} else taxa_index[id][site] = seq;
		}
	}
	num_states = 2; // binary type because the super alignment presents the presence/absence of taxa in the partitions
	site_pattern.resize(nsite, -1);
	clear();
	pattern_index.clear();
	VerboseMode save_mode = verbose_mode; 
	verbose_mode = VB_MIN; // to avoid printing gappy sites in addPattern
	int nseq = getNSeq();
	for (site = 0; site < nsite; site++) {
 		Pattern pat;
 		pat.append(nseq, 0);
		for (seq = 0; seq < nseq; seq++)
			pat[seq] = (taxa_index[seq][site] >= 0)? 1 : 0;
		addPattern(pat, site);
	}
	verbose_mode = save_mode;
	countConstSite();
}

double SuperAlignment::computeObsDist(int seq1, int seq2) {
	int site;
	int diff_pos = 0, total_pos = 0;
	for (site = 0; site < getNSite(); site++) {
		int id1 = taxa_index[seq1][site];
		int id2 = taxa_index[seq2][site];
		if (id1 < 0 || id2 < 0) continue;
		int num_states = partitions[site]->num_states;
		for (Alignment::iterator it = partitions[site]->begin(); it != partitions[site]->end(); it++) 
			if  ((*it)[id1] < num_states && (*it)[id2] < num_states) {
				total_pos += (*it).frequency;
				if ((*it)[id1] != (*it)[id2] )
					diff_pos += (*it).frequency;
			}
	}
	if (!total_pos) 
		return MAX_GENETIC_DIST; // return +INF if no overlap between two sequences
	return ((double)diff_pos) / total_pos;
}

double SuperAlignment::computeDist(int seq1, int seq2) {
	return computeObsDist(seq1, seq2);
}

SuperAlignment::~SuperAlignment()
{
	for (vector<Alignment*>::reverse_iterator it = partitions.rbegin(); it != partitions.rend(); it++)
		delete (*it);
	partitions.clear();
}


void SuperAlignment::printCombinedAlignment(const char *file_name, bool append) {
	vector<Alignment*>::iterator pit;
	int final_length = 0;
	for (pit = partitions.begin(); pit != partitions.end(); pit++)
		final_length += (*pit)->getNSite();
	try {
		ofstream out;
		out.exceptions(ios::failbit | ios::badbit);

		if (append)
			out.open(file_name, ios_base::out | ios_base::app);
		else
			out.open(file_name);
		out << getNSeq() << " " << final_length << endl;
		StrVector::iterator it;
		int max_len = getMaxSeqNameLength();
		if (max_len < 10) max_len = 10;
		int seq_id = 0;
		for (it = seq_names.begin(); it != seq_names.end(); it++, seq_id++) {
			out.width(max_len);
			out << left << (*it) << " ";
			int part = 0;
			for (pit = partitions.begin(); pit != partitions.end(); pit++, part++) {
				int part_seq_id = taxa_index[seq_id][part];
				int nsite = (*pit)->getNSite();
				if (part_seq_id >= 0) {
					for (int i = 0; i < nsite; i++)
						out << (*pit)->convertStateBack((*pit)->getPattern(i) [part_seq_id]);
				} else {
					string str(nsite, '?');
					out << str;
				}
			}
			out << endl;
		}
		out.close();
		cout << "Concatenated alignment was printed to " << file_name << endl;
	} catch (ios::failure) {
		outError(ERR_WRITE_OUTPUT, file_name);
	}	
}

double SuperAlignment::computeUnconstrainedLogL() {
	double logl = 0.0;
	vector<Alignment*>::iterator pit;
	for (pit = partitions.begin(); pit != partitions.end(); pit++)
		logl += (*pit)->computeUnconstrainedLogL();
	return logl;
}
