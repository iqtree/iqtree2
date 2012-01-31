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
#include "phylosupertree.h"

PhyloSuperTree::PhyloSuperTree()
 : PhyloTree()
{
}

PhyloSuperTree::PhyloSuperTree(Params &params) :  PhyloTree() {
	cout << "Reading partition model file " << params.partition_file << " ..." << endl;
	const int MAX_SIZE = 200;
	char buf[MAX_SIZE];
	try {
		ifstream in;
		in.exceptions(ios::failbit | ios::badbit);
		in.open(params.partition_file);

		while (!in.eof()) {
			in.getline(buf, MAX_SIZE, ',');
			if (in.eof()) break;
			part_names.push_back(buf);
			in.getline(buf, MAX_SIZE, ',');
			string model_name = buf;
			in.getline(buf, MAX_SIZE, ',');
			string aln_file = buf;
			in.getline(buf, MAX_SIZE, ',');
			string seq_type = buf;
			in.getline(buf, MAX_SIZE);
			string position_spec = buf;
			cout << "Reading partition " << part_names.back() << " (" << model_name << "," << 
				aln_file << "," << seq_type << "," << position_spec << ") ..." << endl;
		}

		in.clear();
		// set the failbit again
		in.exceptions(ios::failbit | ios::badbit);
		in.close();
	} catch(ios::failure) {
		outError(ERR_READ_INPUT);
	} catch (string str) {
		outError(str);
	}
}


PhyloSuperTree::~PhyloSuperTree()
{
	freePhyloTree();
}


