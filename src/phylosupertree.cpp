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
	try {
		ifstream in;
		in.exceptions(ios::failbit | ios::badbit);
		in.open(params.partition_file);
		in.exceptions(ios::badbit);

		Params origin_params = params;
		//memcpy(&origin_params, &params, sizeof(params));
		string part_name, model_name, aln_file, sequence_type, position_spec;

		while (!in.eof()) {
			getline(in, part_name, ',');
			if (in.eof()) break;
			params = origin_params;
			part_names.push_back(part_name);
			getline(in, model_name, ',');
			if (model_name == "*") model_name = params.model_name;
			getline(in, aln_file, ',');
			if (aln_file == "*") aln_file = params.aln_file;
			getline(in, sequence_type, ',');
			if (sequence_type=="*") sequence_type = params.sequence_type;
			getline(in, position_spec);
			cout << "Reading partition " << part_name << " (" << params.model_name << "," << 
				aln_file << "," << sequence_type << "," << position_spec << ") ..." << endl;
			Alignment *aln = new Alignment((char*)aln_file.c_str(), (char*)sequence_type.c_str(), params.intype);
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


