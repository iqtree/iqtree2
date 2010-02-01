/***************************************************************************
 *   Copyright (C) 2006 by BUI Quang Minh, Steffen Klaere, Arndt von Haeseler   *
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
#include "mtreeset.h"

MTreeSet::MTreeSet()
{
}

MTreeSet::MTreeSet(const char *userTreeFile, bool &is_rooted, int burnin) {
	init(userTreeFile, is_rooted, burnin);
}

void MTreeSet::init(const char *userTreeFile, bool &is_rooted, int burnin) {
	readTrees(userTreeFile, is_rooted, burnin);
	checkConsistency();
}

void MTreeSet::readTrees(const char *infile, bool &is_rooted, int burnin) {
	cout << "Reading tree(s) file " << infile << " ..." << endl;
	ifstream in;
	int count = 1;
	try {
		in.exceptions(ios::failbit | ios::badbit);
		in.open(infile);
		if (burnin > 0) {
			int cnt = 0;
			while (cnt < burnin && !in.eof()) {
				char ch;
				in >> ch;
				if (ch == ';') cnt++;
			}
			cout << cnt << " beginning tree(s) discarded" << endl;
			if (in.eof())
				throw "Burnin value is too large.";
		}
		for (count = 1; !in.eof(); count++) {
			//cout << "Reading tree " << count << " ..." << endl;
			MTree *tree = newTree();
			bool myrooted = is_rooted;
			//tree->userFile = (char*) infile;
			tree->readTree(in, myrooted);
			push_back(tree);
			//cout << "Tree contains " << tree->leafNum - tree->rooted << 
			//" taxa and " << tree->nodeNum-1-tree->rooted << " branches" << endl;
			char ch;
			in.exceptions(ios::goodbit);
			in >> ch;
			if (in.eof()) break;
			in.unget();
			in.exceptions(ios::failbit | ios::badbit);
		}
		cout << count << " tree(s) loaded" << endl;
		//in.exceptions(ios::failbit | ios::badbit);
		in.close();
	} catch (ios::failure) {
		outError(ERR_READ_INPUT, infile);		
	} catch (const char* str) {
		outError(str);
	}
}

void MTreeSet::checkConsistency() {
	if (empty()) 
		return;
	iterator it;
	bool rooted = front()->rooted;
	for (it = begin()+1; it != end(); it++)
		if ((*it)->rooted != rooted) 
			outError("Rooted and unrooted trees are mixed up");

	NodeVector taxa1;
	NodeVector::iterator it2;
	int i;

	for (it = begin(); it != end(); it++) {
		MTree *tree = *it;
		NodeVector taxa;
		tree->getTaxa(taxa);
		sort(taxa.begin(), taxa.end(), nodenamecmp);
		for (it2 = taxa.begin(), i = 0; it2 != taxa.end(); it2++, i++)
			(*it2)->id = i;

		if (it == begin()) 
			taxa1 = taxa;
		else {
			// now check this tree with the first tree	
			if (tree->leafNum != taxa1.size())
				outError("Tree has different number of taxa!");
	
			for (it2 = taxa.begin(), i = 0; it2 != taxa.end(); it2++, i++) {
				if ((*it2)->name != taxa1[i]->name) 
					outError("Tree has different taxa names!");
			}
		}
	}
}

bool MTreeSet::isRooted() {
	if (empty()) return false;
	return (front()->rooted);
}


void MTreeSet::printTrees(const char *ofile, int  brtype)
{
	try {
		ofstream out;
		out.exceptions(ios::failbit | ios::badbit);
		out.open(ofile);
		printTrees(out, brtype);
		out.close();
		cout << "Tree(s) were printed to " << ofile << endl;
	} catch (ios::failure) {
		outError(ERR_WRITE_OUTPUT, ofile);
	}
}

void MTreeSet::printTrees(ostream & out, int brtype) {
	for (iterator  it = begin(); it != end(); it++) {
		(*it)->printTree(out, brtype);
		out << endl;
	}
}

void MTreeSet::convertSplits(SplitGraph &sg, double split_threshold, bool lensum) {
	SplitIntMap hash_ss;
	if (split_threshold == 0.0) {
		convertSplits(sg, hash_ss, lensum);
		return;
	}
	//SplitGraph temp;
	convertSplits(sg, hash_ss, lensum);
	int nsplits = sg.getNSplits();

	double threshold = split_threshold * size();
	
	for (SplitGraph::iterator it = sg.begin(); it != sg.end(); ) {
		//SplitIntMap::iterator ass_it = hash_ss.find(*it);
		int freq_value;
		Split *sp = hash_ss.findSplit(*it, freq_value);
		assert(sp != NULL);
		//Split *sp = ass_it->first;
		if (freq_value <= threshold) {
			if (verbose_mode == VB_DEBUG) {
				sp->report(cout);
			}
			int num = hash_ss.getValue(sg.back());
			//hash_ss.erase(ass_it);
			hash_ss.eraseSplit(sp);
			if (it != sg.end()-1) {
				hash_ss.eraseSplit(sg.back());
				*(*it) = (*sg.back());
			}
			delete sg.back();
			sg.pop_back();
			if (it == sg.end()) break;
			hash_ss.insertSplit(*it, num);
		} else {
			//sg.push_back(new Split(*sp));
			it++;
		}
	}
	/*
	sg.taxa = temp.taxa;
	sg.splits = temp.splits;
	sg.pda = temp.pda;
	sg.sets = temp.sets;
	sg.trees = temp.trees;
	temp.taxa = NULL;
	temp.splits = NULL;
	temp.pda = NULL;
	temp.sets = NULL;
	temp.trees = NULL;
	*/
	cout << nsplits - sg.getNSplits() << " split(s) discarded due to frequency <= " << split_threshold << endl;
}


void MTreeSet::convertSplits(SplitGraph &sg, SplitIntMap &hash_ss, bool lensum) {
	vector<NxsString> taxname(front()->leafNum);
	// make sure that the split system contains at least 1 split
	if (size() == 0)
		return;
	
	front()->getTaxaName(taxname);
	convertSplits(taxname, sg, hash_ss, lensum);
}

void MTreeSet::convertSplits(vector<NxsString> &taxname, SplitGraph &sg, SplitIntMap &hash_ss, bool lensum) {

#ifdef USE_HASH_MAP
	cout << "Using hash_map" << endl;
#else
	cout << "Using map" << endl;
#endif
	cout << "Converting collection of tree(s) into split system..." << endl;
	SplitGraph::iterator itg;
	vector<NxsString>::iterator its;
/*
	for (its = taxname.begin(); its != taxname.end(); its++)
		if (*its == ROOT_NAME) {	
			taxname.erase(its);
			break;
		}*/
	sort(taxname.begin(), taxname.end());
	sg.createBlocks();
	if (size() == 1 && lensum) {
		front()->convertSplits(sg);
		return;
	}

	for (its = taxname.begin(); its != taxname.end(); its++)
		sg.getTaxa()->AddTaxonLabel(*its);

	SplitGraph *isg;
	for (iterator it = begin(); it != end(); it++) {
		MTree *tree = *it;
		if (tree->leafNum != taxname.size())
			outError("Tree has different number of taxa!");
		NodeVector taxa;
		tree->getTaxa(taxa);
		sort(taxa.begin(), taxa.end(), nodenamecmp);
		int i = 0;
		for (NodeVector::iterator it2 = taxa.begin(); it2 != taxa.end(); it2++) {
			if ((*it2)->name != taxname[i]) 
				outError("Tree has different taxa names!");
			(*it2)->id = i++;
		}
		isg = new SplitGraph();
		tree->convertSplits(taxname, *isg);
		//isg->getTaxa()->Report(cout);
		//isg->report(cout);
		for (itg = isg->begin(); itg != isg->end(); itg++) {
			//SplitIntMap::iterator ass_it = hash_ss.find(*itg);
			int value;
			Split *sp = hash_ss.findSplit(*itg, value);
			if (sp != NULL) {
				//Split *sp = ass_it->first;
				if (lensum)
					sp->setWeight(sp->getWeight() + (*itg)->getWeight());
				else
					sp->setWeight(sp->getWeight() + 1);
				hash_ss.setValue(sp, value + 1);
			}
			else {
				sp = new Split(*(*itg));
				if (lensum)
					sp->setWeight((*itg)->getWeight());
				else				
					sp->setWeight(1);
				sg.push_back(sp);
				//SplitIntMap::value_type spair(sp, 1);
				//hash_ss.insert(spair);
				
				hash_ss.insertSplit(sp, 1);
 			}
		}
		delete isg;
	}

	//sg.report(cout);
}


MTreeSet::~MTreeSet()
{
	for (reverse_iterator it = rbegin(); it != rend(); it++) {
		MTree *tree = *it;
		delete tree;
	}
	clear();
}


