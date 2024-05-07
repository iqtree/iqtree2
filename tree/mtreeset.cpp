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
#include "alignment/alignment.h"
#include "utils/gzstream.h"

MTreeSet::MTreeSet()
{
    equal_taxon_set = false;
}

MTreeSet::MTreeSet(const char *userTreeFile, bool &is_rooted, 
	int burnin, int max_count, const char *tree_weight_file) {
	init(userTreeFile, is_rooted, burnin, max_count, tree_weight_file);
}

void readIntVector(const char *file_name, int burnin, int max_count, IntVector &vec) {
	cout << "Reading integer vector file " << file_name << " ..." << endl;
	vec.clear();
	try {
		ifstream in;
		in.exceptions(ios::failbit | ios::badbit);
		in.open(file_name);
		// remove the failbit
		in.exceptions(ios::badbit);

		for (; !in.eof();) {
			int i;
			if(!(in >> i)) break;
			if (burnin > 0) 
				burnin--;
			else if (max_count > 0) {
				vec.push_back(i);
				max_count--;
			}
		}
		in.clear();
		// set the failbit again
		in.exceptions(ios::failbit | ios::badbit);
		in.close();
	} catch(ios::failure) {
		outError(ERR_READ_INPUT);
	}
}

void MTreeSet::init(const char *userTreeFile, bool &is_rooted, int burnin, int max_count,
	const char *tree_weight_file, IntVector *weights, bool compressed) 
{
	readTrees(userTreeFile, is_rooted, burnin, max_count, weights, compressed);
	checkConsistency();

	if (tree_weight_file) 
		readIntVector(tree_weight_file, burnin, max_count, tree_weights);
/*	else if (!weights)
		tree_weights.resize(size(), 1);*/

	if (size() != tree_weights.size())
		outError("Tree file and tree weight file have different number of entries");

}

void MTreeSet::init(StringIntMap &treels, bool &is_rooted, IntVector &weights) {
	//resize(treels.size(), NULL);
	int count = 0;
	//IntVector ok_trees;
	//ok_trees.resize(treels.size(), 0);
	//for (i = 0; i < trees_id.size(); i++) ok_trees[trees_id[i]] = 1;

	for (StringIntMap::iterator it = treels.begin(); it != treels.end(); it++) 
	if (weights[it->second]) {
		count++;
		MTree *tree = newTree();
		stringstream ss(it->first);
		bool myrooted = is_rooted;
		tree->readTree(ss, myrooted);
		NodeVector taxa;
		tree->getTaxa(taxa);
		for (NodeVector::iterator taxit = taxa.begin(); taxit != taxa.end(); taxit++)
			(*taxit)->id = atoi((*taxit)->name.c_str());
		//at(it->second) = tree;
		push_back(tree);
		tree_weights.push_back(weights[it->second]);
		//cout << "Tree " << it->second << ": ";
		//tree->printTree(cout, WT_NEWLINE);
	}
	if (verbose_mode >= VB_MED)
		cout << count << " tree(s) converted" << endl;
	//tree_weights.resize(size(), 1);
}

void MTreeSet::init(StrVector &treels, bool &is_rooted) {
	//resize(treels.size(), NULL);
	int count = 0;
	//IntVector ok_trees;
	//ok_trees.resize(treels.size(), 0);
	//for (i = 0; i < trees_id.size(); i++) ok_trees[trees_id[i]] = 1;

	for (StrVector::iterator it = treels.begin(); it != treels.end(); it++)
    if (!it->empty())
	{
		count++;
		MTree *tree = newTree();
		stringstream ss(*it);
		bool myrooted = is_rooted;
		tree->readTree(ss, myrooted);
		NodeVector taxa;
		tree->getTaxa(taxa);
		for (NodeVector::iterator taxit = taxa.begin(); taxit != taxa.end(); taxit++) {
			if ((*taxit)->name == ROOT_NAME) {
				(*taxit)->id = taxa.size() - 1;
			}
			else {
				(*taxit)->id = atoi((*taxit)->name.c_str());
			}
		}
		//at(it->second) = tree;
		push_back(tree);
		tree_weights.push_back(1);
		//cout << "Tree " << it->second << ": ";
		//tree->printTree(cout, WT_NEWLINE);
	}
	if (verbose_mode >= VB_MED)
		cout << count << " tree(s) converted" << endl;
	//tree_weights.resize(size(), 1);
}

void MTreeSet::init(vector<string> &trees, vector<string> &taxonNames, bool &is_rooted) {
	int count = 0;
	for (vector<string>::iterator it = trees.begin(); it != trees.end(); it++) {
		MTree *tree = newTree();
		stringstream ss(*it);
		tree->readTree(ss, is_rooted);
	    int nseq = taxonNames.size();
	    ASSERT(tree->getNumTaxa() == nseq);
	    for (int seq = 0; seq < nseq; seq++) {
	        string seq_name = taxonNames[seq];
	        Node *node = tree->findLeafName(seq_name);
	        ASSERT(node);
	        ASSERT(node->isLeaf());
	        node->id = seq;
	    }
		push_back(tree);
		tree_weights.push_back(1);
		count++;
	}
    
    //cout << count << " tree(s) converted" << endl;

}

void MTreeSet::readTrees(const char *infile, bool &is_rooted, int burnin, int max_count,
	IntVector *weights, bool compressed) 
{
	cout << "Reading tree(s) file " << infile << " ..." << endl;
	int count, omitted;
/*	IntVector ok_trees;
	if (trees_id) {
		int max_id = *max_element(trees_id->begin(), trees_id->end());
		ok_trees.resize(max_id+1, 0);
		for (IntVector::iterator it = trees_id->begin(); it != trees_id->end(); it++)
			ok_trees[*it] = 1;
		cout << "Restricting to " << trees_id->size() << " trees" << endl;
	}*/
	try {
		istream *in;
		if (compressed) in = new igzstream; else in = new ifstream;
		in->exceptions(ios::failbit | ios::badbit);
		
		if (compressed) ((igzstream*)in)->open(infile); else ((ifstream*)in)->open(infile);
		if (burnin > 0) {
			int cnt = 0;
			while (cnt < burnin && !in->eof()) {
				char ch;
				(*in) >> ch;
				if (ch == ';') cnt++;
			}
			cout << cnt << " beginning tree(s) discarded" << endl;
			if (in->eof())
				throw "Burnin value is too large.";
		}
		for (count = 1, omitted = 0; !in->eof() && count <= max_count; count++) {
			if (!weights || weights->at(count-1)) {
				//cout << "Reading tree " << count << " ..." << endl;
				MTree *tree = newTree();
				bool myrooted = is_rooted;
				//tree->userFile = (char*) infile;
				tree->readTree(*in, myrooted);
				push_back(tree);
				if (weights) 
					tree_weights.push_back(weights->at(count-1)); 
				else tree_weights.push_back(1);
				//cout << "Tree contains " << tree->leafNum - tree->rooted << 
				//" taxa and " << tree->nodeNum-1-tree->rooted << " branches" << endl;
			} else {
				// omit the tree
				//push_back(NULL);
				//in->exceptions(ios::badbit);
				while (!in->eof()) {
					char ch;
					if (!((*in) >> ch)) break;
					if (ch == ';') break;
				}
				omitted++;
			} 
			char ch;
			in->exceptions(ios::goodbit);
			(*in) >> ch;
			if (in->eof()) break;
			in->unget();
			in->exceptions(ios::failbit | ios::badbit);

		}
		cout << size() << " tree(s) loaded (" << countRooted() << " rooted and " << countUnrooted() << " unrooted)" << endl;
		if (omitted) cout << omitted << " tree(s) omitted" << endl;
		//in->exceptions(ios::failbit | ios::badbit);
		if (compressed) ((igzstream*)in)->close(); else ((ifstream*)in)->close();
		// following line was missing which caused small memory leak
		delete in;
	} catch (ios::failure) {
		outError(ERR_READ_INPUT, infile);		
	} catch (const char* str) {
		outError(str);
	}
}

void MTreeSet::checkConsistency() {
    equal_taxon_set = true;
	if (empty()) 
		return;
	iterator it;
	int i;
	bool first = true;
    // it's OK to have rooed and unrooted trees in the set
    /*
    bool rooted = false;
	for (it = begin(), i = 0; it != end(); it++, i++)
	if ((*it)) {
		if (!first && (*it)->rooted != rooted) {
			cout << i+1 << " " << (*it)->rooted << " " << rooted << endl;
			outError("Rooted and unrooted trees are mixed up");
		}
		first = false;
        rooted = (*it)->rooted;
	}
     */

	NodeVector taxa1;
	NodeVector::iterator it2;

	first = true;
	for (it = begin(); it != end(); it++) if (*it) {
		MTree *tree = *it;
		NodeVector taxa;
		tree->getTaxa(taxa);
		sort(taxa.begin(), taxa.end(), nodenamecmp);
		for (it2 = taxa.begin(), i = 0; it2 != taxa.end(); it2++, i++)
			(*it2)->id = i;

		if (first ) {
			taxa1 = taxa;
			first = false;
		} else {
			// now check this tree with the first tree	
            if (tree->leafNum != taxa1.size()) {
                equal_taxon_set = false;
				cout << "Trees have different number of taxa" << endl;
                break;
            }
	
			for (it2 = taxa.begin(), i = 0; it2 != taxa.end(); it2++, i++) {
                if ((*it2)->name != taxa1[i]->name) {
                    equal_taxon_set = false;
					cout << "Trees have different taxa sets" << endl;
                    break;
                }
			}
            if (!equal_taxon_set) break;
		}
	}
}

bool MTreeSet::isRooted() {
	if (empty()) return false;
	return (front()->rooted);
}

int MTreeSet::countRooted() {
    int count = 0;
    for (auto tree = begin(); tree != end(); tree++)
        if ((*tree)->rooted)
            count++;
    return count;
}

/**
 @return number of unrooted trees in the set
 */
int MTreeSet::countUnrooted() {
    return size() - countRooted();
}


void MTreeSet::assignLeafID() {
	for (iterator it = begin(); it != end(); it++)
		(*it)->assignLeafID();
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

void MTreeSet::convertSplits(SplitGraph &sg, double split_threshold, int weighting_type, 
	double weight_threshold) 
{
	SplitIntMap hash_ss;
/*
	if (split_threshold == 0.0) {
		convertSplits(sg, hash_ss, weighting_type, weight_threshold);
		return;
	}*/
	//SplitGraph temp;
	convertSplits(sg, hash_ss, weighting_type, weight_threshold);
	int nsplits = sg.getNSplits();

	double threshold = split_threshold * size();
//	cout << "threshold = " << threshold << endl;
	int count=0;
	for (SplitGraph::iterator it = sg.begin(); it != sg.end(); ) {
		count++;
		//SplitIntMap::iterator ass_it = hash_ss.find(*it);
		int freq_value;
		Split *sp = hash_ss.findSplit(*it, freq_value);
		ASSERT(sp != NULL);
		ASSERT(*sp == *(*it));
		//Split *sp = ass_it->first;
		if (freq_value <= threshold) {
			if (verbose_mode == VB_DEBUG) {
				sp->report(cout);
			}
			int num = hash_ss.getValue(sg.back());
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
	cout << nsplits - sg.getNSplits() << " split(s) discarded because frequency <= " << split_threshold << endl;
}


void MTreeSet::convertSplits(SplitGraph &sg, SplitIntMap &hash_ss, int weighting_type, double weight_threshold) {
	vector<string> taxname(front()->leafNum);
	// make sure that the split system contains at least 1 split
	if (size() == 0)
		return;
	
	front()->getTaxaName(taxname);
	convertSplits(taxname, sg, hash_ss, weighting_type, weight_threshold, NULL);
}

void MTreeSet::convertSplits(vector<string> &taxname, SplitGraph &sg, SplitIntMap &hash_ss, 
	int weighting_type, double weight_threshold, char *tag_str, bool sort_taxa) {

	if (verbose_mode >= VB_MED) {
	#ifdef USE_HASH_MAP
		cout << "Using hash_map" << endl;
	#else
		cout << "Using map" << endl;
	#endif

		cout << "Converting collection of tree(s) into split system..." << endl;
	}
	SplitGraph::iterator itg;
	vector<string>::iterator its;
/*
	for (its = taxname.begin(); its != taxname.end(); its++)
		if (*its == ROOT_NAME) {	
			taxname.erase(its);
			break;
		}*/
	if (sort_taxa) sort(taxname.begin(), taxname.end());
	sg.createBlocks();
	for (its = taxname.begin(); its != taxname.end(); its++)
		sg.getTaxa()->AddTaxonLabel(NxsString(its->c_str()));
/*
	if (size() == 1 && weighting_type != SW_COUNT) {
		front()->convertSplits(taxname, sg);
		return;
	}*/


	SplitGraph *isg;
	int tree_id = 0;
//	cout << "Number of trees: " << size() << endl;
//	cout << "Number of weight: " << tree_weights.size() << endl;
	for (iterator it = begin(); it != end(); it++, tree_id++) {
		if (tree_weights[tree_id] == 0) continue;
		MTree *tree = *it;

		if (tree->leafNum != taxname.size())
			outError("Tree has different number of taxa!");
		if (sort_taxa) {
			NodeVector taxa;
			tree->getTaxa(taxa);
			sort(taxa.begin(), taxa.end(), nodenamecmp);
			int i = 0;
			for (NodeVector::iterator it2 = taxa.begin(); it2 != taxa.end(); it2++) {
				if ((*it2)->name != taxname[i]) {
					cout << "Name 1: " <<  (*it2)->name << endl;
					cout << "Name 2: " <<  taxname[i] << endl;
					outError("Tree has different taxa names!");
				}
				(*it2)->id = i++;
			}
		}
		isg = new SplitGraph();
		tree->convertSplits(taxname, *isg);
		//isg->getTaxa()->Report(cout);
		//isg->report(cout);
		for (itg = isg->begin(); itg != isg->end(); itg++) {
			//SplitIntMap::iterator ass_it = hash_ss.find(*itg);
			int value;
			//if ((*itg)->getWeight()==0.0) cout << "zero weight!" << endl;
			Split *sp = hash_ss.findSplit(*itg, value);
			if (sp != NULL) {
				//Split *sp = ass_it->first;
				if (weighting_type != SW_COUNT)
					sp->setWeight(sp->getWeight() + (*itg)->getWeight() * tree_weights[tree_id]);
				else
					sp->setWeight(sp->getWeight() + tree_weights[tree_id]);
				hash_ss.setValue(sp, value + tree_weights[tree_id]);
			}
			else {
				sp = new Split(*(*itg));
				if (weighting_type != SW_COUNT)
					sp->setWeight((*itg)->getWeight() * tree_weights[tree_id]);
				else				
					sp->setWeight(tree_weights[tree_id]);
				sg.push_back(sp);
				//SplitIntMap::value_type spair(sp, 1);
				//hash_ss.insert(spair);
				
				hash_ss.insertSplit(sp, tree_weights[tree_id]);
 			}
            if (tag_str)
                sp->name += "@" + convertIntToString(tree_id+1);
		}
		delete isg;
	}

	if (weighting_type == SW_AVG_PRESENT) {
		for (itg = sg.begin(); itg != sg.end(); itg++) {
			int value = 0;
			if (!hash_ss.findSplit(*itg, value))
				outError("Internal error ", __func__);
			(*itg)->setWeight((*itg)->getWeight() / value);
		}
	} else if (weighting_type == SW_AVG_ALL) {
		for (itg = sg.begin(); itg != sg.end(); itg++) {
			(*itg)->setWeight((*itg)->getWeight() / tree_weights.size());
		}
	}

	int discarded = 0;	
	for (itg = sg.begin(); itg != sg.end(); )  {
		if ((*itg)->getWeight() <= weight_threshold) {
			discarded++;
			delete (*itg);
			(*itg) = sg.back();
			sg.pop_back(); 
		} else itg++;
	}
	if (discarded)
		cout << discarded << " split(s) discarded because weight <= " << weight_threshold << endl;
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


void MTreeSet::computeRFDist(double *rfdist, int mode, double weight_threshold) {
	// exit if less than 2 trees
	if (size() < 2)
		return;
    if (verbose_mode >= VB_MED) {
#ifdef USE_HASH_MAP
	cout << "Using hash_map" << endl;
#else
	cout << "Using map" << endl;
#endif
    }
	cout << "Computing Robinson-Foulds distance..." << endl;

	vector<string> taxname(front()->leafNum);
	vector<SplitIntMap*> hs_vec;
	vector<SplitGraph*> sg_vec;

	front()->getTaxaName(taxname);


	// converting trees into split system then stored in SplitIntMap for efficiency
	for (iterator it = begin(); it != end(); it++) {
		SplitGraph *sg = new SplitGraph();
		SplitIntMap *hs = new SplitIntMap();

		(*it)->convertSplits(taxname, *sg);
		// make sure that taxon 0 is included
		for (SplitGraph::iterator sit = sg->begin(); sit != sg->end(); sit++) {
			if (!(*sit)->containTaxon(0)) (*sit)->invert();
			hs->insertSplit((*sit), 1);
		}
		hs_vec.push_back(hs);
		sg_vec.push_back(sg);
	}

	// now start the RF computation
	int id = 0;
	for (vector<SplitIntMap*>::iterator hsit = hs_vec.begin(); hsit+1 != hs_vec.end(); hsit++, id++) {
		vector<SplitIntMap*>::iterator end_it = hs_vec.end();
		if (mode == RF_ADJACENT_PAIR) end_it = hsit+2;
		int id2 = id+1;
		for (vector<SplitIntMap*>::iterator hsit2 = hsit+1; hsit2 != end_it; hsit2++, id2++) {
			int diff_splits = 0;
			SplitIntMap::iterator spit;
			for (spit = (*hsit2)->begin(); spit != (*hsit2)->end(); spit++) {
				if (spit->first->getWeight() >= weight_threshold && !(*hsit)->findSplit(spit->first)) diff_splits++;
			}
			for (spit = (*hsit)->begin(); spit != (*hsit)->end(); spit++) {
				if (spit->first->getWeight() >= weight_threshold && !(*hsit2)->findSplit(spit->first)) diff_splits++;
			}
			//int rf_val = (*hsit)->size() + (*hsit2)->size() - 2*common_splits;
			int rf_val = diff_splits;
			if (mode == RF_ADJACENT_PAIR) 
				rfdist[id] = rf_val;
			else {
				rfdist[id*size() + id2] = rfdist[id2*size() + id] = rf_val;
			}
		}
	}
	// delete memory 
	for (id = size()-1; id >= 0; id--) {
		delete hs_vec[id];
		delete sg_vec[id];
	}
}


void MTreeSet::computeRFDist(double *rfdist, MTreeSet *treeset2, bool k_by_k,
	const char *info_file, const char *tree_file, double *incomp_splits)
{
    if (verbose_mode >= VB_MED) {
	// exit if less than 2 trees
#ifdef USE_HASH_MAP
	cout << "Using hash_map" << endl;
#else
	cout << "Using map" << endl;
#endif
    }
	ofstream oinfo;
	ofstream otree;
	if (info_file) oinfo.open(info_file);
	if (tree_file) otree.open(tree_file);
	if (incomp_splits) memset(incomp_splits, 0, size()*treeset2->size()*sizeof(double));

	vector<string> taxname(front()->leafNum);
	vector<SplitIntMap*> hs_vec;
	vector<SplitGraph*> sg_vec;
	vector<NodeVector> nodes_vec;

	front()->getTaxaName(taxname);


	iterator it;
	// converting trees into split system then stored in SplitIntMap for efficiency
	for (iterator it = begin(); it != end(); it++) {
		SplitGraph *sg = new SplitGraph();
		SplitIntMap *hs = new SplitIntMap();
		NodeVector nodes;

		(*it)->convertSplits(taxname, *sg, &nodes);
		// make sure that taxon 0 is included
		int i = 0;
		for (SplitGraph::iterator sit = sg->begin(); sit != sg->end(); sit++, i++) {
			if (!(*sit)->containTaxon(0)) (*sit)->invert();
			hs->insertSplit((*sit), i);
		}
		hs_vec.push_back(hs);
		sg_vec.push_back(sg);
		nodes_vec.push_back(nodes);
	}

	// converting trees into split system then stored in SplitIntMap for efficiency
	for (it = treeset2->begin(); it != treeset2->end(); it++) {
		SplitGraph *sg = new SplitGraph();
		SplitIntMap *hs = new SplitIntMap();
		NodeVector nodes;

		(*it)->convertSplits(taxname, *sg, &nodes);
		// make sure that taxon 0 is included
		int i = 0;
		for (SplitGraph::iterator sit = sg->begin(); sit != sg->end(); sit++, i++) {
			if (!(*sit)->containTaxon(0)) (*sit)->invert();
			hs->insertSplit((*sit), i);
		}
		hs_vec.push_back(hs);
		sg_vec.push_back(sg);
		nodes_vec.push_back(nodes);
	}

	// now start the RF computation
	int id = 0;
	int col_size = hs_vec.size() - size();
	for (vector<SplitGraph*>::iterator hsit = sg_vec.begin(); id < size(); hsit++, id++) {
		int id2 = 0;
        vector<SplitIntMap*>::iterator start_it = (hs_vec.begin() + size());
        vector<SplitIntMap*>::iterator end_it = (hs_vec.end());
        if (k_by_k) {
            // only distance between k-th tree
            start_it = start_it + id;
            end_it = start_it + 1;
        }

		for (vector<SplitIntMap*>::iterator hsit2 = start_it; hsit2 != end_it; hsit2++, id2++) {
			int common_splits = 0;
			int i = 0;
			for (SplitGraph::iterator spit = (*hsit)->begin(); spit != (*hsit)->end(); spit++, i++) {
				if ((*hsit2)->findSplit(*spit)) {
					common_splits++;
					if (info_file && (*spit)->trivial()<0) oinfo << " " << nodes_vec[id][i]->name;
				} else {
					if (info_file && (*spit)->trivial()<0) oinfo << " -" << nodes_vec[id][i]->name;
					nodes_vec[id][i]->name = "-" + nodes_vec[id][i]->name;
					/*if (incomp_splits && !sg_vec[id2+size()]->compatible(*spit))
						nodes_vec[id][i]->name = "-" + nodes_vec[id][i]->name;*/
				} 
			}
			double rf_val = (*hsit)->size() + (*hsit2)->size() - 2*common_splits;
            if (Params::getInstance().normalize_tree_dist) {
                int non_trivial = (*hsit)->size() - (*hsit)->getNTrivialSplits();
                non_trivial += (*hsit2)->size() - ((*hsit2)->begin())->first->getNTaxa();
                rf_val /= non_trivial;
            }
            if (k_by_k)
                rfdist[id] = rf_val;
            else
                rfdist[id*col_size + id2] = rf_val;
			if (info_file) oinfo << endl;
			if (tree_file) { at(id)->printTree(otree); otree << endl; }
			for (i = 0; i < nodes_vec[id].size(); i++)
				if (nodes_vec[id][i]->name[0] == '-') nodes_vec[id][i]->name.erase(0,1);
		}
		if (!incomp_splits || k_by_k) continue;
		id2 = 0;
		// count incompatible splits
		for (vector<SplitGraph*>::iterator hsit3 = sg_vec.begin()+size(); hsit3 != sg_vec.end(); hsit3++, id2++) {
			int num_incomp = 0;
			SplitGraph::iterator spit;
			for (spit = (*hsit)->begin(); spit != (*hsit)->end(); spit++) 
				if (!(*hsit3)->compatible(*spit)) num_incomp++;
			for (spit = (*hsit3)->begin(); spit != (*hsit3)->end(); spit++) 
				if (!(*hsit)->compatible(*spit)) num_incomp++;
					
			incomp_splits[id*col_size + id2] = num_incomp;
		}
	}
	// delete memory 
	for (id = hs_vec.size()-1; id >= 0; id--) {
		delete hs_vec[id];
		delete sg_vec[id];
	}

	if (info_file) {
		oinfo.close();
		cout << "Detailed split occurrences printed to " << info_file << endl;
	}
	if (tree_file) {
		otree.close();
		cout << "Detailed split occurrences on tree printed to " << tree_file << endl;
	}
}

int MTreeSet::sumTreeWeights() {
	int sum = 0;
	for (IntVector::iterator it = tree_weights.begin(); it != tree_weights.end(); it++)
		sum += (*it);
	return sum;
}

int MTreeSet::categorizeDistinctTrees(IntVector &category) {
	if (empty()) return 0;
	if (size() == 1) {
		category.resize(1,0);
		return 1;
	}
	StringIntMap tree_cat_map;
	string root_name = front()->root->name;
	int ncat = 0;
	category.resize(size(),-1);

	int id = 0;
	for (iterator it = begin(); it != end(); it++, id++) {
		(*it)->root = (*it)->findNodeName(root_name);
		if (!(*it)->root || !(*it)->root->isLeaf()) 
			outError("Internal error ", __func__);
		stringstream ostr;
		(*it)->printTree(ostr, WT_TAXON_ID | WT_SORT_TAXA);
		string str = ostr.str();
		//cout << str << endl;
		StringIntMap::iterator map_it = tree_cat_map.find(str);
		if (map_it == tree_cat_map.end()) { // not found
			category[id] = ncat;
			tree_cat_map[str] = ncat;
			ncat++;
		} else {
			category[id] = map_it->second;
		}
	}
	return ncat;
}


/*int MTreeSet::categorizeDistinctTrees(IntVector &category) {
	// exit if less than 2 trees
	if (empty()) return 0;
	if (size() == 1) {
		category.resize(1,0);
		return 1;
	}
#ifdef USE_HASH_MAP
	cout << "Using hash_map" << endl;
#else
	cout << "Using map" << endl;
#endif
	cout << "Checking duplicated trees..." << endl;

	vector<string> taxname(front()->leafNum);
	vector<SplitIntMap*> hs_vec;
	vector<SplitGraph*> sg_vec;

	front()->getTaxaName(taxname);


	// converting trees into split system then stored in SplitIntMap for efficiency
	for (iterator it = begin(); it != end(); it++) {
		SplitGraph *sg = new SplitGraph();
		SplitIntMap *hs = new SplitIntMap();

		(*it)->convertSplits(taxname, *sg);
		// make sure that taxon 0 is included
		for (SplitGraph::iterator sit = sg->begin(); sit != sg->end(); sit++) {
			if (!(*sit)->containTaxon(0)) (*sit)->invert();
			hs->insertSplit((*sit), 1);
		}
		hs_vec.push_back(hs);
		sg_vec.push_back(sg);
	}

	// now start the RF computation
	int id = 0, ncat = 0;
	category.resize(size(),-1);
	for (vector<SplitIntMap*>::iterator hsit = hs_vec.begin(); hsit != hs_vec.end(); hsit++, id++) 
	if (category[id] < 0) {
		category[id] = ncat;
		int id2 = id+1;
		for (vector<SplitIntMap*>::iterator hsit2 = hsit+1; hsit2 != hs_vec.end(); hsit2++, id2++) 
		if (category[id2] < 0) {
			bool equal = true;
			for (SplitIntMap::iterator spit = (*hsit2)->begin(); spit != (*hsit2)->end(); spit++) {
				if (!(*hsit)->findSplit(spit->first)) { equal = false; break; }
			}
			if (equal) category[id2] = ncat;
		}
		ncat++;
	}
	// delete memory 
	for (id = size()-1; id >= 0; id--) {
		delete hs_vec[id];
		delete sg_vec[id];
	}
	return ncat;
}

*/
