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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "phylotree.h"
#include "phylosupertree.h"
#include "phyloanalysis.h"
#include "alignment.h"
#include "superalignment.h"
#include "iqptree.h"
#include "gtrmodel.h"
#include "modeldna.h"
#include "myreader.h"
#include "rateheterogeneity.h"
#include "rategamma.h"
#include "rateinvar.h"
#include "rategammainvar.h"
#include "modeltest_wrapper.h"
#include "modelprotein.h"
#include "stoprule.h"

#include "mtreeset.h"
#include "mexttree.h"
#include "ratemeyerhaeseler.h"
#include "whtest_wrapper.h"
#include "partitionmodel.h"

void readLogLL(Alignment* aln, char *fileName, vector<DoubleVector> &logLLs)
{
	//First read the values from inFile to a DoubleVector
	//int siteNum;
	string currentString;
	cout << "\nReading file containing site's loglikelihood: " << fileName << "...." << endl;
    ifstream inFile;
	try{
		inFile.exceptions (ios::failbit | ios::badbit);
		inFile.open(fileName);
		/**really start reading*/
		//read number of sites
		getline(inFile,currentString);
		//siteNum = convert_int(currentString.c_str());
		//ignore "Site_Lh"		
		inFile.exceptions (ios::badbit);
		while (!inFile.eof())
		{
			DoubleVector _logllVec;
			if ( !(inFile >> currentString) ) break;
			//reading each line of the file
			//remove the badbit
			//set the failbit again
			for (int i = 0; i < aln->getNSite(); i++) {
				double ll;
				if (!(inFile >> ll)) throw "Wrong logLL entry";
				_logllVec.push_back(ll);
			}
			DoubleVector logLL;
			logLL.resize(aln->getNPattern(),0.0);
			for (int i = 0; i < _logllVec.size(); i++)
			{
				int patIndex = aln->getPatternID(i);
				if ( logLL[patIndex] == 0 )
					logLL[patIndex] = _logllVec[i];
				else
					if ( logLL[patIndex] != _logllVec[i] )
						outError("Conflicting between the likelihoods reported for pattern", aln->at(i));
			}
			logLLs.push_back(logLL);
		}/**finish reading*/
		inFile.clear();
		inFile.exceptions (ios::failbit | ios::badbit);
		inFile.close();
	} catch(bad_alloc){
			outError(ERR_NO_MEMORY);
	} catch (const char *str){
			outError(str);
	} catch (char *str){
			outError(str);
	} catch (string str){
			outError(str);
	} catch (ios::failure){
			outError(ERR_READ_INPUT);
	} catch (...){
			outError(ERR_READ_ANY);
	}

}

void computeTreeWeights(DoubleVector &reProb, IntVector &reW) {
	int nDiff = reProb.size();
	reW.resize(nDiff,-1);
	DoubleVector ratio(nDiff,-1.0);
	double sumRatio = 0;
	int i;
	for ( i = 0; i < nDiff; i++ )
	{
		ratio[i] = exp(reProb[i]-reProb[0]);
		sumRatio += ratio[i];
	}
	for ( i = 0; i < nDiff; i++ )
	{
		double temp = (ratio[i]/sumRatio)*1000000;
		reW[i] = (int) floor(temp+0.5);
	}	
}

void runGuidedBootstrap(Params &params, string &original_model, Alignment *alignment, IQPTree &tree) {
    if (!params.user_file) {
		outError("You have to specify user tree file");
    }
    if (!params.siteLL_file) {
		outError("Please provide site log-likelihood file via -gbo option");
    }
    if (!params.second_tree) {
		outError("Please provide target tree file via -sup option");
    }

	// read tree file
	cout << "Reading tree file " << params.second_tree << endl;
	tree.readTree(params.second_tree, params.is_rooted);
    // reindex the taxa in the tree to aphabetical names
    NodeVector taxa;
    tree.getTaxa(taxa);
    sort(taxa.begin(), taxa.end(), nodenamecmp);
    int i = 0;
    for (NodeVector::iterator it = taxa.begin(); it != taxa.end(); it++) {
        (*it)->id = i++;
    }

	// read in trees file
	MTreeSet trees(params.user_file, params.is_rooted, params.tree_burnin);
	vector<DoubleVector> site_lhs;

	// read in corresponding site-log-likelihood for all trees
	readLogLL(alignment, params.siteLL_file, site_lhs);
	cout << site_lhs.size() << " log-likelihood vectors loaded" << endl;
	if (site_lhs.size() != trees.size()) outError("Different number of sitelh vectors");

	// compute RF distance
	int ntrees = trees.size(), j;
	IntVector::iterator it;
	int *rfdist;
	rfdist = new int [ntrees*ntrees];
	memset(rfdist, 0, ntrees*ntrees* sizeof(int));
	trees.computeRFDist(rfdist, RF_ALL_PAIR);
	IntVector diff_tree_ids;
	IntVector checked;
	checked.resize(ntrees, 0);
	
	// identify distinct trees
	for (i = 0; i < ntrees; i++) if (!checked[i]) {
		for (j = i; j < ntrees; j++) if (!checked[j]) {
			if (rfdist[i*ntrees + j] == 0) checked[j] = 1;
		}
		diff_tree_ids.push_back(i);
	}
	cout << diff_tree_ids.size() << " distinct trees detected" << endl;

	// compute multinomial probability for every distinct tree
	DoubleVector prob_vec;
	for (it = diff_tree_ids.begin(); it != diff_tree_ids.end(); it++) {
		double prob;
		alignment->multinomialProb(site_lhs[*it], prob);
		prob_vec.push_back(prob);
	}

	for (i = 0; i < params.gbo_replicates; i++) {
		Alignment* bootstrap_alignment;
		if (alignment->isSuperAlignment())
			bootstrap_alignment = new SuperAlignment;
		else
			bootstrap_alignment = new Alignment;
		bootstrap_alignment->createBootstrapAlignment(alignment);
	}

	// compute tree weights from the log-probability
	IntVector diff_tree_weights;
	computeTreeWeights(prob_vec, diff_tree_weights);
	for (i = 0; i < ntrees; i++) trees.tree_weights[i] = 0;
	for (it = diff_tree_ids.begin(), i = 0; it != diff_tree_ids.end(); it++, i++) {
		trees.tree_weights[*it] = diff_tree_weights[i];
	}
	cout << "sum weights: " << trees.sumTreeWeights() << endl;

	// assign bootstrap support
    SplitGraph sg;
    SplitIntMap hash_ss;
    // make the taxa name
    vector<string> taxname;
    taxname.resize(tree.leafNum);
    tree.getTaxaName(taxname);

    trees.convertSplits(taxname, sg, hash_ss, SW_COUNT, -1);
    cout << sg.size() << " splits found" << endl;
    // compute the percentage of appearance
    sg.scaleWeight(100.0 / trees.sumTreeWeights(), true);
    //	printSplitSet(sg, hash_ss);
    //sg.report(cout);
    cout << "Creating bootstrap support values..." << endl;
	MExtTree mytree(tree);
    mytree.createBootstrapSupport(taxname, trees, sg, hash_ss);
	tree.init(mytree);

	tree.setAlignment(alignment);

    string out_file;
    if (params.out_file)
        out_file = params.out_file;
    else {
        out_file = params.user_file;
        out_file += ".suptree";
    }

    tree.printTree(out_file.c_str());
    cout << "Tree with assigned bootstrap support written to " << out_file << endl;


	delete [] rfdist;
}
