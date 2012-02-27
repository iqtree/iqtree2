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

void readPatternLogLL(Alignment* aln, char *fileName, vector<DoubleVector> &logLLs)
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

void computeExpectedNorFre(Alignment *aln, DoubleVector &logLL, IntVector &expectedNorFre)
{
	//IntVector expectedNorFre;
	if ( logLL.empty()) 
		outError("Error: log likelihood of patterns are not given!");

	int patNum = aln->getNPattern();
	int alignLen = aln->getNSite();		
	//resize the expectedNorFre vector
	expectedNorFre.resize(patNum,-1);

	//Vector containing the likelihood of the pattern p_i
	DoubleVector LL(patNum,-1.0);
	double sumLL = 0; //sum of the likelihood of the patterns in the alignment

	//Compute the likelihood from the logLL
	for ( int i = 0; i < patNum; i++ )
	{
		LL[i] = exp(logLL[i]);
		sumLL += LL[i];
	}

	//Vector containing l_i = p_i*ell/sum_i(p_i)
	DoubleVector ell(patNum, -1.0);
	//Compute l_i
	for ( int i = 0; i < patNum; i++ )
	{
		ell[i] = (double)alignLen * LL[i] / sumLL;
	}


	//Vector containing r_i where r_0 = ell_0; r_{i+1} = ell_{i+1} + r_i - ordinaryRounding(r_i)
	DoubleVector r(patNum, -1.0);
	//Compute r_i and the expected normalized frequencies
	r[0] = ell[0];
	expectedNorFre[0] = (int)floor(ell[0]+0.5); //note that floor(_number+0.5) returns the ordinary rounding of _number
	int sum = expectedNorFre[0];
	for (int j = 1; j < patNum; j++ )
	{
		r[j] = ell[j] + r[j-1] - floor(r[j-1]+0.5);
		expectedNorFre[j] = (int)floor(r[j]+0.5);
		sum += expectedNorFre[j];
	}
	
	//cout << "Number of patterns: " << patNum << ", sum of expected sites: " << sum << endl;
	//return expectedNorFre;
}

void computeTreeWeights(DoubleVector &reProb, IntVector &reW) {
	int nDiff = reProb.size();
	reW.resize(nDiff,-1);
	DoubleVector ratio(nDiff,-1.0);
	double sumRatio = 0;
	int i;
	double max_prob = reProb[0];
	for ( i = 0; i < nDiff; i++ ) 
		if (reProb[i] > max_prob) max_prob = reProb[i];
	
	for ( i = 0; i < nDiff; i++ )
	{
		ratio[i] = exp(reProb[i]-max_prob);
		sumRatio += ratio[i];
	}
	for ( i = 0; i < nDiff; i++ )
	{
		double temp = (ratio[i]/sumRatio)*1000000;
		reW[i] = (int) floor(temp+0.5);
	}	
}

double euclideanDist(IntVector &vec1, IntVector &vec2) {
	if (vec1.size() != vec2.size()) outError("Different vector size ", __func__);
	double dist = 0.0;
	for (int i = 0; i < vec1.size(); i++) 
		dist += (vec1[i]-vec2[i])*(vec1[i]-vec2[i]);
	return sqrt(dist);
}

double computeRELLLogL(DoubleVector &pattern_lh, IntVector &pattern_freq) {
	double lh = 0.0;
	int npat = pattern_lh.size();
	if (npat != pattern_freq.size()) outError("Wrong vector size ", __func__);
	for (int i = 0; i < npat; i++) lh += pattern_freq[i] * pattern_lh[i];
	return lh;
}

/**
	computing Expected Likelihood Weights (ELW) of trees by Strimmer & Rambaut (2002)
*/
void computeExpectedLhWeights(Alignment *aln, vector<DoubleVector> &pattern_lhs, 
	IntVector &treeids, int num_replicates, DoubleVector &elw, DoubleVector *sh_pval = NULL) {
	cout << "Computing Expected Likelihood Weights (ELW) with " << num_replicates << " replicates ..." << endl;
	int i, j, ntrees = treeids.size();
	elw.resize(treeids.size(), 0.0);
	vector<DoubleVector> all_logl;
	// general RELL logl
	for (i = 0; i < num_replicates; i++) {
		IntVector pattern_freq;
		aln->resamplePatternFreq(pattern_freq);
		DoubleVector logl;
		logl.resize(treeids.size(), 0.0);
		j = 0;
		for (IntVector::iterator it = treeids.begin(); it != treeids.end(); it++, j++) {
			logl[j] = computeRELLLogL(pattern_lhs[*it], pattern_freq);
		}
		if (sh_pval) all_logl.push_back(logl);
		double max_logl = logl[0];
		for (j = 0; j < logl.size(); j++) 
			if (max_logl < logl[j]) max_logl = logl[j];
		double sum = 0.0;
		for (j = 0; j < logl.size(); j++) {
			logl[j] = exp(logl[j] - max_logl);
			sum += logl[j];
		}
		for (j = 0; j < logl.size(); j++) 
			elw[j] += (logl[j]/sum);
	}
	// normalize ELW weights to sum of 1
	for (j = 0; j < elw.size(); j++) 
		elw[j] /= num_replicates;

	if (!sh_pval) return;


	// centering step in SH test
	DoubleVector mean_logl;
	mean_logl.resize(ntrees, 0);
	for (i = 0; i < num_replicates; i++) 
		for (j = 0; j < ntrees; j++) {
			mean_logl[j] += all_logl[i][j];
		}
	for (j = 0; j < ntrees; j++) 
		mean_logl[j] /= num_replicates;
	for (i = 0; i < num_replicates; i++) 
		for (j = 0; j < ntrees; j++) {
			all_logl[i][j] -= mean_logl[j];
		}

	// computing delta
	for (i = 0; i < num_replicates; i++) {
		double max_logl = *max_element(all_logl[i].begin(), all_logl[i].end());
		for (j = 0; j < ntrees; j++) all_logl[i][j] = max_logl - all_logl[i][j];
	}

	// computing original delta
	DoubleVector orig_logl;
	orig_logl.resize(ntrees, 0);
	for (j = 0; j < ntrees; j++) {
		int tree_id = treeids[j];
		i = 0;
		for (Alignment::iterator it = aln->begin(); it != aln->end(); it++, i++)
			orig_logl[j] += pattern_lhs[tree_id][i] * it->frequency;
	}
	double max_logl = *max_element(orig_logl.begin(), orig_logl.end());
	for (j = 0; j < ntrees; j++) orig_logl[j] = max_logl - orig_logl[j];
	sh_pval->resize(ntrees, 0);
	for (i = 0; i < num_replicates; i++) 
	for (j = 0; j < ntrees; j++) {
			if (orig_logl[j] < all_logl[i][j]) (*sh_pval)[j] += 1.0;
	}
	for (j = 0; j < ntrees; j++) 
		(*sh_pval)[j] /= num_replicates;
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
    int i = 0, j;
    for (NodeVector::iterator it = taxa.begin(); it != taxa.end(); it++) {
        (*it)->id = i++;
    }

	if (verbose_mode >= VB_DEBUG) {
		cout << "Original pattern freq: ";
		for (i = 0; i < alignment->getNPattern(); i++)
			cout << alignment->at(i).frequency << " ";
		cout << endl;
	}

	// read in trees file
	MTreeSet trees(params.user_file, params.is_rooted, params.tree_burnin);
	vector<DoubleVector> pattern_lhs;
	vector<IntVector> expected_freqs;

	// read in corresponding site-log-likelihood for all trees
	readPatternLogLL(alignment, params.siteLL_file, pattern_lhs);
	cout << pattern_lhs.size() << " log-likelihood vectors loaded" << endl;
	if (pattern_lhs.size() != trees.size()) outError("Different number of sitelh vectors");

	// get distinct trees
	int ntrees = trees.size();
	IntVector::iterator it;
	IntVector diff_tree_ids;
	IntVector tree_category;
	trees.categorizeDistinctTrees(tree_category);
	for (i = 0; i < ntrees; i++) {
		int cat = tree_category[i];
		if (diff_tree_ids.empty() || tree_category[diff_tree_ids.back()] < cat)
			diff_tree_ids.push_back(i);
	}
/*
	int *rfdist;
	rfdist = new int [ntrees*ntrees];
	memset(rfdist, 0, ntrees*ntrees* sizeof(int));
	trees.computeRFDist(rfdist, RF_ALL_PAIR);*/
	// identify distinct trees
	/*
	IntVector checked;
	checked.resize(ntrees, 0);
	for (i = 0; i < ntrees; i++) if (!checked[i]) {
		for (j = i; j < ntrees; j++) if (!checked[j]) {
			if (rfdist[i*ntrees + j] == 0) checked[j] = 1;
		}
		diff_tree_ids.push_back(i);
	}*/
	int ndiff = diff_tree_ids.size();
	cout << diff_tree_ids.size() << " distinct trees detected" << endl;

	// compute multinomial probability for every distinct tree
	DoubleVector prob_vec;
	for (it = diff_tree_ids.begin(); it != diff_tree_ids.end(); it++) {
		double prob;
		alignment->multinomialProb(pattern_lhs[*it], prob);
		prob_vec.push_back(prob);
		IntVector expected_freq;
		computeExpectedNorFre(alignment, pattern_lhs[*it], expected_freq);
		expected_freqs.push_back(expected_freq);
		if (verbose_mode >= VB_DEBUG) {
			for (i = 0; i < expected_freq.size(); i++)
				cout << expected_freq[i] << " ";
			cout << endl;
		}
	}

	IntVector diff_tree_weights;


	if (params.use_elw_method) { 	// compute ELW weights
		DoubleVector elw, sh_pval;
		computeExpectedLhWeights(alignment, pattern_lhs, diff_tree_ids, params.gbo_replicates, elw, &sh_pval);
		string elw_file_name = params.out_prefix;
		elw_file_name += ".elw";
		ofstream elw_file(elw_file_name.c_str());
		elw_file << "Treeid\tELW\tSH-pval" << endl;
		for (i = 0; i < elw.size(); i++) 
			elw_file << diff_tree_ids[i]+1 << "\t" << elw[i] << "\t" << sh_pval[i] << endl;
		elw_file.close();
		cout << "ELW printed to " << elw_file_name << endl;
		diff_tree_weights.resize(diff_tree_ids.size(), 0);
		for (i = 0; i < diff_tree_ids.size(); i++)
			diff_tree_weights[i] = round(elw[i]*1000000);
	} else {
		double own_prob;
		alignment->multinomialProb(*alignment, own_prob);
		cout << "Own prob: " << own_prob << endl;
	
		cout << "Conducting " << params.gbo_replicates << " non-parametric resampling..." << endl;
		if (params.use_rell_method) cout << "Use RELL method (resampling estimated log-likelihoods)" << endl;
		// generate bootstrap samples
		for (i = 0; i < params.gbo_replicates; i++) {
			Alignment* bootstrap_alignment;
			if (alignment->isSuperAlignment())
				bootstrap_alignment = new SuperAlignment;
			else
				bootstrap_alignment = new Alignment;
			IntVector pattern_freq;
			bootstrap_alignment->createBootstrapAlignment(alignment, &pattern_freq);
			double prob;
			bootstrap_alignment->multinomialProb(*alignment, prob);
			double min_dist = -1.0;
			int chosen_id = -1;
			// select best-fit tree by euclidean distance
			for (j = 0; j < expected_freqs.size(); j++) {
				double dist = euclideanDist(pattern_freq, expected_freqs[j]);
				//cout << dist << " ";
				if (dist < min_dist || min_dist < 0) {
					min_dist = dist;
					chosen_id = j;
				}
			}
	
			// select best-fit tree by RELL method 
			double max_logl = 1, max_logl_saved = 0;
			int chosen_id2 = -1, chosen_id2_saved = -1;
			for (j = 0; j < ndiff; j++) {
				int tree_id = diff_tree_ids[j];
				double logl = computeRELLLogL(pattern_lhs[tree_id], pattern_freq);
				//if (verbose_mode >= VB_MAX) cout << logl << endl;
				if (max_logl > 0 || max_logl < logl) {
					max_logl_saved = max_logl;
					max_logl = logl;
					chosen_id2_saved = chosen_id2;
					chosen_id2 = j;
				}
				
			}
			if (verbose_mode >= VB_MED) {
				cout << "Bootstrap " << i+1 << " choose id=" << diff_tree_ids[chosen_id]+1 // <<" dist=" << min_dist 
					<< " id2=" << diff_tree_ids[chosen_id2]+1 << " lh2=" << max_logl 
					<< " id2s=" << diff_tree_ids[chosen_id2_saved]+1 << " lh2s=" << max_logl_saved
					<< " prob=" << prob << endl;
			}
			if (params.use_rell_method)
				diff_tree_ids.push_back(diff_tree_ids[chosen_id2]);
			else
				diff_tree_ids.push_back(diff_tree_ids[chosen_id]);
			prob_vec.push_back(prob);
	
			if (verbose_mode >= VB_DEBUG) {
				for (j = 0; j < pattern_freq.size(); j++) 
					cout << pattern_freq[j] << " ";
				cout << endl;
			}
			delete bootstrap_alignment;
		}
	
		// compute tree weights from the log-probability
		computeTreeWeights(prob_vec, diff_tree_weights);

	} // end of Arndt's method


	int max_weight = 0;
	for (i = 0; i < diff_tree_weights.size(); i++)
		if (diff_tree_weights[i] > max_weight) max_weight = diff_tree_weights[i];
	cout << "Max weight: " << max_weight << endl;
	if (verbose_mode >= VB_MAX) {
		for (i = 0; i < prob_vec.size(); i++) cout << prob_vec[i] << " " << diff_tree_weights[i] << " " << diff_tree_ids[i]+1 << endl;
		cout << endl;
	}

	for (i = 0; i < ntrees; i++) trees.tree_weights[i] = 0;
	for (it = diff_tree_ids.begin(), i = 0; it != diff_tree_ids.end(); it++, i++) {
		trees.tree_weights[*it] += diff_tree_weights[i];
	}
	cout << "Sum weight: " << trees.sumTreeWeights() << endl;

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
	out_file = params.out_prefix;
	out_file += ".suptree";

    tree.printTree(out_file.c_str());
    cout << "Tree with assigned bootstrap support written to " << out_file << endl;

	out_file = params.out_prefix;
	out_file += ".supval";
	tree.writeInternalNodeNames(out_file);

    cout << "Support values written to " << out_file << endl;
	//delete [] rfdist;
}

