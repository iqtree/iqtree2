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
#include <numeric>
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

//#include "zpipe.h"
#include "gzstream.h"

void readPatternLogLL(Alignment* aln, char *fileName, vector<double*> &logLLs, DoubleVector &trees_logl)
{
	//First read the values from inFile to a DoubleVector
	//int siteNum;
	string currentString;
	cout << "\nReading file containing site's loglikelihood: " << fileName << "...." << endl;
    ifstream inFile;
    int i;
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
			double logl = 0.0;
			for (i = 0; i < aln->getNSite(); i++) {
				double ll;
				if (!(inFile >> ll)) throw "Wrong logLL entry";
				_logllVec.push_back(ll);
				logl += ll;
			}
			double *logLL = new double[aln->getNPattern()];
			memset(logLL, 0, sizeof(double) * aln->getNPattern());
			//logLL.resize(aln->getNPattern(),0.0);
			for (i = 0; i < _logllVec.size(); i++)
			{
				int patIndex = aln->getPatternID(i);
				if ( logLL[patIndex] == 0 )
					logLL[patIndex] = _logllVec[i];
				else
					if ( logLL[patIndex] != _logllVec[i] )
						outError("Conflicting between the likelihoods reported for pattern", aln->at(i));
			}
			logLLs.push_back(logLL);
			trees_logl.push_back(logl);
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

void computeExpectedNorFre(Alignment *aln, double *logLL, IntVector &expectedNorFre)
{
	//IntVector expectedNorFre;
/*	if ( logLL.empty()) 
		outError("Error: log likelihood of patterns are not given!");
*/

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

inline double computeRELL(double *pattern_lh, IntVector &pattern_freq) {
	double lh = 0.0;
	int npat = pattern_freq.size();
	//if (npat != pattern_freq.size()) outError("Wrong vector size ", __func__);
	for (int i = 0; i < npat; i++) lh += pattern_freq[i] * pattern_lh[i];
	return lh;
}

/**
	computing Expected Likelihood Weights (ELW) of trees by Strimmer & Rambaut (2002)
*/
void computeExpectedLhWeights(Alignment *aln, vector<double*> &pattern_lhs, 
	IntVector &treeids, int num_replicates, DoubleVector &elw, DoubleVector *sh_pval = NULL) {
	cout << "Computing Expected Likelihood Weights (ELW) with " << num_replicates << " replicates ..." << endl;
	int i, j, ntrees = treeids.size();
	elw.resize(treeids.size(), 0.0);
	vector<DoubleVector> all_logl;
	// general RELL logl
	for (i = 0; i < num_replicates; i++) {
		IntVector pattern_freq;
		aln->createBootstrapAlignment(pattern_freq);
		DoubleVector logl;
		logl.resize(treeids.size(), 0.0);
		j = 0;
		for (IntVector::iterator it = treeids.begin(); it != treeids.end(); it++, j++) {
			logl[j] = computeRELL(pattern_lhs[*it], pattern_freq);
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

void printTrees(const char *ofile, IQPTree &tree, IntVector *weights, bool compression) 
{
	int count = 0;
	try {
		ostream *out;
		if (compression) out = new ogzstream; else out = new ofstream;
		out->exceptions(ios::failbit | ios::badbit);
		if (compression) 
			((ogzstream*)out)->open(ofile);
		else
			((ofstream*)out)->open(ofile);
		(*out) << "[ scale=" << tree.len_scale << " ]" << endl;
		for (StringIntMap::iterator it = tree.treels.begin(); it != tree.treels.end(); it++) 
		if (!weights || weights->at(it->second)) {
			int id = it->second;
			out->precision(10);
			(*out) << "[ lh=" << tree.treels_logl[id];
			if (weights) (*out) << " w=" << weights->at(id);
			(*out) << " ] ";
			(*out) << tree.treels_newick[id] << endl;
			count++;
		}
		if (compression) 
			((ogzstream*)out)->close();
		else
			((ofstream*)out)->close();
		delete out;
		cout << count << " tree(s) printed to " << ofile << endl;
	} catch (ios::failure) {
		outError(ERR_WRITE_OUTPUT, ofile);
	}
}

void printPatternLh(const char *ofile, IQPTree *tree, bool compression) {
	int count = 0, i;
	int scale = 1000;
	try {
		ostream *out;
		if (compression) out = new ogzstream; else out = new ofstream;
		out->exceptions(ios::failbit | ios::badbit);
		if (compression) 
			((ogzstream*)out)->open(ofile/*, ios::out | ios::binary*/);
		else
			((ofstream*)out)->open(ofile/*, ios::out | ios::binary*/);
		int idfirst = tree->treels.begin()->second;
		(*out) << tree->treels.size() << " " << tree->aln->getNSite() <<
				" " << tree->aln->getNPattern() << " " << scale << endl;
		for (i = 0; i < tree->aln->getNSite(); i++)
			(*out) << " " << tree->aln->getPatternID(i);
		(*out) << endl;
		// DO NOT CHANGE 
		for (StringIntMap::iterator it = tree->treels.begin(); it != tree->treels.end(); it++) 
		{
			int id = it->second;
			assert(id < tree->treels_ptnlh.size());
			//out->write((char*)tree->treels_ptnlh[id], sizeof(double)*tree->aln->size());
			out->precision(10);
			(*out) << -tree->treels_logl[id];
			if (id == idfirst) {
				out->precision(6);
				for (int i = 0; i < tree->aln->size(); i++)
					(*out) << " " << -tree->treels_ptnlh[id][i];
			} else {
				for (int i = 0; i < tree->aln->size(); i++) {
					int diff = round((tree->treels_ptnlh[id][i]-tree->treels_ptnlh[idfirst][i])*scale);
					(*out) << " " << diff;
				}
			}
			(*out) << endl;
			count++;
		}
		if (compression) 
			((ogzstream*)out)->close();
		else
			((ofstream*)out)->close();
		delete out;
		cout << count << " pattern log-likelihood vector(s) printed to " << ofile << endl;
	} catch (ios::failure) {
		outError(ERR_WRITE_OUTPUT, ofile);
	}
}

void readPatternLh(const char *infile, IQPTree *tree, bool compression) {
	int count = 0, i;
	int ntrees, nsite, nptn, scale;
	double max_tol = 0.0;
	try {
		istream *in;
		if (compression) in = new igzstream; else in = new ifstream;
		in->exceptions(ios::failbit | ios::badbit);
		if (compression) 
			((igzstream*)in)->open(infile/*, ios::out | ios::binary*/);
		else
			((ifstream*)in)->open(infile/*, ios::out | ios::binary*/);
		(*in) >> ntrees >> nsite >> nptn >> scale;
		if (nsite != tree->aln->getNSite()) outError("Number of sites does not match");
		if (nptn !=  tree->aln->getNPattern()) outError("Number of patterns does not match");
		for (i = 0; i < nsite; i++) {
			int id;
			(*in) >> id;
			if (id != tree->aln->getPatternID(i)) outError("Pattern ID does not match");
		}
		tree->treels_logl.resize(ntrees, 0.0);
		tree->treels_ptnlh.resize(ntrees, NULL);
		for (int id = 0; id < ntrees; id++) 
		{
			double logl;
			(*in) >> logl;
			logl = -logl;
			tree->treels_logl[id] = logl;
			double *pattern_lh = new double[nptn];
			if (id == 0) {
				for (i = 0; i < nptn; i++) {
					(*in) >> pattern_lh[i];
					pattern_lh[i] = -pattern_lh[i];
				}
			} else {
				double sum = 0.0;
				for (i = 0; i < nptn; i++) {
					int diff;
					(*in) >> diff;
					pattern_lh[i] = tree->treels_ptnlh[0][i]+(double)diff/scale;
					sum += pattern_lh[i] * tree->aln->at(i).frequency;
				}
				max_tol = max(max_tol, fabs(sum-logl));
			}
			tree->treels_ptnlh[id] = pattern_lh;
			count++;
		}
		cout << "max tolerance = " << max_tol << endl;
		if (compression) 
			((igzstream*)in)->close();
		else
			((ifstream*)in)->close();
		delete in;
		cout << count << " pattern log-likelihood vector(s) read from " << infile << endl;
	} catch (ios::failure) {
		outError(ERR_READ_INPUT, infile);
	}
}

void computeAllPatternLh(Params &params, IQPTree &tree) {
	/* this part copied from phyloanalysis.cpp */
	tree.optimize_by_newton = params.optimize_by_newton;
	tree.sse = params.SSE;
    try {
		if (!tree.getModelFactory()) {
			if (tree.isSuperTree())
				tree.setModelFactory(new PartitionModel(params, (PhyloSuperTree*)&tree));
			else
				tree.setModelFactory(new ModelFactory(params, &tree));
		}
    } catch (string str) {
    	outError(str);
    } 
    tree.setModel(tree.getModelFactory()->model);
    tree.setRate(tree.getModelFactory()->site_rate);
    tree.setStartLambda(params.lambda);
	if (tree.isSuperTree()) ((PhyloSuperTree*)&tree)->mapTrees();

    int model_df = tree.getModel()->getNDim() + tree.getRate()->getNDim();
    cout << endl;
    cout << "Estimating model parameters for: " << tree.getModelName() << " (" << model_df << " free parameters)" << endl;
    cout << "Fixed branch lengths: " << ((params.fixed_branch_length) ? "Yes" : "No") << endl;
    /* optimize model parameters */
    cout << endl;
    cout << "Optimizing model parameters" << endl;
    double bestTreeScore = tree.getModelFactory()->optimizeParameters(params.fixed_branch_length);
    cout << "Log-likelihood of the current tree: " << bestTreeScore << endl;

    //Update tree score
    tree.curScore = bestTreeScore;
	if (tree.isSuperTree()) ((PhyloSuperTree*)&tree)->computeBranchLengths();
    stringstream best_tree_string;
    tree.printTree(best_tree_string, WT_TAXON_ID + WT_BR_LEN);

	cout << "Computing pattern log-likelihoods of trees in " << params.user_file << " ..." << endl;
	/* now compute the treels_ptnlh */
	try {
		istream *in;
		if (params.do_compression) in = new igzstream; else in = new ifstream;
		in->exceptions(ios::failbit | ios::badbit);
		if (params.do_compression) 
			((igzstream*)in)->open(params.user_file); 
		else 
			((ifstream*)in)->open(params.user_file);
		double max_logl_diff = 0.0;
		char ch;
		(*in) >> ch;
		if (ch == '[') {
			string str;
			(*in) >> str;
			if (str.substr(0,6) == "scale=") {
				tree.len_scale = convert_double(str.substr(6).c_str());
			}
			do {
				(*in) >> ch;
			} while (!in->eof() && ch != ']');
		} else in->unget();
		cout << "Applying branch length scaling: " << tree.len_scale << endl;

		while (!in->eof()) {
			in->exceptions(ios::goodbit);
			(*in) >> ch;
			if (in->eof()) break;
			in->exceptions(ios::failbit | ios::badbit);
			double expected_lh = 0.0;
			if (ch == '[') {
				string str;
				(*in) >> str;
				if (str.substr(0,3) == "lh=") {
					expected_lh = convert_double(str.substr(3).c_str());
				}
				do  {
					(*in) >> ch;
				} while (!in->eof() && ch != ']');
			} else in->unget();

			tree.freeNode();
			tree.readTree(*in, tree.rooted);
			tree.scaleLength(1.0/tree.len_scale); // scale the branch length
			tree.assignLeafNames();
			tree.initializeAllPartialLh();
			tree.clearAllPartialLH();
			if (tree.isSuperTree()) ((PhyloSuperTree*)&tree)->mapTrees();
			double *pattern_lh = new double [tree.aln->getNPattern()];
			if (!params.fixed_branch_length) {
				tree.curScore = tree.optimizeAllBranches();
				tree.computePatternLikelihood(pattern_lh);
			} else {
				tree.curScore = tree.computeLikelihood(pattern_lh);
			}
			if (expected_lh != 0.0)
				max_logl_diff = max(max_logl_diff, fabs(tree.curScore-expected_lh));
			tree.treels_ptnlh.push_back(pattern_lh);
			tree.treels_logl.push_back(tree.curScore);
			if (tree.treels_ptnlh.size() % 500 == 0) 
				cout << tree.treels_ptnlh.size() << " trees evaluated" << endl;
		}

		cout << "Maximal log-likelihood error is " << max_logl_diff << endl << endl; 

		if (params.do_compression) ((igzstream*)in)->close(); else ((ifstream*)in)->close();
		delete in;
	} catch (ios::failure) {
		outError(ERR_READ_INPUT, params.user_file);		
	}

	/* take back the current best tree */
	best_tree_string.seekg(0, ios::beg);
	tree.freeNode();
	tree.readTree(best_tree_string, tree.rooted);
	tree.assignLeafNames();
	tree.initializeAllPartialLh();
	tree.clearAllPartialLH();
	}

void runGuidedBootstrap(Params &params, string &original_model, Alignment *alignment, IQPTree &tree) {

	int i, j;

	clock_t begin_time = clock();

	MTreeSet trees;
	vector<double*> *pattern_lhs = NULL;
	vector<IntVector> expected_freqs;
	DoubleVector *trees_logl = NULL;
	IntVector diff_tree_ids;
	int ntrees = 0;
	IntVector::iterator it;

	if (!tree.save_all_trees) {
		if (!params.user_file) {
			outError("You have to specify user tree file");
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
		i = 0;
		for (NodeVector::iterator it = taxa.begin(); it != taxa.end(); it++) {
			(*it)->id = i++;
		}
		// read in corresponding site-log-likelihood for all trees
		/*trees_logl = new DoubleVector;
		pattern_lhs = new vector<double*>;
		readPatternLogLL(alignment, params.siteLL_file, *pattern_lhs, *trees_logl);*/

		if (params.siteLL_file) {
			// read pattern loglikelihoods from file
			readPatternLh(params.siteLL_file, &tree, params.do_compression);
		} else {
			// compute all pattern log-likelihoods
	        tree.setAlignment(alignment);
			computeAllPatternLh(params, tree);
		} 
		pattern_lhs = &tree.treels_ptnlh;
		trees_logl = &tree.treels_logl;
		if (!params.distinct_trees) {
			// read in trees file
			trees.init(params.user_file, params.is_rooted, params.tree_burnin);
	
			if (pattern_lhs->size() != trees.size()) 
				outError("Different number of sitelh vectors");
			// get distinct trees
			ntrees = trees.size();
			IntVector tree_category;
			trees.categorizeDistinctTrees(tree_category);
			for (i = 0; i < ntrees; i++) {
				int cat = tree_category[i];
				if (diff_tree_ids.empty() || tree_category[diff_tree_ids.back()] < cat)
					diff_tree_ids.push_back(i);
			}
			cout << diff_tree_ids.size() << " distinct trees detected" << endl;
		}

	} else {
		if (tree.treels_ptnlh.empty()) {
			cout << "New bootstrap is not applicable due to no candiate trees" << endl;
			return;
		}
		pattern_lhs = &tree.treels_ptnlh;
		trees_logl = &tree.treels_logl;
		cout << "logl_cutoff = " << tree.logl_cutoff << " after " << tree.max_candidate_trees <<" trees" << endl;
	} 

	if (diff_tree_ids.empty()) {
		diff_tree_ids.resize(pattern_lhs->size());
		ntrees = pattern_lhs->size();
		for (i = 0; i < ntrees; i++) diff_tree_ids[i] = i;
	}

	IntVector origin_freq;
	for (i = 0; i < alignment->getNPattern(); i++) 
		origin_freq.push_back(alignment->at(i).frequency);
	

	if (verbose_mode >= VB_DEBUG) {
		cout << "Original pattern freq: ";
		for (i = 0; i < alignment->getNPattern(); i++)
			cout << alignment->at(i).frequency << " ";
		cout << endl;
	}

	cout << pattern_lhs->size() << " log-likelihood vectors loaded" << endl;

	int ndiff = diff_tree_ids.size();

	// consider only 10,000 trees with highest likelihoods
	if (params.max_candidate_trees > 0 && ndiff > params.max_candidate_trees) {
		DoubleVector neg_logl;
		neg_logl.resize(ndiff);
		for (i = 0; i < ndiff; i++) 
			neg_logl[i] = -trees_logl->at(diff_tree_ids[i]);
		nth_element(neg_logl.begin(), neg_logl.begin() + params.max_candidate_trees, neg_logl.end());
		double logl_cutoff = -neg_logl[params.max_candidate_trees];
		IntVector diff_tree_ids_new;
		diff_tree_ids_new.reserve(params.max_candidate_trees);
		for (i = 0; i < ndiff; i++)
			if (trees_logl->at(diff_tree_ids[i]) > logl_cutoff)
				diff_tree_ids_new.push_back(diff_tree_ids[i]);
		diff_tree_ids = diff_tree_ids_new;
		ndiff = diff_tree_ids.size();
		cout << "Reduce to " << ndiff << " highest likelihood trees with cutoff " << logl_cutoff << endl;
	}

	IntVector orig_diff_tree_ids = diff_tree_ids;

	// compute multinomial probability for every distinct tree
	DoubleVector prob_vec;
	for (it = diff_tree_ids.begin(); it != diff_tree_ids.end(); it++) {
		double prob;
		alignment->multinomialProb((*pattern_lhs)[*it], prob);
		prob_vec.push_back(prob);
		IntVector expected_freq;
		computeExpectedNorFre(alignment, (*pattern_lhs)[*it], expected_freq);
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
		computeExpectedLhWeights(alignment, (*pattern_lhs), diff_tree_ids, params.gbo_replicates, elw, &sh_pval);
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
		//cout << "Own prob: " << own_prob << endl;
	
		cout << "Conducting " << params.gbo_replicates << " non-parametric resampling ";
		if (params.use_rell_method) 
			cout << "using RELL" << endl;
		else 
			cout << "using Euclidean distance" << endl;
		if (params.use_weighted_bootstrap)
			cout << "Multinomial weighting for bootstrap sample ";
		else
			cout << "Equal weighting for bootstrap sample ";

		if (params.use_max_tree_per_bootstrap) 
			cout << "and selecting one tree per bootstrap" << endl;
		else
			cout << "and selecting multiple trees per bootstrap" << endl;

		double accepted_diff = 0.5;
		cout << "Accepted logl difference: " << accepted_diff << endl;

		// generate bootstrap samples
		for (i = 0; i < params.gbo_replicates; i++) {
			IntVector pattern_freq;
			alignment->createBootstrapAlignment(pattern_freq);
			double prob;
			if (params.use_weighted_bootstrap)
				prob = alignment->multinomialProb(pattern_freq);
			else
				prob = 0;
			if (params.use_rell_method) {
				// select best-fit tree by RELL method 
				DoubleVector logl;
				logl.resize(ndiff);
				for (j = 0; j < ndiff; j++) {
					int tree_id = diff_tree_ids[j];
					logl[j] = computeRELL((*pattern_lhs)[tree_id], pattern_freq);
					//if (verbose_mode >= VB_MAX) cout << logl << endl;
				}
				DoubleVector::iterator max_logl = max_element(logl.begin(), logl.end());
				int k = 0;
				if (params.use_max_tree_per_bootstrap) {
					double logl_cutoff = *max_logl - accepted_diff;
					int num_max = 0;
					for (j = 0; j < ndiff; j++) 
						if (logl[j] >= logl_cutoff) num_max++;
					if (num_max == 1) {
						diff_tree_ids.push_back(diff_tree_ids[max_logl - logl.begin()]);
						prob_vec.push_back(prob);
					} else {
						int max_rand = floor(((double)rand() / RAND_MAX)*num_max);
						if (max_rand >= num_max) max_rand = num_max - 1;
						for (j = 0; j < ndiff && max_rand >= 0; j++) 
						if (logl[j] >= logl_cutoff) {
							max_rand--;
							if (max_rand < 0) {
								diff_tree_ids.push_back(diff_tree_ids[j]);
								prob_vec.push_back(prob);
								break;
							}
						}
					}
					if (verbose_mode >= VB_MAX) {
						cout << "Bootstrap " << i+1 <<  " lprob=" << prob << " max_logl=" << 
							*max_logl << " select " << diff_tree_ids[j]+1;
						if (num_max > 1)
						 	cout << "  tie broken " << num_max << endl;
						 else
						 	cout << endl;
					}
				} else {
					DoubleVector weights;
					weights.resize(ndiff);
					for (j = 0; j < ndiff; j++) weights[j] = exp(logl[j] - *max_logl);
					double sum = accumulate(weights.begin(), weights.end(), 0.0);
					for (j = 0; j < ndiff; j++) weights[j] /= sum;
					int max_id = max_element(weights.begin(), weights.end()) - weights.begin();

					double weight_cutoff = weights[max_id] * 0.001;
					for (j = 0; j < ndiff; j++) {
						if (weights[j] >= weight_cutoff) {
							diff_tree_ids.push_back(diff_tree_ids[j]);
							prob_vec.push_back(prob + log(weights[j]));
							k++;
						}
					}
					if (verbose_mode >= VB_MAX) 
						cout << "Bootstrap " << i+1 <<  " lprob=" << prob << " max_id=" << max_id << " max_w=" << 
							weights[max_id] << " " << k << " trees" << endl;
				}
			}
			else {
			// select best-fit tree by euclidean distance
				double min_dist = -1.0;
				int chosen_id = -1;
				for (j = 0; j < expected_freqs.size(); j++) {
					double dist = euclideanDist(pattern_freq, expected_freqs[j]);
					//cout << dist << " ";
					if (dist < min_dist || min_dist < 0) {
						min_dist = dist;
						chosen_id = j;
					}
				}
				diff_tree_ids.push_back(diff_tree_ids[chosen_id]);
				prob_vec.push_back(prob);
				if (verbose_mode >= VB_MAX) {
					cout << "Bootstrap " << i+1 << " choose id=" << diff_tree_ids[chosen_id]+1 // <<" dist=" << min_dist 
						<< " lprob=" << prob << endl;
				}
			}
	
			if (verbose_mode >= VB_DEBUG) {
				for (j = 0; j < pattern_freq.size(); j++) 
					cout << pattern_freq[j] << " ";
				cout << endl;
			}
		}
	
		// compute tree weights from the log-probability
		computeTreeWeights(prob_vec, diff_tree_weights);

	} // end of Arndt's method

	IntVector final_tree_weights;
	final_tree_weights.resize(ntrees, 0);
	//for (i = 0; i < ntrees; i++) trees.tree_weights[i] = 0;
	for (it = diff_tree_ids.begin(), i = 0; it != diff_tree_ids.end(); it++, i++) {
		final_tree_weights[*it] += diff_tree_weights[i];
	}

	// now load in the trees
	if (tree.save_all_trees) {
		trees.init(tree.treels, tree.rooted, final_tree_weights);
		string out_file = params.out_prefix;
		if (params.do_compression) {
			out_file += ".btrees.gz";
			printTrees(out_file.c_str(), tree, &final_tree_weights, params.do_compression);
			out_file = params.out_prefix;
			out_file += ".alltrees.gz";
			printTrees(out_file.c_str(), tree, NULL, params.do_compression);
			if (params.print_site_lh) {
				out_file = params.out_prefix;
				out_file += ".ptnlh.gz";
				printPatternLh(out_file.c_str(), &tree, params.do_compression);
			}
		}
	} else if (params.distinct_trees) {
		trees.init(params.user_file, params.is_rooted, params.tree_burnin, NULL, &final_tree_weights, params.do_compression);
		//trees.init(params.user_file, params.is_rooted, params.tree_burnin, NULL);
/*		if (pattern_lhs->size() != trees.size()) 
			outError("Different number of sitelh vectors");*/
	} 

	int sum_weights = trees.sumTreeWeights();
	if (verbose_mode >= VB_MED) {
		for (i = 0; i < trees.size(); i++)
			if (trees.tree_weights[i] > 0) 
				cout << "Tree " << i+1 << " weight= " << trees.tree_weights[i] * 100 / sum_weights << endl;
	}
	int max_tree_id = max_element(trees.tree_weights.begin(), trees.tree_weights.end()) - trees.tree_weights.begin();
	cout << "max_tree_id = " << max_tree_id+1 << "   max_weight = " << trees.tree_weights[max_tree_id];
	cout << " (" << trees.tree_weights[max_tree_id] * 100 / sum_weights << "%)"<< endl;
	// assign bootstrap support
    SplitGraph sg;
    SplitIntMap hash_ss;
    // make the taxa name
    vector<string> taxname;
    taxname.resize(tree.leafNum);
    tree.getTaxaName(taxname);

	if (!tree.save_all_trees)
    	trees.convertSplits(taxname, sg, hash_ss, SW_COUNT, -1);
    else
    	trees.convertSplits(taxname, sg, hash_ss, SW_COUNT, -1, false);

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

	if (!tree.save_all_trees) {
		out_file = params.out_prefix;
		out_file += ".suptree";
	
		tree.printTree(out_file.c_str());
		cout << "Tree with assigned bootstrap support written to " << out_file << endl;
	}

	out_file = params.out_prefix;
	out_file += ".splits";

	sg.saveFile(out_file.c_str(), true);
	cout << "Split supports printed to NEXUS file " << out_file << endl;

	out_file = params.out_prefix;
	out_file += ".supval";
	tree.writeInternalNodeNames(out_file);

    cout << "Support values written to " << out_file << endl;

	/*
	if (!tree.save_all_trees) {
		for (vector<double* >::reverse_iterator it = pattern_lhs->rbegin(); it != pattern_lhs->rend(); it++)
			delete [] (*it);
		delete pattern_lhs;
		delete trees_logl;
	}*/

	clock_t end_time = clock();
	
	cout << "Time for guided bootstrap: " << (end_time-begin_time)/CLOCKS_PER_SEC << " seconds" << endl << endl;
	//delete [] rfdist;
}

