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

#include "phyloanalysis.h"
#include "alignment.h"
#include "phylotree.h"
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

const int DNA_MODEL_NUM = 14;

string dna_model_names[DNA_MODEL_NUM] =
	{"JC", "F81", "K80", "HKY", "TNef", "TN", "K81", "K81uf",
	"TIMef", "TIM", "TVMef", "TVM", "SYM", "GTR"};

const int AA_MODEL_NUM = 11;

string aa_model_names[AA_MODEL_NUM] =
	{"Dayhoff", "mtMAM", "JTT", "WAG", "cpREV", "mtREV", "rtREV",
	"mtART", "mtZOA", "VT", "LG"};

/**
	testing the best-fit model
	return in params.freq_type and params.rate_type
*/
string modelTest(Params &params, PhyloTree *in_tree)
{
	int nstates = in_tree->aln->num_states;
	if (nstates != 4 && nstates != 20)
		outError("Modeltest only works for DNA or Protein");
	char LRT_model[10];
	char IC_model[10];
	multiset<string> model_list;
	string fscore_name = params.aln_file;
	string fmodel_name = params.aln_file;
	char model_arg[40] = "";
	fscore_name += ".modelscore";
	fmodel_name += ".modeltest";
	ofstream fmodel_test(fmodel_name.c_str());
	fmodel_test.close();
	int model, rate_type;

	string best_model;
	double best_lh = -1000000000.0;

	bool fscore_ok = false;

	if (nstates == 4) {
		ifstream fscore_test(fscore_name.c_str());
		if (fscore_test.is_open()) {
			fscore_test.close();
			cout << fscore_name << " exists, checking this file" << endl;
			/* check if this works */
			if (modeltest(model_arg, fscore_name.c_str(), fmodel_name.c_str(), LRT_model, IC_model) == EXIT_FAILURE)
			{
				cout << "Score file is invalid, assessing all models again" << endl;
				fmodel_test.open(fmodel_name.c_str());
				fmodel_test.close();
			} else fscore_ok = true;
		}
	}
	if (!fscore_ok) {
		/* first compute the likelihood score for all available models */
		ofstream fscore(fscore_name.c_str());
		if (!fscore.is_open())
			outError("cannot write to .modelscore file!");

		fscore << "Tree ";

		PhyloTree *tree_homo = new PhyloTree();
		tree_homo->optimize_by_newton = params.optimize_by_newton;
		tree_homo->copyPhyloTree(in_tree);

		PhyloTree *tree_hetero = new PhyloTree();
		tree_hetero->optimize_by_newton = params.optimize_by_newton;
		tree_hetero->copyPhyloTree(in_tree);

		RateHeterogeneity *rate_class[4];
		rate_class[0] = new RateHeterogeneity();
		rate_class[1] = new RateInvar(NULL);
		rate_class[2] = new RateGamma(4, NULL);
		rate_class[3] = new RateGammaInvar(4, NULL);
		GTRModel *subst_model;
		if (nstates == 4)
			subst_model = new ModelDNA("JC", FREQ_UNKNOWN, in_tree);
		else
			subst_model = new ModelProtein("WAG", FREQ_UNKNOWN, in_tree);

		int num_models = (nstates == 4) ? DNA_MODEL_NUM : AA_MODEL_NUM;

		cout << "Tesing " << num_models*4 << ((nstates==4)?" DNA":" protein") << " models..." << endl;
		for (model = 0; model < num_models; model++) {
			for (rate_type = 0; rate_type <= 3; rate_type+=1) {
				// initialize tree
				PhyloTree *tree;
				if (rate_type == 0) {
					tree = tree_homo;
				} else if (rate_type == 1) {
					tree = tree_homo;
				} else if (rate_type == 2) {
					tree = tree_hetero;
				} else {
					tree = tree_hetero;
				}
				// initialize model
				subst_model->init(((nstates==4)?dna_model_names[model].c_str():aa_model_names[model].c_str()), FREQ_UNKNOWN);
				subst_model->setTree(tree);
				tree->setModel(subst_model);

				// initialize rate
				tree->setRate(rate_class[rate_type]);
				rate_class[rate_type]->setTree(tree);

				// print some infos
				cout << "===> Testing ";
				cout.width(12);
				string str;
				str = subst_model->name;
				str += rate_class[rate_type]->name;
				cout << left << str << " (df=" << subst_model->getNDim()+rate_class[rate_type]->getNDim() << ")";

				// clear all likelihood values
				tree->clearAllPartialLh();

				// optimize model parameters
				double cur_lh = tree->optimizeModel();
				cout.precision(10);
				cout << ":  Log-likelihood " << cur_lh << endl;
				fscore << str << endl;
				fscore.precision(10);
				fscore << "1\t" << -cur_lh;
				fscore.precision(6);
				subst_model->writeParameters(fscore);
				rate_class[rate_type]->writeParameters(fscore);
				fscore << endl;

				if (cur_lh > best_lh) {
					best_lh = cur_lh;
					best_model = subst_model->name + rate_class[rate_type]->name;
				}
			}
		}

		delete tree_hetero;
		delete tree_homo;

		for (rate_type = 3; rate_type >= 0; rate_type--)
			delete rate_class[rate_type];

		fscore.close();

		/* now do the modeltest */
		if (nstates == 4)
		modeltest(model_arg, fscore_name.c_str(), fmodel_name.c_str(), LRT_model, IC_model);
	}

	if (nstates == 4) {
		cout << "Performing ModelTest 3.7 from David Posada..." << endl;
		cout << "  LRT:     " << LRT_model << endl;
		cout << "  AIC:     " << IC_model << endl;
		model_list.insert(LRT_model);
		model_list.insert(IC_model);

		/* applying other tests */
		sprintf(model_arg, " -n%d", in_tree->aln->getNSite());
		modeltest(model_arg, fscore_name.c_str(), fmodel_name.c_str(), LRT_model, IC_model);
		cout << "  AICc:    " << IC_model << endl;
		model_list.insert(IC_model);

		sprintf(model_arg, " -n%d -t%d", in_tree->aln->getNSite(), in_tree->aln->getNSeq());
		modeltest(model_arg, fscore_name.c_str(), fmodel_name.c_str(), LRT_model, IC_model);
		cout << "  AICc_t:  " << IC_model << endl;
		model_list.insert(IC_model);


		sprintf(model_arg, "-n%d -b ", in_tree->aln->getNSite());
		modeltest(model_arg, fscore_name.c_str(), fmodel_name.c_str(), LRT_model, IC_model);
		cout << "  BIC:     " << IC_model << endl;
		model_list.insert(IC_model);

		sprintf(model_arg, "-n%d -t%d -b", in_tree->aln->getNSite(), in_tree->aln->getNSeq());
		modeltest(model_arg, fscore_name.c_str(), fmodel_name.c_str(), LRT_model, IC_model);
		cout << "  BIC_t:   " << IC_model << endl;
		model_list.insert(IC_model);
	} else {
		// FOR protein: no modeltest exists yet
		model_list.insert(best_model);
	}
	/* use the model that is most frequently selected */
	if (model_list.size() > 1) {
		cout << "Tests do not agree on one single model" << endl;
		int max = 0;

		for (multiset<string>::iterator it = model_list.begin(); it != model_list.end(); it++)
			if (model_list.count(*it) > max) {
				max = model_list.count(*it);
				best_model = (*it);
			}
		cout << "Use the model that is most frequently selected: " << best_model << endl;
	} else {
		cout << "Best model: " << (*model_list.begin()) << endl;
	}

	return best_model;
	//return "GTR";
}

void runPhyloAnalysis(Params &params, /*TreesBlock *trees_block, */ Alignment *alignment) {

	clock_t t_begin, t_end;

/*
	cout << "Computing parsimony score..." << endl;
	for (int i = 0; i < trees_block->GetNumTrees(); i++) {
		stringstream strs(trees_block->GetTranslatedTreeDescription(i), ios::in | ios::out | ios::app);
		strs << ";";
		PhyloTree tree;
		bool myrooted = trees_block->IsRootedTree(i);
		tree.readTree(strs, myrooted);
		tree.setAlignment(alignment);
		int score = tree.computeParsimonyScore();
		cout << "Tree " << trees_block->GetTreeName(i) << " has parsimony score of " << score << endl;
	}
*/

	/* initialize tree, either by user tree or BioNJ tree */
	IQPTree tree;
	if (params.user_file) {
		// start the search with user-defined tree
		bool myrooted = params.is_rooted;
		tree.readTree(params.user_file, myrooted);
		tree.setAlignment(alignment);
		if (!tree.dist_matrix) {
			tree.dist_matrix = new double[alignment->getNSeq() * alignment->getNSeq()];
		}
		alignment->computeDist(tree.dist_matrix);
	} else {
		// stepwise addition
		//tree.growTreeML(alignment);
		// create BioNJ tree
		tree.computeBioNJ(params, alignment, tree.dist_matrix);
	}

	if (verbose_mode >= VB_MAX)
		tree.drawTree(cout, WT_INTNODE+WT_BR_LEN);

	/* Fix if negative branch lengths detected */
	double fixed_length = 0.01;
	int fixed_number;
	if (fixed_number = tree.fixNegativeBranch(fixed_length)) {
		cout << "WARNING: " << fixed_number << " branches have no/negative lengths and initialized to " << fixed_length << endl;
		if (verbose_mode >= VB_DEBUG) {
			tree.printTree(cout);
			cout << endl;
		}
	}
	t_begin=clock();

	bool test_only = params.model_name == "TESTONLY";
	/* initialize substitution model */
	if (params.model_name == "TEST" || params.model_name == "TESTONLY") {
		params.model_name = modelTest(params, &tree);
		if (test_only) {
			return;
			t_end=clock();
			params.run_time = (t_end-t_begin);
			printf("Time used: %8.6f seconds.\n", (double)params.run_time / CLOCKS_PER_SEC);
		}
	}
	tree.createModel(params);

	cout << "Model of evolution: " << tree.getModelName() << endl;
	cout << "Fixed branch lengths: " << ((params.fixed_branch_length) ? "Yes" : "No") << endl;
	cout << "Random seed: " << params.ran_seed << endl;

	if (params.parsimony) {
		int score = tree.computeParsimonyScore();
		cout << "User tree has parsimony score of " << score << endl;

		/* NNI with parsimony function */
		tree.searchNNI();
	}
	/* optimize model parameters */
	cout.precision(10);

	cout << "User tree has likelihood score of " << tree.computeLikelihood() << endl;

	cout << "Optimizing model parameters" << endl;
	double score2 = tree.optimizeModel(params.fixed_branch_length);
	cout << "Log-likelihood of the current tree: " << score2 << endl;
	double bestTreeScore = score2;
	/* Optimize branch lengths with likelihood function */
	//cout << "Optimizing branch lengths..." << endl;
	//cout << "Log-likelihood: " << tree.optimizeAllBranches() << endl;

	/* do NNI with likelihood function */
	if (params.iqpnni_iterations > 0) {
		cout << "Performing Nearest Neighbor Interchange..." << endl;
		//cout << "Current tree likelihood: " << tree.optimizeNNIBranches() << endl;
		//tree.optimizeAllBranches();

		//clock_t nniBeginClock, nniEndClock;
		//nniBeginClock = clock();
		double newScore = tree.optimizeNNI();
		//nniEndClock = clock();
		//printf("Time used for first NNI search: %8.6f seconds.\n", (double)(nniBeginClock - nniEndClock) / CLOCKS_PER_SEC);

		cout << "Tree likelihood after NNI: " << newScore << endl;
		if (newScore > bestTreeScore) {
			bestTreeScore = newScore;
			cout << "Found new best tree log-likelihood : " << bestTreeScore << endl;
		}
		else {
			cout << "Tree didn't improve after NNI :( " << endl;
		}

	}
	

	/* do the IQP */
	if (params.k_representative > 0 && params.p_delete > 0.0 && params.iqpnni_iterations > 1) {
		cout << "Start " << params.iqpnni_iterations-1 << " iterations of IQP and NNI" << endl;
		cout << "Number of representative leaves   : " << params.k_representative << endl;
		cout << "Probability of deleting sequences : " << params.p_delete << endl << endl;
		tree.setRepresentNum(params.k_representative);
		tree.setProbDelete(params.p_delete);
		tree.setIQPIterations(params.iqpnni_iterations-1);

		string tree_file_name = params.aln_file;
		tree_file_name += ".treefile";
		double bestscore = tree.doIQPNNI(tree_file_name);
		cout << "Best score found : " << bestscore << endl;
	} else {
		/* do SPR with likelihood function */
		if (params.tree_spr)
			tree.optimizeSPRBranches();
	}

	/* root the tree at the first sequence */
	tree.root = tree.findNodeName(alignment->getSeqName(0));
	assert(tree.root);

	t_end=clock();
	params.run_time = (t_end-t_begin);
	printf("Time used: %8.6f seconds.\n", (double)params.run_time / CLOCKS_PER_SEC);
	if (params.out_file)
		tree.printTree(params.out_file);
	else {
		//tree.printTree(cout);
		//cout << endl;
		if (verbose_mode > VB_MED) {
			tree.printTree(cout, WT_BR_LEN);
			cout << endl;
		}
	}
}

void runPhyloAnalysis(Params &params) {

	Alignment alignment(params.aln_file, params.intype);

	runPhyloAnalysis(params, &alignment);
}


