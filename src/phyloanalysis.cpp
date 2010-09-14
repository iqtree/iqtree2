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
#include "stoprule.h"

#include "mtreeset.h"
#include "mexttree.h"
#include "ratemeyerhaeseler.h"

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
	string fscore_name = params.out_prefix;
	string fmodel_name = params.out_prefix;
	
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

		ModelFactory *model_fac = new ModelFactory();

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

				// initialize model factory
				tree->setModelFactory(model_fac);
				model_fac->model = subst_model;
				model_fac->site_rate = rate_class[rate_type]; 


				// print some infos
				// clear all likelihood values
				tree->clearAllPartialLh();

				// optimize model parameters
				double cur_lh = tree->getModelFactory()->optimizeParameters();
				cout << "===> Testing ";
				cout.width(12);
				string str;
				str = subst_model->name;
				str += rate_class[rate_type]->name;
				cout << left << str << " (df=" << subst_model->getNDim()+rate_class[rate_type]->getNDim() << ")";
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

				tree->setModel(NULL);
				tree->setModelFactory(NULL);
				tree->setRate(NULL);

			}
		}

		delete model_fac;
		delete subst_model;
		for (rate_type = 3; rate_type >= 0; rate_type--)
			delete rate_class[rate_type];
		delete tree_hetero;
		delete tree_homo;


		fscore.close();

		/* now do the modeltest */
		if (nstates == 4)
		modeltest(model_arg, fscore_name.c_str(), fmodel_name.c_str(), LRT_model, IC_model);
	}

	if (nstates == 4) {
		cout << "Performing ModelTest 3.7 (Posada and Crandall, 1998) ..." << endl;
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

void reportPhyloAnalysis(Params &params, string &original_model, Alignment &alignment, IQPTree &tree) {
	int i, j;
	string outfile = params.out_prefix;

	outfile += ".iqtree";
	try {
		ofstream out;
		out.exceptions(ios::failbit | ios::badbit);
		out.open(outfile.c_str());
		out << "ANALYSIS RESULTS OF IQTREE " << VERSION << " built " << __DATE__<< endl << endl;
		out << "Random seed number (for debugging purpose): " << params.ran_seed << endl << endl;
		out << "REFERENCES" << endl << "----------" << endl << endl <<
			   "A manuscript describing IQTREE is currently under preparation." << endl << endl <<
			   "Please always cite: " << endl << endl <<
			   "Le Sy Vinh and Arndt von Haeseler (2004) IQPNNI: moving fast through tree space" << endl <<
			   "and stopping in time. Mol. Biol. Evol., 21(8):1565-1571." << endl << endl <<
			   "If you use the parallel version, please cite: " << endl << endl <<
			   "Bui Quang Minh, Le Sy Vinh, Arndt von Haeseler, and Heiko A. Schmidt (2005)" << endl << 
			   "pIQPNNI - parallel reconstruction of large maximum likelihood phylogenies." << endl <<
			   "Bioinformatics, 21:3794-3796." << endl << endl;

		if (original_model == "TEST" || original_model == "TESTONLY")
		out << "Since you used Modeltest please also cite Posada and Crandall (1998)" << endl
			<< "(see CREDITS section for the citation)" << endl << endl;

		out << "INPUT ALIGNMENT" << endl << "---------------" << endl << endl <<
			   "Alignment file name: " << params.aln_file << endl <<
			   "Type of sequences: " << ((alignment.num_states == 2) ? "Binary" : ((alignment.num_states == 4) ? "DNA" : "Protein")) << endl <<
			   "Number of sequences: " << alignment.getNSeq() << endl <<
			   "Number of sites: " << alignment.getNSite() << endl <<
			   "Number of constant sites: " << round(alignment.frac_const_sites * alignment.getNSite()) <<
			   " (= " << alignment.frac_const_sites*100 << "% of all sites)" << endl << 
			   "Number of patterns: " << alignment.size() << endl << endl;

		if (params.num_bootstrap_samples) 
			out << "BOOTSTRAP ANALYSIS" << endl << "------------------" << endl << endl 
				<< "Type of bootstrap: non-parametric" << endl
				<< "Number of replicates: " << params.num_bootstrap_samples << endl
				<< endl;

		out << "SUBSTITUTION PROCESS" << endl << "--------------------"  << endl << endl <<
			   "Model of substitution: " << tree.getModel()->name << endl << endl;
		out << "Rate matrix R:" << endl << endl;

		double *rate_mat = new double[alignment.num_states*(alignment.num_states-1)/2];
		tree.getModel()->getRateMatrix(rate_mat);
		int k;
		if (alignment.num_states > 4) out << fixed;
		for (i = 0, k = 0; i < alignment.num_states-1; i++)
			for (j = i+1; j < alignment.num_states; j++, k++) {
				out << "  " << alignment.convertStateBack(i) << "-" << alignment.convertStateBack(j) << ": " << rate_mat[k];
				if (alignment.num_states <= 4) out << endl; else
					if (k % 5 == 4) out << endl;
			}
		if (alignment.num_states > 4) out << endl;
		out.unsetf(ios_base::fixed);
		delete [] rate_mat;

		out << endl << "State frequencies: ";
		switch (tree.getModel()->getFreqType()) {
		case FREQ_EMPIRICAL:
			out << "(empirical counts from alignment)" << endl;
			break;
		case FREQ_ESTIMATE:
			out << "(estimated with maximum likelihood)" << endl;
			break;
		case FREQ_USER_DEFINED:
			out << "(user-defined)" << endl;
			break;
		case FREQ_EQUAL:
			out << "(equal frequencies)" << endl;
			break;
		default:
			break;
		}
		out << endl;

		double *state_freqs = new double[alignment.num_states];
		tree.getModel()->getStateFrequency(state_freqs);
		for (i = 0; i < alignment.num_states; i++)
			out << "  pi(" << alignment.convertStateBack(i) << ") = " << state_freqs[i] << endl;
		delete [] state_freqs;

		RateHeterogeneity *rate_model = tree.getRate();
		out << endl << "RATE HETEROGENEITY" << endl << "------------------" << endl << endl;
		out << "Model of rate heterogeneity: " << rate_model->full_name << endl;
		rate_model->writeInfo(out);
		
		if (rate_model->getNRate() > 1 || rate_model->getPInvar() > 0.0) {
			out << endl << "Rates and their respective probabilities used in the likelihood function:" << endl << endl
				<< " Category  Relative rate  Probability" << endl;
			if (rate_model->getPInvar() > 0.0)
				out << "  0         0              " << rate_model->getPInvar() << endl;
			for (i = 0; i < rate_model->getNRate(); i++) {
				out << "  " << i+1 << "         ";
				out.width(14);
				out << left << rate_model->getRate(i) << " " << (1.0 - rate_model->getPInvar()) / rate_model->getNRate() << endl;
			}
		}
// Bootstrap analysis: 
//Display as outgroup: a

		out << endl << "TREE SEARCH" << endl << "-----------" << endl << endl <<
			   "Stopping rule: " << ((params.stop_condition == SC_STOP_PREDICT) ? "Yes" : "No") << endl <<
			   "Number of iterations: " << tree.stop_rule.getNumIterations() << endl <<
			   "Probability of deleting sequences: " << params.p_delete << endl <<
			   "Number of representative leaves: " << params.k_representative << endl << endl;

		out << "MAXIMUM LIKELIHOOD TREE" << endl << "-----------------------" << endl << endl;
		
		out << "NOTE: Tree is UNROOTED although outgroup taxon '" << 
			((params.root) ? params.root : tree.aln->getSeqName(0)) <<"' is drawn at root" << endl;
		if (params.aLRT_replicates > 0 || (params.num_bootstrap_samples && params.compute_ml_tree)) { 
			out << "Numbers in parentheses are ";
			if (params.aLRT_replicates > 0) out << "SH-like aLRT supports";
			if (params.num_bootstrap_samples && params.compute_ml_tree) out << " / bootstrap supports";
			out << " (%)" << endl;
		}
		out << endl;
		tree.setRootNode(params.root);
		tree.sortTaxa();
		tree.drawTree(out);

		out << "Log-likehood of the tree: " << fixed << tree.getBestScore() << endl <<
			   "Unconstrained log-likelihood (without tree): " << alignment.computeUnconstrainedLogL() << endl << endl <<
			   "Tree in newick format:" << endl << endl;
		tree.printResultTree(params, out);
		out << endl;
/*
		if (params.write_intermediate_trees) {
			out << endl << "CONSENSUS OF INTERMEDIATE TREES" << endl << "-----------------------" << endl << endl 
				<< "Number of intermediate trees: " << tree.stop_rule.getNumIterations() << endl
				<< "Split threshold: " << params.split_threshold << endl
				<< "Burn-in: " << params.tree_burnin << endl << endl;
		}*/
		time_t cur_time;
		time (&cur_time);
		
		char *date_str;
		date_str = ctime(&cur_time);
		out.unsetf(ios_base::fixed);
		out << endl << "TIME STAMP" << endl << "----------" << endl << endl << "Date and time: " << date_str << 
			   "Running time: " << (double)params.run_time / CLOCKS_PER_SEC << " seconds" << endl << endl;

		out << "CREDITS" << endl << "-------" << endl << endl 
			<< "Some parts of the code were taken from TREE-PUZZLE package:" << endl << endl 
			<< "Heiko A. Schmidt, Korbinian Strimmer, Martin Vingron, and Arndt von Haeseler" << endl 
			<< "(2002) TREE-PUZZLE: maximum likelihood phylogenetic analysis using quartets" << endl
			<< "and parallel computing. Bioinformatics, 18(3):502-504." << endl << endl 

			<< "The source codes to construct the BIONJ tree were taken from BIONJ software:" << endl << endl 
			<< "Oliver Gascuel (1997) BIONJ: an improved version of the NJ algorithm" << endl 
			<< "nased on a simple model of sequence data. Mol. Bio. Evol., 14:685-695." << endl << endl 

			<< "The Nexus file parser was taken from the Nexus Class Library:" << endl << endl 
			<< "Paul O. Lewis (2003) NCL: a C++ class library for interpreting data files in" << endl 
			<< "NEXUS format. Bioinformatics, 19(17):2330-2331." << endl << endl 

			<< "The Modeltest 3.7 source codes were taken from:" << endl << endl 
			<< "David Posada and Keith A. Crandall (1998) MODELTEST: testing the model of" << endl
			<< "DNA substitution. Bioinformatics, 14(9):817-8." << endl << endl;

		out.close();

	} catch (ios::failure) {
		outError(ERR_WRITE_OUTPUT, outfile);
	}

	cout << "Analysis results written to: " << endl 
		<<  "  IQ-TREE report:           " << params.out_prefix << ".iqtree" << endl
		<<  "  Maximum-likelihood tree:  " << params.out_prefix << ".treefile" << endl;
	if (!params.user_file)
		cout << "  BIONJ tree:               " << params.out_prefix << ".bionj" << endl;
	if (!params.dist_file) {
		cout <<  "  Juke-Cantor distances:    " << params.out_prefix << ".jcdist" << endl;
		if (params.compute_ml_dist)
		cout <<  "  Likelihood distances:     " << params.out_prefix << ".mldist" << endl;
	}
	if (params.mvh_site_rate)
		cout << "  Site-rates by MvH model:  " << params.out_prefix << ".mvhrate" << endl;
	
	if (params.num_bootstrap_samples) 
		cout << endl << "Non-parametric bootstrap results written to:" << endl 
			 << "  Bootstrap alignments:     " << params.out_prefix << ".bootaln" << endl
			 << "  Bootstrap trees:          " << params.out_prefix << ".boottrees" << endl
			 << "  Consensus tree:           " << params.out_prefix << ".contree" << endl;
	cout << endl;

}


void checkZeroDist(Alignment *aln, double *dist) {
	int ntaxa = aln->getNSeq();
	IntVector checked;
	checked.resize(ntaxa, 0);
	int i, j;
	for (i = 0; i < ntaxa-1; i++) {
		if (checked[i]) continue;
		string str = ""; 
		bool first = true;
		for (j = i+1; j < ntaxa; j++)
			if (dist[i*ntaxa+j] <= 1e-6) {
				if (first) str = "ZERO distance between sequences " + aln->getSeqName(i);
				str += ", " + aln->getSeqName(j);
				checked[j] = 1;
				first = false;
			}
		checked[i] = 1;
		if (str != "") outWarning(str);
	}
}

void runPhyloAnalysis(Params &params, string &original_model, Alignment *alignment, IQPTree &tree) {

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
	if (!tree.dist_matrix) {
		tree.dist_matrix = new double[alignment->getNSeq() * alignment->getNSeq()];
		memset(tree.dist_matrix, 0, sizeof(double) * alignment->getNSeq() * alignment->getNSeq());
	}
	if (params.dist_file) {
		cout << "Reading distance matrix file " << params.dist_file << " ..." << endl;
	} else {
		cout << "Computing Juke-Cantor distances..." << endl;
	}
	string dist_file;
	tree.computeDist(params, alignment, tree.dist_matrix, dist_file);

	checkZeroDist(alignment, tree.dist_matrix);

	if (params.user_file) {
		// start the search with user-defined tree
		bool myrooted = params.is_rooted;
		tree.readTree(params.user_file, myrooted);
		tree.setAlignment(alignment);
	} else {
		tree.computeBioNJ(params, alignment, dist_file); // create BioNJ tree
	}

	if (params.root) {
		string str = params.root;
		if (!tree.findNodeName(str)) {
			str = "Specified root name " + str + "not found";
			outError(str);
		}
	}

	/* Fix if negative branch lengths detected */
	double fixed_length = 0.01;
	int fixed_number = tree.fixNegativeBranch(fixed_length);
	if (fixed_number) {
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
	assert(tree.aln);
	tree.optimize_by_newton = params.optimize_by_newton;
	if (!tree.getModelFactory()) tree.setModelFactory(new ModelFactory(params, &tree));
	tree.setModel(tree.getModelFactory()->model);
	tree.setRate(tree.getModelFactory()->site_rate);
	int model_df = tree.getModel()->getNDim() + tree.getRate()->getNDim();
	cout << "Model of evolution: " << tree.getModelName() << " (" << model_df << " free parameters)" << endl;
	cout << "Fixed branch lengths: " << ((params.fixed_branch_length) ? "Yes" : "No") << endl;
	cout << "Random seed: " << params.ran_seed << endl;
	cout << "Lambda used in NNI: " << cmdLambda << endl;

/*
	if (params.parsimony) {
		int score = tree.computeParsimonyScore();
		cout << "User tree has parsimony score of " << score << endl;

		// NNI with parsimony function 
		tree.searchNNI();
	}*/
	/* optimize model parameters */
	cout.precision(10);


	//cout << "User tree has likelihood score of " << tree.computeLikelihood() << endl;

	if (params.parsimony) {
            tree.enable_parsimony = true;
            tree.pars_scores = new double[3000];
            tree.lh_scores = new double[3000];

            for (int i=0; i<3000; i++) {
                tree.pars_scores[i] = 0;
                tree.lh_scores[i] = 0;
            }
        }
	//cout << "Parsimony score: " << tree.computeParsimonyScore() << endl;
	tree.cur_pars_score = tree.computeParsimony();
	//cout << "Fast parsimony score: " << tree.cur_pars_score << endl;

	cout << "Optimizing model parameters" << endl;
	double bestTreeScore = tree.getModelFactory()->optimizeParameters(params.fixed_branch_length);
	cout << "Log-likelihood of the current tree: " << bestTreeScore << endl;
	//Update tree score
	tree.curScore = bestTreeScore;
	
	if ((tree.getModel()->name == "JC") && tree.getRate()->getNDim() == 0)
		params.compute_ml_dist = false;
	if (!params.dist_file && params.compute_ml_dist) {
		stringstream best_tree_string;
		tree.printTree(best_tree_string, WT_BR_LEN + WT_TAXON_ID);
		cout << "Computing ML distances based on estimated model parameters..." << endl;
		tree.computeDist(params, alignment, tree.dist_matrix, dist_file); 
		if (!params.user_file) {
			tree.computeBioNJ(params, alignment, dist_file); // create BioNJ tree
			tree.fixNegativeBranch(fixed_length);
			if (!params.fixed_branch_length)
				tree.curScore = tree.optimizeAllBranches();
			else
				tree.curScore = tree.computeLikelihood();
			cout << "Log-likelihood of the new BIONJ tree: " << tree.curScore << endl;
			if (tree.curScore < bestTreeScore - 1e-6) {
				cout << "The new tree is worse, rolling back the first BIONJ tree..." << endl;
				tree.rollBack(best_tree_string);
				tree.curScore = tree.computeLikelihood();
				cout << "Backup log-likelihood: " << tree.curScore << endl;
			}
		}
	}

        
//        tree.drawTree(cout, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE);
//        string alnfile = params.aln_file;
//        string tree1 = alnfile + ".treefileBEFORE";
//        string tree2 = alnfile + ".treefileAFTER";
//        tree.printTree(tree1.c_str(), WT_TAXON_ID | WT_SORT_TAXA | WT_NEWLINE);

//        tree.perturb();

//        tree.printTree(tree2.c_str(), WT_TAXON_ID | WT_SORT_TAXA | WT_NEWLINE);
//        tree.drawTree(cout, WT_BR_SCALE | WT_INT_NODE | WT_TAXON_ID | WT_NEWLINE);
        
        
	/* Optimize branch lengths with likelihood function */
	//cout << "Optimizing branch lengths..." << endl;
	//cout << "Log-likelihood: " << tree.optimizeAllBranches() << endl;

	/* do NNI with likelihood function */

        
	if (params.min_iterations > 0) {
		cout << "Performing Nearest Neighbor Interchange..." << endl;
		//cout << "Current tree likelihood: " << tree.optimizeNNIBranches() << endl;
		//tree.optimizeAllBranches();

		clock_t nniBeginClock, nniEndClock;
		nniBeginClock = clock();
		//double newScore = tree.optimizeNNI(true);
                double newScore = tree.optimizeNNI(true);
		nniEndClock = clock();
		printf("Time used for first NNI search: %8.6f seconds.\n", (double)(nniEndClock - nniBeginClock) / CLOCKS_PER_SEC);

		//cout << "Tree likelihood after NNI: " << newScore << endl;
		if (newScore > bestTreeScore) {
			bestTreeScore = newScore;
			cout << "Found new best tree log-likelihood : " << bestTreeScore << endl;
		}
		else {
			cout << "Tree didn't improve after NNI :( " << endl;
		}
	}
	/*
	double sum_scaling = 1.0;
	if (!tree.checkEqualScalingFactor(sum_scaling))
		cout << "Scaling factor not equal along the tree" << endl;
	*/
	NodeVector pruned_taxa;
	StrVector linked_name;
	double *saved_dist_mat = tree.dist_matrix;
	double *pattern_lh;
	int num_low_support;
	clock_t mytime;
	if (params.min_iterations > 1 && params.aLRT_threshold <= 100 && params.aLRT_replicates > 0) {
		mytime = clock();
		cout <<"Testing tree branches by SH-like aLRT with " << params.aLRT_replicates << " replicates..." << endl;
		pattern_lh = new double[tree.aln->getNPattern()];
		tree.setRootNode(params.root);
		double score = tree.computeLikelihood(pattern_lh);
		num_low_support = tree.testAllBranches(params.aLRT_threshold, score, pattern_lh, params.aLRT_replicates);
	    tree.printResultTree(params);
		cout << "  " << (((double)clock()) - mytime) / CLOCKS_PER_SEC << " sec." << endl;
		cout << num_low_support << " branches show low support values (<= " << params.aLRT_threshold << "%)" << endl;
		delete [] pattern_lh;
		
		//tree.drawTree(cout);
		cout << "Collapsing stable clades..." << endl;
		tree.collapseStableClade(params.aLRT_threshold, pruned_taxa, linked_name, tree.dist_matrix);
		cout << pruned_taxa.size() << " taxa were pruned from stable clades" << endl;
	}

	if (!pruned_taxa.empty()) {
		cout << "Pruned alignment contains " << tree.aln->getNSeq() << " sequences and " << 
			tree.aln->getNSite() << " sites and " << tree.aln->getNPattern() << " patterns" << endl;
		//tree.clearAllPartialLh();
		tree.initializeAllPartialLh();
		tree.clearAllPartialLh();
		tree.curScore = tree.optimizeAllBranches();
		//cout << "Log-likelihood	after reoptimizing model parameters: " << tree.curScore << endl;
		tree.curScore = tree.optimizeNNI();
		cout << "Log-likelihood after optimizing partial tree: " << tree.curScore << endl;
		/*
		pattern_lh = new double[tree.aln->getNPattern()];
		double score = tree.computeLikelihood(pattern_lh);
		num_low_support = tree.testAllBranches(params.aLRT_threshold, score, pattern_lh, params.aLRT_replicates);
		tree.drawTree(cout);
		delete [] pattern_lh;*/
	}
	// set p delete if ZERO
/*
	if (params.p_delete == 0) {
		int num_high_support = tree.leafNum - 3 - num_low_support; 
		params.p_delete = (2.0 + num_high_support)*2.0 / tree.leafNum;
		if (params.p_delete > 0.5) params.p_delete = 0.5;
	}*/

	tree.setRepresentNum(params.k_representative);
	tree.setProbDelete(params.p_delete);
	tree.setIQPIterations(params.stop_condition, params.stop_confidence, params.min_iterations, params.max_iterations);
	tree.setIQPAssessQuartet(params.iqp_assess_quartet);

        /* do iterated local search */
        if ( ils && params.min_iterations > 1 ) {
            cout << "Star doing iterated local search " << endl;
            tree.doILS( params , perLevel );
            cout << "Optimizing model parameters" << endl;
            double endScore = tree.getModelFactory()->optimizeParameters(params.fixed_branch_length);
            cout << "Best score found : " << endScore << endl;
        }
	/* do the IQPNNI */
        else if (params.k_representative > 0 && params.min_iterations > 1) {
		cout << endl << "START IQPNNI SEARCH WITH THE FOLLOWING PARAMETERS" << endl;
		cout << "Number of representative leaves   : " << params.k_representative << endl;
		cout << "Probability of deleting sequences : " << tree.getProbDelete() << endl;
		cout << "Number of iterations              : ";
		if (params.stop_condition == SC_FIXED_ITERATION) 
			cout << params.min_iterations << endl;
		else cout << "auto-predicted in range [" << params.min_iterations <<"," << 
			params.max_iterations <<"], confidence " << params.stop_confidence << endl;
		cout << "Important quartet assessed on     : " << ((params.iqp_assess_quartet == IQP_DISTANCE) ? "Distance" : "Parsimony") << endl;
		cout << endl;

		tree.doIQPNNI(params);
	} else {
		/* do SPR with likelihood function */
		if (params.tree_spr)
			tree.optimizeSPRBranches();
	}


	if (!pruned_taxa.empty()) {
		tree.disableHeuristic();
		cout << "Restoring full tree..." << endl;
		tree.restoreStableClade(alignment, pruned_taxa, linked_name);
		delete [] tree.dist_matrix;
		tree.dist_matrix = saved_dist_mat;
		tree.initializeAllPartialLh();
		tree.clearAllPartialLh();
		tree.curScore = tree.optimizeAllBranches();
		//cout << "Log-likelihood	after reoptimizing model parameters: " << tree.curScore << endl;
		tree.curScore = tree.optimizeNNI();
		cout << "Log-likelihood	after reoptimizing full tree: " << tree.curScore << endl;
	}

	cout << "Optimizing model parameters" << endl;
	tree.setBestScore(tree.getModelFactory()->optimizeParameters(params.fixed_branch_length));
	cout << "Best score found : " << tree.getBestScore() << endl;

	/* root the tree at the first sequence */
	tree.root = tree.findNodeName(alignment->getSeqName(0));
	assert(tree.root);

	if (params.mvh_site_rate) {
		RateMeyerHaeseler rate_mvh;
		cout << "Computing site-specific rates by " << rate_mvh.full_name << "..." << endl;
		rate_mvh.setTree(&tree);
		rate_mvh.optimizeParameters();
		string rate_file = params.out_prefix;
		rate_file += ".mvhrate";
		rate_mvh.writeSiteRates(rate_file.c_str());
	}

	if (params.aLRT_replicates > 0) {
		mytime = clock();
		cout <<"Testing tree branches by SH-like aLRT with " << params.aLRT_replicates << " replicates..." << endl;
		pattern_lh = new double[tree.aln->getNPattern()];
		tree.setRootNode(params.root);
		double score = tree.computeLikelihood(pattern_lh);
		num_low_support = tree.testAllBranches(params.aLRT_threshold, score, pattern_lh, params.aLRT_replicates);
		cout << "  " << (((double)clock()) - mytime) / CLOCKS_PER_SEC << " sec." << endl;
		cout << num_low_support << " branches show low support values (<= " << params.aLRT_threshold << "%)" << endl;
		delete [] pattern_lh;
	}
	t_end=clock();
	params.run_time = (t_end-t_begin);
	printf("Time used: %8.6f seconds.\n", (double)params.run_time / CLOCKS_PER_SEC);

    tree.printResultTree(params);

	if (params.out_file)
		tree.printTree(params.out_file);
	else {
		//tree.printTree(cout);
		//cout << endl;
		/*
		if (verbose_mode > VB_MED) {
			if (verbose_mode >= VB_DEBUG)
				tree.drawTree(cout, WT_BR_SCALE + WT_INT_NODE + WT_BR_LEN);
			else
				tree.drawTree(cout);
		}*/
	}

}

void runPhyloAnalysis(Params &params) {

	Alignment alignment(params.aln_file, params.sequence_type, params.intype);
	IQPTree tree;
	string original_model = params.model_name;
	if (params.aln_output)  {
		alignment.printPhylip(params.aln_output);
	} else if (params.num_bootstrap_samples == 0) {
		runPhyloAnalysis(params, original_model, &alignment, tree);
		reportPhyloAnalysis(params, original_model, alignment, tree);
	} else {
		// turn off aLRT test
		int saved_aLRT_replicates = params.aLRT_replicates;
		params.aLRT_replicates = 0;
		string treefile_name = params.out_prefix;
		treefile_name += ".treefile";
		string boottrees_name = params.out_prefix;
		boottrees_name += ".boottrees";
		string bootaln_name = params.out_prefix;
		bootaln_name += ".bootaln";
		// first empty the boottrees file
		try {
			// write the tree into .boottrees file
			ofstream tree_out;
			tree_out.exceptions(ios::failbit | ios::badbit);
			tree_out.open(boottrees_name.c_str());
			tree_out.close();
		} catch (ios::failure) {
			outError(ERR_WRITE_OUTPUT, boottrees_name);
		}
		
		// empty the bootaln file
		try {
			// write the tree into .boottrees file
			ofstream tree_out;
			tree_out.exceptions(ios::failbit | ios::badbit);
			tree_out.open(bootaln_name.c_str());
			tree_out.close();
		} catch (ios::failure) {
			outError(ERR_WRITE_OUTPUT, bootaln_name);
		}

		// do bootstrap analysis
		for (int sample = 0; sample < params.num_bootstrap_samples; sample++) {
			cout << endl << "===> START BOOTSTRAP REPLICATE NUMBER " << sample+1 << endl << endl;

			Alignment bootstrap_alignment;
			IQPTree boot_tree;
			cout << "Creating bootstrap alignment..." << endl;
			bootstrap_alignment.createBootstrapAlignment(&alignment);
			bootstrap_alignment.printPhylip(bootaln_name.c_str(), true);
			runPhyloAnalysis(params, original_model, &bootstrap_alignment, boot_tree);
			// read in the output tree file
			string tree_str;
			try {
				ifstream tree_in;
				tree_in.exceptions(ios::failbit | ios::badbit);
				tree_in.open(treefile_name.c_str());
				tree_in >> tree_str;
				tree_in.close();
			} catch (ios::failure) {
				outError(ERR_READ_INPUT, treefile_name);
			}
			// write the tree into .boottrees file
			try {
				ofstream tree_out;
				tree_out.exceptions(ios::failbit | ios::badbit);
				tree_out.open(boottrees_name.c_str(), ios_base::out | ios_base::app);
				tree_out << tree_str << endl;
				tree_out.close();
			} catch (ios::failure) {
				outError(ERR_WRITE_OUTPUT, boottrees_name);
			}
			if (params.num_bootstrap_samples == 1) 
				reportPhyloAnalysis(params, original_model, bootstrap_alignment, boot_tree);
		}

		cout << endl << "===> COMPUTE CONSENSUS TREE FROM " << 
			params.num_bootstrap_samples << " BOOTSTRAP TREES" << endl << endl;
		computeConsensusTree(boottrees_name.c_str(), 0, 
			params.split_threshold, NULL, params.out_prefix); 

		if (params.compute_ml_tree) {
			cout << endl << "===> START ANALYSIS ON THE ORIGINAL ALIGNMENT" << endl << endl;
			params.aLRT_replicates = saved_aLRT_replicates;
			runPhyloAnalysis(params, original_model, &alignment, tree);

			cout << endl << "===> ASSIGN BOOTSTRAP SUPPORTS TO THE TREE FROM ORIGINAL ALIGNMENT" << endl << endl;
			MExtTree ext_tree;
			assignBootstrapSupport(boottrees_name.c_str(), 0, 
				treefile_name.c_str(), false, treefile_name.c_str(), params.out_prefix, ext_tree); 
			tree.copyTree(&ext_tree);
			reportPhyloAnalysis(params, original_model, alignment, tree);
		}
	}	
}



void assignBootstrapSupport(const char *input_trees, int burnin, const char *target_tree, bool rooted, 
	const char *output_tree, const char *out_prefix, MExtTree &mytree) {
	//bool rooted = false;
	// read the tree file

	mytree.init(target_tree, rooted);
	// reindex the taxa in the tree to aphabetical names
	NodeVector taxa;
	mytree.getTaxa(taxa);
	sort(taxa.begin(), taxa.end(), nodenamecmp);
	int i = 0;
	for (NodeVector::iterator it = taxa.begin(); it != taxa.end(); it++) {
		(*it)->id = i++;
	}

	// read the bootstrap tree file
	MTreeSet boot_trees(input_trees, rooted, burnin);
	/*
	string filename = params.boot_trees;
	filename += ".nolen";
	boot_trees.printTrees(filename.c_str(), false);
	return;
	*/
	SplitGraph sg;
	SplitIntMap hash_ss;
	// make the taxa name
	vector<string> taxname;
	taxname.resize(mytree.leafNum);
	mytree.getTaxaName(taxname);

	boot_trees.convertSplits(taxname, sg, hash_ss, SW_COUNT);
	// compute the percentage of appearance
	sg.scaleWeight(100.0/boot_trees.size(), true);
//	printSplitSet(sg, hash_ss);
	//sg.report(cout);
	cout << "Creating bootstrap support values..." << endl;
	mytree.createBootstrapSupport(taxname, boot_trees, sg, hash_ss);
	//mytree.scaleLength(100.0/boot_trees.size(), true);
	string out_file;
	if (output_tree) 
		out_file = output_tree;
	else {
		if (out_prefix)
			out_file = out_prefix;
		else
			out_file = target_tree;
		out_file += ".suptree";
	}
	
	mytree.printTree(out_file.c_str());
	cout << "Tree with assigned bootstrap support written to " << out_file << endl;
}

void computeConsensusTree(const char *input_trees, int burnin, double cutoff, 
	const char *output_tree, const char *out_prefix) {
	bool rooted = false;

	// read the bootstrap tree file
	MTreeSet boot_trees(input_trees, rooted, burnin);
	string first_taxname = boot_trees.front()->root->name;
	//if (params.root) first_taxname = params.root;	

	SplitGraph sg;

	boot_trees.convertSplits(sg, cutoff, SW_COUNT);
	//sg.report(cout);

	cout << "Creating greedy consensus tree..." << endl;
	MTree mytree;
	SplitGraph maxsg;
	sg.findMaxCompatibleSplits(maxsg);
	//maxsg.saveFile(cout);
	cout << "convert compatible split system into tree..." << endl;
	mytree.convertToTree(maxsg);
	Node *node = mytree.findNodeName(first_taxname);
	if (node) mytree.root = node;
	mytree.scaleLength(100.0/boot_trees.size(), true);

	mytree.getTaxaID(maxsg.getSplitsBlock()->getCycle());
	//maxsg.saveFile(cout);

	string out_file;

	if (output_tree) 
		out_file = output_tree;
	else {
		if (out_prefix)
			out_file = out_prefix;
		else
			out_file = input_trees;
		out_file += ".contree";
	}

	
	mytree.printTree(out_file.c_str(), WT_BR_CLADE);
	cout << "Consensus tree written to " << out_file << endl;
}

void computeConsensusNetwork(const char *input_trees, int burnin, double cutoff, 
	const char *output_tree, const char *out_prefix) {
	bool rooted = false;

	// read the bootstrap tree file
	MTreeSet boot_trees(input_trees, rooted, burnin);

	SplitGraph sg;
	//SplitIntMap hash_ss;

	boot_trees.convertSplits(sg, cutoff, SW_SUM);

	string out_file;

	if (output_tree) 
		out_file = output_tree;
	else {
		if (out_prefix)
			out_file = out_prefix;
		else
			out_file = input_trees;
		out_file += ".nex";
	}


	try {
		ofstream out;
		out.exceptions(ios::failbit | ios::badbit);
		out.open(out_file.c_str());
		sg.saveFile(out);
		out.close();
		cout << "Consensus network printed to " << out_file << endl;
	} catch (ios::failure) {
		outError(ERR_WRITE_OUTPUT, out_file);
	}

}
