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
#include <iqtree_config.h>
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
//#include "modeltest_wrapper.h"
#include "modelprotein.h"
#include "stoprule.h"

#include "mtreeset.h"
#include "mexttree.h"
#include "ratemeyerhaeseler.h"
#include "whtest_wrapper.h"
#include "partitionmodel.h"
#include "guidedbootstrap.h"
#include "modelset.h"

//const int DNA_MODEL_NUM = 14;
clock_t t_begin, t_end;

const int BIN_MODEL_NUM = 2;
string bin_model_names[BIN_MODEL_NUM] = {"JC2","GTR2"};

const int DNA_MODEL_NUM = 22;
string dna_model_names[DNA_MODEL_NUM] ={"JC", "F81", "K80", "HKY", "TNe", "TN", "K81", "K81u", "TPM2", "TPM2u",
    "TPM3", "TPM3u", "TIMe", "TIM", "TIM2e", "TIM2", "TIM3e", "TIM3", "TVMe", "TVM", "SYM", "GTR"};
/*string dna_model_names[DNA_MODEL_NUM] ={"JC", "F81", "K80", "HKY", "TNe", "TN", "K81", "K81u", 
	"TIMe", "TIM", "TVMe", "TVM", "SYM", "GTR"};*/

const int AA_MODEL_NUM = 18;
string aa_model_names[AA_MODEL_NUM] ={"Dayhoff", "mtMAM", "JTT", "WAG", "cpREV", "mtREV", "rtREV",
    "mtART", "mtZOA", "VT", "LG", "DCMut", "PMB", "HIVb", "HIVw", "JTTDCMut", "FLU", "Blosum62"};

/**
 * check if the model file contains correct information
 * @param model_file model file names
 * @param model_name (OUT) vector of model names
 * @param lh_scores (OUT) vector of tree log-likelihoods
 * @param df_vec (OUT) vector of degrees of freedom (or K)
 * @return TRUE if success, FALSE failed.
 */
bool checkModelFile(string model_file, StrVector &model_names, DoubleVector &lh_scores, IntVector &df_vec) {
	if (!fileExists(model_file)) return false;
    cout << model_file << " exists, checking this file" << endl;
	ifstream in;
	try {
		in.exceptions(ios::failbit | ios::badbit);
		in.open(model_file.c_str());
		in.exceptions(ios::badbit);
		string str;
		in >> str;
		if (str != "Model") throw false;
		in >> str;
		if (str != "df") throw false;
		in >> str;
		if (str != "LnL") throw false;
		while (!in.eof()) {
			in >> str;
			if (in.eof()) break;
			model_names.push_back(str);
			int df;
			double logl;
			in >> df >> logl;
			df_vec.push_back(df);
			lh_scores.push_back(logl);
			//cout << str << " " << df << " " << logl << endl;
		}
		in.clear();
		// set the failbit again
		in.exceptions(ios::failbit | ios::badbit);
		in.close();
	} catch (bool ret) {
		in.close();
		return ret;
	} catch (ios::failure) {
		outError("Cannot read file ", model_file);
	}
	return true;
}
	
/**
        testing the best-fit model
        return in params.freq_type and params.rate_type
 */
string modelTest(Params &params, PhyloTree *in_tree) {
    int nstates = in_tree->aln->num_states;
    if (nstates != 2 && nstates != 4 && nstates != 20)
        outError("Test of best-fit models only works for DNA or Protein");
    string fmodel_str = params.out_prefix;
    fmodel_str += ".model";
   
    int num_models = (nstates == 2) ? BIN_MODEL_NUM : ( (nstates == 4) ? DNA_MODEL_NUM : AA_MODEL_NUM);
    int model, rate_type;

    string best_model;
	StrVector model_names;
	DoubleVector lh_scores;
	IntVector df_vec;
	/* first check the model file */
	bool ok_model_file = checkModelFile(fmodel_str, model_names, lh_scores, df_vec);
	ok_model_file &= (model_names.size() == num_models*4);
	ofstream fmodel;
	if (!ok_model_file) {
		model_names.clear();
		lh_scores.clear();
		df_vec.clear();
		fmodel.open(fmodel_str.c_str());
		if (!fmodel.is_open())
			outError("cannot write to ", fmodel_str);
		fmodel << "Model\tdf\tLnL" << endl;
		fmodel.precision(4);
		fmodel << fixed;
	} else {
		cout << fmodel_str << " seems to be a correct model file" << endl;
	} 

	PhyloTree *tree_homo = new PhyloTree();
	tree_homo->optimize_by_newton = params.optimize_by_newton;
	tree_homo->sse = params.SSE;
	tree_homo->copyPhyloTree(in_tree);

	PhyloTree *tree_hetero = new PhyloTree();
	tree_hetero->optimize_by_newton = params.optimize_by_newton;
	tree_hetero->sse = params.SSE;
	tree_hetero->copyPhyloTree(in_tree);

	RateHeterogeneity * rate_class[4];
	rate_class[0] = new RateHeterogeneity();
	rate_class[1] = new RateInvar(-1, NULL);
	rate_class[2] = new RateGamma(params.num_rate_cats, -1, params.gamma_median, NULL);
	rate_class[3] = new RateGammaInvar(params.num_rate_cats, -1, params.gamma_median, -1, NULL);
	GTRModel *subst_model;
	if (nstates == 4)
		subst_model = new ModelDNA("JC", FREQ_UNKNOWN, in_tree);
	else
		subst_model = new ModelProtein("WAG", FREQ_UNKNOWN, in_tree);

	ModelFactory *model_fac = new ModelFactory();

	int ssize = in_tree->aln->getNSite(); // sample size
	if (params.model_test_sample_size) ssize = params.model_test_sample_size;
	cout << "Testing " << num_models * 4 << ((nstates == 4) ? " DNA" : " protein") << " models (sample size: " << ssize << ") ..." << endl;
	cout << "Model         -LnL         df AIC          AICc         BIC" << endl;
	DoubleVector AIC_scores;
	DoubleVector AICc_scores;
	DoubleVector BIC_scores;
	
	for (model = 0; model < num_models; model++) {
		for (rate_type = 0; rate_type <= 3; rate_type += 1) {
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
			subst_model->init((nstates == 2) ? bin_model_names[model].c_str() : ((nstates == 4) ? dna_model_names[model].c_str() : aa_model_names[model].c_str()), FREQ_UNKNOWN);
			subst_model->setTree(tree);
			tree->setModel(subst_model);
			// initialize rate
			tree->setRate(rate_class[rate_type]);
			rate_class[rate_type]->setTree(tree);

			// initialize model factory
			tree->setModelFactory(model_fac);
			model_fac->model = subst_model;
			model_fac->site_rate = rate_class[rate_type];

			string str;
			str = subst_model->name;
			str += rate_class[rate_type]->name;

			// print some infos
			// clear all likelihood values
			tree->clearAllPartialLH();

			// optimize model parameters
			int df = subst_model->getNDim() + rate_class[rate_type]->getNDim();
			double cur_lh;
			if (!ok_model_file) {
				cur_lh = tree->getModelFactory()->optimizeParameters(false, false);
				model_names.push_back(str);
				fmodel << str << "\t" << df << "\t" << cur_lh << endl;
			} else {
				// sanity check
				if (str != model_names[model*4+rate_type] || df != df_vec[model*4+rate_type])
					outError("Incorrect model file, please delete it and rerun again: ", fmodel_str); 
				cur_lh = lh_scores[model*4 + rate_type];
			}
			double AIC_score = -2*cur_lh + 2 * df;
			double AICc_score = AIC_score + 2.0*df*(df+1)/(ssize -df-1);
			double BIC_score = -2*cur_lh + df * log(ssize);
			AIC_scores.push_back(AIC_score);
			AICc_scores.push_back(AICc_score);
			BIC_scores.push_back(BIC_score);
			cout.width(13);
			cout << left << str << " ";
			cout.precision(3);
			cout.width(12);
			cout << fixed;
			cout << -cur_lh << " ";
			cout.width(2);
			cout << df << " ";
			cout.width(12);
			cout << AIC_score << " ";
			cout.width(12);
			cout << AICc_score << " " << BIC_score;
			cout << endl;
			tree->setModel(NULL);
			tree->setModelFactory(NULL);
			tree->setRate(NULL);

		}
	}
	cout.unsetf(ios::fixed);
	int model_aic = min_element(AIC_scores.begin(), AIC_scores.end()) - AIC_scores.begin();
	cout << "Akaike Information Criterion:           " << model_names[model_aic] << endl;
	int model_aicc = min_element(AICc_scores.begin(), AICc_scores.end()) - AICc_scores.begin();
	cout << "Corrected Akaike Information Criterion: " << model_names[model_aicc] << endl;
	int model_bic = min_element(BIC_scores.begin(), BIC_scores.end()) - BIC_scores.begin();
	cout << "Bayesian Information Criterion:         " << model_names[model_bic] << endl;
	switch (params.model_test_criterion) {
		case MTC_AIC: best_model = model_names[model_aic]; break;
		case MTC_AICC: best_model = model_names[model_aicc]; break;
		case MTC_BIC: best_model = model_names[model_bic]; break;
	}
	delete model_fac;
	delete subst_model;
	for (rate_type = 3; rate_type >= 0; rate_type--)
		delete rate_class[rate_type];
	delete tree_hetero;
	delete tree_homo;

	if (!ok_model_file) fmodel.close();
	cout << "Best-fit model: " << best_model << endl;
    return best_model;
}

void reportReferences(ofstream &out, string &original_model) {
	out << "Two manuscripts are currently under preparation:" << endl << endl <<
		"Bui Quang Minh, Minh Anh Thi Nguyen, and Arndt von Haeseler (2012) Ultra-fast" << endl <<
		"approximation for phylogenetic bootstrap. Submitted." << endl << endl <<
		"Lam-Tung Nguyen, Heiko A. Schmidt, Bui Quang Minh, and Arndt von Haeseler (2012)" << endl <<
		"IQ-TREE: Efficient algorithm for phylogenetic inference by maximum likelihood" << endl <<
		"and important quartet puzzling. In prep." << endl << endl <<
		"For the original IQPNNI algorithm please cite: " << endl << endl <<
		"Le Sy Vinh and Arndt von Haeseler (2004) IQPNNI: moving fast through tree space" << endl <<
		"and stopping in time. Mol. Biol. Evol., 21(8):1565-1571." << endl << endl;
/*		"*** If you use the parallel version, please cite: " << endl << endl <<
		"Bui Quang Minh, Le Sy Vinh, Arndt von Haeseler, and Heiko A. Schmidt (2005)" << endl <<
		"pIQPNNI - parallel reconstruction of large maximum likelihood phylogenies." << endl <<
		"Bioinformatics, 21:3794-3796." << endl << endl;*/

// 	if (original_model == "TEST" || original_model == "TESTONLY")
// 		out << "Since you used Modeltest please also cite Posada and Crandall (1998)" << endl << endl;
}

void reportAlignment(ofstream &out, Alignment &alignment) {
	out <<  "Input data: "  << alignment.getNSeq() << " sequences with " <<
			alignment.getNSite() << " " << ((alignment.num_states == 2) ? "binary" : ((alignment.num_states == 4) ? "nucleotide" : "amino-acid")) <<
			" sites" << endl <<
			"Number of constant sites: " << round(alignment.frac_const_sites * alignment.getNSite()) <<
			" (= " << alignment.frac_const_sites * 100 << "% of all sites)" << endl <<
			"Number of site patterns: " << alignment.size() << endl << endl;
}

void reportModel(ofstream &out, PhyloTree &tree) {
	int i, j;
	out << "Model of substitution: " << tree.getModel()->name << endl << endl;
	out << "Rate parameter R:" << endl << endl;

	double *rate_mat = new double[tree.aln->num_states * tree.aln->num_states];
	if (!tree.getModel()->isSiteSpecificModel())
		tree.getModel()->getRateMatrix(rate_mat);
	else
		((ModelSet*)tree.getModel())->front()->getRateMatrix(rate_mat);
	int k;
	if (tree.aln->num_states > 4) out << fixed;
	if (tree.getModel()->isReversible()) {
		for (i = 0, k = 0; i < tree.aln->num_states - 1; i++)
			for (j = i + 1; j < tree.aln->num_states; j++, k++) {
				out << "  " << tree.aln->convertStateBack(i) << "-" << tree.aln->convertStateBack(j) << ": " << rate_mat[k];
				if (tree.aln->num_states <= 4) out << endl;
				else
					if (k % 5 == 4) out << endl;
			}
	} else { // non-reversible model
		for (i = 0, k = 0; i < tree.aln->num_states; i++)
			for (j = 0; j < tree.aln->num_states; j++)
				if (i!=j) {
				out << "  " << tree.aln->convertStateBack(i) << "-" << tree.aln->convertStateBack(j) << ": " << rate_mat[k];
				if (tree.aln->num_states <= 4) out << endl;
				else
					if (k % 5 == 4) out << endl;
				k++;
			}

	}


	if (tree.aln->num_states > 4) out << endl;
	out.unsetf(ios_base::fixed);
	delete [] rate_mat;

	out << endl << "State frequencies: ";
	if (tree.getModel()->isSiteSpecificModel())
		out << "(site specific frequencies)" << endl << endl;
	else {
		if (!tree.getModel()->isReversible()) 
			out << "(inferred from Q matrix)" << endl; 
		else
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

		double *state_freqs = new double[tree.aln->num_states];
		tree.getModel()->getStateFrequency(state_freqs);
		for (i = 0; i < tree.aln->num_states; i++)
			out << "  pi(" << tree.aln->convertStateBack(i) << ") = " << state_freqs[i] << endl;
		delete [] state_freqs;
		out << endl;
		// report Q matrix
		double *q_mat = new double[tree.aln->num_states * tree.aln->num_states];
		tree.getModel()->getQMatrix(q_mat);

		out << "Rate matrix Q:" << endl << endl;
		for (i = 0, k = 0; i < tree.aln->num_states; i++) {
			out << "  " << tree.aln->convertStateBack(i);
			for (j = 0; j < tree.aln->num_states; j++, k++) {
				out << "  ";
				out.width(8);
				out << q_mat[k];
			}
			out << endl;
		}
		out << endl;
		delete [] q_mat;
	}
}

void reportRate(ofstream &out, PhyloTree &tree) {
	int i;
	RateHeterogeneity *rate_model = tree.getRate();
	out << "Model of rate heterogeneity: " << rate_model->full_name << endl;
	rate_model->writeInfo(out);

	if (rate_model->getNDiscreteRate() > 1 || rate_model->getPInvar() > 0.0) {
		out << endl << " Category  Relative_rate  Proportion" << endl;
		if (rate_model->getPInvar() > 0.0)
			out << "  0         0              " << rate_model->getPInvar() << endl;
		int cats = rate_model->getNDiscreteRate();
		DoubleVector prop;
		if (rate_model->getGammaShape() > 0)
			prop.resize(cats, (1.0 - rate_model->getPInvar()) / rate_model->getNRate());
		else {
			prop.resize(cats, 0.0);
			for (i = 0; i < tree.aln->getNPattern(); i++)
				prop[rate_model->getPtnCat(i)] += tree.aln->at(i).frequency;
			for (i = 0; i < cats; i++)
				prop[i] /= tree.aln->getNSite();
		}
		for (i = 0; i < cats; i++) {
			out << "  " << i + 1 << "         ";
			out.width(14);
			out << left << rate_model->getRate(i) << " " << prop[i];
			out << endl;
		}
	}
/*
	if (rate_model->getNDiscreteRate() > 1 || rate_model->isSiteSpecificRate())
		out << endl << "See file " << rate_file << " for site-specific rates and categories" << endl;*/
	out << endl;
}

void reportTree(ofstream &out, Params &params, PhyloTree &tree, double tree_lh, double lh_variance) {

	int zero_branches = tree.countZeroBranches();
	if (zero_branches > 0) {
		out << "WARNING: " << zero_branches << " branches of ZERO lengths and should be treated with caution!" << endl;
		cout << endl << "WARNING: " << zero_branches  << " branches of ZERO lengths and should be treated with caution!" << endl;
		out << "         Such branches are denoted by '***' in the figure" << endl << endl;
	}
	tree.sortTaxa();
	//tree.setExtendedFigChar();
	tree.drawTree(out);

	out << "Log-likehood of the tree: " << fixed << tree_lh << " (s.e. " << sqrt(lh_variance) << ")" << endl
		<< "Unconstrained log-likelihood (without tree): " << tree.aln->computeUnconstrainedLogL() << endl
		<< "Total tree length: " << tree.treeLength() << endl << endl
		<< "Tree in newick format:" << endl << endl;

    tree.printTree(out, WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA);

	out << endl << endl;
}

void reportCredits(ofstream &out) {
	out << "CREDITS" << endl << "-------" << endl << endl
		<< "Some parts of the code were taken from TREE-PUZZLE package:" << endl << endl
		<< "Heiko A. Schmidt, Korbinian Strimmer, Martin Vingron, and Arndt von Haeseler" << endl
		<< "(2002) TREE-PUZZLE: maximum likelihood phylogenetic analysis using quartets" << endl
		<< "and parallel computing. Bioinformatics, 18(3):502-504." << endl << endl

		<< "The source codes to construct the BIONJ tree were taken from BIONJ software:" << endl << endl
		<< "Oliver Gascuel (1997) BIONJ: an improved version of the NJ algorithm" << endl
		<< "based on a simple model of sequence data. Mol. Bio. Evol., 14:685-695." << endl << endl

		<< "The Nexus file parser was taken from the Nexus Class Library:" << endl << endl
		<< "Paul O. Lewis (2003) NCL: a C++ class library for interpreting data files in" << endl
		<< "NEXUS format. Bioinformatics, 19(17):2330-2331." << endl << endl

		<< "The Modeltest 3.7 source codes were taken from:" << endl << endl
		<< "David Posada and Keith A. Crandall (1998) MODELTEST: testing the model of" << endl
		<< "DNA substitution. Bioinformatics, 14(9):817-8." << endl << endl;
}

void reportPhyloAnalysis(Params &params, string &original_model, Alignment &alignment, IQPTree &tree) {
    string outfile = params.out_prefix;

    outfile += ".iqtree";
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(outfile.c_str());
        out << "IQ-TREE " << iqtree_VERSION_MAJOR << "." << iqtree_VERSION_MINOR << "." << iqtree_VERSION_PATCH << " built " << __DATE__ << endl << endl;
        if (params.partition_file)
               out << "Partition file name: " << params.partition_file << endl;
        if (params.aln_file)
               out << "Input file name: " << params.aln_file << endl;

		if (params.user_file)
			out << "User tree file name: " << params.user_file << endl;
		out << "Type of analysis: ";
		if (params.compute_ml_tree)
			out << "tree reconstruction";
		if (params.num_bootstrap_samples > 0) {
			if (params.compute_ml_tree) out << " + ";
			out << "non-parametric bootstrap (" << params.num_bootstrap_samples << " replicates)";
		}
		out << endl;
    	out << "Random seed number: " << params.ran_seed << endl << endl;
        out << "REFERENCES" << endl << "----------" << endl << endl;
        reportReferences(out, original_model);

        out << "SEQUENCE ALIGNMENT" << endl << "------------------" << endl << endl;
        if (tree.isSuperTree()) {
			out << "Input data: " << alignment.getNSeq() << " taxa with " <<
					alignment.getNSite() << " partitions" << endl << endl;
        	PhyloSuperTree *stree = (PhyloSuperTree*)&tree;
        	int part = 0;
        	for (PhyloSuperTree::iterator it = stree->begin(); it != stree->end(); it++, part++) {
        		out << "FOR PARTITION " << stree->part_info[part].name <<":" << endl << endl;
        		reportAlignment(out, *((*it)->aln));
        	}
        } else reportAlignment(out, alignment);


        out << "SUBSTITUTION PROCESS" << endl << "--------------------" << endl << endl;
        if (tree.isSuperTree()) {
        	out << "Full partition model with separate branch lengths and models between partitions" << endl << endl;
        	PhyloSuperTree *stree = (PhyloSuperTree*)&tree;
        	int part = 0;
        	for (PhyloSuperTree::iterator it = stree->begin(); it != stree->end(); it++, part++) {
        		out << "FOR PARTITION " << stree->part_info[part].name <<":" << endl << endl;
        		reportModel(out, *(*it));
        	}
        } else
        	reportModel(out, tree);

        out << "RATE HETEROGENEITY" << endl << "------------------" << endl << endl;
        if (tree.isSuperTree()) {
        	PhyloSuperTree *stree = (PhyloSuperTree*)&tree;
        	int part = 0;
        	for (PhyloSuperTree::iterator it = stree->begin(); it != stree->end(); it++, part++) {
        		out << "FOR PARTITION " << stree->part_info[part].name <<":" << endl << endl;
        		reportRate(out, *(*it));
        	}
        } else
        	reportRate(out, tree);

        // Bootstrap analysis:
        //Display as outgroup: a

		if (original_model == "WHTEST") {
			out << "TEST OF MODEL HOMOGENEITY" << endl << "-------------------------" << endl << endl;
			out << "Delta of input data:                 " << params.whtest_delta << endl;
			out << ".95 quantile of Delta distribution:  " << params.whtest_delta_quantile << endl;
			out << "Number of simulations performed:     " << params.whtest_simulations << endl;
			out << "P-value:                             " << params.whtest_p_value << endl;
			if (params.whtest_p_value < 0.05) {
				out << "RESULT: Homogeneity assumption is rejected (p-value cutoff 0.05)" << endl;
			} else {
				out << "RESULT: Homogeneity assumption is NOT rejected (p-value cutoff 0.05)" << endl;
			}
			out << endl << "*** For this result please cite:" << endl << endl;
			out << "G. Weiss and A. von Haeseler (2003) Testing substitution models" << endl
				<< "within a phylogenetic tree. Mol. Biol. Evol, 20(4):572-578" << endl << endl;
		}

        out << "TREE SEARCH" << endl << "-----------" << endl << endl <<
                "Stopping rule: " << ((params.stop_condition == SC_STOP_PREDICT) ? "Yes" : "No") << endl <<
                "Number of iterations: " << tree.stop_rule.getNumIterations() << endl <<
                "Probability of deleting sequences: " << params.p_delete << endl <<
                "Number of representative leaves: " << params.k_representative << endl <<
                "NNI log-likelihood cutoff: " << tree.getNNICutoff() << endl << endl;

		if (params.compute_ml_tree) {
			out << "MAXIMUM LIKELIHOOD TREE" << endl << "-----------------------" << endl << endl;

			tree.setRootNode(params.root);
			out << "NOTE: Tree is UNROOTED although outgroup taxon '" <<
					tree.root->name << "' is drawn at root" << endl;
			if (params.partition_file) 
				out << "NOTE: Branch lengths are weighted average over all partitions" << endl <<
					   "      (weighted by the number of sites in the partitions)" << endl;
			if (params.aLRT_replicates > 0 || params.gbo_replicates || (params.num_bootstrap_samples && params.compute_ml_tree)) {
				out << "Numbers in parentheses are ";
				if (params.aLRT_replicates > 0) out << "SH-aLRT supports";
				if (params.num_bootstrap_samples && params.compute_ml_tree) out << " standard bootstrap supports";
				if (params.gbo_replicates) out << " ultra-fast bootstrap supports";
				out << " (%)" << endl;
			}
			out << endl;
			reportTree(out, params, tree, tree.getBestScore(), tree.logl_variance);

			if (tree.isSuperTree()) {
				PhyloSuperTree *stree = (PhyloSuperTree*)&tree;
				int part = 0;
				for (PhyloSuperTree::iterator it = stree->begin(); it != stree->end(); it++, part++) {
					out << "FOR PARTITION " << stree->part_info[part].name <<":" << endl << endl;
					string root_name;
					if (params.root) root_name = params.root;
					else root_name = (*it)->aln->getSeqName(0);
					(*it)->root = (*it)->findNodeName(root_name);
					assert((*it)->root);
					reportTree(out, params, *(*it),(*it)->computeLikelihood(),
						(*it)->computeLogLVariance());
				}
			}


		}
        /*
                        if (params.write_intermediate_trees) {
                                out << endl << "CONSENSUS OF INTERMEDIATE TREES" << endl << "-----------------------" << endl << endl
                                        << "Number of intermediate trees: " << tree.stop_rule.getNumIterations() << endl
                                        << "Split threshold: " << params.split_threshold << endl
                                        << "Burn-in: " << params.tree_burnin << endl << endl;
                        }*/


		if (params.consensus_type == CT_CONSENSUS_TREE) {
        	out << "CONSENSUS TREE" << endl << "--------------" << endl << endl;
        	out << "Consensus tree is constructed from " << (params.num_bootstrap_samples ? params.num_bootstrap_samples : params. gbo_replicates) << " bootstrap trees" << endl
				<< "Branches with bootstrap support >" << floor(params.split_threshold * 1000)/10 << "% are kept";
			if (params.split_threshold == 0.0) out << " (extended consensus)";
			if (params.split_threshold == 0.5) out << " (majority-rule consensus)";
			if (params.split_threshold >= 0.99) out << " (strict consensus)";

        	out << endl << "Branch lengths are optimized by maximum likelihood on original alignment" << endl;
            out << "Numbers in parentheses are bootstrap supports (%)" << endl << endl;

        	string con_file = params.out_prefix;
        	con_file += ".contree";
        	bool rooted = false;

			tree.freeNode();
			tree.readTree(con_file.c_str(), rooted);
			double fixed_length = 0.01;
			tree.fixNegativeBranch(fixed_length);
			tree.setAlignment(tree.aln);
		    tree.initializeAllPartialLh();
			if (tree.isSuperTree()) ((PhyloSuperTree*)&tree)->mapTrees();
			tree.optimizeAllBranches();
			tree.printTree(con_file.c_str(), WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA);
			tree.sortTaxa();
        	tree.drawTree(out);
        	out << endl << "Consensus tree in newick format: " << endl << endl;
	        tree.printResultTree(out);
        	out << endl << endl;
		}

        time_t cur_time;
        time(&cur_time);

        char *date_str;
        date_str = ctime(&cur_time);
        out.unsetf(ios_base::fixed);
        out << "TIME STAMP" << endl << "----------" << endl << endl << "Date and time: " << date_str <<
                "Running time: " << (double) params.run_time / CLOCKS_PER_SEC << " seconds" << endl << endl;

		//reportCredits(out); // not needed, now in the manual
        out.close();

    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, outfile);
    }

    cout << endl << "Analysis results written to: " << endl
            << "  IQ-TREE report:           " << params.out_prefix << ".iqtree" << endl;
    if (params.compute_ml_tree)
        cout << "  Maximum-likelihood tree:  " << params.out_prefix << ".treefile" << endl;
    if (!params.user_file)
        cout << "  BIONJ tree:               " << params.out_prefix << ".bionj" << endl;
    if (!params.dist_file) {
        cout << "  Juke-Cantor distances:    " << params.out_prefix << ".jcdist" << endl;
        if (params.compute_ml_dist)
            cout << "  Likelihood distances:     " << params.out_prefix << ".mldist" << endl;
        if (params.partition_file)
            cout << "  Concatenated alignment:   " << params.out_prefix << ".concat" << endl;
    }
	if (tree.getRate()->getGammaShape() > 0)
        cout << "  Gamma-distributed rates:  " << params.out_prefix << ".rate" << endl;


    if (tree.getRate()->isSiteSpecificRate() || tree.getRate()->getPtnCat(0) >= 0)
        cout << "  Site-rates by MH model:   " << params.out_prefix << ".rate" << endl;

    if (params.print_site_lh)
		cout << "  Site log-likelihoods:     " << params.out_prefix << ".sitelh" << endl;

    if (params.write_intermediate_trees)
		cout << "  All intermediate trees:   " << params.out_prefix << ".treels" << endl;

	if (params.gbo_replicates) {
		cout << endl 
			 <<"Ultra-fast bootstrap approximation results written to:" << endl
			 << "  Split support values:     " << params.out_prefix << ".splits" << endl
			 << "  Consensus tree:           " << params.out_prefix << ".contree" << endl;
	}

	cout <<     "  Screen log file:          " << params.out_prefix << ".log" << endl;
/*	if (original_model == "WHTEST")
		cout <<"  WH-TEST report:           " << params.out_prefix << ".whtest" << endl;*/
    cout << endl;

}

void checkZeroDist(Alignment *aln, double *dist) {
    int ntaxa = aln->getNSeq();
    IntVector checked;
    checked.resize(ntaxa, 0);
    int i, j;
    for (i = 0; i < ntaxa - 1; i++) {
        if (checked[i]) continue;
        string str = "";
        bool first = true;
        for (j = i + 1; j < ntaxa; j++)
            if (dist[i * ntaxa + j] <= 1e-6) {
                if (first) str = "ZERO distance between sequences " + aln->getSeqName(i);
                str += ", " + aln->getSeqName(j);
                checked[j] = 1;
                first = false;
            }
        checked[i] = 1;
        if (str != "") outWarning(str);
    }
}

void printSiteLh(const char*filename, IQPTree &tree, double *ptn_lh = NULL, bool append = false) {
	int i;
	double *pattern_lh;
	if (!ptn_lh) {
		pattern_lh = new double[tree.getAlnNPattern()];
		tree.computePatternLikelihood(pattern_lh);
	} else pattern_lh = ptn_lh;

	try {
		ofstream out;
		out.exceptions(ios::failbit | ios::badbit);
		if (append) {
			out.open(filename, ios::out | ios::app);
		}
		else {
			out.open(filename);
			out << tree.getAlnNSite() << endl;
		}
		IntVector pattern_index;
		tree.aln->getSitePatternIndex(pattern_index);
		out << "Site_Lh   ";
		for (i = 0; i < tree.getAlnNSite(); i++)
			out << " " << pattern_lh[pattern_index[i]];
		out << endl;
		out.close();
		if (!append) cout << "Site log-likelihoods printed to " << filename << endl;
	} catch (ios::failure) {
		outError(ERR_WRITE_OUTPUT, filename);
	}

	if (!ptn_lh) delete [] pattern_lh;
}

void runPhyloAnalysis(Params &params, string &original_model, Alignment *alignment, IQPTree &tree) {
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
/*
    if (!tree.dist_matrix) {
        tree.dist_matrix = new double[alignment->getNSeq() * alignment->getNSeq()];
        memset(tree.dist_matrix, 0, sizeof (double) * alignment->getNSeq() * alignment->getNSeq());
    }*/
    clock_t begin_time = clock();
    params.startTime = begin_time;
    if (params.dist_file) {
        cout << "Reading distance matrix file " << params.dist_file << " ..." << endl;
    } else {
        cout << "Computing Juke-Cantor distances..." << endl;
    }
    string dist_file;
    double longest_dist = tree.computeDist(params, alignment, tree.dist_matrix, dist_file);
    checkZeroDist(alignment, tree.dist_matrix);
    if (longest_dist > MAX_GENETIC_DIST * 0.99) {
		cout << "Some distances are too long, computing observed distances..." << endl;
		longest_dist = tree.computeObsDist(params, alignment, tree.dist_matrix, dist_file);
		assert(longest_dist <= 1.0);
    }
    // start the search with user-defined tree
    if (params.user_file) {
        cout << "Reading user tree file " << params.user_file << " ..." << endl;
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
    double fixed_length = 0.001;
    int fixed_number = tree.fixNegativeBranch(fixed_length);
    if (fixed_number) {
        cout << "WARNING: " << fixed_number << " branches have non-positive lengths and initialized to " << fixed_length << endl;
        if (verbose_mode >= VB_DEBUG) {
            tree.printTree(cout);
            cout << endl;
        }
    }
    t_begin = clock();
    bool test_only = params.model_name == "TESTONLY";
    /* initialize substitution model */
    if (params.model_name == "TEST" || params.model_name == "TESTONLY") {
        params.model_name = modelTest(params, &tree);
        if (test_only) {
            return;
            t_end = clock();
            params.run_time = (t_end - t_begin);
            printf("Time used: %8.6f seconds.\n", (double) params.run_time / CLOCKS_PER_SEC);
        }
    }

    if (params.model_name == "WHTEST") {
    	if (alignment->num_states != 4)
    		outError("Weiss & von Haeseler test of model homogeneity only works for DNA");
    	params.model_name="GTR+G";
    }

    assert(tree.aln);
    tree.optimize_by_newton = params.optimize_by_newton;
    tree.sse = params.SSE;
	if (params.gbo_replicates) params.speed_conf = 1.0;
    if (params.speed_conf == 1.0)
        tree.disableHeuristic();
    else
        tree.setSpeed_conf(params.speed_conf);
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
    clock_t t_tree_search_start, t_tree_search_end;
    t_tree_search_start = clock();
    cout << endl;
	tree.setParams(params);
    cout << "ML-TREE SEARCH START WITH THE FOLLOWING PARAMETERS:" << endl;
    cout << "Model of evolution: " << tree.getModelName() << " with ";
    switch (tree.getModel()->getFreqType()) {
		case FREQ_EQUAL: cout << "equal"; break;
		case FREQ_EMPIRICAL: cout << "counted"; break;
		case FREQ_USER_DEFINED: cout << "user-defined"; break;
		case FREQ_ESTIMATE: cout << "optimized"; break;
		default: outError("Wrong specified state frequencies");
	}
	cout << " frequencies (" << model_df << " free parameters)" << endl;
    cout << "Fixed branch lengths: " << ((params.fixed_branch_length) ? "Yes" : "No") << endl;
    cout << "Lambda for local search: " << params.lambda << endl;

    if (params.speed_conf != 1.0) {
        cout <<"Confidence value for speed up NNI: ";
		if (params.new_heuristic)
			cout << "Using 50%*"<< params.speed_conf <<endl;
		else cout << "N" << params.speed_conf << " * delta" << params.speed_conf <<endl;
    } else {
        cout <<"Speed up NNI: disabled " << endl;
    }
	cout << "NNI cutoff: " << params.nni_cutoff << endl;


//    if (params.parsimony) {
//		int score = tree.computeParsimonyScore();
//		cout << "User tree has parsimony score of " << score << endl;
//
//		// NNI with parsimony function
//		tree.searchNNI();
//	}

    cout.precision(10);
    //cout << "User tree has likelihood score of " << tree.computeLikelihood() << endl;
    if (params.parsimony) {
        tree.enable_parsimony = true;
        tree.pars_scores = new double[3000];
        tree.lh_scores = new double[3000];
        for (int i = 0; i < 3000; i++) {
            tree.pars_scores[i] = 0;
            tree.lh_scores[i] = 0;
        }
        //cout << "Parsimony score: " << tree.computeParsimonyScore() << endl;
        tree.cur_pars_score = tree.computeParsimony();
        //cout << "Fast parsimony score: " << tree.cur_pars_score << endl;
    }

    /* optimize model parameters */
    cout << endl;
    cout << "Optimizing model parameters and branch lengths" << endl;
    double bestTreeScore = tree.getModelFactory()->optimizeParameters(params.fixed_branch_length);
    cout.precision(10);
    cout << "Log-likelihood of the current tree: " << bestTreeScore << endl;

    //Update tree score
    tree.curScore = bestTreeScore;
	if (tree.isSuperTree()) ((PhyloSuperTree*)&tree)->computeBranchLengths();
	/*
    if ((tree.getModel()->name == "JC") && tree.getRate()->getNDim() == 0)
        params.compute_ml_dist = false;*/

    if (!params.dist_file && params.compute_ml_dist) {
        stringstream best_tree_string;
        tree.printTree(best_tree_string, WT_BR_LEN + WT_TAXON_ID);
        cout << "Computing ML distances based on estimated model parameters...";
        double *ml_dist = NULL;
        longest_dist = tree.computeDist(params, alignment, ml_dist, dist_file);
        cout << " " << double(clock() - begin_time) / CLOCKS_PER_SEC << " sec" << endl;

	    if (longest_dist > MAX_GENETIC_DIST * 0.99) {
			cout << "Some ML distances are too long, using old distances..." << endl;
	    } else {
			memmove(tree.dist_matrix, ml_dist, sizeof(double)*alignment->getNSeq()*alignment->getNSeq());

			if (!params.user_file) {
				tree.computeBioNJ(params, alignment, dist_file); // create BioNJ tree
				tree.fixNegativeBranch(fixed_length);
				if (tree.isSuperTree()) ((PhyloSuperTree*)&tree)->mapTrees();
				if (!params.fixed_branch_length)
					tree.curScore = tree.optimizeAllBranches();
				else
					tree.curScore = tree.computeLikelihood();
				cout << "Log-likelihood of the new BIONJ tree: " << tree.curScore << endl;
				if (tree.curScore < bestTreeScore - 1e-5) {
					cout << "The new tree is worse, rolling back the first BIONJ tree..." << endl;
					tree.rollBack(best_tree_string);
					if (tree.isSuperTree()) {
						((PhyloSuperTree*)&tree)->mapTrees();
						tree.optimizeAllBranches();
					}
					tree.curScore = tree.computeLikelihood();
					cout << "Backup log-likelihood: " << tree.curScore << endl;
				}
				double elapsedTime = getCPUTime(params.startTime);
				cout << "Time elapsed: " << elapsedTime << endl;
			}
        }
        delete [] ml_dist;
    }

    /* do NNI with likelihood function */
	//bool saved_estimate_nni = estimate_nni_cutoff;
	//estimate_nni_cutoff = false; // do not estimate NNI cutoff based on initial BIONJ tree

    if (params.min_iterations > 0) {
        cout << endl;
        cout << "Performing local search with NNI moves ... " << endl;
        //cout << "Current tree likelihood: " << tree.optimizeNNIBranches() << endl;
        //tree.optimizeAllBranches();
        clock_t nniBeginClock, nniEndClock;
        nniBeginClock = clock();
        tree.optimizeNNI();
		tree.curScore = tree.optimizeAllBranches(100, 0.0001);
        //tree.optimizeNNINew();
        nniEndClock = clock();
        cout << "First NNI search required :" << (double) (nniEndClock - nniBeginClock) / CLOCKS_PER_SEC << "s" << endl;
        if (tree.curScore > bestTreeScore) {
            bestTreeScore = tree.curScore ;
            cout << "Found new best tree log-likelihood : " << bestTreeScore << endl;
        } else {
            cout << "The local search cannot improve the tree likelihood :( " << endl;
        }
		double elapsedTime = getCPUTime(params.startTime);
		cout << "CPU time elapsed: " << elapsedTime << endl;
    }

	//estimate_nni_cutoff = saved_estimate_nni;

	if (original_model == "WHTEST") {
		cout << endl << "Testing model homogeneity by Weiss & von Haeseler (2003)..." << endl;
		WHTest(params, tree);
	}


    /*double sum_scaling = 1.0;
    if (!tree.checkEqualScalingFactor(sum_scaling))
            cout << "Scaling factor not equal along the tree" << endl;*/

    NodeVector pruned_taxa;
    StrVector linked_name;
    double *saved_dist_mat = tree.dist_matrix;
    double *pattern_lh;
    int num_low_support;
    clock_t mytime;

    pattern_lh = new double[tree.getAlnNPattern()];

    if (params.aLRT_threshold <= 100 && (params.aLRT_replicates > 0 || params.localbp_replicates > 0)) {
        mytime = clock();
        cout << "Testing tree branches by SH-like aLRT with " << params.aLRT_replicates << " replicates..." << endl;
        tree.setRootNode(params.root);
        tree.computePatternLikelihood(pattern_lh, &tree.curScore);
        num_low_support = tree.testAllBranches(params.aLRT_threshold, tree.curScore, pattern_lh, params.aLRT_replicates, params.localbp_replicates);
        tree.printResultTree();
        cout << "  " << (((double) clock()) - mytime) / CLOCKS_PER_SEC << " sec." << endl;
        cout << num_low_support << " branches show low support values (<= " << params.aLRT_threshold << "%)" << endl;

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
        tree.clearAllPartialLH();
        tree.curScore = tree.optimizeAllBranches();
        //cout << "Log-likelihood	after reoptimizing model parameters: " << tree.curScore << endl;
        tree.curScore = tree.optimizeNNI();
        cout << "Log-likelihood after optimizing partial tree: " << tree.curScore << endl;
        /*
        pattern_lh = new double[tree.getAlnSize()];
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

	//tree.setParams(params);


    /* do iterated local search */
//    if (ils && params.min_iterations > 1) {
//        cout << "Start doing iterated local search " << endl;
//        tree.doILS(params, perLevel);
//        cout << "Optimizing model parameters" << endl;
//        double endScore = tree.getModelFactory()->optimizeParameters(params.fixed_branch_length);
//        cout << "Best score found : " << endScore << endl;
//    }

	/* evaluating all trees in user tree file */

    /* do the IQPNNI */
    if (params.k_representative > 0/* && params.min_iterations >= 1*/) {
        cout << endl << "START IQPNNI SEARCH WITH THE FOLLOWING PARAMETERS" << endl;
        cout << "Number of representative leaves   : " << params.k_representative << endl;
        cout << "Probability of deleting sequences : " << tree.getProbDelete() << endl;
        cout << "Number of iterations              : ";
        if (params.stop_condition == SC_FIXED_ITERATION)
            cout << params.min_iterations << endl;
        else cout << "predicted in [" << params.min_iterations << "," <<
                params.max_iterations << "] (confidence " << params.stop_confidence << ")" << endl;
        cout << "Important quartet assessed on     : " << ((params.iqp_assess_quartet == IQP_DISTANCE) ?
        	"Distance" : ((params.iqp_assess_quartet == IQP_PARSIMONY) ? "Parsimony" : "Bootstrap")) << endl;
       	cout << "SSE instructions                  : " << ((tree.sse)? "Yes" : "No") << endl;
		cout << "Branch length optimization method : " << ((tree.optimize_by_newton) ? "Newton" : "Brent") << endl;
        cout << endl;
        tree.doIQPNNI();
    } else {
        /* do SPR with likelihood function */
		if (params.tree_spr) {
			//tree.optimizeSPRBranches();
			cout << "Doing SPR Search" << endl;
			cout << "Start tree.optimizeSPR()" << endl;
			double spr_score = tree.optimizeSPR();
			cout << "Finish tree.optimizeSPR()" << endl;
			//double spr_score = tree.optimizeSPR(tree.curScore, (PhyloNode*) tree.root->neighbors[0]->node);
			if (spr_score <= tree.curScore) {
				cout << "SPR search did not found any better tree" << endl;
			} else {
				tree.curScore = spr_score;
				cout << "Found new BETTER SCORE by SPR: " << spr_score << endl;
				double nni_score = tree.optimizeNNI();
				cout << "Score by NNI: " << nni_score << endl;
			}
		}

    }

    if (!pruned_taxa.empty()) {
        tree.disableHeuristic();
        cout << "Restoring full tree..." << endl;
        tree.restoreStableClade(alignment, pruned_taxa, linked_name);
        delete [] tree.dist_matrix;
        tree.dist_matrix = saved_dist_mat;
        tree.initializeAllPartialLh();
        tree.clearAllPartialLH();
        tree.curScore = tree.optimizeAllBranches();
        //cout << "Log-likelihood	after reoptimizing model parameters: " << tree.curScore << endl;
        tree.curScore = tree.optimizeNNI();
        cout << "Log-likelihood	after reoptimizing full tree: " << tree.curScore << endl;
    }

    if (params.min_iterations) {
        cout << endl;
        cout << "Optimizing model parameters" << endl;
        tree.setBestScore(tree.getModelFactory()->optimizeParameters(params.fixed_branch_length));
    } else
        tree.setBestScore(tree.curScore);
    cout << endl;
    cout.precision(10);
    cout << "BEST SCORE FOUND : " << tree.getBestScore() << endl;
    t_tree_search_end = clock();
    double treeSearchTime = (double) (t_tree_search_end - t_tree_search_start) / CLOCKS_PER_SEC;

    /* root the tree at the first sequence */
    tree.root = tree.findLeafName(alignment->getSeqName(0));
    assert(tree.root);

    double myscore = tree.getBestScore();
    tree.computePatternLikelihood(pattern_lh, &myscore);
    //tree.computeLikelihood(pattern_lh);

	// compute logl variance
	tree.logl_variance = tree.computeLogLVariance();

	if (params.print_site_lh) {
		string site_lh_file = params.out_prefix;
		site_lh_file += ".sitelh";
		printSiteLh(site_lh_file.c_str(), tree, pattern_lh);
	}


    if (params.mvh_site_rate) {
	    RateMeyerHaeseler *rate_mvh = new RateMeyerHaeseler(params.rate_file, &tree, params.rate_mh_type);
        cout << endl << "Computing site-specific rates by " << rate_mvh->full_name << "..." << endl;
        rate_mvh->runIterativeProc(params, tree);
    	cout << endl << "BEST SCORE FOUND : " << tree.getBestScore() << endl;
		string mhrate_file = params.out_prefix;
		mhrate_file += ".mhrate";
		tree.getRate()->writeSiteRates(mhrate_file.c_str());

		if (params.print_site_lh) {
			string site_lh_file = params.out_prefix;
			site_lh_file += ".mhsitelh";
			printSiteLh(site_lh_file.c_str(), tree);
		}
    }

    if ((params.aLRT_replicates > 0 || params.localbp_replicates > 0)) {
        mytime = clock();
        cout << endl;
        cout << "Testing tree branches by SH-like aLRT with " << params.aLRT_replicates << " replicates..." << endl;
        tree.setRootNode(params.root);
		//if (tree.isSuperTree()) ((PhyloSuperTree*)&tree)->mapTrees();
        num_low_support = tree.testAllBranches(params.aLRT_threshold, myscore, pattern_lh, params.aLRT_replicates, params.localbp_replicates);
        //cout << num_low_support << " branches show low support values (<= " << params.aLRT_threshold << "%)" << endl;
        cout << "CPU Time used:  " << (((double) clock()) - mytime) / CLOCKS_PER_SEC << " sec." << endl;
        //delete [] pattern_lh;
/*
		string out_file = params.out_prefix;
		out_file += ".alrt";
		tree.writeInternalNodeNames(out_file);

		cout << "Support values written to " << out_file << endl;*/
    }


	string rate_file = params.out_prefix;
	rate_file += ".rate";
	tree.getRate()->writeSiteRates(rate_file.c_str());

	if (tree.isSuperTree()) {
		PhyloSuperTree *stree = (PhyloSuperTree*) &tree;
		int part = 0;
		for (PhyloSuperTree::iterator it = stree->begin(); it != stree->end(); it++, part++) {
			rate_file = params.out_prefix;
			rate_file = rate_file + "." + stree->part_info[part].name + ".rate";
			(*it)->getRate()->writeSiteRates(rate_file.c_str());
		}
	}

    if (params.gbo_replicates > 0) {
		if (!params.online_bootstrap)
			runGuidedBootstrap(params, alignment, tree);
		else
			tree.summarizeBootstrap(params);
	}

	cout << "Total tree length: " << tree.treeLength() << endl;

    t_end = clock();
    params.run_time = (t_end - t_begin);
    cout << endl;
    cout << "CPU time used for tree reconstruction: " << treeSearchTime << "s (" << convert_time(treeSearchTime) << ") "<< endl;
    cout << "CPU total time used: " << (double) params.run_time / CLOCKS_PER_SEC << " (" << convert_time((double) params.run_time / CLOCKS_PER_SEC) << " )"<< "s" << endl;
    //printf( "Total time used: %8.6f seconds.\n", (double) params.run_time / CLOCKS_PER_SEC );

    tree.printResultTree();
    if (params.out_file)
        tree.printTree(params.out_file);
        //tree.printTree(params.out_file,WT_BR_LEN_FIXED_WIDTH);

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

    delete [] pattern_lh;

/*	if (tree.getRate()->isSiteSpecificRate() || tree.getRate()->getPtnCat(0) >= 0) {
		string rate_file = params.out_prefix;
		rate_file += ".mhrate";
		tree.getRate()->writeSiteRates(rate_file.c_str());
	}*/
}


void evaluateTrees(Params &params, IQPTree *tree) {
	if (!params.treeset_file) return;
	MTreeSet trees(params.treeset_file, params.is_rooted, params.tree_burnin, params.tree_max_count);
	//if (trees.size() == 1) return;
	string tree_file = params.treeset_file;
	tree_file += ".eval";
	ofstream treeout;
	if (!params.fixed_branch_length) {
		treeout.open(tree_file.c_str());
	}
	string score_file = params.treeset_file;
	score_file += ".treelh";
	ofstream scoreout(score_file.c_str());
	for (MTreeSet::iterator it = trees.begin(); it != trees.end(); it++) {
		cout << "Tree " << (it-trees.begin())+1;
		tree->copyTree(*it);
		double fixed_length = 0.01;
		//int fixed_number = tree->fixNegativeBranch(fixed_length);
		tree->fixNegativeBranch(fixed_length);
		tree->initializeAllPartialLh();
		if (tree->isSuperTree()) ((PhyloSuperTree*)tree)->mapTrees();
		if (!params.fixed_branch_length) {
			tree->curScore = tree->optimizeAllBranches();
			treeout << "[ lh=" << tree->curScore << " ]";
			tree->printTree(treeout);
			treeout << endl;
		}
		else
			tree->curScore = tree->computeLikelihood();
		cout << " / LogL: " << tree->curScore << endl;
		if (params.print_site_lh) {
			string site_lh_file = params.treeset_file;
			site_lh_file += ".sitelh";
			printSiteLh(site_lh_file.c_str(), *tree, NULL, (it != trees.begin()));
		}
		scoreout << tree->curScore << endl;
	}
	if (!params.fixed_branch_length) {
		treeout.close();
		cout << "Trees with optimized branch lengths printed to " << tree_file << endl;
	}
	scoreout.close();
	cout << "Tree log-likelihoods printed to " << score_file << endl;
}

void runPhyloAnalysis(Params &params) {
	Alignment *alignment;
	IQPTree *tree;
	if (params.partition_file) {
		tree = new PhyloSuperTree(params);
		alignment = tree->aln;
	} else {
		alignment = new Alignment(params.aln_file, params.sequence_type, params.intype);
		tree = new IQPTree(alignment);
    }
    string original_model = params.model_name;
	if (params.concatenate_aln) {
		Alignment aln(params.concatenate_aln, params.sequence_type, params.intype);
		cout << "Concatenating " << params.aln_file << " with " << params.concatenate_aln << " ..." << endl;
		alignment->concatenateAlignment(&aln);
	}

    if (params.aln_output) {
		if (params.gap_masked_aln) {
			Alignment out_aln;
		    Alignment masked_aln(params.gap_masked_aln, params.sequence_type, params.intype);
			out_aln.createGapMaskedAlignment(&masked_aln, alignment);
        	out_aln.printPhylip(params.aln_output, false, params.aln_site_list,
        		params.aln_nogaps, params.ref_seq_name);
			string str = params.gap_masked_aln;
			str += ".sitegaps";
			out_aln.printSiteGaps(str.c_str());
		} else
		if (params.aln_output_format == ALN_PHYLIP)
        	alignment->printPhylip(params.aln_output, false, params.aln_site_list,
        		params.aln_nogaps, params.ref_seq_name);
        else if (params.aln_output_format == ALN_FASTA)
        	alignment->printFasta(params.aln_output, false, params.aln_site_list,
        	params.aln_nogaps, params.ref_seq_name);
    } else if (params.gbo_replicates > 0 && params.user_file && params.second_tree) {
    	runGuidedBootstrap(params, alignment, *tree);
	} else if (params.avh_test) {
		runAvHTest(params, alignment, *tree);
    } else if (params.num_bootstrap_samples == 0) {
		alignment->checkGappySeq();
        runPhyloAnalysis(params, original_model, alignment, *tree);
		if (params.gbo_replicates && params.online_bootstrap) {

			cout << endl << "Computing consensus tree..." << endl;
			string splitsfile = params.out_prefix;
			splitsfile += ".splits.nex";
			//cout << splitsfile << endl;
			computeConsensusTree(splitsfile.c_str(), 0, 1e6, -1,
					params.split_threshold, NULL, params.out_prefix, NULL, &params);
		}
        if (original_model != "TESTONLY")
            reportPhyloAnalysis(params, original_model, *alignment, *tree);
        if (params.treeset_file) {
        	evaluateTrees(params, tree);
        }
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
        string bootlh_name = params.out_prefix;
        bootlh_name += ".bootlh";
        // first empty the boottrees file
        try {
            ofstream tree_out;
            tree_out.exceptions(ios::failbit | ios::badbit);
            tree_out.open(boottrees_name.c_str());
            tree_out.close();
        } catch (ios::failure) {
            outError(ERR_WRITE_OUTPUT, boottrees_name);
        }

        // empty the bootaln file
        try {
            ofstream tree_out;
            tree_out.exceptions(ios::failbit | ios::badbit);
            tree_out.open(bootaln_name.c_str());
            tree_out.close();
        } catch (ios::failure) {
            outError(ERR_WRITE_OUTPUT, bootaln_name);
        }

		clock_t start_time = clock();

        // do bootstrap analysis
        for (int sample = 0; sample < params.num_bootstrap_samples; sample++) {
            cout << endl << "===> START BOOTSTRAP REPLICATE NUMBER " << sample + 1 << endl << endl;

            Alignment* bootstrap_alignment;
            cout << "Creating bootstrap alignment..." << endl;
            if (alignment->isSuperAlignment())
            	bootstrap_alignment = new SuperAlignment;
            else
            	bootstrap_alignment = new Alignment;
            bootstrap_alignment->createBootstrapAlignment(alignment);
			if (params.print_tree_lh) {
				double prob;
				bootstrap_alignment->multinomialProb(*alignment, prob);
				ofstream boot_lh;
				if (sample == 0)
					boot_lh.open(bootlh_name.c_str());
				else
					boot_lh.open(bootlh_name.c_str(), ios_base::out | ios_base::app);
				boot_lh << "0\t" << prob << endl;
				boot_lh.close();
			}
            IQPTree *boot_tree;
            if (alignment->isSuperAlignment())
            	boot_tree = new PhyloSuperTree((SuperAlignment*)bootstrap_alignment, (PhyloSuperTree*)tree);
           	else
            	boot_tree = new IQPTree(bootstrap_alignment);
            bootstrap_alignment->printPhylip(bootaln_name.c_str(), true);
            runPhyloAnalysis(params, original_model, bootstrap_alignment, *boot_tree);
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
                reportPhyloAnalysis(params, original_model, *bootstrap_alignment, *boot_tree);
            delete bootstrap_alignment;
        }

		if (params.consensus_type == CT_CONSENSUS_TREE) {

			cout << endl << "===> COMPUTE CONSENSUS TREE FROM " <<
					params.num_bootstrap_samples << " BOOTSTRAP TREES" << endl << endl;
			computeConsensusTree(boottrees_name.c_str(), 0, 1e6, -1,
					params.split_threshold, NULL, params.out_prefix, NULL, &params);
		}

        if (params.compute_ml_tree) {
            cout << endl << "===> START ANALYSIS ON THE ORIGINAL ALIGNMENT" << endl << endl;
            params.aLRT_replicates = saved_aLRT_replicates;
            runPhyloAnalysis(params, original_model, alignment, *tree);

            cout << endl << "===> ASSIGN BOOTSTRAP SUPPORTS TO THE TREE FROM ORIGINAL ALIGNMENT" << endl << endl;
            MExtTree ext_tree;
            assignBootstrapSupport(boottrees_name.c_str(), 0, 1e6,
                    treefile_name.c_str(), false, treefile_name.c_str(),
                    params.out_prefix, ext_tree, NULL, &params);
            tree->copyTree(&ext_tree);
            reportPhyloAnalysis(params, original_model, *alignment, *tree);
        }  else if (params.consensus_type == CT_CONSENSUS_TREE) {
        	int mi = params.min_iterations;
        	STOP_CONDITION sc = params.stop_condition;
        	params.min_iterations = 0;
        	params.stop_condition = SC_FIXED_ITERATION;
            runPhyloAnalysis(params, original_model, alignment, *tree);
        	params.min_iterations = mi;
        	params.stop_condition = sc;
    		tree->setIQPIterations(params.stop_condition, params.stop_confidence, params.min_iterations, params.max_iterations);
            reportPhyloAnalysis(params, original_model, *alignment, *tree);
        } else    cout << endl;

		cout << "CPU total time for bootstrap: " << (clock() - start_time) / CLOCKS_PER_SEC << " seconds." << endl << endl;
		cout << "Non-parametric bootstrap results written to:" << endl
			<< "  Bootstrap alignments:     " << params.out_prefix << ".bootaln" << endl
			<< "  Bootstrap trees:          " << params.out_prefix << ".boottrees" << endl;
		if (params.consensus_type == CT_CONSENSUS_TREE)
			cout << "  Consensus tree:           " << params.out_prefix << ".contree" << endl;
		cout << endl;
    }

	delete tree;
	delete alignment;
}

void assignBootstrapSupport(const char *input_trees, int burnin, int max_count, const char *target_tree, bool rooted,
        const char *output_tree, const char *out_prefix, MExtTree &mytree,
        const char* tree_weight_file, Params *params)
{
    //bool rooted = false;
    // read the tree file
	cout << "Reading tree " << target_tree << " ..." << endl;
    mytree.init(target_tree, rooted);
    // reindex the taxa in the tree to aphabetical names
    NodeVector taxa;
    mytree.getTaxa(taxa);
    sort(taxa.begin(), taxa.end(), nodenamecmp);
    int i = 0;
    for (NodeVector::iterator it = taxa.begin(); it != taxa.end(); it++) {
        (*it)->id = i++;
    }

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

    // read the bootstrap tree file
	double scale=100.0;
	if (params->scaling_factor > 0) scale = params->scaling_factor;

    MTreeSet boot_trees;
    if (params && detectInputFile((char*)input_trees) == IN_NEXUS) {
		sg.init(*params);
		for (SplitGraph::iterator it = sg.begin(); it != sg.end(); it++)
			hash_ss.insertSplit((*it), (*it)->getWeight());
		StrVector sgtaxname;
		sg.getTaxaName(sgtaxname);
		i = 0;
		for (StrVector::iterator sit = sgtaxname.begin(); sit != sgtaxname.end(); sit++, i++) {
			Node *leaf = mytree.findLeafName(*sit);
			if (!leaf) outError("Tree does not contain taxon ", *sit);
			leaf->id = i;
		}
		scale /= sg.maxWeight();
    } else {
    	boot_trees.init(input_trees, rooted, burnin, max_count, tree_weight_file);
    	boot_trees.convertSplits(taxname, sg, hash_ss, SW_COUNT, -1);
    	scale /= boot_trees.sumTreeWeights();
    }
    //sg.report(cout);
	cout << "Rescaling split weights by " << scale << endl;
	if (params->scaling_factor < 0)
	  sg.scaleWeight(scale, true);
	else {
	  sg.scaleWeight(scale, false, params->numeric_precision);
	}

     cout << sg.size() << " splits found" << endl;
   // compute the percentage of appearance
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

	if (out_prefix)
		out_file = out_prefix;
	else
		out_file = target_tree;
	out_file += ".supval";
	mytree.writeInternalNodeNames(out_file);

    cout << "Support values written to " << out_file << endl;

}

void computeConsensusTree(const char *input_trees, int burnin, int max_count, double cutoff, double weight_threshold,
        const char *output_tree, const char *out_prefix, const char *tree_weight_file, Params *params) {
    bool rooted = false;

    // read the bootstrap tree file
	/*
    MTreeSet boot_trees(input_trees, rooted, burnin, tree_weight_file);
    string first_taxname = boot_trees.front()->root->name;
    //if (params.root) first_taxname = params.root;

    SplitGraph sg;

    boot_trees.convertSplits(sg, cutoff, SW_COUNT, weight_threshold);*/
	
	
    //sg.report(cout);

    SplitGraph sg;
    SplitIntMap hash_ss;
    // make the taxa name
    //vector<string> taxname;
    //taxname.resize(mytree.leafNum);
    //mytree.getTaxaName(taxname);

    // read the bootstrap tree file
	double scale=100.0;
	if (params->scaling_factor > 0) scale = params->scaling_factor;

    MTreeSet boot_trees;
    if (params && detectInputFile((char*)input_trees) == IN_NEXUS) {
		 char *user_file = params->user_file;
		params->user_file = (char*)input_trees;
		sg.init(*params);
		params->user_file = user_file;
		for (SplitGraph::iterator it = sg.begin(); it != sg.end(); it++)
			hash_ss.insertSplit((*it), (*it)->getWeight());
/*		StrVector sgtaxname;
		sg.getTaxaName(sgtaxname);
		i = 0;
		for (StrVector::iterator sit = sgtaxname.begin(); sit != sgtaxname.end(); sit++, i++) {
			Node *leaf = mytree.findLeafName(*sit);
			if (!leaf) outError("Tree does not contain taxon ", *sit);
			leaf->id = i;
		}*/
		scale /= sg.maxWeight();
    } else {
    	boot_trees.init(input_trees, rooted, burnin, max_count, tree_weight_file);
    	boot_trees.convertSplits(sg, cutoff, SW_COUNT, weight_threshold);
    	scale /= boot_trees.sumTreeWeights();
		cout << sg.size() << " splits found" << endl;
    }
    //sg.report(cout);
	cout << "Rescaling split weights by " << scale << endl;
	if (params->scaling_factor < 0)
	  sg.scaleWeight(scale, true);
	else {
	  sg.scaleWeight(scale, false, params->numeric_precision);
	}


    cout << "Creating greedy consensus tree..." << endl;
    MTree mytree;
    SplitGraph maxsg;
    sg.findMaxCompatibleSplits(maxsg);
    
    if (verbose_mode >= VB_MED) maxsg.saveFileStarDot(cout);
    cout << "convert compatible split system into tree..." << endl;
    mytree.convertToTree(maxsg);
    //cout << "done" << endl;
    string taxname = sg.getTaxa()->GetTaxonLabel(0);
    Node *node = mytree.findLeafName(taxname);
    if (node) mytree.root = node;
   // mytree.scaleLength(100.0 / boot_trees.sumTreeWeights(), true);

   // mytree.getTaxaID(maxsg.getSplitsBlock()->getCycle());
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

void computeConsensusNetwork(const char *input_trees, int burnin, int max_count, double cutoff, double weight_threshold,
        const char *output_tree, const char *out_prefix, const char* tree_weight_file) {
    bool rooted = false;

    // read the bootstrap tree file
    MTreeSet boot_trees(input_trees, rooted, burnin, max_count, tree_weight_file);

    SplitGraph sg;
    //SplitIntMap hash_ss;

    boot_trees.convertSplits(sg, cutoff, SW_SUM, weight_threshold);

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


    sg.saveFile(out_file.c_str(), IN_NEXUS);
    cout << "Consensus network printed to " << out_file << endl;

}
