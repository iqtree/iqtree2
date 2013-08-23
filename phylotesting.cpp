/*
 * phylotesting.cpp
 *
 *  Created on: Aug 23, 2013
 *      Author: minh
 */



#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iqtree_config.h>
#include "phylotree.h"
#include "iqtree.h"
#include "phylosupertree.h"
#include "phylotesting.h"

#include "gtrmodel.h"
#include "modeldna.h"
#include "myreader.h"
#include "rateheterogeneity.h"
#include "rategamma.h"
#include "rateinvar.h"
#include "rategammainvar.h"
//#include "modeltest_wrapper.h"
#include "modelprotein.h"
#include "modelbin.h"
#include "modelcodon.h"
#include "timeutil.h"


//const int DNA_MODEL_NUM = 14;

const int BIN_MODEL_NUM = 2;
string bin_model_names[BIN_MODEL_NUM] = { "JC2", "GTR2" };

const int DNA_MODEL_NUM = 22;
string dna_model_names[DNA_MODEL_NUM] = { "JC", "F81", "K80", "HKY", "TNe",
		"TN", "K81", "K81u", "TPM2", "TPM2u", "TPM3", "TPM3u", "TIMe", "TIM",
		"TIM2e", "TIM2", "TIM3e", "TIM3", "TVMe", "TVM", "SYM", "GTR" };
/*string dna_model_names[DNA_MODEL_NUM] ={"JC", "F81", "K80", "HKY", "TNe", "TN", "K81", "K81u",
 "TIMe", "TIM", "TVMe", "TVM", "SYM", "GTR"};*/

const int AA_MODEL_NUM = 18;
string aa_model_names[AA_MODEL_NUM] = { "Dayhoff", "mtMAM", "JTT", "WAG",
		"cpREV", "mtREV", "rtREV", "mtART", "mtZOA", "VT", "LG", "DCMut", "PMB",
		"HIVb", "HIVw", "JTTDCMut", "FLU", "Blosum62" };

void printSiteLh(const char*filename, PhyloTree *tree, double *ptn_lh,
		bool append, const char *linename) {
	int i;
	double *pattern_lh;
	if (!ptn_lh) {
		pattern_lh = new double[tree->getAlnNPattern()];
		tree->computePatternLikelihood(pattern_lh);
	} else
		pattern_lh = ptn_lh;

	try {
		ofstream out;
		out.exceptions(ios::failbit | ios::badbit);
		if (append) {
			out.open(filename, ios::out | ios::app);
		} else {
			out.open(filename);
			out << 1 << " " << tree->getAlnNSite() << endl;
		}
		IntVector pattern_index;
		tree->aln->getSitePatternIndex(pattern_index);
		if (!linename)
			out << "Site_Lh   ";
		else {
			out.width(10);
			out << left << linename;
		}
		for (i = 0; i < tree->getAlnNSite(); i++)
			out << " " << pattern_lh[pattern_index[i]];
		out << endl;
		out.close();
		if (!append)
			cout << "Site log-likelihoods printed to " << filename << endl;
	} catch (ios::failure) {
		outError(ERR_WRITE_OUTPUT, filename);
	}

	if (!ptn_lh)
		delete[] pattern_lh;
}

/**
 * check if the model file contains correct information
 * @param model_file model file names
 * @param model_name (OUT) vector of model names
 * @param lh_scores (OUT) vector of tree log-likelihoods
 * @param df_vec (OUT) vector of degrees of freedom (or K)
 * @return TRUE if success, FALSE failed.
 */
bool checkModelFile(string model_file, StrVector &model_names,
		DoubleVector &lh_scores, IntVector &df_vec) {
	if (!fileExists(model_file))
		return false;
	cout << model_file << " exists, checking this file" << endl;
	ifstream in;
	try {
		in.exceptions(ios::failbit | ios::badbit);
		in.open(model_file.c_str());
		in.exceptions(ios::badbit);
		string str;
		in >> str;
		if (str != "Model")
			throw false;
		in >> str;
		if (str != "df")
			throw false;
		in >> str;
		if (str != "LnL")
			throw false;
		getline(in, str);
		while (!in.eof()) {
			in >> str;
			if (in.eof())
				break;
			model_names.push_back(str);
			int df;
			double logl;
			in >> df >> logl;
			df_vec.push_back(df);
			lh_scores.push_back(logl);
			getline(in, str);
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

string modelTest(Params &params, PhyloTree *in_tree, vector<ModelInfo> &model_info) {
	if (in_tree->isSuperTree()) {
		// select model for each partition
		return (string)"JC";
	}
	int nstates = in_tree->aln->num_states;
	if (nstates != 2 && nstates != 4 && nstates != 20)
		outError("Test of best-fit models only works for Binary/DNA/Protein");
	string fmodel_str = params.out_prefix;
	fmodel_str += ".model";
	string sitelh_file = params.out_prefix;
	sitelh_file += ".sitelh";
	in_tree->params = &params;

	int num_models = (nstates == 2) ? BIN_MODEL_NUM :
					((nstates == 4) ? DNA_MODEL_NUM : AA_MODEL_NUM);
	int model, rate_type;

	string best_model;
	StrVector model_names;
	DoubleVector lh_scores;
	IntVector df_vec;
	/* first check the model file */
	bool ok_model_file = (params.print_site_lh) ? false :
		checkModelFile(fmodel_str, model_names, lh_scores, df_vec);
	ok_model_file &= (model_names.size() == num_models * 4);
	ofstream fmodel;
	if (!ok_model_file) {
		model_names.clear();
		lh_scores.clear();
		df_vec.clear();
		fmodel.open(fmodel_str.c_str());
		if (!fmodel.is_open())
			outError("cannot write to ", fmodel_str);
		fmodel << "Model\tdf\tLnL";
		if (nstates == 2)
			fmodel << "\t0\t1";
		else if (nstates == 4)
			fmodel << "\tA-C\tA-G\tA-T\tC-G\tC-T\tG-T\tA\tC\tG\tT";
		fmodel << "\talpha\tpinv" << endl;
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
	rate_class[2] = new RateGamma(params.num_rate_cats, -1, params.gamma_median,
			NULL);
	rate_class[3] = new RateGammaInvar(params.num_rate_cats, -1,
			params.gamma_median, -1, NULL);
	GTRModel *subst_model = NULL;
	if (nstates == 2)
		subst_model = new ModelBIN("JC2", "", FREQ_UNKNOWN, "", in_tree);
	else if (nstates == 4)
		subst_model = new ModelDNA("JC", "", FREQ_UNKNOWN, "", in_tree);
	else if (nstates == 20)
		subst_model = new ModelProtein("WAG", "", FREQ_UNKNOWN, "", in_tree);
	else if (in_tree->aln->codon_table)
		subst_model = new ModelCodon("GY", "", FREQ_UNKNOWN, "", in_tree);

	assert(subst_model);

	ModelFactory *model_fac = new ModelFactory();

	int ssize = in_tree->aln->getNSite(); // sample size
	if (params.model_test_sample_size)
		ssize = params.model_test_sample_size;
	cout << "Testing " << num_models * 4
			<< ((nstates == 2) ? "binary" : ((nstates == 4) ? " DNA" : " protein"))
			<< " models (sample size: " << ssize << ") ..." << endl;
	cout << "Model         -LnL         df  AIC          AICc         BIC" << endl;

	if (params.print_site_lh) {
		ofstream sitelh_out(sitelh_file.c_str());
		if (!sitelh_out.is_open())
			outError("Cannot write to file ", sitelh_file);
		sitelh_out << num_models*4 << " " << in_tree->getAlnNSite() << endl;
		sitelh_out.close();
	}

	double* AIC_scores = new double[num_models*4];
	double* AICc_scores = new double[num_models*4];
	double* BIC_scores = new double[num_models*4];
	double* LogL_scores = new double[num_models*4];
	int *model_rank = new int[num_models*4];

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
			subst_model->init((nstates == 2) ? bin_model_names[model].c_str() :
							((nstates == 4) ? dna_model_names[model].c_str() : aa_model_names[model].c_str()),
					"", FREQ_UNKNOWN, "");
			subst_model->setTree(tree);
			tree->params = &params;

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
			int df = subst_model->getNDim() + rate_class[rate_type]->getNDim() + tree->branchNum;
			double cur_lh;
			if (!ok_model_file) {
				cur_lh = tree->getModelFactory()->optimizeParameters(false, false);
				model_names.push_back(str);
				// print information to .model file
				fmodel << str << "\t" << df << "\t" << cur_lh;
				if (nstates == 4) {
					int nrates = tree->getModel()->getNumRateEntries();
					double *rate_mat = new double[nrates];
					tree->getModel()->getRateMatrix(rate_mat);
					for (int rate = 0; rate < nrates; rate++)
						fmodel << "\t" << rate_mat[rate];
					delete [] rate_mat;
				}
				if (nstates <= 4) {
					double *freqs = new double[nstates];
					tree->getModel()->getStateFrequency(freqs);
					for (int freq = 0; freq < nstates; freq++)
						fmodel << "\t" << freqs[freq];
					delete [] freqs;
				}
				double alpha = tree->getRate()->getGammaShape();
				fmodel << "\t";
				if (alpha > 0) fmodel << alpha; else fmodel << "NA";
				fmodel << "\t";
				double pinvar = tree->getRate()->getPInvar();
				if (pinvar > 0) fmodel << pinvar << endl; else fmodel << "NA" << endl;
				const char *model_name = (params.print_site_lh) ? str.c_str() : NULL;
				if (params.print_site_lh)
					printSiteLh(sitelh_file.c_str(), tree, NULL, true, model_name);
			} else {
				// sanity check
				if (str != model_names[model * 4 + rate_type] || df != df_vec[model * 4 + rate_type])
					outError("Incorrect model file, please delete it and rerun again: ", fmodel_str);
				cur_lh = lh_scores[model * 4 + rate_type];
			}
			double AIC_score = -2 * cur_lh + 2 * df;
			double AICc_score = AIC_score + 2.0 * df * (df + 1) / (ssize - df - 1);
			double BIC_score = -2 * cur_lh + df * log(ssize);
			LogL_scores[model*4 + rate_type] = cur_lh;
			AIC_scores[model*4 + rate_type] = AIC_score;
			AICc_scores[model*4 + rate_type] = AICc_score;
			BIC_scores[model*4 + rate_type] = BIC_score;
			cout.width(13);
			cout << left << str << " ";
			cout.precision(3);
			cout << fixed;
			cout.width(12);
			cout << -cur_lh << " ";
			cout.width(3);
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
	//cout.unsetf(ios::fixed);
	int model_aic = min_element(AIC_scores, AIC_scores + (num_models*4)) - AIC_scores;
	cout << "Akaike Information Criterion:           " << model_names[model_aic]
			<< endl;
	int model_aicc = min_element(AICc_scores, AICc_scores + (num_models*4)) - AICc_scores;
	cout << "Corrected Akaike Information Criterion: "
			<< model_names[model_aicc] << endl;
	int model_bic = min_element(BIC_scores, BIC_scores + (num_models*4)) - BIC_scores;
	cout << "Bayesian Information Criterion:         " << model_names[model_bic]
			<< endl;

	/* computing model weights */
	double AIC_sum = 0.0, AICc_sum = 0.0, BIC_sum = 0.0;
	if (params.model_test_criterion == MTC_BIC) {
		sort_index(BIC_scores, BIC_scores + (num_models*4), model_rank);
	} else if (params.model_test_criterion == MTC_AIC) {
		sort_index(AIC_scores, AIC_scores + (num_models*4), model_rank);
	} else {
		sort_index(AICc_scores, AICc_scores + (num_models*4), model_rank);
	}

	for (model = 0; model < num_models*4; model++) {
		ModelInfo info;
		info.name = model_names[model_rank[model]];
		info.logl = LogL_scores[model_rank[model]];
		info.AIC_score = AIC_scores[model_rank[model]];
		info.AICc_score = AICc_scores[model_rank[model]];
		info.BIC_score = BIC_scores[model_rank[model]];
		info.AIC_weight = exp(-0.5*(info.AIC_score-AIC_scores[model_aic]));
		info.AICc_weight = exp(-0.5*(info.AICc_score-AICc_scores[model_aicc]));
		info.BIC_weight = exp(-0.5*(info.BIC_score-BIC_scores[model_bic]));
		info.AIC_conf = false;
		info.AICc_conf = false;
		info.BIC_conf = false;
		model_info.push_back(info);
		AIC_sum += info.AIC_weight;
		AICc_sum += info.AICc_weight;
		BIC_sum += info.BIC_weight;
	}

	vector<ModelInfo>::iterator it;
	for (it = model_info.begin(); it != model_info.end(); it++) {
		it->AIC_weight /= AIC_sum;
		it->AICc_weight /= AICc_sum;
		it->BIC_weight /= BIC_sum;
	}

	/* compute confidence set for BIC */
    AIC_sum = 0.0;
    AICc_sum = 0.0;
    BIC_sum = 0.0;
	for (it = model_info.begin(), model = 0; it != model_info.end(); it++, model++) {
		BIC_scores[model] = it->BIC_score;
		AIC_scores[model] = it->AIC_score;
		AICc_scores[model] = it->AICc_score;
	}
	sort_index(BIC_scores, BIC_scores+(num_models*4), model_rank);
	for (model = 0; model < num_models*4; model++) {
		model_info[model_rank[model]].BIC_conf = true;
		BIC_sum += model_info[model_rank[model]].BIC_weight;
		if (BIC_sum > 0.95) break;
	}
	/* compute confidence set for AIC */
	sort_index(AIC_scores, AIC_scores+(num_models*4), model_rank);
	for (model = 0; model < num_models*4; model++) {
		model_info[model_rank[model]].AIC_conf = true;
		AIC_sum += model_info[model_rank[model]].AIC_weight;
		if (AIC_sum > 0.95) break;
	}

	/* compute confidence set for AICc */
	sort_index(AICc_scores, AICc_scores+(num_models*4), model_rank);
	for (model = 0; model < num_models*4; model++) {
		model_info[model_rank[model]].AICc_conf = true;
		AICc_sum += model_info[model_rank[model]].AICc_weight;
		if (AICc_sum > 0.95) break;
	}

	delete [] model_rank;
	delete [] LogL_scores;
	delete [] BIC_scores;
	delete [] AICc_scores;
	delete [] AIC_scores;
	switch (params.model_test_criterion) {
	case MTC_AIC:
		best_model = model_names[model_aic];
		break;
	case MTC_AICC:
		best_model = model_names[model_aicc];
		break;
	case MTC_BIC:
		best_model = model_names[model_bic];
		break;
	}
	delete model_fac;
	delete subst_model;
	for (rate_type = 3; rate_type >= 0; rate_type--)
		delete rate_class[rate_type];
	delete tree_hetero;
	delete tree_homo;

	if (!ok_model_file)
		fmodel.close();
	cout << "Best-fit model: " << best_model << endl;
	if (params.print_site_lh)
		cout << "Site log-likelihoods per model printed to " << sitelh_file << endl;
	return best_model;
}

int countDistinctTrees(const char *filename, bool rooted, IQTree *tree, IntVector &distinct_ids) {
	StringIntMap treels;
	try {
		ifstream in;
		in.exceptions(ios::failbit | ios::badbit);
		in.open(filename);
		// remove the failbit
		in.exceptions(ios::badbit);
		int tree_id;
		for (tree_id = 0; !in.eof(); tree_id++) {
			tree->freeNode();
			tree->readTree(in, rooted);
			tree->setAlignment(tree->aln);
			tree->setRootNode((char*)tree->aln->getSeqName(0).c_str());
		    StringIntMap::iterator it = treels.end();
		    ostringstream ostr;
		    tree->printTree(ostr, WT_TAXON_ID | WT_SORT_TAXA);
		    it = treels.find(ostr.str());
		    if (it != treels.end()) { // already in treels
		    	distinct_ids.push_back(it->second);
		    } else {
		    	distinct_ids.push_back(-1);
		    	treels[ostr.str()] = tree_id;
		    }
			char ch;
			in.exceptions(ios::goodbit);
			(in) >> ch;
			if (in.eof()) break;
			in.unget();
			in.exceptions(ios::failbit | ios::badbit);

		}
		in.close();
	} catch (ios::failure) {
		outError("Cannot read file ", filename);
	}
	return treels.size();
}

//const double TOL_RELL_SCORE = 0.01;

void evaluateTrees(Params &params, IQTree *tree, vector<TreeInfo> &info, IntVector &distinct_ids)
{
	if (!params.treeset_file)
		return;
	cout << endl;
	//MTreeSet trees(params.treeset_file, params.is_rooted, params.tree_burnin, params.tree_max_count);
	cout << "Reading trees in " << params.treeset_file << " ..." << endl;
	int ntrees = countDistinctTrees(params.treeset_file, params.is_rooted, tree, distinct_ids);
	if (ntrees < distinct_ids.size()) {
		cout << "WARNING: " << distinct_ids.size() << " trees detected but only " << ntrees << " distinct trees will be evaluated" << endl;
	} else {
		cout << ntrees << " distinct trees detected" << endl;
	}
	if (ntrees == 0) return;
	ifstream in(params.treeset_file);

	//if (trees.size() == 1) return;
	string tree_file = params.treeset_file;
	tree_file += ".trees";
	ofstream treeout;
	//if (!params.fixed_branch_length) {
		treeout.open(tree_file.c_str());
	//}
	string score_file = params.treeset_file;
	score_file += ".treelh";
	ofstream scoreout;
	if (params.print_tree_lh)
		scoreout.open(score_file.c_str());
	string site_lh_file = params.treeset_file;
	site_lh_file += ".sitelh";
	if (params.print_site_lh) {
		ofstream site_lh_out(site_lh_file.c_str());
		site_lh_out << ntrees << " " << tree->getAlnNSite() << endl;
		site_lh_out.close();
	}

	double time_start = getCPUTime();

	int *boot_samples = NULL;
	int boot;
	//double *saved_tree_lhs = NULL;
	double *tree_lhs = NULL;
	double *pattern_lh = NULL;
	double *pattern_lhs = NULL;
	double *orig_tree_lh = NULL;
	double *max_lh = NULL;
	double *lhdiff_weights = NULL;
	int nptn = tree->getAlnNPattern();
	if (params.topotest_replicates && ntrees > 1) {
		size_t mem_size = (size_t)params.topotest_replicates*nptn*sizeof(int) +
				ntrees*params.topotest_replicates*sizeof(double) +
				(nptn + ntrees*3 + params.topotest_replicates*2)*sizeof(double) +
				ntrees*sizeof(TreeInfo) +
				params.do_weighted_test*(ntrees * nptn * sizeof(double) + ntrees*ntrees*sizeof(double));
		cout << "Note: " << ((double)mem_size/1024)/1024 << " MB of RAM required!" << endl;
		if (mem_size > getMemorySize()-100000)
			outWarning("The required memory does not fit in RAM!");
		cout << "Creating " << params.topotest_replicates << " bootstrap replicates..." << endl;
		if (!(boot_samples = new int [params.topotest_replicates*nptn]))
			outError(ERR_NO_MEMORY);
		for (boot = 0; boot < params.topotest_replicates; boot++)
			tree->aln->createBootstrapAlignment(boot_samples + (boot*nptn), params.bootstrap_spec);
		//if (!(saved_tree_lhs = new double [ntrees * params.topotest_replicates]))
		//	outError(ERR_NO_MEMORY);
		if (!(tree_lhs = new double [ntrees * params.topotest_replicates]))
			outError(ERR_NO_MEMORY);
		if (params.do_weighted_test) {
			if (!(lhdiff_weights = new double [ntrees * ntrees]))
				outError(ERR_NO_MEMORY);
			if (!(pattern_lhs = new double[ntrees* nptn]))
				outError(ERR_NO_MEMORY);
		}
		if (!(pattern_lh = new double[nptn]))
			outError(ERR_NO_MEMORY);
		if (!(orig_tree_lh = new double[ntrees]))
			outError(ERR_NO_MEMORY);
		if (!(max_lh = new double[params.topotest_replicates]))
			outError(ERR_NO_MEMORY);
	}
	int tree_index, tid, tid2;
	info.resize(ntrees);
	//for (MTreeSet::iterator it = trees.begin(); it != trees.end(); it++, tree_index++) {
	for (tree_index = 0, tid = 0; tree_index < distinct_ids.size(); tree_index++) {

		cout << "Tree " << tree_index + 1;
		if (distinct_ids[tree_index] >= 0) {
			cout << " / identical to tree " << distinct_ids[tree_index]+1 << endl;
			// ignore tree
			char ch;
			do {
				in >> ch;
			} while (!in.eof() && ch != ';');
			continue;
		}
		tree->freeNode();
		tree->readTree(in, params.is_rooted);
		tree->setAlignment(tree->aln);
		tree->initializeAllPartialLh();
		tree->fixNegativeBranch(false);
		if (tree->isSuperTree())
			((PhyloSuperTree*) tree)->mapTrees();
		if (!params.fixed_branch_length) {
			tree->curScore = tree->optimizeAllBranches(100, 0.001);
		} else {
			tree->curScore = tree->computeLikelihood();
		}
		treeout << "[ tree " << tree_index+1 << " lh=" << tree->curScore << " ]";
		tree->printTree(treeout);
		treeout << endl;
		if (params.print_tree_lh)
			scoreout << tree->curScore << endl;

		cout << " / LogL: " << tree->curScore << endl;

		if (pattern_lh) {
			tree->computePatternLikelihood(pattern_lh, &(tree->curScore));
			if (params.do_weighted_test)
				memcpy(pattern_lhs + tid*nptn, pattern_lh, nptn*sizeof(double));
		}
		if (params.print_site_lh) {
			string tree_name = "Tree" + convertIntToString(tree_index+1);
			printSiteLh(site_lh_file.c_str(), tree, pattern_lh, true, tree_name.c_str());
		}
		info[tid].logl = tree->curScore;

		if (!params.topotest_replicates || ntrees <= 1) {
			tid++;
			continue;
		}
		// now compute RELL scores
		orig_tree_lh[tid] = tree->curScore;
		double *tree_lhs_offset = tree_lhs + (tid*params.topotest_replicates);
		for (boot = 0; boot < params.topotest_replicates; boot++) {
			double lh = 0.0;
			int *this_boot_sample = boot_samples + (boot*nptn);
			for (int ptn = 0; ptn < nptn; ptn++)
				lh += pattern_lh[ptn] * this_boot_sample[ptn];
			tree_lhs_offset[boot] = lh;
		}
		tid++;
	}

	assert(tid == ntrees);

	if (params.topotest_replicates && ntrees > 1) {
		double *tree_probs = new double[ntrees];
		memset(tree_probs, 0, ntrees*sizeof(double));
		int *tree_ranks = new int[ntrees];

		/* perform RELL BP method */
		cout << "Performing RELL test..." << endl;
		int *maxtid = new int[params.topotest_replicates];
		double *maxL = new double[params.topotest_replicates];
		int *maxcount = new int[params.topotest_replicates];
		memset(maxtid, 0, params.topotest_replicates*sizeof(int));
		memcpy(maxL, tree_lhs, params.topotest_replicates*sizeof(double));
		for (boot = 0; boot < params.topotest_replicates; boot++)
			maxcount[boot] = 1;
		for (tid = 1; tid < ntrees; tid++) {
			double *tree_lhs_offset = tree_lhs + (tid * params.topotest_replicates);
			for (boot = 0; boot < params.topotest_replicates; boot++)
				if (tree_lhs_offset[boot] > maxL[boot] + params.ufboot_epsilon) {
					maxL[boot] = tree_lhs_offset[boot];
					maxtid[boot] = tid;
					maxcount[boot] = 1;
				} else if (tree_lhs_offset[boot] > maxL[boot] - params.ufboot_epsilon &&
						random_double() <= 1.0/(maxcount[boot]+1)) {
					maxL[boot] = max(maxL[boot],tree_lhs_offset[boot]);
					maxtid[boot] = tid;
					maxcount[boot]++;
				}
		}
		for (boot = 0; boot < params.topotest_replicates; boot++)
			tree_probs[maxtid[boot]] += 1.0;
		for (tid = 0; tid < ntrees; tid++) {
			tree_probs[tid] /= params.topotest_replicates;
			info[tid].rell_confident = false;
			info[tid].rell_bp = tree_probs[tid];
		}
		sort_index(tree_probs, tree_probs + ntrees, tree_ranks);
		double prob_sum = 0.0;
		// obtain the confidence set
		for (tid = ntrees-1; tid >= 0; tid--) {
			info[tree_ranks[tid]].rell_confident = true;
			prob_sum += tree_probs[tree_ranks[tid]];
			if (prob_sum > 0.95) break;
		}

		// sanity check
		for (tid = 0, prob_sum = 0.0; tid < ntrees; tid++)
			prob_sum += tree_probs[tid];
		if (fabs(prob_sum-1.0) > 0.01)
			outError("Internal error: Wrong ", __func__);

		delete [] maxcount;
		delete [] maxL;
		delete [] maxtid;

		/* now do the SH test */
		cout << "Performing KH and SH test..." << endl;
		// SH centering step
		for (boot = 0; boot < params.topotest_replicates; boot++)
			max_lh[boot] = -DBL_MAX;
		double *avg_lh = new double[ntrees];
		for (tid = 0; tid < ntrees; tid++) {
			avg_lh[tid] = 0.0;
			double *tree_lhs_offset = tree_lhs + (tid * params.topotest_replicates);
			for (boot = 0; boot < params.topotest_replicates; boot++)
				avg_lh[tid] += tree_lhs_offset[boot];
			avg_lh[tid] /= params.topotest_replicates;
			for (boot = 0; boot < params.topotest_replicates; boot++) {
				max_lh[boot] = max(max_lh[boot], tree_lhs_offset[boot] - avg_lh[tid]);
			}
		}

		double orig_max_lh = orig_tree_lh[0];
		int orig_max_id = 0;
		double orig_2ndmax_lh = -DBL_MAX;
		int orig_2ndmax_id = -1;
		// find the max tree ID
		for (tid = 1; tid < ntrees; tid++)
			if (orig_max_lh < orig_tree_lh[tid]) {
				orig_max_lh = orig_tree_lh[tid];
				orig_max_id = tid;
			}
		// find the 2nd max tree ID
		for (tid = 0; tid < ntrees; tid++)
			if (tid != orig_max_id && orig_2ndmax_lh < orig_tree_lh[tid]) {
				orig_2ndmax_lh = orig_tree_lh[tid];
				orig_2ndmax_id = tid;
			}


		// SH compute p-value
		for (tid = 0; tid < ntrees; tid++) {
			double *tree_lhs_offset = tree_lhs + (tid * params.topotest_replicates);
			// SH compute original deviation from max_lh
			info[tid].kh_pvalue = 0.0;
			info[tid].sh_pvalue = 0.0;
			int max_id = (tid != orig_max_id) ? orig_max_id : orig_2ndmax_id;
			double orig_diff = orig_tree_lh[max_id] - orig_tree_lh[tid] - avg_lh[tid];
			double *max_kh = tree_lhs + (max_id * params.topotest_replicates);
			for (boot = 0; boot < params.topotest_replicates; boot++) {
				if (max_lh[boot] - tree_lhs_offset[boot] > orig_diff)
					info[tid].sh_pvalue += 1.0;
				//double max_kh_here = max(max_kh[boot]-avg_lh[max_id], tree_lhs_offset[boot]-avg_lh[tid]);
				double max_kh_here = (max_kh[boot]-avg_lh[max_id]);
				if (max_kh_here - tree_lhs_offset[boot] > orig_diff)
					info[tid].kh_pvalue += 1.0;
			}
			info[tid].sh_pvalue /= params.topotest_replicates;
			info[tid].kh_pvalue /= params.topotest_replicates;
		}

		if (params.do_weighted_test) {

			cout << "Computing pairwise logl difference variance ..." << endl;
			/* computing lhdiff_weights as 1/sqrt(lhdiff_variance) */
			for (tid = 0; tid < ntrees; tid++) {
				double *pattern_lh1 = pattern_lhs + (tid * nptn);
				lhdiff_weights[tid*ntrees+tid] = 0.0;
				for (tid2 = tid+1; tid2 < ntrees; tid2++) {
					double lhdiff_variance = tree->computeLogLDiffVariance(pattern_lh1, pattern_lhs + (tid2*nptn));
					lhdiff_weights[tid*ntrees+tid2] = 1.0/sqrt(lhdiff_variance);
					lhdiff_weights[tid2*ntrees+tid] = lhdiff_weights[tid*ntrees+tid2];
				}
			}

			// Weighted KH and SH test
			cout << "Performing WKH and WSH test..." << endl;
			for (tid = 0; tid < ntrees; tid++) {
				double *tree_lhs_offset = tree_lhs + (tid * params.topotest_replicates);
				info[tid].wkh_pvalue = 0.0;
				info[tid].wsh_pvalue = 0.0;
				double worig_diff = -DBL_MAX;
				int max_id = -1;
				for (tid2 = 0; tid2 < ntrees; tid2++)
					if (tid2 != tid) {
						double wdiff = (orig_tree_lh[tid2] - orig_tree_lh[tid])*lhdiff_weights[tid*ntrees+tid2];
						if (wdiff > worig_diff) {
							worig_diff = wdiff;
							max_id = tid2;
						}
					}
				for (boot = 0; boot < params.topotest_replicates; boot++) {
					double wmax_diff = -DBL_MAX;
					for (tid2 = 0; tid2 < ntrees; tid2++)
						if (tid2 != tid)
							wmax_diff = max(wmax_diff,
									(tree_lhs[tid2*params.topotest_replicates+boot] - avg_lh[tid2] -
									tree_lhs_offset[boot] + avg_lh[tid]) * lhdiff_weights[tid*ntrees+tid2]);
					if (wmax_diff > worig_diff)
						info[tid].wsh_pvalue += 1.0;
					wmax_diff = (tree_lhs[max_id*params.topotest_replicates+boot] - avg_lh[max_id] -
							tree_lhs_offset[boot] + avg_lh[tid]);
					if (wmax_diff >  orig_tree_lh[max_id] - orig_tree_lh[tid])
						info[tid].wkh_pvalue += 1.0;
				}
				info[tid].wsh_pvalue /= params.topotest_replicates;
				info[tid].wkh_pvalue /= params.topotest_replicates;
			}
		}
		/* now to ELW - Expected Likelihood Weight method */
		cout << "Performing ELW test..." << endl;

		for (boot = 0; boot < params.topotest_replicates; boot++)
			max_lh[boot] = -DBL_MAX;
		for (tid = 0; tid < ntrees; tid++) {
			double *tree_lhs_offset = tree_lhs + (tid * params.topotest_replicates);
			for (boot = 0; boot < params.topotest_replicates; boot++)
				max_lh[boot] = max(max_lh[boot], tree_lhs_offset[boot]);
		}
		double *sumL = new double[params.topotest_replicates];
		memset(sumL, 0, sizeof(double) * params.topotest_replicates);
		for (tid = 0; tid < ntrees; tid++) {
			double *tree_lhs_offset = tree_lhs + (tid * params.topotest_replicates);
			for (boot = 0; boot < params.topotest_replicates; boot++) {
				tree_lhs_offset[boot] = exp(tree_lhs_offset[boot] - max_lh[boot]);
				sumL[boot] += tree_lhs_offset[boot];
			}
		}
		for (tid = 0; tid < ntrees; tid++) {
			double *tree_lhs_offset = tree_lhs + (tid * params.topotest_replicates);
			tree_probs[tid] = 0.0;
			for (boot = 0; boot < params.topotest_replicates; boot++) {
				tree_probs[tid] += (tree_lhs_offset[boot] / sumL[boot]);
			}
			tree_probs[tid] /= params.topotest_replicates;
			info[tid].elw_confident = false;
			info[tid].elw_value = tree_probs[tid];
		}

		sort_index(tree_probs, tree_probs + ntrees, tree_ranks);
		prob_sum = 0.0;
		// obtain the confidence set
		for (tid = ntrees-1; tid >= 0; tid--) {
			info[tree_ranks[tid]].elw_confident = true;
			prob_sum += tree_probs[tree_ranks[tid]];
			if (prob_sum > 0.95) break;
		}

		// sanity check
		for (tid = 0, prob_sum = 0.0; tid < ntrees; tid++)
			prob_sum += tree_probs[tid];
		if (fabs(prob_sum-1.0) > 0.01)
			outError("Internal error: Wrong ", __func__);
		delete [] sumL;

		delete [] tree_ranks;
		delete [] tree_probs;

	}
	if (max_lh)
		delete [] max_lh;
	if (orig_tree_lh)
		delete [] orig_tree_lh;
	if (pattern_lh)
		delete [] pattern_lh;
	if (pattern_lhs)
		delete [] pattern_lhs;
	if (lhdiff_weights)
		delete [] lhdiff_weights;
	if (tree_lhs)
		delete [] tree_lhs;
	//if (saved_tree_lhs)
	//	delete [] saved_tree_lhs;
	if (boot_samples)
		delete [] boot_samples;

	if (params.print_tree_lh) {
		scoreout.close();
	}

	treeout.close();
	in.close();

	cout << "Time for evaluating all trees: " << getCPUTime() - time_start << " sec." << endl;

}


void evaluateTrees(Params &params, IQTree *tree) {
	vector<TreeInfo> info;
	IntVector distinct_ids;
	evaluateTrees(params, tree, info, distinct_ids);
}

