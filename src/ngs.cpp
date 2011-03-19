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

/*
	collection of classes for Next-generation sequencing 
*/

#include "ngs.h"
#include "modeltest_wrapper.h"

/****************************************************************************
        NGSAlignment
 ****************************************************************************/

NGSAlignment::NGSAlignment(const char *filename) : AlignmentPairwise() {
	readFritzFile(filename);
}

NGSAlignment::NGSAlignment(int nstate, int ncat, int *freq) : AlignmentPairwise() {
	num_states = nstate;
	ncategory = ncat;
	int total_size = ncategory*num_states*num_states;
	pair_freq = new int[total_size];
	memcpy(pair_freq, freq, total_size * sizeof(int));
}

void NGSAlignment::readFritzFile(const char *filename) {
	cout << "Reading Fritz file " << filename << " ..." << endl;
	try {
		ifstream in;
		in.exceptions(ios::failbit | ios::badbit);
		in.open(filename);
		in.clear();
		int i, total_size;
		string tmp;
		in >> tmp;
		ncategory = convert_int(tmp.c_str());
		if (ncategory < 1) throw "Wrong number of positions";
		in >> tmp;
		num_states = convert_int(tmp.c_str());
		total_size = ncategory*num_states*num_states;
		if (num_states < 1) throw "Wrong number of states";
		pair_freq = new int[total_size];
		for (i=0; i < total_size; i++) {
			in >> tmp;
			int count = convert_int(tmp.c_str());
			if (count < 0) throw "Wrong count";
			pair_freq[i] = count;
		}
		// set the failbit again
		in.exceptions(ios::failbit | ios::badbit);
		in.close();
	} catch (const char *str) {
		outError(str);
	} catch (string str) {
		outError(str);
	} catch(ios::failure) {
		outError(ERR_READ_INPUT);
	}

	cout << ncategory << " matrices of size " << num_states << endl;
}

void NGSAlignment::computeStateFreq (double *stateFrqArr) {
	int cat, i, j, id = 0;
	int state_count[num_states];
	memset(state_count, 0, sizeof(int)*num_states);
	for (cat = 0, id = 0; cat < ncategory; cat++) {
		for (i = 0; i < num_states; i++)
			for (j = 0; j < num_states; j++, id++) {
				state_count[i] += pair_freq[id];
				state_count[j] += pair_freq[id];
			}
	}

	int sum_count = 0;
	for (i = 0; i < num_states; i++) sum_count += state_count[i];
	if (sum_count == 0) throw "Empty data observed";
	for (i = 0; i < num_states; i++) stateFrqArr[i] = double(state_count[i]) / sum_count;
	if (verbose_mode >= VB_MIN) {
		cout << "Empirical state frequencies: ";
		for (i = 0; i < num_states; i++) 
			cout << stateFrqArr[i] << " ";
		cout << endl;
	}
}

void NGSAlignment::computeSumPairFreq (int *sum_pair_freq) {
	int cat, id, i, j;
	memset(sum_pair_freq, 0, sizeof(int)*num_states*num_states);
	for (cat = 0, id = 0; cat < ncategory; cat++) {
		for (i = 0; i < num_states; i++)
			for (j = 0; j < num_states; j++, id++) {
				sum_pair_freq[i*num_states+j] += pair_freq[id];
			}
	}
}

void NGSAlignment::computeEmpiricalRate (double *rates) {
	int i, j, k, cat, id;
	assert(rates);
	double **pair_rates = (double**) new double[num_states];
	for (i = 0; i < num_states; i++) {
		pair_rates[i] = new double[num_states];
		memset(pair_rates[i], 0, sizeof(double)*num_states);
	}

	for (cat = 0, id = 0; cat < ncategory; cat++) {
		for (i = 0; i < num_states; i++)
			for (j = 0; j < num_states; j++, id++) {
				pair_rates[i][j] += pair_freq[id];
			}
	}

	k = 0;
	double last_rate = pair_rates[num_states-2][num_states-1] + pair_rates[num_states-1][num_states-2];
	if (last_rate == 0.0) throw "Last rate entry is ZERO";
	for (i = 0; i < num_states-1; i++)
		for (j = i+1; j < num_states; j++)
			rates[k++] = (pair_rates[i][j] + pair_rates[j][i]) / last_rate;
	if (verbose_mode >= VB_MIN) {
		cout << "Empirical rates: ";
		for (k = 0; k < num_states*(num_states-1)/2; k++)
			cout << rates[k] << " ";
		cout << endl;
	}
	for (i = num_states-1; i >= 0; i--) {
		delete [] pair_rates[i];
	}
	delete [] pair_rates;
}

double NGSAlignment::computeEmpiricalDist(int cat) {
	int i;
	int trans_size = num_states*num_states;
	int *pair_pos = pair_freq + (cat*trans_size);
	int match_pos = 0, total_pos = 0;
	for (i = 0; i < num_states; i++) 
		match_pos += pair_pos[i*num_states+i];
	for (i = 0; i < trans_size; i++) 
		total_pos += pair_pos[i];
	if (total_pos == 0) total_pos = 1;
	return (double)(total_pos - match_pos) / total_pos;
}


double NGSAlignment::computeFunctionCat(int cat, double value) {
	int trans_size = num_states*num_states;
	double lh = 0.0;
	double trans_mat[trans_size];
	int i;

	tree->getModelFactory()->computeTransMatrix(value, trans_mat);
	int *pair_pos = pair_freq + cat*trans_size;

	for (i = 0; i < trans_size; i++) if (pair_pos[i]) {
		if (trans_mat[i] <= 0) throw "Negative transition probability";
		lh -= pair_pos[i] * log(trans_mat[i]);
	}
	return lh;
}


double NGSAlignment::computeFuncDervCat(int cat, double value, double &df, double &ddf) {
	int trans_size = num_states*num_states;
	double lh = 0.0;
	df = 0.0; ddf = 0.0;
	int i;
	double derv1 = 0.0, derv2 = 0.0;
	double trans_mat[trans_size], trans_derv1[trans_size], trans_derv2[trans_size];
	

	tree->getModelFactory()->computeTransDerv(value, trans_mat, trans_derv1, trans_derv2);
	int *pair_pos = pair_freq + cat*trans_size;
	for (i = 0; i < trans_size; i++) if (pair_pos[i] > 0) {
		if (trans_mat[i] <= 0) throw "Negative transition probability";
		double d1 = trans_derv1[i] / trans_mat[i];
		derv1 += pair_pos[i] * d1;
		derv2 += pair_pos[i] * (trans_derv2[i]/trans_mat[i] - d1 * d1);
		lh -= pair_pos[i] * log(trans_mat[i]);
	}
	//df -= derv1 * rate_val;
	//ddf -= derv2 * rate_val * rate_val;
	df -= derv1;
	ddf -= derv2;
	return lh;
}

/****************************************************************************
        NGSRate
 ****************************************************************************/
NGSRate::NGSRate(PhyloTree *tree) {
	phylo_tree = tree;
	ncategory = ((NGSAlignment*)tree->aln)->ncategory;
	rates = new double[ncategory];
	int i;
	for (i = 0; i < ncategory; i++) {
		rates[i] = ((NGSAlignment*)tree->aln)->computeEmpiricalDist(i);
		if (rates[i] < 1e-6) rates[i] = 1e-6;
	}
	
	name = "+F";
	name += convertIntToString(ncategory);
	full_name = name;
	is_categorized = true;

}

double NGSRate::optimizeParameters() {
	int cat;
	double negative_lh;
	for (cat = 0; cat < ncategory; cat++) {
		optimizing_cat = cat;
		if (phylo_tree->optimize_by_newton) 
			rates[cat] = minimizeNewtonSafeMode(1e-6, rates[cat], 10.0, 1e-6, negative_lh);
		else
			rates[cat] = minimizeOneDimenSafeMode(1e-6, rates[cat], 10.0, 1e-6, &negative_lh);
	}
	return phylo_tree->computeLikelihood();
}


double NGSRate::computeFunction(double value) {
	return ((NGSAlignment*)phylo_tree->aln)->computeFunctionCat(optimizing_cat, value);
}
double NGSRate::computeFuncDerv(double value, double &df, double &ddf) {
	return ((NGSAlignment*)phylo_tree->aln)->computeFuncDervCat(optimizing_cat, value, df, ddf);
}

void NGSRate::writeInfo(ostream &out) {
}

/****************************************************************************
        NGSTree
 ****************************************************************************/

NGSTree::NGSTree(Params &params, NGSAlignment *alignment) {
	aln = alignment;
    model = NULL;
    site_rate = NULL;
    model_factory = NULL;
    optimize_by_newton = params.optimize_by_newton;
    //tree.sse = params.SSE;
    sse = false;
}

double NGSTree::computeLikelihood(double *pattern_lh) {
	return -((NGSAlignment*)aln)->computeFunction(1.0);
}

double NGSTree::optimizeAllBranches(int iterations) {
	return computeLikelihood();
}

/****************************************************************************
        main function
 ****************************************************************************/

void reportNGSAnalysis(const char *file_name, Params &params, NGSAlignment &aln, NGSTree &tree, 
	DoubleMatrix &rate_info, StrVector &rate_name) {
	ofstream out(file_name);
	out.setf(ios::fixed,ios::floatfield);

	int i, j, k;


	double rate_param[aln.num_states * aln.num_states];
	double rate_matrix[aln.num_states * aln.num_states];

	out << "Input file: " << params.ngs_file << endl;
	out << "Model of evolution: " << tree.getModel()->name << endl << endl;

	out << "Substitution process assuming one homogeneous model among all positions:" << endl;

	out << "Rate parameters: " << endl;

	tree.getModel()->getRateMatrix(rate_param);

	if (tree.getModel()->name == "UNREST") {
		for (i = 0, k=0; i < aln.num_states; i++)
			for (j = 0; j < aln.num_states; j++)
				if (i != j)
					rate_matrix[i*aln.num_states+j] = rate_param[k++];
	} else {
		for (i = 0, k=0; i < aln.num_states-1; i++)
			for (j = i+1; j < aln.num_states; j++, k++)
				rate_matrix[i*aln.num_states+j] = rate_matrix[j*aln.num_states+i] = rate_param[k];
	}

	for (i = 0; i < aln.num_states; i++) {
		for (j = 0; j < aln.num_states; j++) {
			if (j > 0) out << " \t";
			if (j != i) out << rate_matrix[i*aln.num_states+j]; else out << "-";
		}
		out << endl;
	}
	out << endl;
	out << "State frequencies: " << endl;

	double state_freq[aln.num_states];
	tree.getModel()->getStateFrequency(state_freq);

	for (i = 0; i < aln.num_states; i++) out << state_freq[i] << " \t";
	out << endl << endl;

	out << "Q matrix can be obtained by multiplying rate parameters with state frequencies" << endl << endl;

	out << "Log-likelihood: " << tree.computeLikelihood() << endl << endl;

	out << "Inferred posisiton-specific rates under one model or position-specific model: " << endl;

	out << "Position\tOne_model_rate";
	for (StrVector::iterator it = rate_name.begin(); it != rate_name.end(); it++)
		out << "\t" << (*it);
	out << endl;
	for (i = 0; i < aln.ncategory; i++) {
		out << i+1 << '\t' << tree.getRate()->getRate(i);
		DoubleVector *rate_vec = &rate_info[i];
		for (DoubleVector::iterator dit = rate_vec->begin(); dit != rate_vec->end(); dit++)
			out << "\t" << *dit;
		out << endl;
	}
	out.close();
	cout << endl << "Results written to: " << file_name << endl << endl;
}

bool checkFreq(int *pair_freq, int n) {
	int i, count = 0;
	for (i=0; i < n*n; i++)
		if (pair_freq[i] != 0) count++;
	if (count <= n) return false;
	return true;
}

void testSingleRateModel(Params &params, NGSAlignment &aln, NGSTree &tree, string model, 
	int *freq, DoubleVector &rate_info, StrVector &rate_name, bool write_info) {

	char model_name[20];
	NGSAlignment sum_aln(aln.num_states, 1, freq);

	NGSTree sum_tree(params, &sum_aln);
	sum_aln.tree = &sum_tree;

	if (model == "") 
		sprintf(model_name, "GTR+F1");
	else
		sprintf(model_name, "%s+F1", model.c_str());
	try {
		params.model_name = model_name;
		sum_tree.setModelFactory(new ModelFactory(params, &sum_tree));
		sum_tree.setModel(sum_tree.getModelFactory()->model);
		sum_tree.setRate(sum_tree.getModelFactory()->site_rate);
    	double bestTreeScore = sum_tree.getModelFactory()->optimizeParameters(false, write_info);
		cout << "LogL: " << bestTreeScore;
		cout << " / Rate: " << sum_tree.getRate()->getRate(0) << endl;
    } catch (...) {
		cout << "Skipped due to sparse matrix" << endl;
		//rate_info.push_back(MIN_SITE_RATE);
		rate_info.insert(rate_info.end(), rate_name.size(), MIN_SITE_RATE);
		return;
    }
    //return sum_tree.getRate()->getRate(0);
	rate_info.push_back(sum_tree.getRate()->getRate(0));

    double rate_mat[aln.num_states*aln.num_states];
    memset(rate_mat, 0, aln.num_states*aln.num_states*sizeof(double));
    sum_tree.getModel()->getRateMatrix(rate_mat);
    rate_info.insert(rate_info.end(), rate_mat, rate_mat+sum_tree.getModel()->getNumRateEntries());

	if (tree.getModel()->isReversible()) {
		sum_tree.getModel()->getStateFrequency(rate_mat);
		rate_info.insert(rate_info.end(), rate_mat, rate_mat+aln.num_states);
    }
}


/*

void testSingleRateModel(Params &params, NGSAlignment &aln, NGSTree &tree, string model, int *sum_freq) {
	char model_name[20];

	NGSAlignment sum_aln(aln.num_states, 1, sum_freq);

	NGSTree sum_tree(params, &sum_aln);
	sum_aln.tree = &sum_tree;

	if (model == "") 
		sprintf(model_name, "GTR+F1");
	else
		sprintf(model_name, "%s+F1", model.c_str());
	params.model_name = model_name;
    sum_tree.setModelFactory(new ModelFactory(params, &sum_tree));
    sum_tree.setModel(sum_tree.getModelFactory()->model);
    sum_tree.setRate(sum_tree.getModelFactory()->site_rate);

    double bestTreeScore = sum_tree.getModelFactory()->optimizeParameters(false, false);
    cout << "Log-likelihood of null model: " << bestTreeScore << endl;
    cout << "Rate (or distance) of null model: " << sum_tree.getRate()->getRate(0) << endl;
    double lh_diff = 2*(tree.computeLikelihood() - bestTreeScore);
    cout << "2(lnL1 - lnL0) = " << lh_diff << endl;
    cout << "p-value (chi-square test, df = " << aln.ncategory-1 << "): " << computePValueChiSquare(lh_diff, aln.ncategory-1) << endl;

	string out_file = params.out_prefix;
	out_file += ".ngs_e";
	DoubleVector tmp;
	reportNGSAnalysis(out_file.c_str(), params, sum_aln, sum_tree, tmp);

}*/


void runNGSAnalysis(Params &params) {
	char model_name[20];
	// read input file, initialize NGSAlignment
	NGSAlignment aln(params.ngs_file);
	cout.setf(ios::fixed,ios::floatfield);

	params.freq_type = FREQ_ESTIMATE;

	// initialize NGSTree
	NGSTree tree(params, &aln);
	aln.tree = &tree;

	// initialize Model 
	string original_model = params.model_name;
	if (params.model_name == "") 
		sprintf(model_name, "GTR+F%d", aln.ncategory);
	else
		sprintf(model_name, "%s+F%d", params.model_name.c_str(), aln.ncategory);
	params.model_name = model_name;
    tree.setModelFactory(new ModelFactory(params, &tree));
    tree.setModel(tree.getModelFactory()->model);
    tree.setRate(tree.getModelFactory()->site_rate);

    int model_df = tree.getModel()->getNDim() + tree.getRate()->getNDim();
    cout << endl;
    cout << "Model of evolution: " << tree.getModelName() << " (" << model_df << " free parameters)" << endl;
    cout << endl;

	// optimize model parameters and rate scaling factors
    cout << "Optimizing model parameters" << endl;
    double bestTreeScore = tree.getModelFactory()->optimizeParameters(false, true);
    cout << "Log-likelihood: " << bestTreeScore << endl;


	DoubleMatrix part_rate(aln.ncategory);
	StrVector rate_name;


	int i, j;

	rate_name.push_back("Varying_model_rate");

	if (tree.getModel()->isReversible()) {
		for (i = 0; i < aln.num_states-1; i++) 
			for (j = i+1; j < aln.num_states; j++) {
				stringstream x;
				x << aln.convertStateBack(i) << "->" << aln.convertStateBack(j);
				rate_name.push_back(x.str());
			}
		for (i = 0; i < aln.num_states; i++) {
			stringstream x;
			x << aln.convertStateBack(i);
			rate_name.push_back(x.str());
		}
	} else {
		for (i = 0; i < aln.num_states; i++) 
			for (j = 0; j < aln.num_states; j++) if (j != i) {
				stringstream x;
				x << aln.convertStateBack(i) << "->" << aln.convertStateBack(j);
				rate_name.push_back(x.str());
			}
	}


	VerboseMode vb_saved = verbose_mode;
	verbose_mode = VB_QUIET;

	cout << endl << "--> INFERING RATE ASSUMING POSITION-SPECIFIC MODEL..." << endl << endl;
	for (int pos = 0; pos < aln.ncategory; pos++) {
		cout << "Position " << pos << " / ";
		int *pair_pos = aln.pair_freq + (pos*aln.num_states*aln.num_states);
		testSingleRateModel(params, aln, tree, original_model, pair_pos, part_rate[pos], rate_name, false);
	}


	verbose_mode = vb_saved;

	int sum_freq[aln.num_states*aln.num_states];
	cout << endl << "-->INFERING RATE UNDER EQUAL-RATE NULL MODEL..." << endl << endl;
	aln.computeSumPairFreq(sum_freq);
	DoubleVector null_rate;
	testSingleRateModel(params, aln, tree, original_model, sum_freq, null_rate, rate_name, true);

	// report running results
	string out_file = params.out_prefix;
	out_file += ".ngs";
	reportNGSAnalysis(out_file.c_str(), params, aln, tree, part_rate, rate_name);

}