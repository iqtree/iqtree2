/*
 * modelmixture.cpp
 *
 *  Created on: Nov 29, 2014
 *      Author: minh
 */

#include "modelgtr.h"
#include "modelnonrev.h"
#include "modeldna.h"
#include "modelprotein.h"
#include "modelbin.h"
#include "modelcodon.h"
#include "modelmorphology.h"
#include "modelset.h"
#include "modelmixture.h"

const double MIN_MIXTURE_PROP = 0.0;
const double MAX_MIXTURE_PROP = 1.0;
//const double MIN_MIXTURE_RATE = 0.01;
//const double MAX_MIXTURE_RATE = 100.0;

ModelSubst* createModel(string model_str, string model_desc, StateFreqType freq_type, string freq_params,
		PhyloTree* tree, bool count_rates)
{
	ModelSubst *model = NULL;
	//cout << "Numstates: " << tree->aln->num_states << endl;
	string model_params = model_desc;
	size_t pos = model_str.find(OPEN_BRACKET);
	if (pos != string::npos) {
		if (model_str.rfind(CLOSE_BRACKET) != model_str.length()-1)
			outError("Close bracket not found at the end of ", model_str);
		model_params = model_str.substr(pos+1, model_str.length()-pos-2);
		model_str = model_str.substr(0, pos);
	}
	/*
	if ((model_str == "JC" && tree->aln->seq_type == SEQ_DNA) ||
		(model_str == "POISSON" && tree->aln->seq_type == SEQ_PROTEIN) ||
		(model_str == "JC2" && tree->aln->seq_type == SEQ_BINARY) ||
		(model_str == "JCC" && tree->aln->seq_type == SEQ_CODON) ||
		(model_str == "MK" && tree->aln->seq_type == SEQ_MORPH))
	{
		model = new ModelSubst(tree->aln->num_states);
	} else */
	if ((model_str == "GTR" && tree->aln->seq_type == SEQ_DNA) ||
		(model_str == "GTR2" && tree->aln->seq_type == SEQ_BINARY) ||
		(model_str == "GTR20" && tree->aln->seq_type == SEQ_PROTEIN)) {
		model = new ModelGTR(tree, count_rates);
		if (freq_params != "")
			((ModelGTR*)model)->readStateFreq(freq_params);
		if (model_params != "")
			((ModelGTR*)model)->readRates(model_params);
		((ModelGTR*)model)->init(freq_type);
	} else if (model_str == "UNREST") {
		freq_type = FREQ_EQUAL;
		//params.optimize_by_newton = false;
		tree->optimize_by_newton = false;
		model = new ModelNonRev(tree, count_rates);
		((ModelNonRev*)model)->init(freq_type);
	} else if (tree->aln->seq_type == SEQ_BINARY) {
		model = new ModelBIN(model_str.c_str(), model_params, freq_type, freq_params, tree, count_rates);
	} else if (tree->aln->seq_type == SEQ_DNA) {
		model = new ModelDNA(model_str.c_str(), model_params, freq_type, freq_params, tree, count_rates);
	} else if (tree->aln->seq_type == SEQ_PROTEIN) {
		model = new ModelProtein(model_str.c_str(), model_params, freq_type, freq_params, tree, count_rates);
	} else if (tree->aln->seq_type == SEQ_CODON) {
		model = new ModelCodon(model_str.c_str(), model_params, freq_type, freq_params, tree, count_rates);
	} else if (tree->aln->seq_type == SEQ_MORPH) {
		model = new ModelMorphology(model_str.c_str(), model_params, freq_type, freq_params, tree);
	} else {
		outError("Unsupported model type");
	}

	return model;
}

ModelMixture::ModelMixture(string model_name, string model_list, ModelsBlock *models_block, StateFreqType freq, string freq_params, PhyloTree *tree, bool count_rates)
	: ModelGTR(tree, count_rates)
{
	if (freq_params != "")
		readStateFreq(freq_params);
	init(freq);

	const int MAX_MODELS = 64;
	size_t cur_pos = 0;
	int m;
	name = full_name = (string)"MIX" + OPEN_BRACKET;
	for (m = 0; m < MAX_MODELS && cur_pos < model_list.length(); m++) {
		size_t pos = model_list.find(',', cur_pos);
		if (pos == string::npos)
			pos = model_list.length();
		if (pos <= cur_pos)
			outError("One model name in the mixture is empty.");
		string this_name = model_list.substr(cur_pos, pos-cur_pos);
		double rate = 1.0;
		if (this_name.find(':') != string::npos) {
			rate = convert_double(this_name.substr(this_name.find(':')+1).c_str());
			this_name = this_name.substr(0, this_name.find(':'));
		}
		cur_pos = pos+1;
		string model_desc;
		NxsModel *nxsmodel = models_block->findModel(this_name);
		if (nxsmodel) model_desc = nxsmodel->description;
		push_back((ModelGTR*)createModel(this_name, model_desc, freq, freq_params, tree, count_rates));
		back()->total_num_subst = rate;
		if (m > 0) {
			name += ',';
			full_name += ',';
		}
		name += back()->name;
		full_name += back()->full_name;
	}
	name += CLOSE_BRACKET;
	full_name += CLOSE_BRACKET;

	int nmixtures = size();
	prop = aligned_alloc<double>(nmixtures);

	double sum = 0.0;
	int i;
	// initialize rates as increasing
	for (i = 0; i < nmixtures; i++) {
		prop[i] = (double)(nmixtures-i);
//		sum += prop[i]*at(i)->total_num_subst;
		sum += prop[i];
	}
	for (i = 0; i < nmixtures; i++) {
		 prop[i] /= sum;
	}

	// rescale total_num_subst such that the global rate is 1
	for (i = 0, sum = 0.0; i < nmixtures; i++)
		sum += prop[i]*at(i)->total_num_subst;
	for (i = 0; i < nmixtures; i++)
		at(i)->total_num_subst /= sum;

	fix_prop = (nmixtures == 1);
	// use central eigen etc. stufffs

	if (eigenvalues) delete [] eigenvalues;
	if (eigenvectors) delete [] eigenvectors;
	if (inv_eigenvectors) delete [] inv_eigenvectors;
	if (eigen_coeff) delete [] eigen_coeff;

	eigenvalues = new double[num_states*nmixtures];
	eigenvectors = new double[num_states*num_states*nmixtures];
	inv_eigenvectors = new double[num_states*num_states*nmixtures];
	int ncoeff = num_states*num_states*num_states;
	eigen_coeff = new double[ncoeff*nmixtures];

	// assigning memory for individual models
	m = 0;
	for (iterator it = begin(); it != end(); it++, m++) {
        // first copy memory for eigen stuffs
        memcpy(&eigenvalues[m*num_states], (*it)->eigenvalues, num_states*sizeof(double));
        memcpy(&eigenvectors[m*num_states*num_states], (*it)->eigenvectors, num_states*num_states*sizeof(double));
        memcpy(&inv_eigenvectors[m*num_states*num_states], (*it)->inv_eigenvectors, num_states*num_states*sizeof(double));
        memcpy(&eigen_coeff[m*ncoeff], (*it)->eigen_coeff, ncoeff*sizeof(double));
        // then delete
		if ((*it)->eigenvalues) delete [] (*it)->eigenvalues;
		if ((*it)->eigenvectors) delete [] (*it)->eigenvectors;
		if ((*it)->inv_eigenvectors) delete [] (*it)->inv_eigenvectors;
		if ((*it)->eigen_coeff) delete [] (*it)->eigen_coeff;

        // and assign new memory
		(*it)->eigenvalues = &eigenvalues[m*num_states];
		(*it)->eigenvectors = &eigenvectors[m*num_states*num_states];
		(*it)->inv_eigenvectors = &inv_eigenvectors[m*num_states*num_states];
		(*it)->eigen_coeff = &eigen_coeff[m*ncoeff];
	}
	decomposeRateMatrix();
}

ModelMixture::~ModelMixture() {
	if (prop)
		aligned_free(prop);
	for (reverse_iterator rit = rbegin(); rit != rend(); rit++) {
		(*rit)->eigen_coeff = NULL;
		(*rit)->eigenvalues = NULL;
		(*rit)->eigenvectors = NULL;
		(*rit)->inv_eigenvectors = NULL;
		delete (*rit);
	}
}

int ModelMixture::getNDim() {
	int dim = (fix_prop) ? 0: (size()-1);
	for (iterator it = begin(); it != end(); it++)
		dim += (*it)->getNDim();
	return dim;
}

double ModelMixture::targetFunk(double x[]) {
	getVariables(x);
	// always decompose rate matrix since total_num_subst might change
	decomposeRateMatrix();
//	for (iterator it = begin(); it != end(); it++)
//		if ((*it)->getNDim() > 0)
//			(*it)->decomposeRateMatrix();
	assert(phylo_tree);
	phylo_tree->clearAllPartialLH();
	return -phylo_tree->computeLikelihood();
}

void ModelMixture::decomposeRateMatrix() {
	for (iterator it = begin(); it != end(); it++)
		(*it)->decomposeRateMatrix();
}

void ModelMixture::setVariables(double *variables) {
	int dim = 0;
	for (iterator it = begin(); it != end(); it++) {
		(*it)->setVariables(&variables[dim]);
		dim += (*it)->getNDim();
	}
	if (fix_prop) return;
	int i, ncategory = size();
//	variables[dim+1] = prop[0]*at(0)->total_num_subst;
//	for (i = 2; i < ncategory; i++)
//		variables[dim+i] = variables[dim+i-1] + prop[i-1]*at(i-1)->total_num_subst;
	variables[dim+1] = prop[0];
	for (i = 2; i < ncategory; i++)
		variables[dim+i] = variables[dim+i-1] + prop[i-1];
}

void ModelMixture::getVariables(double *variables) {
	int dim = 0;
	for (iterator it = begin(); it != end(); it++) {
		(*it)->getVariables(&variables[dim]);
		dim += (*it)->getNDim();
	}
	if (fix_prop) return;
	int i, ncategory = size();
	double *y = new double[ncategory+1];
	y[0] = 0; y[ncategory] = 1.0;
	memcpy(y+1, variables+dim+1, (ncategory-1) * sizeof(double));
	std::sort(y+1, y+ncategory);
	double sum = 0.0;
	for (i = 0; i < ncategory; i++) {
		prop[i] = (y[i+1]-y[i]);
	}
	for (i = 0, sum = 0.0; i < ncategory; i++)
		sum += prop[i]*at(i)->total_num_subst;
	for (i = 0; i < ncategory; i++)
		at(i)->total_num_subst /= sum;

	if (verbose_mode >= VB_MAX) {
		for (i = 0; i < ncategory; i++)
			cout << "Component " << i << " prop=" << prop[i] << endl;
	}
	delete [] y;

}

void ModelMixture::setBounds(double *lower_bound, double *upper_bound, bool *bound_check) {
	int dim = 0;
	for (iterator it = begin(); it != end(); it++) {
		(*it)->setBounds(&lower_bound[dim], &upper_bound[dim], &bound_check[dim]);
		dim += (*it)->getNDim();
	}
	if (fix_prop) return;
	int i, ncategory = size();
	for (i = 1; i < ncategory; i++) {
		lower_bound[dim+i] = MIN_MIXTURE_PROP;
		upper_bound[dim+i] = MAX_MIXTURE_PROP;
		bound_check[dim+i] = false;
	}
}

void ModelMixture::writeInfo(ostream &out) {
	for (int i = 0; i < size(); i++) {
		if (!fix_prop)
			out << "Weight of mixture component " << i << " (" << at(i)->name << "): "
			<< prop[i] << endl;
		at(i)->writeInfo(out);
	}
}

void ModelMixture::writeParameters(ostream &out) {
	for (iterator it = begin(); it != end(); it++) {
		(*it)->writeParameters(out);
	}
}
