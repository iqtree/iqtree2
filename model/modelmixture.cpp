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

ModelSubst* createModel(string model_str, ModelsBlock *models_block, StateFreqType freq_type, string freq_params,
		PhyloTree* tree, bool count_rates)
{
	ModelSubst *model = NULL;
	//cout << "Numstates: " << tree->aln->num_states << endl;
	string model_params;
	NxsModel *nxsmodel = models_block->findModel(model_str);
	if (nxsmodel) model_params = nxsmodel->description;
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

ModelMixture::ModelMixture(string model_name, string model_list, ModelsBlock *models_block,
		StateFreqType freq, string freq_params, PhyloTree *tree, bool optimize_weights, bool count_rates)
	: ModelGTR(tree, count_rates)
{
	const int MAX_MODELS = 64;
	size_t cur_pos;
	int m;

	vector<NxsModel*> freq_vec;
	DoubleVector freq_rates;
	DoubleVector freq_weights;
	fix_prop = false;

	if (freq == FREQ_MIXTURE) {
		for (m = 0, cur_pos = 0; m < MAX_MODELS && cur_pos < freq_params.length(); m++) {
			size_t pos = freq_params.find(',', cur_pos);
			if (pos == string::npos)
				pos = freq_params.length();
			if (pos <= cur_pos)
				outError("One frequency name in the mixture is empty.");
			string this_name = freq_params.substr(cur_pos, pos-cur_pos);
			double rate = 1.0, weight = 1.0;
			size_t pos_rate = this_name.find(':');
			if (pos_rate != string::npos) {
				size_t pos_weight = this_name.find(':', pos_rate+1);
				if (pos_weight == string::npos) {
					rate = convert_double(this_name.substr(pos_rate+1).c_str());
				} else {
					rate = convert_double(this_name.substr(pos_rate+1, pos_weight-pos_rate-1).c_str());
					weight = convert_double(this_name.substr(pos_weight+1).c_str());
					fix_prop = true;
					if (weight <= 0.0)
						outError("Mixture component weight is negative!");
				}
				this_name = this_name.substr(0, pos_rate);
			}
			freq_rates.push_back(rate);
			freq_weights.push_back(weight);
			cur_pos = pos+1;
			if (this_name == "empirical") {
				freq_vec.push_back(NULL);
			} else {
				NxsModel *freq_mod = models_block->findModel(this_name);
				if (!freq_mod)
					outError("Frequency mixture name not found ", this_name);
				if (!(freq_mod->flag & NM_FREQ)) {
					cout << freq_mod->flag << endl;
					outError("Frequency mixture name does not corresponding to frequency model ", this_name);
				}
				freq_vec.push_back(freq_mod);
			}
		}
		init(FREQ_USER_DEFINED);
	} else {
		if (freq_params != "")
			readStateFreq(freq_params);
		init(freq);
	}

	DoubleVector weights;
	name = full_name = (string)"MIX" + OPEN_BRACKET;
	if (model_list == "") model_list = model_name;
	for (m = 0, cur_pos = 0; m < MAX_MODELS && cur_pos < model_list.length(); m++) {
		size_t pos = model_list.find(',', cur_pos);
		if (pos == string::npos)
			pos = model_list.length();
		if (pos <= cur_pos)
			outError("One model name in the mixture is empty.");
		string this_name = model_list.substr(cur_pos, pos-cur_pos);
		double rate = 1.0, weight = 1.0;
		size_t pos_rate = this_name.find(':');
		if (pos_rate != string::npos) {
			size_t pos_weight = this_name.find(':', pos_rate+1);
			if (pos_weight == string::npos) {
				rate = convert_double(this_name.substr(pos_rate+1).c_str());
			} else {
				rate = convert_double(this_name.substr(pos_rate+1, pos_weight-pos_rate-1).c_str());
				weight = convert_double(this_name.substr(pos_weight+1).c_str());
				fix_prop = true;
				if (weight <= 0.0)
					outError("Mixture component weight is negative!");
			}
			this_name = this_name.substr(0, pos_rate);
		}
		cur_pos = pos+1;
		ModelGTR* model;
		if (freq == FREQ_MIXTURE) {
			for(int f = 0; f != freq_vec.size(); f++) {
				if (freq_vec[f])
					model = (ModelGTR*)createModel(this_name, models_block, FREQ_USER_DEFINED, freq_vec[f]->description, tree, count_rates);
				else
					model = (ModelGTR*)createModel(this_name, models_block, FREQ_EMPIRICAL, "", tree, count_rates);
				model->total_num_subst = rate * freq_rates[f];
				push_back(model);
				weights.push_back(weight * freq_weights[f]);
				if (m+f > 0) {
					name += ',';
					full_name += ',';
				}
				if (freq_vec[f]) {
					name += model->name + "+F{" +freq_vec[f]->name + "}";
					full_name += model->full_name + "+F{" +freq_vec[f]->name + "}";
				} else {
					name += model->name + "+F";
					full_name += model->full_name + "+F";
				}
			}
		} else {
			model = (ModelGTR*)createModel(this_name, models_block, freq, freq_params, tree, count_rates);
			model->total_num_subst = rate;
			push_back(model);
			weights.push_back(weight);
			if (m > 0) {
				name += ',';
				full_name += ',';
			}
			name += model->name;
			full_name += model->full_name;
		}
	}

	name += CLOSE_BRACKET;
	full_name += CLOSE_BRACKET;

	int nmixtures = size();
	prop = aligned_alloc<double>(nmixtures);

	double sum = 0.0;
	int i;
	if (fix_prop) {
		for (i = 0, sum = 0.0; i < nmixtures; i++) {
			prop[i] = weights[i];
			sum += prop[i];
		}
	} else {
		// initialize rates as increasing
		for (i = 0, sum = 0.0; i < nmixtures; i++) {
			prop[i] = random_double();
			sum += prop[i];
		}
	}
	// normalize weights to 1.0
	for (i = 0; i < nmixtures; i++)
		 prop[i] /= sum;

	// rescale total_num_subst such that the global rate is 1
	for (i = 0, sum = 0.0; i < nmixtures; i++)
		sum += prop[i]*at(i)->total_num_subst;
	for (i = 0; i < nmixtures; i++)
		at(i)->total_num_subst /= sum;

	if (optimize_weights) fix_prop = false;
	fix_prop |= (nmixtures == 1);
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
//	decomposeRateMatrix();
	for (iterator it = begin(); it != end(); it++)
		if ((*it)->getNDim() > 0)
			(*it)->decomposeRateMatrix();
	assert(phylo_tree);
	phylo_tree->clearAllPartialLH();
	if (prop[size()-1] < 0.0) return 1.0e+12;
	return -phylo_tree->computeLikelihood();
}

double ModelMixture::optimizeParameters(double epsilon) {
	double score = ModelGTR::optimizeParameters(epsilon);
	if (getNDim() == 0) return score;
	// now rescale Q matrices to have proper interpretation of branch lengths
	double sum;
	int i, ncategory = size();
	for (i = 0, sum = 0.0; i < ncategory; i++)
		sum += prop[i]*at(i)->total_num_subst;
	for (i = 0; i < ncategory; i++)
		at(i)->total_num_subst /= sum;
	decomposeRateMatrix();
	return score;
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
//	double sum = 0.0;
	for (i = 0; i < ncategory; i++) {
		prop[i] = (y[i+1]-y[i]);
	}
//	for (i = 0, sum = 0.0; i < ncategory; i++)
//		sum += prop[i]*at(i)->total_num_subst;
//	for (i = 0; i < ncategory; i++)
//		at(i)->total_num_subst /= sum;

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
	int i;
	for (i = 0; i < size(); i++) {
		at(i)->writeInfo(out);
	}
//	if (fix_prop) return;
	cout << "Mixture weights:";
	for (i = 0; i < size(); i++)
		cout << " " << prop[i];
	cout << endl;
}

void ModelMixture::writeParameters(ostream &out) {
	for (iterator it = begin(); it != end(); it++) {
		(*it)->writeParameters(out);
	}
}
