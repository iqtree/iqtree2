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

ModelSubst* createModel(string model_str, StateFreqType freq_type, string freq_params,
		PhyloTree* tree, bool count_rates)
{
	ModelSubst *model = NULL;
	//cout << "Numstates: " << tree->aln->num_states << endl;
	string model_params;
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
	} else if (model_str == "MIX" || model_str == "MIXTURE") {
		model = new ModelMixture(model_str, model_params, freq_type, freq_params, tree, count_rates);
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

ModelMixture::ModelMixture(string model_name, string model_list, StateFreqType freq, string freq_params, PhyloTree *tree, bool count_rates)
	: ModelGTR(tree, count_rates)
{
//	cout << model_name << " " << model_list << endl;
	if (freq_params != "")
		readStateFreq(freq_params);
	init(freq);

	const int MAX_MODELS = 64;
	size_t cur_pos = 0;
	name = full_name = (string)"MIXTURE" + OPEN_BRACKET;
	for (int m = 0; m < MAX_MODELS && cur_pos < model_list.length(); m++) {
		size_t pos = model_list.find(',', cur_pos);
		if (pos == string::npos)
			pos = model_list.length();
		if (pos <= cur_pos)
			outError("One model name in the mixture is empty.");
		string this_name = model_list.substr(cur_pos, pos-cur_pos);
		cur_pos = pos+1;
		push_back(createModel(this_name, freq, freq_params, tree, count_rates));
		if (m > 0) {
			name += ',';
			full_name += ',';
		}
		name += back()->name;
		full_name += back()->full_name;
	}
	name += CLOSE_BRACKET;
	full_name += CLOSE_BRACKET;
	cout << "mixture model: " << name << endl;
}

ModelMixture::~ModelMixture() {
	for (reverse_iterator rit = rbegin(); rit != rend(); rit++)
		delete (*rit);
}

