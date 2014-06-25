#include "modelpomo.h"
#include <stdlib.h>
#include <assert.h>
#include <string.h>

ModelPoMo::ModelPoMo(PhyloTree *tree, bool count_rate) : GTRModel(tree, count_rate) {
	init("PoMo", "", FREQ_EQUAL, "");
}

void ModelPoMo::init(const char *model_name, string model_params, StateFreqType freq, string freq_params) {
	initMoranWithMutation();
	GTRModel::init(freq);
}

void ModelPoMo::initMoranWithMutation() {
	int i=0;

	for (i=0; i<num_states*(num_states-1)/2; i++) {
		rates[i] = 1.0;
	}
	for (i=0; i<num_states; i++) {
		state_freq[i] = 1.0/num_states;
	}
	/* TODO Initialize rate matrix here. */
	/* TODO Future: Initialize state freqs here. */
}

int ModelPoMo::getNDim() {
	return 0;
}

