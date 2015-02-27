#include "modelpomo.h"
#include <stdlib.h>
#include <assert.h>
#include <string.h>

ModelPoMo::ModelPoMo(PhyloTree *tree, bool count_rate) : ModelNonRev(tree, count_rate) {
	init("PoMo", "", FREQ_EQUAL, "");
}

void ModelPoMo::init(const char *model_name, string model_params, StateFreqType freq, string freq_params) {
	mutation_prob = new double[16];
	freq_fixed_states = new double[4];
	// TODO: need to reoptimize later
	int i;
	for (i = 0; i < 16; i++) mutation_prob[i] = 1.0;
	// TODO: check these starting values
	for (i = 0; i < 4; i++) freq_fixed_states[i] = 0.2;
	initMoranWithBoundaryMutation();
	ModelGTR::init(freq);
}

ModelPoMo::~ModelPoMo() {
	// TODO: Check with ModelPoMo::init!
	delete [] freq_fixed_states;
	delete [] mutation_prob;
}

double ModelPoMo::computeNormConst() {
	int i, j;
	int N = phylo_tree->aln->virtual_pop_size;
	double harmonic = 0.0;
	for (i = 1; i < N; i++)
		harmonic += 1.0/(double)i;

	double norm_fixed = 0.0, norm_polymorphic = 0.0;
	for (i = 0; i < 4; i++)
		norm_fixed += freq_fixed_states[i];
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++)
			if (i != j)
				norm_polymorphic += freq_fixed_states[i] * freq_fixed_states[j] * mutCoeff(i, j);
	}
	norm_polymorphic *= N * harmonic;

	return 1.0/(norm_fixed + norm_polymorphic);
}

void ModelPoMo::computeStateFreq() {
	double norm = computeNormConst();
	int state;
	int N = phylo_tree->aln->virtual_pop_size;

	for (state = 0; state < num_states; state++)
		if (isFixed(state))
			state_freq[state] = freq_fixed_states[state]*norm;
		else {
			int k, X, Y;
			decomposeState(state, k, X, Y);
			state_freq[state] = norm*freq_fixed_states[X]*freq_fixed_states[Y]*mutCoeff(X, Y)*N*N / (k*(N-k));
		}
}

void ModelPoMo::initMoranWithBoundaryMutation() {
	int state1, state2;
	int N = phylo_tree->aln->virtual_pop_size;

	if (verbose_mode >= VB_MED) cout << "PoMo rate matrix:" << endl;

	computeStateFreq();

	// Loop over rows (transition starting from state1)
	for (state1 = 0; state1 < num_states; state1++) {
		double row_sum = 0.0;
		// Loop over columns (transition to state2)
		for (state2 = 0; state2 < num_states; state2++)
			if (state2 != state1) {
				row_sum += (rate_matrix[state1*num_states+state2] = computeProbBoundaryMutation(state1, state2));
			}
		rate_matrix[state1*num_states+state1] = -(row_sum);
	}
}

void ModelPoMo::initMoranWithMutation() {

	// // This code was used to run a dummy JC69 model with 58 states.
	// int i=0;

	// for (i=0; i<num_states*(num_states-1)/2; i++) {
	// 	rates[i] = 1.0;
	// }
	// for (i=0; i<num_states; i++) {
	// 	state_freq[i] = 1.0/num_states;
	// }

	// TODO: Recheck this.
	// Initialize rate matrix Q[state1,state2] = transition rate from
	// state1 to state2.
	int state1, state2;
	if (verbose_mode >= VB_MED) cout << "PoMo rate matrix:" << endl;
	// Loop over rows (transition starting from state1)
	for (state1 = 0; state1 < num_states; state1++) {
		double row_sum = 0.0;
		// Loop over columns (transition to state2)
		for (state2 = 0; state2 < num_states; state2++) {
			if (state1 == state2) {
				// Q = P - I
				rate_matrix[state1*num_states+state2] = computeProb(state1, state2) - 1.0;
			} else {
				rate_matrix[state1*num_states+state2] = computeProb(state1, state2);
			}
			// Compute row sum.
			row_sum += rate_matrix[state1*num_states+state2];
			if (verbose_mode >= VB_MED)
				cout << rate_matrix[state1*num_states+state2] << "\t";
		}
		if (verbose_mode >= VB_MED) cout << endl;
		// Check row sum.
		if (fabs(row_sum) > 0.000001) outError("Row sum not equal 0");
	}

	/* TODO: Initialize state freqs here. */
}

double ModelPoMo::computeP(int i, int major, int minor) {
	// Cf. Moran model with mutation.
	int N = phylo_tree->aln->virtual_pop_size;
	return (1.0 - mutation_prob[major*4+minor])*(i)*(N-i)/(N*N) +
			mutation_prob[minor*4+major]*(N-i)*(N-i)/(N*N);
}

double ModelPoMo::computeR(int i, int major, int minor) {
	// Cf. Moran model with mutation.
	int N = phylo_tree->aln->virtual_pop_size;
	return (mutation_prob[major*4+minor]+mutation_prob[minor*4+major])*i*(N-i)/(N*N) +
			(1-mutation_prob[minor*4+major])*(N-i)*(N-i)/(N*N) +
			(1-mutation_prob[major*4+minor])*i*i/(N*N);
}

void ModelPoMo::decomposeState(int state, int &i, int &nt1, int &nt2) {
	int N = phylo_tree->aln->virtual_pop_size;
	if (state < 4) {
		// Fixed A, C, G or T
		i = N;
		nt1 = state;
		nt2 = -1; // -1 for unknown nt
	} else if (state < 4+(N-1)) {
		// (iA,N-iC)
		i = state-3;
		nt1 = 0; // A
		nt2 = 1; // C
	} else if (state < 4+2*(N-1)) {
		// (iA,N-iG)
		i = state-3-(N-1);
		nt1 = 0; // A
		nt2 = 2; // G
	} else if (state < 4+3*(N-1)) {
		// (iA,N-iT)
		i = state-3-2*(N-1);
		nt1 = 0; // A
		nt2 = 3; // T
	} else if (state < 4+4*(N-1)) {
		// (iC,N-iG)
		i = state-3-3*(N-1);
		nt1 = 1; // C
		nt2 = 2; // G
	} else if (state < 4+5*(N-1)) {
		// (iC,N-iT)
		i = state-3-4*(N-1);
		nt1 = 1; // C
		nt2 = 3; // T
	} else if (state < 4+6*(N-1)) {
		// (iG,N-iT)
		i = state-3-5*(N-1);
		nt1 = 2; // G
		nt2 = 3; // T
	} else {
		outError("State exceeds limit");
	}
}

double ModelPoMo::computeProb(int state1, int state2) {
	int N = phylo_tree->aln->virtual_pop_size;

	// Both states are decomposed into the abundance of the first
	// allele as well as the nucleotide of the first and the second
	// allele.
	int i1=0, i2=0, nt1=-1, nt2=-1, nt3=-1, nt4=-1;
	decomposeState(state1, i1, nt1, nt2);
	decomposeState(state2, i2, nt3, nt4);
	// Either the first nucleotides match or the first of state 1 with
	// the second of state 2 or the first of state 2 with the second
	// of state 1.  Additionally, we have to consider fixed states as
	// special cases.
	if (nt1 == nt3 && (nt2==nt4 || nt2==-1 || nt4 == -1)) {
		if (i1==i2) {
			if (i1==N) {
				// e.g.: 10A -> 10A
				double sum = 0;
				for (int nt=0; nt < 4; nt++)
					if (nt != nt1) sum += computeR(N, nt1, nt);
				return sum-2.0;
			} else
				// e.g.: nA(N-n)C -> nA(N-n)C where 0<n<N
				return computeR(i1, nt1, nt2);
		} else if (i1+1==i2)
			// e.g.: 2A8C -> 3A7C or 9A1C -> 10A
			if (nt2 == -1)
				// e.g. 10A ->
				return computeP(i1, nt1, nt4);
			else
				return computeP(i1, nt1, nt2);
		else if (i1-1 == i2)
			// e.g.: 3A7C -> 2A8C or 10A -> 9A1C
			if (nt2 == -1)
				// e.g. 10A -> 9A1C
				return computeP(N-i1, nt4, nt1);
			else
				// e.g. 9A1C -> 8A2C
				return computeP(N-i1, nt2, nt1);
		else
			// 0 for all others
			return 0.0;
	} else if (nt1 == nt4 && nt2 == -1 && i2 == 1)  {
		// e.g.: 10G -> 1A9G
		return computeP(0, nt3, nt1);
	} else if (nt2 == nt3  && i1 == 1 && nt4 == -1) {
		// E.g.: 1A9G -> 10G
		return computeP(N-1, nt2, nt1);
	} else
		// 0 for all other transitions
		return 0.0;
}

bool ModelPoMo::isFixed(int state) {
	return (state < 4);
}

bool ModelPoMo::isPolymorphic(int state) {
	return (!isFixed(state));
}

double ModelPoMo::mutCoeff(int nt1, int nt2) {
	assert(nt1!=nt2 && nt1<4 && nt2<4);
	if (nt2 < nt1) {
		int tmp=nt1;
		nt1=nt2;
		nt2=tmp;
	}
	if (nt1==0) return mutation_prob[nt2-1];
	if (nt1==1) return mutation_prob[nt2+1];
	if (nt1==2) return mutation_prob[5];
	assert(0);
}

double ModelPoMo::computeProbBoundaryMutation(int state1, int state2) {
	int N = phylo_tree->aln->virtual_pop_size;

	assert(state1 != state2);

	// Both states are decomposed into the abundance of the first
	// allele as well as the nucleotide of the first and the second
	// allele.
	int i1=0, i2=0, nt1=-1, nt2=-1, nt3=-1, nt4=-1;
	decomposeState(state1, i1, nt1, nt2);
	decomposeState(state2, i2, nt3, nt4);
	// Either the first nucleotides match or the first of state 1 with
	// the second of state 2 or the first of state 2 with the second
	// of state 1.  Additionally, we have to consider fixed states as
	// special cases.
	if (nt1 == nt3 && (nt2==nt4 || nt2==-1 || nt4 == -1)) {
		assert(i1 != i2); // because state1 != state2
		if (i1+1==i2)
			// e.g.: 2A8C -> 3A7C or 9A1C -> 10A
			if (nt2 == -1)
				// e.g. 10A ->
				assert(0);
			else
				return double(i1*(N-i1)) / double(N*N);
		else if (i1-1 == i2)
			// e.g.: 3A7C -> 2A8C or 10A -> 9A1C
			if (nt2 == -1)
				// e.g. 10A -> 9A1C
				return mutCoeff(nt1,nt4) * state_freq[nt4];
			else
				// e.g. 9A1C -> 8A2C
				return double(i1*(N-i1)) / double(N*N);
		else
			// 0 for all others
			return 0.0;
	} else if (nt1 == nt4 && nt2 == -1 && i2 == 1)  {
		// e.g.: 10G -> 1A9G
		return mutCoeff(nt1,nt3) * state_freq[nt3];
	} else if (nt2 == nt3  && i1 == 1 && nt4 == -1) {
		// E.g.: 1A9G -> 10G
		return double(i1*(N-i1)) / double(N*N);
	} else
		// 0 for all other transitions
		return 0.0;
}

int ModelPoMo::getNDim() {
	return 9;
}

void ModelPoMo::setBounds(double *lower_bound, double *upper_bound, bool *bound_check) {
	int i;
	for (i = 1; i <= 3; i++) {
		lower_bound[i] = 0.0001;
		// TODO: check this bound!
		upper_bound[i] = 100;
		bound_check[i] = false;
	}
	for (i = 4; i <= 9; i++) {
		//TODO: check this bound!
		lower_bound[i] = MIN_RATE;
		lower_bound[i] = MAX_RATE;
		bound_check[i] = false;
	}
}

void ModelPoMo::setVariables(double *variables) {
	int i;
	for (i = 1; i <= 3; i++) {
		variables[i] = freq_fixed_states[i-1];
	}
	for (i = 4; i <= 9; i++) {
		variables[i] = mutation_prob[i-4];
	}

}

void ModelPoMo::getVariables(double *variables) {
	int i;
	for (i = 1; i <= 3; i++) {
		freq_fixed_states[i-1] = variables[i];
	}
	for (i = 4; i <= 9; i++) {
		mutation_prob[i-4] = variables[i];
	}
	initMoranWithBoundaryMutation();
}

