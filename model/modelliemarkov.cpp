/*
 * modelliemarkov.cpp
 *
 *  Created on: 24/05/2016
 *      Author: Michael Woodhams
 */

/*
 * TO DO: 
 * Currently all LM models are treated as non time-reversible.
 * In fact, a few are time reversible (TIME_REVERSIBLE[model_num]).
 * Allowing this to be recognized by rest of program would speed up
 * analysis of these models significantly.
 */
/*
 * TO DO:
 * Currently symmetry permutation is applied every time setRates is called.
 * Would be more efficient to apply it just once to basis in constructor.
 */
#include "modelliemarkov.h"
#include <float.h>
#undef NDEBUG
#include <assert.h>

/*
 * TO DO: It is inconvenient to have all of these const declarations
 * here, as it is so much to scroll past to find the actual code.
 * Find some way to be able to put these at end of .cpp file,
 * or maybe in .h file (not preferred, as nothing external to this
 * file needs or is allowed to see the definitions.)
 */

/*
 * Basis matrices as displayed are columns sum to zero convention.
 * For RY pairing, columns/rows are in order A, G, C, T
 * For WS pairing, columns/rows are in order A, T, C, G
 * For MK pairing, columns/rows are in order A, C, G, T
 * Only the off-diagonal entries are stored. rates array is in
 * rows sum to zero convention, so we read down the columns
 * to enter them into the array (skipping diagonal entries.)
 *
 * From table 1 of Woodhams et al Syst Biol 64 p638-650 (2015) DOI:10.1093/sysbio/syv021
 */
/* A2:
 *  0 +2 -1 -1
 * +2  0 -1 -1
 * -1 -1  0 +2
 * -1 -1 +2  0
 */
const static double A2[] = {+2,-1,-1,+2,-1,-1,-1,-1,+2,-1,-1,+2};
/* B:
 * 0 0 + -
 * 0 0 - +
 * + - 0 0
 * - + 0 0
 */
const static double B[]  = { 0,+1,-1, 0,-1,+1,+1,-1, 0,-1,+1, 0};
/* C:
 * 0 0 + -
 * 0 0 - +
 * - + 0 0
 * + - 0 0
 */
const static double C[]  = { 0,-1,+1, 0,+1,-1,+1,-1, 0,-1,+1, 0};
/* D1:
 * - + 0 0
 * + - 0 0
 * 0 0 + -
 * 0 0 - +
 */
const static double D1[] = {+1, 0, 0,+1, 0, 0, 0, 0,-1, 0, 0,-1};
/* D:
 * + + + +
 * + + + +
 * - - - -
 * - - - -
 */
const static double D[]  = {+1,-1,-1,+1,-1,-1,+1,+1,-1,+1,+1,-1};
/* E1:
 * + + + +
 * - - - -
 * 0 0 0 0
 * 0 0 0 0
 */
const static double E1[] = {-1, 0, 0,+1, 0, 0,+1,-1, 0,+1,-1, 0};
/* E2:
 * 0 0 0 0
 * 0 0 0 0
 * + + + +
 * - - - -
 */
const static double E2[] = { 0,+1,-1, 0,+1,-1, 0, 0,-1, 0, 0,+1};
/* F1:
 * + + - -
 * - - + +
 * 0 0 0 0
 * 0 0 0 0
 */
const static double F1[] = {-1, 0, 0,+1, 0, 0,-1,+1, 0,-1,+1, 0};
/* F2:
 * 0 0 0 0
 * 0 0 0 0
 * + + - -
 * - - + +
 */
const static double F2[] = { 0,+1,-1, 0,+1,-1, 0, 0,+1, 0, 0,-1};
/* G1:
 * + - 0 0
 * + - 0 0
 * - + 0 0
 * - + 0 0
 */
const static double G1[] = {+1,-1,-1,-1,+1,+1, 0, 0, 0, 0, 0, 0};
/* G2:
 * 0 0 + -
 * 0 0 + -
 * 0 0 - +
 * 0 0 - +
 */
const static double G2[] = { 0, 0, 0, 0, 0, 0,+1,+1,-1,-1,-1,+1};

// Lengths of these arrays stored in MODEL_PARAMS
const static double *BASIS_11[]   = {};
const static double *BASIS_22B[]  = {A2};
const static double *BASIS_33A[]  = {A2,B };
const static double *BASIS_33B[]  = {A2,C };
const static double *BASIS_33C[]  = {A2,D1};
const static double *BASIS_34[]   = {A2,D };
const static double *BASIS_44A[]  = {D, E1,E2};
const static double *BASIS_44B[]  = {A2,D, D1};
const static double *BASIS_45A[]  = {A2,B, D };
const static double *BASIS_45B[]  = {A2,C, D };
const static double *BASIS_56A[]  = {A2,B, C, D1};
const static double *BASIS_56B[]  = {A2,D, E1,E2};
const static double *BASIS_57A[]  = {A2,B, E1,E2};
const static double *BASIS_57B[]  = {A2,B, F1,F2};
const static double *BASIS_57C[]  = {A2,B, G1,G2};
const static double *BASIS_511A[] = {A2,D1,E1,E2};
const static double *BASIS_511B[] = {A2,D1,F1,F2};
const static double *BASIS_511C[] = {A2,D1,G1,G2};
const static double *BASIS_516[]  = {A2,D, G1,G2};
const static double *BASIS_66[]   = {A2,B, C, D, D1};
const static double *BASIS_67A[]  = {A2,B, D, E1,E2};
const static double *BASIS_67B[]  = {A2,C, D, E1,E2};
const static double *BASIS_68A[]  = {A2,D, D1,E1,E2};
const static double *BASIS_68B[]  = {A2,D, D1,G1,G2};
const static double *BASIS_617A[] = {A2,B, D, G1,G2};
const static double *BASIS_617B[] = {A2,C, D, G1,G2};
const static double *BASIS_88[]   = {A2,D, D1,E1,E2,F1,F2};
const static double *BASIS_810A[] = {A2,B, C, D, D1,E1,E2};
const static double *BASIS_810B[] = {A2,B, C, D, D1,G1,G2};
const static double *BASIS_816[]  = {A2,D, D1,E1,E2,G1,G2};
const static double *BASIS_817[]  = {A2,B, D, E1,E2,G1,G2};
const static double *BASIS_818[]  = {A2,B, D, E1,E2,F1,F2};
const static double *BASIS_920A[] = {A2,B, C, D1,E1,E2,F1,F2};
const static double *BASIS_920B[] = {A2,B, C, D1,F1,F2,G1,G2};
const static double *BASIS_1012[] = {A2,B, C, D, D1,E1,E2,F1,F2};
const static double *BASIS_1034[] = {A2,B, C, D, D1,E1,E2,G1,G2};
const static double *BASIS_1212[] = {A2,B, C, D, D1,E1,E2,F1,F2,G1,G2};

const static int NUM_LM_MODELS = 37;
const static double **BASES[] = 
            {BASIS_11,  BASIS_22B, BASIS_33A, BASIS_33B, BASIS_33C,
	     BASIS_34,  BASIS_44A, BASIS_44B, BASIS_45A, BASIS_45B,
	     BASIS_56A, BASIS_56B, BASIS_57A, BASIS_57B, BASIS_57C,
	     BASIS_511A,BASIS_511B,BASIS_511C,BASIS_516, BASIS_66,
	     BASIS_67A, BASIS_67B, BASIS_68A, BASIS_68B, BASIS_617A,
	     BASIS_617B,BASIS_88,  BASIS_810A,BASIS_810B,BASIS_816,
	     BASIS_817, BASIS_818, BASIS_920A,BASIS_920B,BASIS_1012,
	     BASIS_1034,BASIS_1212};
const static string MODEL_NAMES[] = 
            { "1.1",  "2.2b", "3.3a", "3.3b",  "3.3c",
	      "3.4",  "4.4a", "4.4b", "4.5a",  "4.5b",
	      "5.6a", "5.6b", "5.7a", "5.7b",  "5.7c",
	      "5.11a","5.11b","5.11c","5.16",  "6.6",
	      "6.7a", "6.7b", "6.8a", "6.8b",  "6.17a",
	      "6.17b","8.8",  "8.10a","8.10b", "8.16",
	      "8.17", "8.18", "9.20a","9.20b","10.12",
	     "10.34","12.12"};
const static int MODEL_PARAMS[] = 
             {0,1,2,2,2,2,3,3,3,3,4,4,4,4,4,4,4,4,4,
              5,5,5,5,5,5,5,7,7,7,7,7,7,8,8,9,9,11};
const static bool TIME_REVERSIBLE[] = 
             {true, true, true, false,true,
	      true, true, true, false,false,
	      false,false,false,false,false,
	      false,false,false,false,false,
	      false,false,false,false,false,
	      false,false,false,false,false,
	      false,false,false,false,false,
	      false,false};
/*
 * symmetry 3 (empty string) is only for models with full symmetry
 * (RY, WS, MK models are isomorphic). FULL_SYMMETRY identifies
 * these models (1.1, 3.3a, 4.4a, 6.7a, 9.20b, 12.12).
 * (this is cosmetic - names of these models don't have RY, WS or MK appended.)
 */
const static string SYMMETRY[] = {"RY","WS","MK",""};
const static int SYMMETRY_PERM[][12] = 
                   {{1,0,2,6,7,8,3,4,5,9,11,10},
		    {7,6,8,4,3,5,1,0,2,11,10,9},
		    {0,1,2,3,4,5,6,7,8,9,10,11},
		    {0,1,2,3,4,5,6,7,8,9,10,11}};
const static bool FULL_SYMMETRY[] = 
             {true, false,true, false,false,
	      false,true, false,false,false,
	      false,false,false,false,false,
	      false,false,false,false,false,
	      true, false,false,false,false,
	      false,false,false,false,false,
	      false,false,false,true, false,
	      false,true};
const static int NUM_RATES = 12;

ModelLieMarkov::ModelLieMarkov(string model_name, PhyloTree *tree, string model_params, bool count_rates)
	: ModelNonRev(tree)
{
    assert(NUM_RATES==getNumRateEntries());
    parseModelName(model_name,&model_num,&symmetry);
    if (model_num<0) {
        // should never happen - model_name should have been accepted 
        // by validModelName before constructor was called.
        cerr << "Bad model name in ModelLieMarkov constructor" << endl;
        exit(1);
    }
    basis = BASES[model_num];
    num_params = MODEL_PARAMS[model_num];
    model_parameters = new double [num_params];
    for (int i=0; i< num_params; i++) model_parameters[i] = 0;
    this->setRates();
	/*
	 * I'm not sure how to correctly handle count_rates, so for now I'm just
	 * avoiding the problem. Actual IQTree programmers can fix this.
	 * Whatever happens should leave model_parameters[] and rates[]
	 * consistent with each other.
	 */
    if (count_rates)
        cerr << "WARNING: count_rates=TRUE not implemented in ModelLieMarkov constructor -- ignored" << endl;
	/* phylo_tree->aln->computeEmpiricalRateNonRev(rates); */
    if (model_params != "") {
        cerr << "WARNING: Supplying model params to constructor not yet properly implemented -- ignored" << endl;
        // TODO: parse model_params into model_parameters, then call setRates().
    }
    name = "LM"+MODEL_NAMES[model_num]+SYMMETRY[symmetry];
    full_name = "Lie Markov model "+MODEL_NAMES[model_num]+SYMMETRY[symmetry]+" (non reversible)";
}

ModelLieMarkov::~ModelLieMarkov() {
  // Do nothing, for now. model_parameters is reclaimed in ~ModelNonRev
}


/* static */ bool ModelLieMarkov::validModelName(string model_name) {
    int model_num, symmetry;
    parseModelName(model_name,&model_num,&symmetry);
    return (model_num!=-1);
}


/*
 * Model names are like LM3.3a or LM6.6WS.
 * All start with LM. They may end with RY, WS, MK or nothing.
 * In between, the name must be on the list MODEL_NAMES
 *
 * Returns number of entry on MODEL_NAMES in model_num (-1 if not found),
 * and symmetry is 0 for RY, 1 for WS, 2 for MK.
 * (If no RY, WS or MK at end of name, assume RY.)
 */

/* static */ void ModelLieMarkov::parseModelName(string model_name, int* model_num, int* symmetry) {
    *model_num = -1; // not found yet
    int len = model_name.length();
    string base_name;
    if (model_name.find("LM")==0) {
        // found "LM" at start of model name
        if (model_name.find("RY")==len-2) {
	    // found "RY" at end
	    *symmetry = 0;
            base_name = model_name.substr(2,len-4);
        } else if (model_name.find("WS")==len-2) {
	    // found "WS" at end
            *symmetry = 1;
            base_name = model_name.substr(2,len-4);
        } else if (model_name.find("MK")==len-2) {
	    // found "MK" at end
            *symmetry = 2;
            base_name = model_name.substr(2,len-4);
	} else {
	    // did not find RY, WS or MK, assume RY symmetry
  	    *symmetry = 0;
            base_name = model_name.substr(2,len-2);
	}
	// search for basename in MODEL_NAMES
	for (int i=0; i<NUM_LM_MODELS; i++) {
	    if (MODEL_NAMES[i].compare(base_name)==0) {
	        *model_num = i;
	        break;
	    }
	}
        // set full symmetry if have a fully symmetric model
        if (*model_num>=0 && FULL_SYMMETRY[*model_num]) *symmetry = 3;
    }
    return;
}

/*
 * Technically bounds are +/- 1, but on the boundaries there will be
 * mutation rates equal to zero, which may cause problems later.
 */
void ModelLieMarkov::setBounds(double *lower_bound, double *upper_bound, bool *bound_check) {
	int i, ndim = getNDim();

	for (i = 1; i <= ndim; i++) {
		lower_bound[i] = -0.9999;
		upper_bound[i] =  0.9999;
		bound_check[i] = false; // I've no idea what this does, so don't know if 'false' is appropriate - MDW.
	}
}

/*
 * Set rates from model_parameters
 */
void ModelLieMarkov::setRates() {
    memset(rates, 0, NUM_RATES*sizeof(double));
    double max_abs = 0;
    for (int param=0; param<num_params; param++) {
        // COMMENT: is this abs() or fabs()? abs is for int type, whereas fabs for double 
        max_abs = (fabs(model_parameters[param])>max_abs ? fabs(model_parameters[param]) : max_abs);
        for (int rate=0; rate<NUM_RATES; rate++) 
            rates[rate] += model_parameters[param]*basis[param][SYMMETRY_PERM[symmetry][rate]];
    }
    double min_unnorm = DBL_MAX;
    for (int rate=0; rate<NUM_RATES; rate++) 
        min_unnorm = (rates[rate]<min_unnorm ? rates[rate] : min_unnorm);
    double norm = (max_abs==0 ? 0 : -max_abs/min_unnorm);
    for (int rate=0; rate<NUM_RATES; rate++) 
        rates[rate]=1+norm*rates[rate];
    if (verbose_mode >= VB_DEBUG) {
      cout << "LM setRates params = (";
      for (int param=0; param<num_params; param++) 
	cout << model_parameters[param] << ",";
      cout << ")\nrates = (";
      for (int rate=0; rate<NUM_RATES; rate++) 
	cout << rates[rate] << ",";
      cout << ")" << endl;
    }
}

