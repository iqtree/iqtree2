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
#ifdef USE_EIGEN3
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
using namespace Eigen;
#endif
#include "modelliemarkov.h"
#include <float.h>
#undef NDEBUG
#include <assert.h>
#include <complex>
/*
 * TO DO: It is inconvenient to have all of these const declarations
 * here, as it is so much to scroll past to find the actual code.
 * Find some way to be able to put these at end of .cpp file,
 * or maybe in .h file (not preferred, as nothing external to this
 * file needs or is allowed to see the definitions.)
 */

/*
 * Enum of the Lie-Markov basis matrices. D, E1 and E2 are the
 * ones which affect base frequencies and so get treated differently
 * when params.freq_type is not FREQ_ESTIMATE, so they are given the 
 * highest numbers.
 */

enum BASIS_MATRIX_TYPE {BM_A=0, BM_A2, BM_B, 
                        BM_C, BM_D1, BM_F1, 
                        BM_F2, BM_G1, BM_G2,
			BM_D, BM_E1, BM_E2};

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

/* A:
 * +3 +1 +1 +1
 * +1 +3 +1 +1
 * +1 +1 +3 +1
 * +1 +1 +1 +3
 */
const static double A[] = {+1,+1,+1,+1,+1,+1,+1,+1,+1,+1,+1,+1};

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
const static double  D[] = {+1,-1,-1,+1,-1,-1,+1,+1,-1,+1,+1,-1};
const static double mD[] = {-1,+1,+1,-1,+1,+1,-1,-1,+1,-1,-1,+1}; // -D
/* E1:
 * + + + +
 * - - - -
 * 0 0 0 0
 * 0 0 0 0
 */
const static double  E1[] = {-1, 0, 0,+1, 0, 0,+1,-1, 0,+1,-1, 0};
const static double mE1[] = {+1, 0, 0,-1, 0, 0,-1,+1, 0,-1,+1, 0}; // -E1
const static double tE1[] = {-2, 0, 0,+2, 0, 0,+2,-2, 0,+2,-2, 0}; // 2 E1

/* E2:
 * 0 0 0 0
 * 0 0 0 0
 * + + + +
 * - - - -
 */
const static double  E2[] = { 0,+1,-1, 0,+1,-1, 0, 0,-1, 0, 0,+1};
const static double mE2[] = { 0,-1,+1, 0,-1,+1, 0, 0,+1, 0, 0,-1}; // -E2
const static double tE2[] = { 0,+2,-2, 0,+2,-2, 0, 0,-2, 0, 0,+2}; // 2 E2

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

// same order as BASIS_MATRIX_TYPE enum.
const static double *LM_BASIS_MATRICES[] = {A,A2,B,C,D1,F1,F2,G1,G2,D,E1,E2};

/*
// never found a use for this. 
enum LM_MODEL = {LM11, LM22b, LM33a, LM33b, LM33c, LM34, LM44a, LM44b,
                 LM45a, LM45b, LM56a, LM56b, LM57a, LM57b, LM57c, 
                 LM511a, LM511b, LM511c, LM516, LM66, LM67a, LM57b,
		 LM68a, LM68b, LM617a, LM617b, LM88, LM810a, LM810b,
		 LM816, LM817, LM818, LM920a, LM920b, LM1012, LM1034, LM1212};
*/

// Lengths of these arrays (minus one) stored in MODEL_PARAMS. Note BM_D, BM_E1, BM_E2 must be at end
// (as these are treated differently in setBasis() for some freq_type values) (spacing emphasises this split)
const static BASIS_MATRIX_TYPE BASIS_11[]   = {BM_A};
const static BASIS_MATRIX_TYPE BASIS_22B[]  = {BM_A,BM_A2};
const static BASIS_MATRIX_TYPE BASIS_33A[]  = {BM_A,BM_A2,BM_B };
const static BASIS_MATRIX_TYPE BASIS_33B[]  = {BM_A,BM_A2,BM_C };
const static BASIS_MATRIX_TYPE BASIS_33C[]  = {BM_A,BM_A2,BM_D1};
const static BASIS_MATRIX_TYPE BASIS_34[]   = {BM_A,BM_A2,                                          BM_D};
const static BASIS_MATRIX_TYPE BASIS_44A[]  = {BM_A,                                                BM_D,BM_E1,BM_E2};
const static BASIS_MATRIX_TYPE BASIS_44B[]  = {BM_A,BM_A2,BM_D1,                                    BM_D};
const static BASIS_MATRIX_TYPE BASIS_45A[]  = {BM_A,BM_A2,BM_B,                                     BM_D};
const static BASIS_MATRIX_TYPE BASIS_45B[]  = {BM_A,BM_A2,BM_C,                                     BM_D};
const static BASIS_MATRIX_TYPE BASIS_56A[]  = {BM_A,BM_A2,BM_B, BM_C, BM_D1};
const static BASIS_MATRIX_TYPE BASIS_56B[]  = {BM_A,BM_A2,                                          BM_D,BM_E1,BM_E2};
const static BASIS_MATRIX_TYPE BASIS_57A[]  = {BM_A,BM_A2,BM_B,                                          BM_E1,BM_E2};
const static BASIS_MATRIX_TYPE BASIS_57B[]  = {BM_A,BM_A2,BM_B, BM_F1,BM_F2};
const static BASIS_MATRIX_TYPE BASIS_57C[]  = {BM_A,BM_A2,BM_B, BM_G1,BM_G2};
const static BASIS_MATRIX_TYPE BASIS_511A[] = {BM_A,BM_A2,BM_D1,                                         BM_E1,BM_E2};
const static BASIS_MATRIX_TYPE BASIS_511B[] = {BM_A,BM_A2,BM_D1,BM_F1,BM_F2};
const static BASIS_MATRIX_TYPE BASIS_511C[] = {BM_A,BM_A2,BM_D1,BM_G1,BM_G2};
const static BASIS_MATRIX_TYPE BASIS_516[]  = {BM_A,BM_A2,BM_G1,BM_G2,                              BM_D};
const static BASIS_MATRIX_TYPE BASIS_66[]   = {BM_A,BM_A2,BM_B, BM_C, BM_D1,                        BM_D};
const static BASIS_MATRIX_TYPE BASIS_67A[]  = {BM_A,BM_A2,BM_B,                                     BM_D,BM_E1,BM_E2};
const static BASIS_MATRIX_TYPE BASIS_67B[]  = {BM_A,BM_A2,BM_C,                                     BM_D,BM_E1,BM_E2};
const static BASIS_MATRIX_TYPE BASIS_68A[]  = {BM_A,BM_A2,BM_D1,                                    BM_D,BM_E1,BM_E2};
const static BASIS_MATRIX_TYPE BASIS_68B[]  = {BM_A,BM_A2,BM_D1,BM_G1,BM_G2,                        BM_D};
const static BASIS_MATRIX_TYPE BASIS_617A[] = {BM_A,BM_A2,BM_B, BM_G1,BM_G2,                        BM_D};
const static BASIS_MATRIX_TYPE BASIS_617B[] = {BM_A,BM_A2,BM_C, BM_G1,BM_G2,                        BM_D};
const static BASIS_MATRIX_TYPE BASIS_88[]   = {BM_A,BM_A2,BM_D1,BM_F1,BM_F2,                        BM_D,BM_E1,BM_E2};
const static BASIS_MATRIX_TYPE BASIS_810A[] = {BM_A,BM_A2,BM_B, BM_C, BM_D1,                        BM_D,BM_E1,BM_E2};
const static BASIS_MATRIX_TYPE BASIS_810B[] = {BM_A,BM_A2,BM_B, BM_C, BM_D1,BM_G1,BM_G2,            BM_D};
const static BASIS_MATRIX_TYPE BASIS_816[]  = {BM_A,BM_A2,BM_D1,BM_G1,BM_G2,                        BM_D,BM_E1,BM_E2};
const static BASIS_MATRIX_TYPE BASIS_817[]  = {BM_A,BM_A2,BM_B, BM_G1,BM_G2,                        BM_D,BM_E1,BM_E2};
const static BASIS_MATRIX_TYPE BASIS_818[]  = {BM_A,BM_A2,BM_B, BM_F1,BM_F2,                        BM_D,BM_E1,BM_E2};
const static BASIS_MATRIX_TYPE BASIS_920A[] = {BM_A,BM_A2,BM_B, BM_C, BM_D1,BM_F1,BM_F2,                 BM_E1,BM_E2};
const static BASIS_MATRIX_TYPE BASIS_920B[] = {BM_A,BM_A2,BM_B, BM_C, BM_D1,BM_F1,BM_F2,BM_G1,BM_G2};
const static BASIS_MATRIX_TYPE BASIS_1012[] = {BM_A,BM_A2,BM_B, BM_C, BM_D1,BM_F1,BM_F2,            BM_D,BM_E1,BM_E2};
const static BASIS_MATRIX_TYPE BASIS_1034[] = {BM_A,BM_A2,BM_B, BM_C, BM_D1,BM_G1,BM_G2,            BM_D,BM_E1,BM_E2};
const static BASIS_MATRIX_TYPE BASIS_1212[] = {BM_A,BM_A2,BM_B, BM_C, BM_D1,BM_F1,BM_F2,BM_G1,BM_G2,BM_D,BM_E1,BM_E2};

const static int NUM_LM_MODELS = 37;
const static BASIS_MATRIX_TYPE *BASES[] = 
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
 * Base frequency Degrees of Freedom, by model. This is the number
 * of matrices out of D, E1, E2 in the model. 
 * BDF=0 => equilibrium base frequencies are pi_A=pi_G=pi_C=pi_T = 1/4
 * BDF=1 => eqbm freqs pi_A=pi_G, pi_C=pi_T (for RY pairing)
 * BDF=2 => eqbm freqs pi_A+pi_G = pi_C+pi_T = 1/2
 * BDF=3 => arbitrary eqbm freqs
 */

const static int BDF[] = 
  {0,0,0,0,0, // 1.1,   2.2b,  3.3a,  3.3b,  3.3c
   1,3,1,1,1, // 3.4,   4.4a,  4.4b,  4.5a,  4.5b
   0,3,2,0,0, // 5.6a,  5.6b,  5.7a,  5.7b,  5.7c
   2,0,0,1,1, // 5.11a, 5.11b, 5.11c, 5.16,  6.6
   3,3,3,1,1, // 6.7a,  6.7b,  6.8a,  6.8b,  6.17a
   1,3,3,1,3, // 6.17b, 8.8,   8.10a, 8.10b, 8.16
   3,3,2,0,3, // 8.17,  8.18,  9.20a, 9.20b, 10.12
   3,3};      // 10.34, 12.12
/*
 * For the TRANSFORM_* arrays:
 * Each shows how to modify a basis matrix to enforce a fixed base
 * frequency vector. The base frequency vector is encoded as
 * tauRY = pi_A+pi_G-pi_C-pi_T
 * tauAG = pi_A-pi_G
 * tauCT = pi_C-pi_T
 * (or for WS symmetry think of them as tauWS, tauAT, tauCG, and
 * for MK symmetry think of them as tauMK, tauAC, tauGT)
 * then transformed basis matrix X = tauRY*TRANSFORM_X[0]+tauAG*TRANSFORM_X[1]+tauCT*TRANSFORM_X[2]
 */
const static double *TRANSFORM_A[]  = {D,    tE1,  tE2};
const static double *TRANSFORM_A2[] = {mD,    E1,   E2};
const static double *TRANSFORM_B[]  = {NULL, mE2,  mE1};
const static double *TRANSFORM_C[]  = {NULL,  E2,  mE1};
const static double *TRANSFORM_D1[] = {NULL,  E1,  mE2};
const static double *TRANSFORM_F1[] = {mE1,  NULL, NULL};
const static double *TRANSFORM_F2[] = {mE2,  NULL, NULL};
const static double *TRANSFORM_G1[] = {NULL,  mD,  NULL};
const static double *TRANSFORM_G2[] = {NULL, NULL,  mD};

const static double **BASIS_TRANSFORM[] = {
  TRANSFORM_A,  TRANSFORM_A2, TRANSFORM_B, 
  TRANSFORM_C,  TRANSFORM_D1, TRANSFORM_F1, 
  TRANSFORM_F2, TRANSFORM_G1, TRANSFORM_G2};

/*
 * symmetry 3 (empty string) is only for models with full symmetry
 * (RY, WS, MK models are isomorphic). FULL_SYMMETRY identifies
 * these models (1.1, 3.3a, 4.4a, 6.7a, 9.20b, 12.12).
 * (this is cosmetic - names of these models don't have RY, WS or MK appended.)
 */
const static string SYMMETRY[] = {"RY","WS","MK",""};
/* 
 * The definitions of tau (piToTau(), below) determine the unpermuted
 * base orderings for each symmetry. For RY, order = (A,G,C,T); for
 * WS order is (A,T,C,G). For MK, (A,C,G,T). Then the unpermuted
 * rates need to be permuted to the (A,C,G,T) order native to iqtree.
 */

const static int SYMMETRY_PERM[][12] = 
                   {{1,0,2,6,7,8,3,4,5,9,11,10}, // RY
		    {1,2,0,6,8,7,9,11,10,3,4,5}, // WS
		    {0,1,2,3,4,5,6,7,8,9,10,11}, // MK
		    {1,0,2,6,7,8,3,4,5,9,11,10}}; // sym=3 uses RY permutation
const static bool FULL_SYMMETRY[] = 
             {true, false,true, false,false, // 1.1,   2.2b,  3.3a,  3.3b,  3.3c
	      false,true, false,false,false, // 3.4,   4.4a,  4.4b,  4.5a,  4.5b
	      false,false,false,false,false, // 5.6a,  5.6b,  5.7a,  5.7b,  5.7c
	      false,false,false,false,false, // 5.11a, 5.11b, 5.11c, 5.16,  6.6
	      true, false,false,false,false, // 6.7a,  6.7b,  6.8a,  6.8b,  6.17a 
	      false,false,false,false,false, // 6.17b, 8.8,   8.10a, 8.10b, 8.16
	      false,false,false,true, false, // 8.17,  8.18,  9.20a, 9.20b, 10.12
	      false,true};                   // 10.34, 12.12
const static int NUM_RATES = 12;

const double MIN_LIE_WEIGHT = -0.9999;
const double MAX_LIE_WEIGHT =  0.9999;

ModelLieMarkov::ModelLieMarkov(string model_name, PhyloTree *tree, string model_params, StateFreqType freq_type, string freq_params, bool count_rates)
	: ModelNonRev(tree) {
        init(model_name.c_str(), model_params, freq_type, freq_params);
}

void ModelLieMarkov::init(const char *model_name, string model_params, StateFreqType freq, string freq_params)
{
    assert(NUM_RATES==getNumRateEntries());
    parseModelName(model_name,&model_num,&symmetry);
    name = "LM"+MODEL_NAMES[model_num]+SYMMETRY[symmetry];
    full_name = "Lie Markov model "+MODEL_NAMES[model_num]+SYMMETRY[symmetry]+" (non reversible)";
    if (model_num<0) {
        // should never happen - model_name should have been accepted 
        // by validModelName before constructor was called.
        cerr << "Bad model name in ModelLieMarkov constructor" << endl;
        abort();
    }
    freq_type = freq;
    if (freq_params != "") {
	readStateFreq(freq_params);
    }
    setBasis(); // sets basis and num_params

    if (model_parameters)
        delete[] model_parameters;
    model_parameters = new double [num_params];
    memset(model_parameters, 0, sizeof(double)*num_params);
    this->setRates();
	/*
	 * I'm not sure how to correctly handle count_rates, so for now I'm just
	 * avoiding the problem. Actual IQTree programmers can fix this.
	 * Whatever happens should leave model_parameters[] and rates[]
	 * consistent with each other.
	 */
     // Minh's answer: count_rates is not used anymore. This behaviour is correct!
//    if (count_rates)
//        cerr << "WARNING: count_rates=TRUE not implemented in ModelLieMarkov constructor -- ignored" << endl;
	/* phylo_tree->aln->computeEmpiricalRateNonRev(rates); */
    if (model_params != "") {
//        cerr << "WARNING: Supplying model params to constructor not yet properly implemented -- ignored" << endl;
        DoubleVector vec;
        convert_double_vec(model_params.c_str(), vec);
        if (vec.size() != num_params) 
            outError("String '"+ model_params + "' does not have exactly " + convertIntToString(num_params) + " parameters");
        for (int i = 0; i < num_params; i++) {
            if (vec[i] <= MIN_LIE_WEIGHT || vec[i] >= MAX_LIE_WEIGHT)
                outError("Weights for Lie Markov model must be between " + convertDoubleToString(MIN_LIE_WEIGHT) + " and " +
                    convertDoubleToString(MAX_LIE_WEIGHT));
            model_parameters[i] = vec[i];
            fixed_parameters = true;
        }
        setRates();
    }
    ModelNonRev::init(freq_type);
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
		lower_bound[i] = MIN_LIE_WEIGHT;
		upper_bound[i] = MAX_LIE_WEIGHT;
        // answer: this is correct, typically no bound check is performed
		bound_check[i] = false; // I've no idea what this does, so don't know if 'false' is appropriate - MDW.
	}
}


/*
 * tau[0] = pi_R-pi_Y = pi_A+pi_G-pi_C-pi_T (for RY pairing)
 *          (or pi_W-pi_S for WS pairing, or pi_M-pi_K for MK pairing)
 * tau[1] = pi_A-pi_G (RY pairing), pi_A-pi_T (WS pairing), pi_A-pi_C (MK pairing)
 * tau[2] = pi_C-pi_T (RY pairing), pi_C-pi_G (WS pairing), pi_G-pi_T (MK pairing)
 */
// Writes into a length 3 tau vector calculated from a given pi, for given pairing/symmetry
static void piToTau(double* pi, double* tau, int sym) {
  switch (sym) {
  case 0: // RY
  case 3: // full symmetry
    tau[0] = pi[0]+pi[2]-pi[1]-pi[3];
    tau[1] = pi[0]-pi[2];
    tau[2] = pi[1]-pi[3];
    break;
  case 1: // WS
    tau[0] = pi[0]+pi[3]-pi[1]-pi[2];
    tau[1] = pi[0]-pi[3];
    tau[2] = pi[1]-pi[2];
    break;
  case 2: // MK
    tau[0] = pi[0]+pi[1]-pi[2]-pi[3];
    tau[1] = pi[0]-pi[1];
    tau[2] = pi[2]-pi[3];
    break;
  default: outError("Can't happen");
  } // switch
}

// Writes into a length 4 pi vector calculated from given tau, for given pairing/symmetry.
static void tauToPi(double* tau, double* pi, int sym) {
  switch (sym) {
  case 0: // RY
  case 3: // Full symmetry
    // tau[0] = A+G-C-T, tau[1]=A-G, tau[2]=C-T
    pi[0] = 0.25 + 0.25*tau[0] + 0.5*tau[1]; // pi_A
    pi[1] = 0.25 - 0.25*tau[0] + 0.5*tau[2]; // pi_C
    pi[2] = 0.25 + 0.25*tau[0] - 0.5*tau[1]; // pi_G
    pi[3] = 0.25 - 0.25*tau[0] - 0.5*tau[2]; // pi_T
    break;
  case 1: // WS
    // tau[0] = A+T-C-G, tau[1]=A-T, tau[2]=C-G
    pi[0] = 0.25 + 0.25*tau[0] + 0.5*tau[1]; // pi_A
    pi[1] = 0.25 - 0.25*tau[0] + 0.5*tau[2]; // pi_C
    pi[2] = 0.25 - 0.25*tau[0] - 0.5*tau[2]; // pi_G
    pi[3] = 0.25 + 0.25*tau[0] - 0.5*tau[1]; // pi_T
    break;
  case 2: // MK
    // tau[0] = A+C-G-T, tau[1]=A-C, tau[2]=G-T
    pi[0] = 0.25 + 0.25*tau[0] + 0.5*tau[1]; // pi_A
    pi[1] = 0.25 + 0.25*tau[0] - 0.5*tau[1]; // pi_C
    pi[2] = 0.25 - 0.25*tau[0] + 0.5*tau[2]; // pi_G
    pi[3] = 0.25 - 0.25*tau[0] - 0.5*tau[2]; // pi_T
    break;
    default: outError("Can't happen");
  }
}

/**
 * Uses model_num, symmetry to populate 'basis' array.
 */

void ModelLieMarkov::setBasis() {
  if (getFreqType() == FREQ_UNKNOWN)
    switch (BDF[model_num]) {
    case 0:
        freq_type = FREQ_EQUAL;
        break;
    default:
        freq_type = FREQ_ESTIMATE;
        break;
    }

  if (getFreqType() == FREQ_EMPIRICAL || 
      getFreqType() == FREQ_USER_DEFINED || 
      getFreqType() == FREQ_EQUAL) {
    int bdf = BDF[model_num];
    // There are no free parameters for base frequencies:
    num_params = MODEL_PARAMS[model_num]-bdf;
    // This populates field state_freq. (TODO: this call might be redundant - check)
    if (bdf>0 && getFreqType() == FREQ_EQUAL) {
      outWarning("You have demanded equal base frequencies (-f q) for a Lie-Markov\nmodel which can produce unequal base frequencies. It is better to just\nselect a Lie-Markov model which can only produce equal base frequencies.\nThese models are 1.1, 2.2b, 3.3a, 3.3b, 3.3c, 5.6a, 5.7b, 5.7c, 5.11b, 5.11c, 9.20b.");
      bdf = 0; // need bdf=0 basis matrices.
    }

    init_state_freq(getFreqType());
    // state_freq is in order {pi_A, pi_C, pi_G, pi_T}
    double* tau = new double[3];
    piToTau(state_freq,tau,symmetry);
    
    // Now zero tau entries which BDF forces to be zero, and print warnings
    bool canMatchFreq = true;
    switch (bdf) {
    case 0:
      canMatchFreq = ((fabs(tau[0])<0.001 && fabs(tau[1])<0.001) || fabs(tau[2])<0.001);
      tau[0] = 0; tau[1] = 0; tau[2] = 0;
      break;
    case 1:
      canMatchFreq = (fabs(tau[1])<0.001 || fabs(tau[2])<0.001);
      tau[1] = 0; tau[2] = 0;
      break;
    case 2:
      canMatchFreq = (fabs(tau[0])<0.001);
      tau[0] = 0;
      break;
    case 3:
      break;
    default: outError("Can't happen");
    } // switch
    if (!canMatchFreq) {
      // MDW to Minh: I suspect there is a better way, please recode if there is.
      double* eqbm = new double[4];
      tauToPi(tau,eqbm,symmetry);
      char* buffer = new char[200];
      snprintf(buffer,200,"Model %s cannot achieve requested equilibrium base frequencies\n(%5.3f,%5.3f,%5.3f,%5.3f).\nInstead it will use equilibrium base frequencies (%5.3f,%5.3f,%5.3f,%5.3f).\n",
	       name.c_str(),state_freq[0],state_freq[1],state_freq[2],state_freq[3],eqbm[0],eqbm[1],eqbm[2],eqbm[3]);
      outWarning(buffer);
      delete[] eqbm;
      delete[] buffer;
    }

    basis = new double*[num_params+1];
    for (int i=0;i<=num_params;i++) {
      int basisIndex = BASES[model_num][i];
      double* unpermuted_rates = new double[NUM_RATES];
      memcpy(unpermuted_rates, LM_BASIS_MATRICES[basisIndex], NUM_RATES* sizeof(double));
      for (int tauIndex=0; tauIndex<3; tauIndex++) {
	const double* transformationMatrix = BASIS_TRANSFORM[basisIndex][tauIndex];
	if (tau[tauIndex]!=0 && transformationMatrix != NULL) {
	  for (int rate=0; rate<NUM_RATES; rate++) {
	    unpermuted_rates[rate] = unpermuted_rates[rate]+tau[tauIndex]*transformationMatrix[rate];
	  } // for rate
	} // if tau && !=NULL
      } // for tauIndex
      double* permuted_rates = new double[NUM_RATES];
      for (int rate=0; rate<NUM_RATES; rate++) {
	permuted_rates[rate] = unpermuted_rates[SYMMETRY_PERM[symmetry][rate]];
      }
      basis[i] = permuted_rates;
    } // for i
    delete[] tau;
  } else {
      assert(getFreqType() == FREQ_ESTIMATE); // only other legal possibility
    num_params = MODEL_PARAMS[model_num];
    basis = new double*[num_params+1];
    for (int i=0;i<=num_params;i++) {
      const double* unpermuted_rates = LM_BASIS_MATRICES[BASES[model_num][i]];
      double* permuted_rates = new double[NUM_RATES];
      for (int rate=0; rate<NUM_RATES; rate++) {
	permuted_rates[rate] = unpermuted_rates[SYMMETRY_PERM[symmetry][rate]];
      } // for rate
    basis[i] = permuted_rates;
    } // for i
  } // if getFreqType() ... else ...
}


/*
 * Set rates from model_parameters
 */
void ModelLieMarkov::setRates() {
    memset(rates, 0, NUM_RATES*sizeof(double));  // rates = 0
    double* aprime = basis[0]; // the only basis matrix with all offdiagonals non-negative, and trace non-zero
    double max_abs = 0;
    for (int param=0; param<num_params; param++) {
        // COMMENT: is this abs() or fabs()? abs is for int type, whereas fabs for double 
        max_abs = (fabs(model_parameters[param])>max_abs ? fabs(model_parameters[param]) : max_abs);
        for (int rate=0; rate<NUM_RATES; rate++) 
            rates[rate] += model_parameters[param]*basis[param+1][rate];
        // basis[0] is 'A' matrix which doesn't get a parameter.
    }
    double min_unnorm = DBL_MAX;
    for (int rate=0; rate<NUM_RATES; rate++) {
        double ratio = rates[rate]/aprime[rate]; 
        min_unnorm = (ratio<min_unnorm ? ratio : min_unnorm);
    }
    double norm = (max_abs==0 ? 0 : -max_abs/min_unnorm);
    for (int rate=0; rate<NUM_RATES; rate++) 
        rates[rate]=aprime[rate]+norm*rates[rate];
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

void ModelLieMarkov::decomposeRateMatrix() {
    ModelNonRev::decomposeRateMatrix();
    if (phylo_tree->params->matrix_exp_technique == MET_SCALING_SQUARING) 
        return;
    if (phylo_tree->params->matrix_exp_technique == MET_EIGEN3LIB_DECOMPOSITION) {
        // using Eigen library
        decomposeRateMatrixEigen3lib();
        return;
    }
    if (phylo_tree->params->matrix_exp_technique == MET_LIE_MARKOV_DECOMPOSITION) {
        decomposeRateMatrixClosedForm();
        return;
    }
}

void ModelLieMarkov::decomposeRateMatrixEigen3lib() {
#ifdef USE_EIGEN3
  nondiagonalizable = false; // until proven otherwise
    Matrix4d mat(rate_matrix);
    mat.transpose();
    EigenSolver<Matrix4d> eigensolver(mat);
    assert (eigensolver.info() == Eigen::Success);
    Map<Vector4cd,Aligned> eval(ceval);
    eval = eigensolver.eigenvalues();
    Map<Matrix4cd,Aligned> evec(cevec);
    evec = eigensolver.eigenvectors();
    if (abs(evec.determinant())<1e-10) {
      // limit of 1e-10 is something of a guess. 1e-12 was too restrictive.
      nondiagonalizable = true; // will use scaled squaring instead of eigendecomposition for matrix exponentiation
      return;
    }
    Map<Matrix4cd,Aligned> inv_evec(cinv_evec);
    inv_evec = evec.inverse();
//    cout << "det(evecs): " << setprecision(25) << evec.determinant() << endl;
//    int i, j;
//    for (i = 0; i < 4; i++) {
//        ceval[i] = eval(i);
//        for (j = 0; j < 4; j++) {
//            cevec[j*4+i] = evec(i, j);
//            cinv_evec[j*4+i] = inv_evec(i, j);
//        }
//    }
   Matrix4cd eval_diag = eval.asDiagonal();
//   cout << "eigenvalues:" << endl << eval_diag << endl;
//   cout << "columns right eigenvectors" << endl << evec << endl;    
//   cout << "row left eigenvectors" << endl << inv_evec << endl;    
//   cout << "rate_matrix: " << endl << mat << endl;
//   cout << "rate_matrix*eigenvectors:" << endl << mat*evec << endl;
//   cout << "eigenvectors*eigenvalues: " << endl << (evec*eval_diag) << endl;
//   cout << "diff: " << endl << (mat*evec - evec*eval_diag) << endl; 
//   cout << "check: " << endl << (inv_evec * mat * evec - eval_diag) << endl;
//    Matrix4cd cmat = mat;
//    cout << "check: " << endl << (evec * eval_diag * inv_evec) << endl;
    Matrix4cd check = inv_evec * mat * evec - eval_diag;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            assert(abs(check(i,j)) < 1e-4);
#else
    outError("Please install Eigen3 library for this option ", __func__);
#endif

}

const static int a2index[] = {-1, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
const static int  bindex[] = {-1,-1, 1,-1,-1,-1,-1,-1, 1,-1, 1,-1, 1, 1, 1,-1,-1,-1,-1, 1, 1,-1,-1,-1, 1,-1,-1, 1, 1 -1, 1, 1, 1, 1, 1, 1, 1};
const static int  cindex[] = {-1,-1,-1, 1,-1,-1,-1,-1,-1, 1, 2,-1,-1,-1,-1,-1,-1,-1,-1, 2,-1, 1,-1,-1,-1, 1,-1, 2, 2,-1,-1,-1, 2, 2, 2, 2, 2};
const static int  dindex[] = {-1,-1,-1,-1,-1, 1, 0, 1, 2, 2,-1, 1,-1,-1,-1,-1,-1,-1, 1, 3, 2, 2, 1, 1, 2, 2, 1, 3, 3, 1, 2, 2,-1,-1, 3, 3, 3};
const static int d1index[] = {-1,-1,-1,-1, 1,-1,-1, 2,-1,-1, 3,-1,-1,-1,-1, 1, 1, 1,-1, 4,-1,-1, 2, 2,-1,-1, 2, 4, 4, 2,-1,-1, 3, 3, 4, 4, 4};
const static int e1index[] = {-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1, 2, 2,-1,-1, 2,-1,-1,-1,-1, 3, 3, 3,-1,-1,-1, 3, 5,-1, 3, 3, 3, 4,-1, 5, 5, 5};
const static int e2index[] = {-1,-1,-1,-1,-1,-1, 2,-1,-1,-1,-1, 3, 3,-1,-1, 3,-1,-1,-1,-1, 4, 4, 4,-1,-1,-1, 4, 6,-1, 4, 4, 4, 5,-1, 6, 6, 6};
const static int f1index[] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 2,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1, 5,-1,-1,-1,-1, 5, 6, 4, 7,-1, 7};
const static int f2index[] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 3,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1, 6,-1,-1,-1,-1, 6, 7, 5, 8,-1, 8};
const static int g1index[] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 2,-1,-1, 2, 2,-1,-1,-1,-1, 3, 3, 3,-1,-1, 5, 5, 5,-1,-1, 6,-1, 7, 9};
const static int g2index[] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 3,-1,-1, 3, 3,-1,-1,-1,-1, 4, 4, 4,-1,-1, 6, 6, 6,-1,-1, 7,-1, 8,10};

void ModelLieMarkov::decomposeRateMatrixClosedForm() {
    // Lie Markov eigen decomposition with closed formula
    int i, j;
    double a = 1.0, a2 = 0, b = 0, c = 0, d = 0, d1 = 0, e1 = 0, e2 = 0, f1 = 0, f2 = 0, g1 = 0, g2 = 0;

    if (a2index[model_num] >= 0) a2 = model_parameters[a2index[model_num]];
    if ( bindex[model_num] >= 0)  b = model_parameters[ bindex[model_num]];
    if ( cindex[model_num] >= 0)  c = model_parameters[ cindex[model_num]];
    if ( dindex[model_num] >= 0)  d = model_parameters[ dindex[model_num]];
    if (d1index[model_num] >= 0) d1 = model_parameters[d1index[model_num]];
    if (e1index[model_num] >= 0) e1 = model_parameters[e1index[model_num]];
    if (e2index[model_num] >= 0) e2 = model_parameters[e2index[model_num]];
    if (f1index[model_num] >= 0) f1 = model_parameters[f1index[model_num]];
    if (f2index[model_num] >= 0) f2 = model_parameters[f2index[model_num]];
    if (g1index[model_num] >= 0) g1 = model_parameters[g1index[model_num]];
    if (g2index[model_num] >= 0) g2 = model_parameters[g2index[model_num]];
    
    // following code is from Christina
    
    
    if (name.substr(0,5) == "LM1.1" || name.substr(0,6) == "LM4.4a" || name.substr(0,6) == "LM5.6b") {

        /******** eigenvalues *********/
        //Eigenvalues = {0, -4*a, -4*a, -4*a}
        ceval[0] = 0.0; ceval[1] = ceval[2] = ceval[3] = -4.0*a;
        
        /******** right eigenvectors *********/
        //v0 = {1, 1, 1, 1}
        cevec[0] = cevec[4] = cevec[8] = cevec[12] = 1.0;
        //v1 = {-((a-a2-d)/(a-a2+d)), 1, -((a-a2-d)/(a-a2+d)), 1}
        cevec[1] = cevec[9] = -((a-a2-d)/(a-a2+d));
        cevec[5] = cevec[13] = 1.0;
        //v2 = {(2*e2)/(2*a+a2+e1+e2), -((2*a+a2+e1-e2)/(2*a+a2+e1+e2)), 0, 1}
        double temp = 1.0/(2*a+a2+e1+e2);
        cevec[2] = 2*e2*temp; cevec[6] = -(2*a+a2+e1-e2)*temp; cevec[10] = 0.0; cevec[14] = 1.0;
        //v3 = {-((2*a+a2-e1+e2)/(2*a+a2+e1+e2)), (2*e1)/(2*a+a2+e1+e2), 1, 0}
        cevec[3] = cevec[6]; cevec[7] = cevec[2]; cevec[11] = 1.0; cevec[15] = 0.0;

        /******** left eigenvectors *********/
        //v0 = {-((-2*std::pow(a,2.)+a*a2+std::pow(a2,2.)-2*a*d-a2*d-2*a*e1+2*a2*e1)/(2*std::pow(a,2.)-a*a2-std::pow(a2,2.)-2*a*d-a2*d-2*a*e2+2*a2*e2)), -1, 0, -1}
        cinv_evec[0] = -((-2*std::pow(a,2.)+a*a2+std::pow(a2,2.)-2*a*d-a2*d-2*a*e1+2*a2*e1)/(2*std::pow(a,2.)-a*a2-std::pow(a2,2.)-2*a*d-a2*d-2*a*e2+2*a2*e2));
        cinv_evec[4] =  -1; cinv_evec[8] = 0; cinv_evec[12] = -1;  
        //v1 = {-((-2*std::pow(a,2.)+a*a2+std::pow(a2,2.)+2*a*d+a2*d-2*a*e2+2*a2*e2)/(2*std::pow(a,2.)-a*a2-std::pow(a2,2.)-2*a*d-a2*d-2*a*e2+2*a2*e2)), 1, -1, 0}
        cinv_evec[1] = -((-2*std::pow(a,2.)+a*a2+std::pow(a2,2.)+2*a*d+a2*d-2*a*e2+2*a2*e2)/(2*std::pow(a,2.)-a*a2-std::pow(a2,2.)-2*a*d-a2*d-2*a*e2+2*a2*e2));
        cinv_evec[5] = 1; cinv_evec[9] = -1; cinv_evec[13] = 0;
        //v2 = {-((-2*std::pow(a,2.)+a*a2+std::pow(a2,2.)-2*a*d-a2*d+2*a*e1-2*a2*e1)/(2*std::pow(a,2.)-a*a2-std::pow(a2,2.)-2*a*d-a2*d-2*a*e2+2*a2*e2)), -1, 0, 1}
        cinv_evec[2] = -((-2*std::pow(a,2.)+a*a2+std::pow(a2,2.)-2*a*d-a2*d+2*a*e1-2*a2*e1)/(2*std::pow(a,2.)-a*a2-std::pow(a2,2.)-2*a*d-a2*d-2*a*e2+2*a2*e2));
        cinv_evec[6] = -1; cinv_evec[10] = 0; cinv_evec[14] = 1;
        //v3 = {1, 1, 1, 0}}
        cinv_evec[3] = 1; cinv_evec[7] = 1; cinv_evec[11] = 1; cinv_evec[15] = 0;
        
    } else if (name.substr(0,6) == "LM2.2b" || name.substr(0,5) == "LM3.4" || name.substr(0,6) == "LM5.16") {
        //5.16, 3.4, 2.2b, 1.1

        /******** eigenvalues *********/
        //Eigenvalues = {0, -4*(a-a2), -2*(2a+a2), -2(2a+a2)}
        ceval[0] = 0.0; ceval[1] = -4.0*(a-a2); ceval[2] = ceval[3] = -2*(2*a+a2);

        /******** right eigenvectors *********/
        //v0 = {1, 1, 1, 1}
        cevec[0] = cevec[4] = cevec[8] = cevec[12] = 1.0;
        //v1 = {-((6*(a*a2-std::pow(a2,2.)-a2*d))/(3*a*a2-3*std::pow(a2,2.)+3*a2*d+2*a*g2-2*a2*g2))+(3*a*a2-3*std::pow(a2,2.)-3*a2*d-2*a*g1+2*a2*g1)/(3*a*a2-3*std::pow(a2,2.)+3*a2*d+2*a*g2-2*a2*g2),
        //  -((-3*a*a2+3*std::pow(a2,2.)-3*a2*d+2*a*g2-2*a2*g2)/(3*a*a2-3*std::pow(a2,2.)+3*a2*d+2*a*g2-2*a2*g2)), 
        //  -((3*a*a2-3*std::pow(a2,2.)-3*a2*d-2*a*g1+2*a2*g1)/(3*a*a2-3*std::pow(a2,2.)+3*a2*d+2*a*g2-2*a2*g2)), 1}
        cevec[1] = -((6*(a*a2-std::pow(a2,2.)-a2*d))/(3*a*a2-3*std::pow(a2,2.)+3*a2*d+2*a*g2-2*a2*g2))+(3*a*a2-3*std::pow(a2,2.)-3*a2*d-2*a*g1+2*a2*g1)/(3*a*a2-3*std::pow(a2,2.)+3*a2*d+2*a*g2-2*a2*g2);
        cevec[5] = -((-3*a*a2+3*std::pow(a2,2.)-3*a2*d+2*a*g2-2*a2*g2)/(3*a*a2-3*std::pow(a2,2.)+3*a2*d+2*a*g2-2*a2*g2));
        cevec[9] = -((3*a*a2-3*std::pow(a2,2.)-3*a2*d-2*a*g1+2*a2*g1)/(3*a*a2-3*std::pow(a2,2.)+3*a2*d+2*a*g2-2*a2*g2));
        cevec[13] = 1.0;       
        //v2 = {0, -1, 0, 1}
        cevec[2] = cevec[10] = 0.0; cevec[6] = -1.0; cevec[14] = 1.0;
        //v3 = {-1, 0, 1, 0}
        cevec[3] = -1.0; cevec[7] = cevec[15] = 0.0; cevec[11] = 1.0;

        /******** left eigenvectors *********/
        //v0 = {-((-a+a2-d)/(a-a2-d)), -1, (2*g2)/(3*a2+g1-g2), -((3*a2-g1-g2)/(3*a2+g1-g2))}
        double temp = 1.0/(3*a2+g1-g2);
        cinv_evec[0] = -((-a+a2-d)/(a-a2-d)); cinv_evec[4] = -1.0; cinv_evec[8] = 2*g2*temp; cinv_evec[12] = -(3*a2-g1-g2)*temp;
        //v1 = {1, 1, -((3*a2+g1+g2)/(3*a2+g1-g2)), -((2*g1)/(3*a2+g1-g2))}
        cinv_evec[1] = cinv_evec[5] = 1.0; cinv_evec[9] = -(3*a2+g1+g2)*temp; cinv_evec[13] = -(2*g1)*temp;
        //v2 = {-((-a+a2-d)/(a-a2-d)), -1, 0, 1}
        cinv_evec[2] = -((-a+a2-d)/(a-a2-d)); cinv_evec[6] = -1.0; cinv_evec[10] = 0.0; cinv_evec[14] = 1.0; 
        //v3 = {1, 1, 1, 0}
        cinv_evec[3] = cinv_evec[7] = cinv_evec[11] = 1.0; cinv_evec[15] = 0.0;
    } else {
    
        double p0 = 4*a*a*a-3*a*a2*a2-a2*a2*a2+a*b*b+a2*b*b+a*c*c-a2*c*c-a*d1*d1+a2*d1*d1+2*a*f1*g1+a2*f1*g1-d1*f1*g1+b*f2*g1+c*f2*g1-b*f1*g2+c*f1*g2-2*a*f2*g2-a2*f2*g2-d1*f2*g2;
        double q0 = 4*a*a*d-4*a*a2*d+a2*a2*d-b*b*d+c*c*d-d*d1*d1-2*a*e1*g1-a2*e1*g1+d1*e1*g1-b*e2*g1-c*e2*g1+b*e1*g2-c*e1*g2+2*a*e2*g2+a2*e2*g2+d1*e1*g2;
        double r0 = 4*a*a*e1-2*a*a2*e1-2*a2*a2*e1-2*a*d1*e1+2*a2*d1*e1+2*a*b*e2-2*a2*b*e2+2*a*c*e2-2*a2*c*e2+4*a*d*f1+2*a2*d*f1-2*d*d1*f1+2*b*d*f2+2*c*d*f2+2*e2*f1*g2-2*e1*f2*g2;
        double s0 = 4*a*a*e2-2*a*a2*e2-2*a2*a2*e2+2*a*d1*e2-2*a2*d1*e2+2*a*b*e1-2*a2*b*e1-2*a*c*e1+2*a2*c*e1+4*a*d*f2+2*a2*d*f2+2*d*d1*f2+2*b*d*f1-2*c*d*f1+2*e2*f1*g1-2*e1*f2*g1;
        
        // eigenvector for 0 eigenvalue (base frequencies)
        ceval[0] = 0.0;
        double inv_p0 = 0.25/p0;
        double v0[] = {(p0+q0+r0)*inv_p0, (p0+q0-r0)*inv_p0, (p0-q0+s0)*inv_p0, (p0-q0-s0)*inv_p0};
        for (i = 0; i < 4; i++) {
            cinv_evec[i] = (state_freq[i] = v0[i]);
        }
        
        std::complex<double> complex1(1, sqrt(3));
        std::complex<double> complex2(1,-sqrt(3));

        std::complex<double> alpha=-36*a2*a2-12*b*b+12*c*c-12*d1*d1+24*f1*g1-24*f2*g2;
        std::complex<double> beta=432*a2*a2*a2-432*a2*b*b-432*a2*d1*d1-432*a2*f1*g1+432*d1*f1*g1-432*b*f2*g1-432*c*f2*g1+432*b*f1*g2-432*c*f1*g2+432*a2*f2*g2+432*d1*f2*g2;

        cout.unsetf(ios::fixed);
        cout.precision(10);

        cout << "alpha: " << alpha << endl;
        cout << "beta:  " << beta << endl;
        
//        assert(alpha != 0.0 || beta != 0.0);
        
        ceval[1] = -4.0*a-(std::pow(2.0,1/3.0))*alpha/(3.0*std::pow(beta+std::pow(4.0*alpha*alpha*alpha+beta*beta,1/2.),1/3.))+std::pow(beta+std::pow(4.0*alpha*alpha*alpha+beta*beta,1/2.),1/3.)/(3.0*std::pow(2.0,1/3.));
        ceval[2] = -4.0*a+complex1*alpha/((3.0*std::pow(2.0,2/3.))*(std::pow(beta+std::pow(4.0*alpha*alpha*alpha+beta*beta,1/2.),1/3.)))-complex2*std::pow(beta+std::pow(4.0*alpha*alpha*alpha+beta*beta,1/2.),1/3.)/(6.0*std::pow(2.0,1/3.));
        ceval[3] =-4.0*a+complex2*alpha/((3.0*std::pow(2.0,2/3.))*(std::pow(beta+std::pow(4.0*alpha*alpha*alpha+beta*beta,1/2.),1/3.)))-complex1*std::pow(beta+std::pow(4.0*alpha*alpha*alpha+beta*beta,1/2.),1/3.)/(6.0*std::pow(2.0,1/3.));


        for (i = 1; i < 4; i++) {
            std::complex<double> p, q, r, s, t;
            p = (8*a*a-2*f2*g2+a*ceval[i])*ceval[i];
            q = (-4*a*a2-4*a2*a2+8*a*d+4*a2*d+2*b*g1+2*c*g1-2*e1*g1+2*f1*g1-4*a*g2-2*a2*g2-2*d1*g2+2*e2*g2-a2*ceval[i]+d*ceval[i]-g2*ceval[i])*ceval[i];
            r = (-4*a*b+4*a2*b-4*a*c+4*a2*c+8*a*e1-2*a2*e1-2*d1*e1+2*b*e2+2*c*e2-4*a*f1-2*a2*f1+4*d*f1+2*d1*f1-2*b*f2-2*c*f2-4*f1*g2-b*ceval[i]-c*ceval[i]+e1*ceval[i]-f1*ceval[i])*ceval[i];
            s = (4*a*a2-2*a2*a2-2*b*b+2*c*c-4*a*d1+4*a2*d1-2*d1*d1+2*b*e1-2*c*e1+8*a*e2-2*a2*e2+2*d1*e2-2*b*f1+2*c*f1-4*a*f2-2*a*f2+4*d*f2-2*d1*f2+a2*ceval[i]-d1*ceval[i]+e2*ceval[i]-f2*ceval[i])*ceval[i];
            t = (16*a*a-12*a2*a2-4*b*b+4*c*c-4*d1*d1+8*f1*g1+8*a*ceval[i]+ceval[i]*ceval[i])*ceval[i];
            
            std::complex<double> vi[] = {4*v0[0]+p+q+r, 4*v0[1]+p+q-r, 4*v0[2]+p-q+s, 4*v0[3]+p-q-s+t};
            for (j = 0; j < 4; j++)
                cinv_evec[i*4+j] = vi[j];
        }
    }

    cout << "complex eigenvalues:";
    for (i = 0; i < 4; i++) 
        cout << " " << ceval[i];
    cout << endl;

    cout << "complex eigenvectors: " << endl;
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++)
            cout << " " << cevec[i*4+j];
        cout << endl;
    }

    cout << "complex inv_eigenvectors: " << endl;
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++)
            cout << " " << cinv_evec[i*4+j];
        cout << endl;
    }
    
	/* check eigenvalue equation */
	std::complex<double> zero;
    int error = 0;
	for (j = 0; j < num_states; j++) {
		for (i = 0, zero = 0.0; i < num_states; i++) {
			for (int k = 0; k < num_states; k++) 
                zero += rate_matrix[i*num_states+k] * cevec[k*num_states+j];
			zero -= ceval[j] * cevec[i*num_states+j];
			if (abs(zero) > 1.0e-5) {
                cout << "too large error[" << i << "," << j << "]: " << zero << endl;
				error = 1;
				break;
			}
		}
	}

	for (i = 0; i < num_states; i++) {
		for (j = 0, zero = 0.0; j < num_states; j++) {
			for (int k = 0; k < num_states; k++) 
                zero += cinv_evec[i*num_states+k] * rate_matrix[k*num_states+j];
			zero -= ceval[i] * cinv_evec[i*num_states+j];
			if (abs(zero) > 1.0e-5) {
                cout << "too large inv_error[" << i << "," << j << "]: " << zero << endl;
				error = 1;
				break;
			}
		}
	}
	if (error) {
		cerr << "\nERROR: Eigensystem doesn't satisfy eigenvalue equation!\n";
		cerr << "Rate matrix Q: " << endl;
		for (i = 0; i < num_states; i++) {
			for (j = 0; j < num_states; j++) cout << rate_matrix[i*num_states+j] << " ";
            cerr << endl;
		}
		cout << "State frequencies: " << endl;
		for (i = 0; i < num_states; i++) cout << state_freq[i] << " ";
		cout << endl;
	}

    
}

void ModelLieMarkov::computeTransMatrix(double time, double *trans_matrix) {
#ifdef USE_EIGEN3
  MatrixExpTechnique technique = phylo_tree->params->matrix_exp_technique;
  if (technique == MET_SCALING_SQUARING || (technique == MET_EIGEN3LIB_DECOMPOSITION && nondiagonalizable)) {
        Matrix4d A = Map<Matrix4d>(rate_matrix);
        A = (A.transpose() * time).exp();
        Map<Matrix4d> P(trans_matrix);
        P = A.transpose();

        int i, j;
        for (i = 0; i < 4; i++) {
            double sum = 0.0;
            for (j = 0; j < 4; j++)
                sum += (trans_matrix[i*4+j]);
            assert(fabs(sum-1.0) < 1e-4);
        }
    } else if (technique == MET_EIGEN3LIB_DECOMPOSITION) {
    // and nondiagonalizable == false, else we used scaled squaring
        int i;
        Vector4cd ceval_exp;
        for (i = 0; i < 4; i++)
            ceval_exp(i) = exp(ceval[i]*time);
        Matrix4cd cevectors(cevec);
        Matrix4cd cinv_evectors(cinv_evec);
        Matrix4cd res = cevectors * ceval_exp.asDiagonal() * cinv_evectors;
	// if assertions fail, it may be due to cevec having near-zero
	// determinant, and a fix could be to relax the test for
	// nondiagonalizable in ModelLieMarkov::decomposeRateMatrixEigen3lib()
        for (i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                trans_matrix[i*4+j] = res(j, i).real();
                assert(fabs(res(j,i).imag()) < 1e-6);
                assert(trans_matrix[i*4+j] >= -0.000001);
                assert(trans_matrix[i*4+j] <=  1.000001);
                if (trans_matrix[i*4+j] < 0)
                    trans_matrix[i*4+j] = 0.0;
                if (trans_matrix[i*4+j] > 1)
                    trans_matrix[i*4+j] = 1.0;
            }
            assert(fabs(trans_matrix[i*4]+trans_matrix[i*4+1]+trans_matrix[i*4+2]+trans_matrix[i*4+3]-1.0) < 1e-4);
        }        
    } else
        ModelNonRev::computeTransMatrix(time, trans_matrix);

#else
    ModelNonRev::computeTransMatrix(time, trans_matrix);
#endif
}

