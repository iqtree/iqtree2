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
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
using namespace Eigen;
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

const static int NUM_LM_MODELS = 38;
// Note: really just 37 models, 38th is to provide "StrSym" as an alias for WS6.6
const static int STR_SYM_INDEX = 37; // entry 37 in BASES, MODEL_NAMES etc is strand symmetric model.
const static BASIS_MATRIX_TYPE *BASES[] = 
            {BASIS_11,  BASIS_22B, BASIS_33A, BASIS_33B, BASIS_33C,
	     BASIS_34,  BASIS_44A, BASIS_44B, BASIS_45A, BASIS_45B,
	     BASIS_56A, BASIS_56B, BASIS_57A, BASIS_57B, BASIS_57C,
	     BASIS_511A,BASIS_511B,BASIS_511C,BASIS_516, BASIS_66,
	     BASIS_67A, BASIS_67B, BASIS_68A, BASIS_68B, BASIS_617A,
	     BASIS_617B,BASIS_88,  BASIS_810A,BASIS_810B,BASIS_816,
	     BASIS_817, BASIS_818, BASIS_920A,BASIS_920B,BASIS_1012,
	     BASIS_1034,BASIS_1212,BASIS_66};
const static string MODEL_NAMES[] = 
            { "1.1",  "2.2b", "3.3a", "3.3b",  "3.3c",
	      "3.4",  "4.4a", "4.4b", "4.5a",  "4.5b",
	      "5.6a", "5.6b", "5.7a", "5.7b",  "5.7c",
	      "5.11a","5.11b","5.11c","5.16",  "6.6",
	      "6.7a", "6.7b", "6.8a", "6.8b",  "6.17a",
	      "6.17b","8.8",  "8.10a","8.10b", "8.16",
	      "8.17", "8.18", "9.20a","9.20b","10.12",
	      "10.34","12.12","strsym"};
const static int MODEL_PARAMS[] = 
             {0,1,2,2,2,2,3,3,3,3,4,4,4,4,4,4,4,4,4,
              5,5,5,5,5,5,5,7,7,7,7,7,7,8,8,9,9,11,5};
const static bool TIME_REVERSIBLE[] = 
             {true, true, true, false,true,
	      true, true, true, false,false,
	      false,false,false,false,false,
	      false,false,false,false,false,
	      false,false,false,false,false,
	      false,false,false,false,false,
	      false,false,false,false,false,
	      false,false,false};
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
   3,3,1};    // 10.34, 12.12, strsym
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
	      false,true, false};            // 10.34, 12.12, StrSym
const static int NUM_RATES = 12;

/*
 * In principle, Lie Markov parameters range in [-1,1]. However values
 * very near the boundary are biologically implausible (will put at least
 * one entry in the rate matrix very close to zero) and can cause 
 * numerical problems, so we restrict the range a bit.
 */

const double MIN_LIE_WEIGHT = -0.98;
const double MAX_LIE_WEIGHT =  0.98;

ModelLieMarkov::ModelLieMarkov(string model_name, PhyloTree *tree, string model_params, StateFreqType freq_type, string freq_params)
	: ModelMarkov(tree, false) {
  init(model_name.c_str(), model_params, freq_type, freq_params);
        
        // show warning if the user is running AliSim without inference mode but has not yet specified model parameters
        if (Params::getInstance().alisim_active && !Params::getInstance().alisim_inference_mode && model_params.length() == 0 && getNParams()>0)
            outWarning("Without Inference Mode, we strongly recommend users to specify model parameters for more accuracy simulations. Users could use <Model_Name>{<param_0>/.../<param_n>} to specify the model parameters. For the model "+model_name+", users should provide "+convertIntToString(getNParams())+" params (see User Manuals).");
}

void ModelLieMarkov::init(const char *model_name, string model_params, StateFreqType freq, string freq_params)
{
    // TODO: why is freq_params not handled here?

    ASSERT(NUM_RATES==getNumRateEntries());
    StateFreqType expected_freq_type; // returned by getLieMarkovModelInfo but not used here
    getLieMarkovModelInfo((string)model_name, name, full_name, model_num, symmetry, expected_freq_type);

    if (model_num<0) {
        // should never happen - model_name should have been accepted 
        // by validModelName before constructor was called.
        cerr << "Bad model name in ModelLieMarkov constructor" << endl;
        abort();
    }

    setBasis(); // sets basis and num_params

//    if (model_parameters)
//        delete[] model_parameters;
    model_parameters = new double [num_params];
    memset(model_parameters, 0, sizeof(double)*num_params);
    this->setRates();
    // param optfromgiven only has effect if model_params != ""
    if (model_params != "") {
        DoubleVector vec;
        
        // detect the seperator
        char separator = ',';
        if (model_params.find('/') != std::string::npos)
            separator = '/';
        
        convert_double_vec_with_distributions(model_params.c_str(), vec, separator);
        if (vec.size() != num_params) 
            outError("String '"+ model_params + "' does not have exactly " + convertIntToString(num_params) + " parameters");
        for (int i = 0; i < num_params; i++) {
            if (vec[i] <= MIN_LIE_WEIGHT || vec[i] >= MAX_LIE_WEIGHT)
                outError("Weights for Lie Markov model must be between " + convertDoubleToString(MIN_LIE_WEIGHT) + " and " +
                    convertDoubleToString(MAX_LIE_WEIGHT));
            model_parameters[i] = vec[i];
            fixed_parameters = !Params::getInstance().optimize_from_given_params;
        }
        setRates();
    }

    if (freq_type == FREQ_UNKNOWN || expected_freq_type == FREQ_EQUAL) freq_type = expected_freq_type;
    
    // read user-specified frequencies
    if (freq_params.length() > 0)
    {
        // ignore user-specified frequencies if the model has equal state frequencies
        if (expected_freq_type == FREQ_EQUAL)
            outWarning("The model "+(string)model_name+" has equal state frequencies. Therefore, user-specified frequencies will be ignored.");
        // otherwise, read user-specified frequencies
        else
            readFreqs(expected_freq_type, freq_params);
    }
    
    ModelMarkov::init(freq_type);
    
    // initialize random state frequencies if AliSim is running without inference mode
    if (Params::getInstance().alisim_active && !Params::getInstance().alisim_inference_mode && (freq_type == FREQ_ESTIMATE || freq_type == FREQ_EMPIRICAL)){
        // initializing state_freqs from expected_freq_type
        initStateFreqsAliSim(expected_freq_type);
    }
}

/**
     initialize random state frequencies when running AliSim without inference mode
*/
void ModelLieMarkov::initStateFreqsAliSim(StateFreqType expected_freq_type)
{
    switch (expected_freq_type) {
        case FREQ_ESTIMATE:
        case FREQ_EMPIRICAL:
        {
            random_frequencies_from_distributions(state_freq);
            break;
        }
        case FREQ_DNA_1212:
        case FREQ_DNA_1221:
        case FREQ_DNA_1122:
        {
            // randomly generate an input frequency which is less than 0.5
            int num_freqs = 1;
            double* freqs = new double[num_freqs];
            freqs[0] = random_number_from_distribution_with_upperbound("uniform", 0.5);
            
            // set state freqs
            mappingFreqs(expected_freq_type, freqs);
            
            // delete freqs
            delete[] freqs;
            
            break;
        }
        case FREQ_DNA_RY:
        case FREQ_DNA_WS:
        case FREQ_DNA_MK:
        {
            // randomly generate two pairs of frequencies
            int num_freqs = 2;
            double* freqs = new double[num_freqs];
            for (int i = 0; i < num_freqs; i++)
                freqs[i] = random_number_from_distribution_with_upperbound("uniform", 0.5);
            
            // set state freqs
            mappingFreqs(expected_freq_type, freqs);
            
            // delete freqs
            delete[] freqs;
            
            break;
        }
        default:
            break;
    }
}

/**
     mapping state frequencies from user-specified/random frequencies
*/
void ModelLieMarkov::mappingFreqs(StateFreqType expected_freq_type, double *freqs)
{
    switch (expected_freq_type) {
        case FREQ_DNA_1212:
        case FREQ_DNA_1221:
        case FREQ_DNA_1122:
        {
            // validate the input freqs[0], it must be less than 0.5
            if (freqs[0] >= 0.5)
                outError("The input base frequency must be less than 0.5. Please check and try again!");
                
            // set state freqs
            if (expected_freq_type == FREQ_DNA_1212)
            {
                state_freq[0] = state_freq[2] = freqs[0];
                state_freq[1] = state_freq[3] = 0.5 - freqs[0];
            }
            else if (expected_freq_type == FREQ_DNA_1221)
            {
                state_freq[0] = state_freq[3] = freqs[0];
                state_freq[1] = state_freq[2] = 0.5 - freqs[0];
            } else if (expected_freq_type == FREQ_DNA_1122)
            {
                state_freq[0] = state_freq[1] = freqs[0];
                state_freq[2] = state_freq[3] = 0.5 - freqs[0];
            }
            
            break;
        }
        case FREQ_DNA_RY:
        case FREQ_DNA_WS:
        case FREQ_DNA_MK:
        {
            // validate the input freqs[0], and freqs[1], they must be less than 0.5
            if (freqs[0] >= 0.5 || freqs[1] >= 0.5)
                outError("The input base frequencies must be less than 0.5. Please check and try again!");
            
            // set state freqs
            if (expected_freq_type == FREQ_DNA_RY)
            {
                state_freq[0] = freqs[0];
                state_freq[1] = freqs[1];
                state_freq[2] = 0.5 - freqs[0];
                state_freq[3] = 0.5 - freqs[1];
            }
            else if (expected_freq_type == FREQ_DNA_WS)
            {
                state_freq[0] = freqs[0];
                state_freq[1] = freqs[1];
                state_freq[2] = 0.5 - freqs[1];
                state_freq[3] = 0.5 - freqs[0];
            } else if (expected_freq_type == FREQ_DNA_MK)
            {
                state_freq[0] = freqs[0];
                state_freq[1] = 0.5 - freqs[0];
                state_freq[2] = freqs[1];
                state_freq[3] = 0.5 - freqs[1];
            }
            break;
        }
        default:
            break;
    }
}

/**
     read user-specified state frequencies
*/
void ModelLieMarkov::readFreqs(StateFreqType expected_freq_type, string freq_params)
{
    // detect the seperator
    char separator = ',';
    if (freq_params.find('/') != std::string::npos)
        separator = '/';
    
    switch (expected_freq_type) {
        case FREQ_ESTIMATE:
        case FREQ_EMPIRICAL:
        {
            // extract/generate freqs one by one
            convert_double_array_with_distributions(freq_params, state_freq, num_states, separator);
            
            // normalize state freqs
            normalize_frequencies(state_freq, num_states, -1, true);
            
            break;
        }
        case FREQ_DNA_1212:
        case FREQ_DNA_1221:
        case FREQ_DNA_1122:
        {
            // extract/generate input base frequency
            int num_freqs = 1;
            double* freqs = new double[num_freqs];
            freqs[0] = convert_double_with_distribution_and_upperbound(freq_params, 0.5);
            
            // set state freqs
            mappingFreqs(expected_freq_type, freqs);
            
            // delete freqs
            delete[] freqs;
            
            break;
        }
        case FREQ_DNA_RY:
        case FREQ_DNA_WS:
        case FREQ_DNA_MK:
        {
            // extract/generate two freqs
            int num_freqs = 2;
            double* freqs = new double[num_freqs];
            // validate the number of items
            size_t num_separators = std::count(freq_params.begin(), freq_params.end(), separator);
            if (num_separators != num_freqs - 1)
                outError("The number of frequencies in "+freq_params+" is "+convertIntToString(num_separators+1)+", which is different from the expected number of frequencies ("+convertIntToString(num_freqs)+") for this model. Please check and try again!");

            // extract/generate double numbers one by one
            for (int i = 0; i < num_freqs; i++) {
                // extract sub-string by separator
                size_t pos = freq_params.find(separator);
                string token = freq_params.substr(0, pos);
                
                // convert/generate a double
                freqs[i] = convert_double_with_distribution_and_upperbound(token, 0.5);
                
                // remove the current double/distribution name from tmp_str
                freq_params.erase(0, pos + 1);
            }
            
            // set state freqs
            mappingFreqs(expected_freq_type, freqs);
            
            // delete freqs
            delete[] freqs;
            
            break;
        }
        default:
            break;
    }
    
    // update freq_type
    freq_type = FREQ_USER_DEFINED;
}

// Note to Minh: I see ModelUnrest also lacks checkpointing.
// I think this code could be copied straight over.
// If modelDNA had a setRates() method,
// we could make virtual setRates in ModelMarkov and
// perhaps move this code all into there. - MDW
void ModelLieMarkov::startCheckpoint() {
    checkpoint->startStruct("ModelLieMarkov" + name);
}

void ModelLieMarkov::saveCheckpoint() {
    // saves model_parameters
    startCheckpoint();
    if (num_params > 0)
        CKP_ARRAY_SAVE(num_params, model_parameters);
    endCheckpoint();
    ModelMarkov::saveCheckpoint();
}

void ModelLieMarkov::restoreCheckpoint() {
    ModelMarkov::restoreCheckpoint();
    // restores model_parameters
    startCheckpoint();
    if (num_params > 0)
        CKP_ARRAY_RESTORE(num_params, model_parameters);
    endCheckpoint();
    setRates();                        // updates rate matrix
    decomposeRateMatrix();             // updates eigen system.
    if (phylo_tree)
        phylo_tree->clearAllPartialLH();
}

void ModelLieMarkov::writeInfo(ostream &out) {
    int i;
    out << "Model parameters: ";
    if (num_params>0) out << model_parameters[0];
    for (i=1; i < num_params; i++) out << "," << model_parameters[i];
    out << endl;
}

/*static*/ void ModelLieMarkov::getLieMarkovModelInfo(string model_name, string &name, string &full_name, int &model_num, int &symmetry, StateFreqType &def_freq) {
    parseModelName(model_name,&model_num,&symmetry);
    // Special case, just because it is confusing
    if (model_name == "2.2a" || model_name == "RY2.2a" ||
        model_name == "WS2.2a" || model_name == "MK2.2a") {
      cerr << "Model 2.2a does not exist, do you mean 2.2b?\n";
    }
    if (model_num<0) {
        // model not found
	name = "";
        full_name = "";
	model_num = -1;
	symmetry = -1;
	def_freq = FREQ_UNKNOWN;
	return;
    }

    // name and full_name:
    // Special case for strand symmetric model.
    if (model_num == STR_SYM_INDEX) {
      name = "StrSym"; // Can't use MODEL_NAMES[STR_SYM_INDEX] as this is all lowercase, as it must be for parseModelName to work.
      full_name = "Strand Symmetric model (alias WS6.6) (non reversible)";
    } else {
      name = SYMMETRY[symmetry]+MODEL_NAMES[model_num];
      full_name = "Lie Markov model "+SYMMETRY[symmetry]+MODEL_NAMES[model_num]
	+ (TIME_REVERSIBLE[model_num] ? "" : " (non reversible)");
    }

    // def_freq
    int bdf = BDF[model_num];
    if (bdf==0) {
      def_freq=FREQ_EQUAL;
    } else if (bdf==1) {
      switch(symmetry) {
      case 0:
	def_freq=FREQ_DNA_1212;
	break;
      case 1:
	def_freq=FREQ_DNA_1221;
	break;
      case 2:
	def_freq=FREQ_DNA_1122;
	break;
      case 3:
      default:
	cerr << "Can't happen" << endl;
        abort();
      }
    } else if (bdf==2) {
      switch(symmetry) {
      case 0:
	def_freq=FREQ_DNA_RY;
	break;
      case 1:
	def_freq=FREQ_DNA_WS;
	break;
      case 2:
	def_freq=FREQ_DNA_MK;
	break;
      case 3:
      default:
	cerr << "Can't happen" << endl;
        abort();
      }
    } else if (bdf==3) {
      def_freq=FREQ_ESTIMATE;
    }

    return;
}


ModelLieMarkov::~ModelLieMarkov() {
  // Do nothing, for now. model_parameters is reclaimed in ~ModelMarkov
  
  // BQM: Do something now
    if (model_parameters)
        delete [] model_parameters;
}

/*
 * Return 'true' if freq type is compatible with this Lie-Markov model.
 * NOTE: Any freq_type exept FREQ_USER_DEFINED, FREQ_EMPIRICAL or
 * FREQ_ESTIMATE is at best redundant, and worst incompatible.
 * The above three are really the only freq types which should be used
 * with an LM model.
 * Actually, the +F1123s are valid with a BDF=3 LM model, but
 * I haven't coded for this possibility so reject it. Could be fixed.
 * 
 * Also for FREQ_USER_DEFINED and FREQ_EMPIRICAL, for LM models with 
 * BDF<3, compatibility depends on the given base freqs. There
 * is code elsewhere which prints a warning if incompatible base freqs
 * (and actual model base freqs will be 'close to' the requested freqs.)
 *
 * (update - this code is now unused, but left in for possible future use.)
 */
bool  ModelLieMarkov::validFreqType() {
  int bdf=BDF[model_num];
  switch(getFreqType()) {
    case FREQ_USER_DEFINED:
    case FREQ_EMPIRICAL:
    case FREQ_ESTIMATE:
        return true;
    case FREQ_UNKNOWN:
    case FREQ_CODON_1x4:
    case FREQ_CODON_3x4:
    case FREQ_CODON_3x4C:
    case FREQ_MIXTURE:
    case FREQ_DNA_1112:
    case FREQ_DNA_1121:
    case FREQ_DNA_1211:
    case FREQ_DNA_2111:
    case FREQ_DNA_1123:
    case FREQ_DNA_1213:
    case FREQ_DNA_1231:
    case FREQ_DNA_2113:
    case FREQ_DNA_2131:
    case FREQ_DNA_2311:
        return false;
    case FREQ_EQUAL:
        return (bdf==0);
    case FREQ_DNA_RY:   return("+FRY");
        return (bdf==2 && symmetry==0);
    case FREQ_DNA_WS:   return("+FWS");
        return (bdf==2 && symmetry==1);
    case FREQ_DNA_MK:   return("+FMK");
        return (bdf==2 && symmetry==2);
    case FREQ_DNA_1122: return("+F1122");
        return (bdf==1 && symmetry==2);
    case FREQ_DNA_1212: return("+F1212");
        return (bdf==1 && symmetry==0);
    case FREQ_DNA_1221: return("+F1221");
        return (bdf==1 && symmetry==1);
    default: throw("Unrecoginzed freq_type in validFreqType - can't happen");
    }
}

/*
 * Overrides ModelMarkov::getNDimFreq().
 * The degrees of freedom in base frequencies are already accounted
 * for in num_param, so no more should be added.
 */

//int ModelLieMarkov::getNDimFreq() { 
//	return 0;
//}

/*
 * Some LM models are time reversible. Currently this is used in 
 * ModelFactory::getNParameters() to adjust the degrees of freedom
 * by one. Should the code ever be changed such that TR LM models
 * are given an unrooted tree and optimized by TR methods,
 * ModelFactory::getNParameters() may need changing.
 */
bool ModelLieMarkov::isReversible() {
    // TODO: crash when setting reversible to true
//    ASSERT(is_reversible == TIME_REVERSIBLE[model_num]);
//    return(TIME_REVERSIBLE[model_num]);
    return false;
}

/* static */ bool ModelLieMarkov::validModelName(string model_name) {
    int model_num, symmetry;
    parseModelName(model_name,&model_num,&symmetry);
    return (model_num!=-1);
}

/*
 * Model names are like 3.3a or WS6.6.
 * The model name is something in the list MODEL_NAMES, optionally
 * prefixed by "RY", "WS" or "MK" to set the distinguished pair. 
 * If the model has full symmetry, prefix is irrelevant and is ignored.
 * If the model does not have full symmetry and has no prefix, "RY"
 * pair is assumed.
 *
 * Returns number of entry on MODEL_NAMES in model_num (-1 if not found),
 * and symmetry is 0 for RY, 1 for WS, 2 for MK, 3 for full symmetry.
 *
 * SPECIAL CASE: "StrSym" (case insensitive) is a synonym for WS6.6 
 * (strand symmetric). A minor misfeature is that RY, WS and MK will
 * be accepted as prefixes to StrSym (e.g. "ryStrsym" is an alias for StrSym) 
 */


/* static */ void ModelLieMarkov::parseModelName(string model_name, int* model_num, int* symmetry) {
    int len = model_name.length();
    string base_name;
    string name_lower = model_name;
    for (string::iterator it = name_lower.begin(); it != name_lower.end(); it++)
	(*it) = tolower(*it);
    if (name_lower.find("ry")==0) {
      // found "RY" at start of model name
      *symmetry = 0;
      base_name = name_lower.substr(2,len-2);
    } else if (name_lower.find("ws")==0) {
      // found "WS" at start of model name
      *symmetry = 1;
      base_name = name_lower.substr(2,len-2);
    } else if (name_lower.find("mk")==0) {
      // found "MK" at start of model name
      *symmetry = 2;
      base_name = name_lower.substr(2,len-2);
    } else { 
      // Found no prefix
      *symmetry = 0;
      base_name = name_lower;
    } 
    // search for basename in MODEL_NAMES
    *model_num = -1; // not found yet
    for (int i=0; i<NUM_LM_MODELS; i++) {
      if (MODEL_NAMES[i].compare(base_name)==0) {
	*model_num = i;
	break;
      }
    }
    // Special case: strand symmetric model has WS symmetry
    if (*model_num == STR_SYM_INDEX) *symmetry = 1;
    // set full symmetry if have a fully symmetric model
    if (*model_num>=0 && FULL_SYMMETRY[*model_num]) *symmetry = 3;
    return;
}

/*
 * Overrides ModelMarkov::getName().
 * Avoids appending +FO to name, as this is implied by how LM models 
 * work.
 * Minh: you might chose to remove this override, if you like "+FO"
 * to be on LM model names.
 */

string ModelLieMarkov::getName() {
    switch(getFreqType()) {
    case FREQ_ESTIMATE:
        return name;
    case FREQ_EMPIRICAL:
        return name+"+F";
    case FREQ_USER_DEFINED:
        return name+"+FU";
    case FREQ_EQUAL:
      return name;
    default:
       	cerr << "Bad freq_type for a Lie-Markov model. Can't happen" << endl;
        abort();
    }
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
		// If we end up with optimum on boundary, try restarting
		// optimization from different start point, because LM models
		// have a tendancy to find local maxima on boundary.
		bound_check[i] = true; 
	}
}

void ModelLieMarkov::setVariables(double *variables) {
	int nrate = getNDim();

    // non-reversible case
    if (!is_reversible) {
        if (nrate > 0)
            memcpy(variables+1, model_parameters, nrate*sizeof(double));
        return;
    }

	if (freq_type == FREQ_ESTIMATE) nrate -= (num_states-1);
	if (nrate > 0)
		memcpy(variables+1, rates, nrate*sizeof(double));
	if (freq_type == FREQ_ESTIMATE) {
        // 2015-09-07: relax the sum of state_freq to be 1, this will be done at the end of optimization
		int ndim = getNDim();
		memcpy(variables+(ndim-num_states+2), state_freq, (num_states-1)*sizeof(double));
    }
}

bool ModelLieMarkov::getVariables(double *variables) {
	int nrate = getNDim();
	int i;
	bool changed = false;

    // non-reversible case
    if (!is_reversible) {
        for (i = 0; i < nrate && !changed; i++)
            changed = (model_parameters[i] != variables[i+1]);
        if (changed) {
            memcpy(model_parameters, variables+1, nrate * sizeof(double));
            setRates();
        }
        return changed;
    }

	if (freq_type == FREQ_ESTIMATE) nrate -= (num_states-1);
	if (nrate > 0) {
		for (i = 0; i < nrate; i++)
			changed |= (rates[i] != variables[i+1]);
		memcpy(rates, variables+1, nrate * sizeof(double));
	}

	if (freq_type == FREQ_ESTIMATE) {
        // 2015-09-07: relax the sum of state_freq to be 1, this will be done at the end of optimization
        // 2015-09-07: relax the sum of state_freq to be 1, this will be done at the end of optimization
		int ndim = getNDim();
		for (i = 0; i < num_states-1; i++)
			changed |= (state_freq[i] != variables[i+ndim-num_states+2]);
		memcpy(state_freq, variables+(ndim-num_states+2), (num_states-1)*sizeof(double));


	}
	return changed;
}

/*
 * Lie Markov model parameter restart strategy:
 * If no parameters (in 'guess') are on the boundary of parameter space,
 * no restart is needed. 
 * If restart is needed, the first attempt is to take parameters which are on
 * boundary, halve them and change sign. This means our restart should be
 * well away from the local optimum on the boundary, so if the local optimum
 * is not the global optimum, this restart will hopefully go somewhere else.
 * 
 * On subsequent restarts (iterations 2 to 5), every parameter starts halfway 
 * between 0 and the boundary. The parameters are split into two sections 
 * of (nearly) equal size, and all restart parameters within a section have 
 * the same sign. The pattern of signs is:
 * iteration group1 group2
 *     2       -      -
 *     3       +      +
 *     4       -      +
 *     5       +      -
 * (this ordering gives maximal difference between restarts 2 and 3,
 * which we hope will be more likely to find the way to a different
 * local optimum.)
 */
const int MAX_ITER = 5;
bool ModelLieMarkov::restartParameters(double guess[], int ndim, double lower[], double upper[], bool bound_check[], int iteration) {
    int i;
    bool restart = false;
    if (iteration <= MAX_ITER) {
        for (i = 1; i <= ndim; i++) {
            if (fabs(guess[i]-lower[i]) < 1e-4 || fabs(guess[i]-upper[i]) < 1e-4) {
	        restart = true; break;
	    } // if (fabs...
        } // for
    } // if iteration <= MAX_ITER
    if (restart) {
        if (iteration == 1) {
            unsigned int signbits = 0;
            for (i = ndim; i>0; i--) {
	         if (fabs(guess[i]-lower[i]) < 1e-4 || fabs(guess[i]-upper[i]) < 1e-4) {
                     guess[i] *= -0.5;
	         } // if (fabs...
		 signbits = (signbits << 1) + ((guess[i]>0) ? 1 : 0);
	    } // for
        } else {
            int halfN = ndim/2;
            double sign1 = (iteration==2 || iteration==4) ? -1 : 1;
            double sign2 = (iteration==2 || iteration==5) ? -1 : 1;
            for (i=1; i<=halfN; i++) {
                guess[i] = sign1 * upper[i]/2;
            }
            for (i=halfN+1; i<=ndim; i++) {
                guess[i] = sign2 * upper[i]/2;
            }
	}
	if (verbose_mode >= VB_MED) {
            cout << "Lie Markov Restart estimation at the boundary, iteration " << iteration;
            if (verbose_mode >= VB_MAX) {
                cout << ", new start point:" << std::endl << guess[1] ;
                for (i = 2; i <= ndim; i++) cout << "," << guess[i]; 
            }
            cout << std::endl;
	}
    } else {
        if (iteration > 1 && verbose_mode >= VB_MAX)
	  cout << "Lie Markov restarts ended at iteration " << iteration-1 << std::endl;
    
    } // if restart else
    return (restart);
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

  // BQM 2017-05-02: set reversibility
  // TODO: crash when setting reversible to true
  //setReversible(TIME_REVERSIBLE[model_num]);
  setReversible(false);

  // if not otherwise specified, use FREQ_ESTIMATE.
  if (getFreqType() == FREQ_UNKNOWN) freq_type = FREQ_ESTIMATE;

  /* 
   * Note I've chosen to be picky here, and reject almost all <model>+F
   * frequency constraints. With some effort, I could be less picky:
   * validFreqType() can detect when the +F... is redundant rather than
   * contradictory. In some cases, a submodel could be used, e.g.
   * RY5.6b+FQ is RY2.2b.  
   */

  if (getFreqType() != FREQ_EMPIRICAL && 
      getFreqType() != FREQ_USER_DEFINED && 
      getFreqType() != FREQ_ESTIMATE) {
      // Note to Minh: this is formatted horribly - one hugely long line - if you know how to tidily output 
      // multiline throw, please fix.
      throw("Lie-Markov models can only have base frequencies specified as\nempirical (-f c, <model>+FC or default), user defined (<model>+F{<freqs>})\nor estimated/optimized (-f o, <model>+FO).\nEach Lie-Markov model has its own base frequency constraints (corresponding\nto one of +FQ, +F1122,+F1212, +F1221, +FRY, +FWS, +FMK or unconstrained).\nImposing extra constraints is either redundant, makes the model no longer\nLie-Markov, or makes it a lower dimensioned Lie-Markov model.\n");
      //throw("Invalid base frequency constraints for a Lie-Markov model");
  }

  if (getFreqType() == FREQ_EMPIRICAL || 
      getFreqType() == FREQ_USER_DEFINED) {
    int bdf = BDF[model_num];
    // There are no free parameters for base frequencies:
    num_params = MODEL_PARAMS[model_num]-bdf;
    // This populates field state_freq. (TODO: this call might be redundant - check)

    init_state_freq(getFreqType());
    // state_freq is in order {pi_A, pi_C, pi_G, pi_T}
    double tau[3];
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
      double eqbm[4];
      tauToPi(tau,eqbm,symmetry);
      char buffer[200];
      snprintf(buffer,200,"Model %s cannot achieve requested equilibrium base frequencies\n(%5.3f,%5.3f,%5.3f,%5.3f).\nInstead it will use equilibrium base frequencies (%5.3f,%5.3f,%5.3f,%5.3f).\n",
	       name.c_str(),state_freq[0],state_freq[1],state_freq[2],state_freq[3],eqbm[0],eqbm[1],eqbm[2],eqbm[3]);
      outWarning(buffer);
    }

    basis = new double*[num_params+1];
    for (int i=0;i<=num_params;i++) {
      int basisIndex = BASES[model_num][i];
      double unpermuted_rates[NUM_RATES];
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
  } else {
      ASSERT(getFreqType() == FREQ_ESTIMATE); // only other legal possibility
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
    return ModelMarkov::decomposeRateMatrix();
    /*
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
     */
}

void ModelLieMarkov::decomposeRateMatrixEigen3lib() {
  nondiagonalizable = false; // until proven otherwise
    Matrix4d mat(rate_matrix);
    mat.transpose();
    EigenSolver<Matrix4d> eigensolver(mat);
    ASSERT (eigensolver.info() == Eigen::Success);
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
            ASSERT(abs(check(i,j)) < 1e-4);
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
    double a = 1., a2 = 0, b = 0, c = 0, d = 0, d1 = 0, e1 = 0, e2 = 0, f1 = 0, f2 = 0, g1 = 0, g2 = 0;

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
    
    // following code is from Cassius
    
    
    if (name.find("1.1") != string::npos) {
    	a = 1./3.;
    	//Eigenvalues = {0, -4*a, -4*a, -4*a}
		ceval[0] = 0.0; ceval[1] = ceval[2] = ceval[3] = -4.0*a;
		//v0 = {1, 1, 1, 1}
		cevec[0] = cevec[1] = cevec[2] = cevec[3] = 1.0;
		//v1 = {-1, 0, 0, 1}
		cevec[4] = -1.0;
		cevec[5] = 0.0;
		cevec[6] = 0.0;
		cevec[7]= 1.0;
		//v2 = {-1, 0, 1, 0}
		cevec[8] =  -1.0;
		cevec[9] =  0.0;
		cevec[10] =  1.0;
		cevec[11] = 0.0;
		//v3 = {-1, 1, 0, 0}
		cevec[12] = -1.0;
		cevec[13] =  1.0;
		cevec[14] = 0.0;
		cevec[15] =  0.0;

		/*Inverses*/
		cinv_evec[1] = cinv_evec[2] = cinv_evec[3] = cinv_evec[5] = cinv_evec[6] = cinv_evec[9] = cinv_evec[11] = cinv_evec[14] = cinv_evec[15] = -0.25;
		cinv_evec[0] = cinv_evec[4] = cinv_evec[8] = cinv_evec[12] = 0.25;
		cinv_evec[7] = cinv_evec[10] = cinv_evec[13] = 0.75;

    } else if (name.find("2.2b") != string::npos) {
        /******** eigenvalues *********/
        //Eigenvalues = {0, -4 (a - a2), -2 (2 a + a2), -2 (2 a + a2)}
    	a = 1./3.;
        a2 = -rate_matrix[1] + a;
        ceval[0] = 0.0; ceval[1] = -4.0*(a - a2); ceval[2] = ceval[3] = -2.0*(2.0*a + a2);
        /*ceval[0] = 0.0; ceval[1] = 0.0; ceval[2] = ceval[3] = 0.0;*/

        /******** right eigenvectors *********/
		// {{1, 1, 1, 1}, {-1, 1, -1, 1}, {0, -1, 0, 1}, {-1, 0, 1, 0}}

		//v0 = {1, 1, 1, 1}
		cevec[0] = cevec[1] = cevec[2] = cevec[3] = 1.0;
		//v1 = {-1, 1, -1, 1}
		cevec[4] = -1.0;
		cevec[5] =  1.0;
		cevec[6] =  -1.0;
		cevec[7] =  1.0;
		//v2 ={0, -1, 0, 1}
		cevec[8] =  0.0;
		cevec[9] =-1.0;
		cevec[10] = 0.0;
		cevec[11] = 1.0;
		//v3 = {-1, 0, 1, 0}
		cevec[12] = -1.0;
		cevec[13]  = 0.0;
		cevec[14] = 1.0;
		cevec[15] = 0.0;

		/*INVERSES*/

		cinv_evec[0] = cinv_evec[4] = cinv_evec[5] = cinv_evec[8] = cinv_evec[12] = cinv_evec[13] = 0.25;
		cinv_evec[1] = cinv_evec[9] = -0.25;
		cinv_evec[2] = cinv_evec[10] = cinv_evec[15] = 0.;
		cinv_evec[3] = cinv_evec[6] = -0.5;
		cinv_evec[14] = cinv_evec[11] = 0.5;

    }  else if (name.find("3.3a") != string::npos) {
        a = -rate_matrix[0]/3.;
        a2 = (rate_matrix[2] - a)/2.;
        b = rate_matrix[1] + a2 - a;
        /******** eigenvalues *********/
		//Eigenvalues = {0, -4 (a - a2), -2 (2 a + a2 - b), -2 (2 a + a2 + b)}
		ceval[0] = 0.0; ceval[1] = -4.0*(a - a2); ceval[2] = -2.0*(2.0*a + a2 - b); ceval[3] = ceval[2] -4.0*b;

		/******** right eigenvectors  *********/
		// {{1, 1, 1, 1}, {-1, 1, -1, 1}, {-1, -1, 1, 1}, {1, -1, -1, 1}}

		//v0 = {1, 1, 1, 1}
		cevec[0] = cevec[1] = cevec[2] = cevec[3] = 1.0;
		//v1 = {-1, 1, -1, 1}
		cevec[4] =  -1.0;
		cevec[5] = 1.0;
		cevec[6] = -1.0;
		cevec[7] = 1.0;
		//v2 ={-1, -1, 1, 1}
		cevec[8] = -1.0;
		cevec[9]  = -1.0;
		cevec[10] = 1.0;
		cevec[11] = 1.0;
		//v3 = {1, -1, -1, 1}
		cevec[12] = 1.0;
		cevec[13] = -1.0;
		cevec[14] = -1.0;
		cevec[15] = 1.0;

		/******** INVERSE *********/
		cinv_evec[0] =cinv_evec[3] =cinv_evec[4] =cinv_evec[5] =cinv_evec[8] =cinv_evec[10] =cinv_evec[12] =cinv_evec[13] =cinv_evec[14] =cinv_evec[15] = 0.25;
		cinv_evec[1] =cinv_evec[2] =cinv_evec[6] =cinv_evec[7] =cinv_evec[9] =cinv_evec[11] = -0.25;

     }else if (name.find("3.3b") != string::npos) {

    	a = -rate_matrix[0]/3.;
		a2 = (rate_matrix[2] - a)/2.;
		c = -rate_matrix[1] - a2 + a;

		/* cout <<"Los parametros son a = " << a << " a2 =  " << a2 << " c = " << c << endl;*/

		/******** eigenvalues *********/
		//{0, -4 (a - a2), -2 (2 a + a2 - I c), -2 (2 a + a2 + I c)} std::sqrt()
		ceval[0] = 0.0; ceval[1] = -4.0*(a - a2);
		ceval[2] = complex<double> (-2.0*(2.0*a + a2), -2.0* c);
		ceval[3] = complex<double> (-2.0*(2.0*a + a2), 2.0* c);
		/******** right eigenvectors *********/
		// {{1, 1, 1, 1}, {-1, 1, -1, 1}, {I, -1, -I, 1}, {-I, -1, I, 1}}

		//v0 = {1, 1, 1, 1}
		cevec[0] = cevec[1] = cevec[2] = cevec[3] = 1.0;
		//v1 = {-1, 1, -1, 1}
		cevec[4] = -1.0;
		cevec[5] = 1.0;
		cevec[6] = -1.0;
		cevec[7] = 1.0;
		//v2 ={I, -1, -I, 1}
		cevec[8] = complex<double> (0., 1.);
		cevec[9] = -1.0;
		cevec[10] = complex<double> (0., -1.);
		cevec[11] = 1.0;
		//v3 = {-I, -1, I, 1}
		cevec[12] = complex<double> (0., -1.);
		cevec[13] = -1.0;
		cevec[14] = complex<double> (0., 1.);
		cevec[15] = 1.0;
		/******** INVERSE *********/
		cinv_evec[0] =cinv_evec[4] =cinv_evec[5] =cinv_evec[8] =cinv_evec[12] =cinv_evec[13] =cinv_evec[14] =cinv_evec[15] = 0.25;
		cinv_evec[1] =cinv_evec[6] =cinv_evec[7] =cinv_evec[9] = -0.25;
		cinv_evec[2] =cinv_evec[11] =complex<double> (0., -0.25);
		cinv_evec[3] =cinv_evec[10] = complex<double> (0., 0.25);

     } else if (name.find("3.3c") != string::npos) {
		a = -(rate_matrix[0] + rate_matrix[5])/6. ;
		d1 = -rate_matrix[0] - 3.*a;
		a2 = -rate_matrix[1] + a;
		/******** eigenvalues *********/
		//Eigenvalues = {0, -4 (a - a2), -2 (2 a + a2), -2 (2 a + a2)}
		ceval[0] = 0.0; ceval[1] = -4.0*(a - a2); ceval[2] = -2.0*(2.0*a + a2 - d1); ceval[3] = -2.0*(2.0*a + a2 + d1);

		/******** right eigenvectors *********/
		// {{1, 1, 1, 1}, {-((a - a2 - d)/(a - a2 + d)), 1, -((a - a2 - d)/(a - a2 + d)), 1}, {0, -1, 0, 1}, {-1, 0, 1, 0}}

		//v0 = {1, 1, 1, 1}
		cevec[0] = cevec[1] = cevec[2] = cevec[3] = 1.0;
		//v1 = {-((a - a2 - d)/(a - a2 + d)), 1, -((a - a2 - d)/(a - a2 + d)), 1}
		cevec[4]  = -1.;
		cevec[5]  = 1.;
		cevec[6]  = -1.;
		cevec[7]  = 1.;
		//v2 ={0, -1, 0, 1}
		cevec[8] =0.0;
		cevec[9] =  -1.0;
		cevec[10] = 0.0;
		cevec[11] =  1.0;
		//v3 = {-1, 0, 1, 0}}
		cevec[12] =  -1.0;
		cevec[13] =  0.0;
		cevec[14] = 1.0;
		cevec[15] =  0.0;

		 /*INVERSE*/

		cinv_evec[0] =cinv_evec[4] =cinv_evec[5] =cinv_evec[8] =cinv_evec[12] =cinv_evec[13] = 0.25;
		cinv_evec[1] =cinv_evec[9] = -0.25;
		cinv_evec[2] =cinv_evec[7] =cinv_evec[10] =cinv_evec[15] = 0.;
		cinv_evec[3] =cinv_evec[6] = -0.5;
		cinv_evec[11] =cinv_evec[14] = 0.5;

     } else if (name.find("3.4") != string::npos) {
		a = -(rate_matrix[0] + rate_matrix[5])/6. ;
		d = rate_matrix[0] + 3.*a;
		a2 = -rate_matrix[1] + a - d;
		double deno1 = a-a2 + d;
		double deno2 = a-a2;
		if (abs(deno1) < 1.0e-7 || abs(deno2) < 1.0e-7) {
		    // call numerical method if denominator == 0
		    nondiagonalizable = true;
		} else {
			nondiagonalizable = false;
			/******** eigenvalues *********/
			//Eigenvalues = {0, -4 (a - a2), -2 (2 a + a2), -2 (2 a + a2)}
			ceval[0] = 0.0; ceval[1] = -4.0*(a - a2); ceval[2] = ceval[3] = -2.0*(2.0*a + a2);

			/******** right eigenvectors *********/
			// {{1, 1, 1, 1}, {-((a - a2 - d)/(a - a2 + d)), 1, -((a - a2 - d)/(a - a2 + d)), 1}, {0, -1, 0, 1}, {-1, 0, 1, 0}}

			//v0 = {1, 1, 1, 1}
			cevec[0] = cevec[1] = cevec[2] = cevec[3] = 1.0;
			//v1
			cevec[4] = -(a - a2 - d)/deno1;
			cevec[5] = 1.0;
			cevec[6] = cevec[4];
			cevec[7] = 1.0;
			//v2
			cevec[8] = 0.0;
			cevec[9] = -1.0;
			cevec[10] = 0.0;
			cevec[11] = 1.0;
			//v3
			cevec[12] = -1.0;
			cevec[13] = 0.0;
			cevec[14] = 1.0;
			cevec[15] = 0.0;
			/*INVERSE*/
			cinv_evec[0] =cinv_evec[5] =cinv_evec[8] =cinv_evec[13] = 0.25 + 0.25*d/deno2;
			cinv_evec[1] =cinv_evec[9] = -cinv_evec[0];
			cinv_evec[4] =cinv_evec[12] = -cinv_evec[0] + 0.5;
			cinv_evec[14] =cinv_evec[11] = 0.5;
			cinv_evec[3] =cinv_evec[6] = -0.5;
			cinv_evec[2] =cinv_evec[7] =cinv_evec[10] =cinv_evec[15] = 0.;
		}

    } else if (name.find("4.4a") != string::npos) {
		e1 = (rate_matrix[0] - rate_matrix[10])*0.5;
		e2 = (rate_matrix[5] - rate_matrix[15])*0.5;
		a = (-rate_matrix[0] + rate_matrix[4])*0.25;
		d = rate_matrix[4] - a - e1;
		double deno = a+d+e1;
		double deno2 = 4*a;

		if (abs(deno) < 1.0e-7|| abs(deno2) < 1.0e-7) {
			// call numerical method if denominator == 0
			nondiagonalizable = true;
		} else {
			nondiagonalizable = false;
			deno = 1/deno;
			deno2 = 1/deno2;
			/******** eigenvalues *********/
			//Eigenvalues = {0, -4 a, -4 a, -4 a}
			ceval[0] = 0.0; ceval[1] = ceval[2] = ceval[3] = -4.0*a;

			/******** right eigenvectors *********/
			// {{1, 1, 1, 1}, {(-a + d + e2)/(a + d + e1), 0, 0, 1}, {-((a + d - e1)/(a + d + e1)), 0, 1, 0},
			// {-((a - d + e2)/(a + d + e1)), 1, 0, 0}}

			//v0
			cevec[0] = cevec[1] = cevec[2] = cevec[3] = 1.0;
			//v1
			cevec[4] = (-a + d + e2)*deno;
			cevec[5] = 0.0;
			cevec[6] = 0.0;
			cevec[7] = 1.0;
			//v2
			cevec[8] = -(a + d - e1)*deno;
			cevec[9] = 0.0;
			cevec[10] = 1.0;
			cevec[11] = 0.0;
			//v3 =
			cevec[12] = -(a - d + e2)*deno;
			cevec[13] = 1.0;
			cevec[14] = 0.0;
			cevec[15] = 0.0;
			/******** INVERSE *********/

			double auxd = d*deno2; double auxe1 = e1*deno2; double auxe2 = e2*deno2;

			cinv_evec[0] = 0.25 + auxd + auxe1;
			cinv_evec[1] = cinv_evec[2] = cinv_evec[3] = -cinv_evec[0];
			cinv_evec[4] =  0.25 - auxd +auxe2;
			cinv_evec[5] = cinv_evec[6] = -cinv_evec[4];
			cinv_evec[7] = cinv_evec[6] + 1.;
			cinv_evec[8] = 0.25 + auxd - auxe1;
			cinv_evec[9] = cinv_evec[11] = -cinv_evec[8] ;
			cinv_evec[10] = cinv_evec[9] + 1.;
			cinv_evec[12] = 0.25 - auxd - auxe2;
			cinv_evec[14] = cinv_evec[15] = -cinv_evec[12];
			cinv_evec[13] = cinv_evec[14] + 1.;
		}
    } else if (name.find("4.4b") != string::npos) {
		d = (rate_matrix[4] - rate_matrix[3])/2.;
		d1 = (-rate_matrix[0] + rate_matrix[5] + 2.*d)/2.;
		a = (-rate_matrix[0] + d - d1)/3.;
		a2 = -rate_matrix[1] + a - d;

		double deno = a - a2 + d;
		double deno2 = a-a2;

		if (abs(deno) < 1.0e-7 || abs(deno2) < 1.0e-7) {
				// call numerical method if denominator == 0
				nondiagonalizable = true;
		} else {
			nondiagonalizable = false;

			/******** eigenvalues *********/
			//Eigenvalues = {0, -4 (a - a2), -2 (2 a + a2 - d1), -2 (2 a + a2 + d1)}
			ceval[0] = 0.0; ceval[1] = -4.0*(a-a2); ceval[2] = -2.0*(2.0*a + a2 -d1); ceval[3] = ceval[2]-4.0*d1;

			/******** right eigenvectors *********/
			// {{1, 1, 1, 1}, {(-a + a2 + d)/(a - a2 + d), 1, (-a + a2 + d)/(a - a2 + d), 1},
			// {0, -1, 0, 1}, {-1, 0, 1, 0}}

			//v0
			cevec[0] = cevec[1] = cevec[2] = cevec[3] = 1.0;
			//v1
			cevec[4] = (-a + a2 + d)/deno;
			cevec[5] = 1.0;
			cevec[6] = cevec[4];
			cevec[7] = 1.0;
			//v2
			cevec[8] = 0.0;
			cevec[9] = -1.0;
			cevec[10] = 0.0;
			cevec[11] = 1.0;
			//v3 =
			cevec[12] = -1.0;
			cevec[13] = 0.0;
			cevec[14] = 1.0;
			cevec[15] = 0.0;
			/*INVERSE*/
			cinv_evec[0] = 0.25 + 0.25*d/deno2;
			cinv_evec[8] = cinv_evec[5] = cinv_evec[13] = cinv_evec[0];
			cinv_evec[1] = cinv_evec[9] = -cinv_evec[0];
			cinv_evec[2] = cinv_evec[10] = cinv_evec[7] = cinv_evec[15] = 0.;
			cinv_evec[3] = cinv_evec[6] = -0.5;
			cinv_evec[11] = cinv_evec[14] = 0.5;
			cinv_evec[4] = -cinv_evec[0] + 0.5;
			cinv_evec[12] = cinv_evec[4];
		}

    } else if (name.find("4.5a") != string::npos) {

		d = (rate_matrix[0] - rate_matrix[5])/2.;
		a = -(rate_matrix[0] - d)/3.;
		a2 = (rate_matrix[2] - a - d)/2.;
		b = rate_matrix[1] - a + a2 + d;
		double deno = a - a2 + d;
        if (abs(deno) < 1.0e-7 || abs(a-a2) < 1.0e-7) {
            // call numerical method if denominator == 0
            nondiagonalizable = true;
        } else {
            nondiagonalizable = false;
            /******** eigenvalues *********/
            //Eigenvalues = {0, -4 (a - a2), -2 (2 a + a2 - b), -2 (2 a + a2 + b)}
            ceval[0] = 0.0; ceval[1] = -4.0*(a-a2); ceval[2] =-2.0*(2.0*a + a2 - b); ceval[3] = ceval[2]- 4.0*b;

            /******** right eigenvectors *********/
            // {{1, 1, 1, 1}, {-((a - a2 - d)/(a - a2 + d)), 1, -((a - a2 - d)/(a - a2 + d)), 1},
            // {-1, -1, 1, 1}, {1, -1, -1, 1}}

            //v0
            cevec[0] = cevec[1] = cevec[2] = cevec[3] = 1.0;
            //v1
            cevec[4] = -(a - a2 - d)/deno;
            cevec[5] = 1.0;
            cevec[6] = cevec[4];
            cevec[7] = 1.0;
            //v2
            cevec[8] = -1.0;
            cevec[9] = -1.0;
            cevec[10] = 1.0;
            cevec[11] = 1.0;
            //v3 =
            cevec[12] = 1.0;
            cevec[13] = -1.0;
            cevec[14] = -1.0;
            cevec[15] = 1.0;
             /*INVERSE */
            double auxd = 0.25*d/(a-a2);
            cinv_evec[0] = cinv_evec[8] = cinv_evec[5] = cinv_evec[13] = 0.25 + auxd;
            cinv_evec[2] = cinv_evec[11] = cinv_evec[6] = cinv_evec[7] = -0.25;
            cinv_evec[3] = cinv_evec[10] = cinv_evec[14] = cinv_evec[15] = 0.25;
            cinv_evec[1] = cinv_evec[9] = -cinv_evec[0];
            cinv_evec[4] = cinv_evec[12] = 0.25 - auxd;
        }
	} else if (name.find("4.5b") != string::npos) {
		d = (rate_matrix[0] - rate_matrix[5])/2.;
		a = -(rate_matrix[0] + rate_matrix[5])/6.;
		c = (rate_matrix[3] - rate_matrix[1])/2.;
		a2 = -rate_matrix[1] + a - d - c;
		double deno = a-a2 + d;
		if (abs(deno) < 1.0e-7 || abs(a-a2) < 1.0e-7) {
			nondiagonalizable = true;
		} else {
			nondiagonalizable = false;
			deno = 1/deno;
			/******** eigenvalues *********/
			//Eigenvalues = {0, -4 (a - a2), -2 (2 a + a2 - I c), -2 (2 a + a2 + I c)}
			ceval[0] = 0.0; ceval[1] = -4.0*(a-a2); ceval[2] = complex<double> (-2. *(2.* a + a2), -2.*c); ceval[3] = complex<double>  (-2.*(2.*a + a2), 2.*c);

			/******** right eigenvectors *********/
			// {{1, 1, 1, 1}, {-((a - a2 - d)/(a - a2 + d)), 1, -((a - a2 - d)/(a - a2 + d)), 1},
			// {I, -1, -I, 1}, {-I, -1, I, 1}}

			//v0
			cevec[0] = cevec[1] = cevec[2] = cevec[3] = 1.0;
			//v1
			cevec[4] = -(a - a2 - d)*deno;
			cevec[5] = 1.0;
			cevec[6] = cevec[4];
			cevec[7] = 1.0;
			//v2
			cevec[8] = complex<double> (0., 1.);
			cevec[9] = -1.0;
			cevec[10] = complex<double> (0., -1.);
			cevec[11] = 1.0;
			//v3 =
			cevec[12] = complex<double> (0., -1.);
			cevec[13] = -1.0;
			cevec[14] = complex<double> (0., 1.);
			cevec[15] = 1.0;
			/*INVERSE*/
			double auxd = 0.25*d/(a-a2);
			cinv_evec[0] =cinv_evec[8] =cinv_evec[5] =cinv_evec[13] = 0.25 + auxd;
			cinv_evec[1] =cinv_evec[9] =-cinv_evec[0];
			cinv_evec[4] =cinv_evec[12] = 0.25 - auxd;
			cinv_evec[10] =cinv_evec[3] = complex<double> (0, 0.25);
			cinv_evec[2] =cinv_evec[11] = complex<double> (0, -0.25);
			cinv_evec[6] =cinv_evec[7] = -0.25;
			cinv_evec[14] =cinv_evec[15] = 0.25;
		}

    } else if (name.find("5.6a") != string::npos) {
    	/*Cassius note: The formulas are complex; I think they are not worth it to be computed*/
    	nondiagonalizable = true;

	} else if (name.find("5.6b") != string::npos) {
        e1 = (rate_matrix[0] - rate_matrix[10])*0.5;
        e2 = (rate_matrix[5] - rate_matrix[15])*0.5;
        a = -(rate_matrix[0] + rate_matrix[5] -e1 - e2)/6.;
        d = rate_matrix[0] + 3.*a - e1;
        a2 = -rate_matrix[1] + a - d + e2;
        double deno56b = (2.0* a + a2 + e1 + e2);
        if (abs(deno56b) < 1.0e-7 || abs(a - a2 + d) < 1.0e-7 || abs(2.*a + a2) < 1.0e-7 || abs(a - a2) < 1.0e-7) {
        	nondiagonalizable = true;
        } else {
			nondiagonalizable = false;
		   /******** eigenvalues *********/
			//Eigenvalues = {0, -4 (a - a2), -2 (2 a + a2), -2 (2 a + a2)}
			ceval[0] = 0.0; ceval[1] = -4.0*(a - a2); ceval[2] = ceval[3] =-2.0*(2.0*a + a2);

			/******** right eigenvectors *********/
			// {{1, 1, 1, 1},
			//{-((a - a2 - d)/(a - a2 + d)), 1, -((a - a2 - d)/(a - a2 + d)), 1},
			// {(2 e2)/(2 a + a2 + e1 + e2),-((2 a + a2 + e1 - e2)/(2 a + a2 + e1 + e2)), 0, 1},
			// {-((2 a + a2 - e1 + e2)/(2 a + a2 + e1 + e2)), (2 e1)/(2 a + a2 + e1 + e2), 1, 0}}

			//v0
			cevec[0] = cevec[1] = cevec[2] = cevec[3] = 1.0;
			//v1
			cevec[4] = -((a - a2 - d)/(a - a2 + d));
			cevec[5] = 1.0;
			cevec[6] = cevec[4];
			cevec[7] = 1.0;
			//v2;   we declare an auxiliar denominator in order to reuse it

			cevec[8] = (2.0* e2)/deno56b;
			cevec[9] = -1.0 + cevec[8];
			cevec[10] = 0.0;
			cevec[11] = 1.0;
			//v3
			cevec[13] = (2.0* e1)/deno56b;
			cevec[12] = -1.0+cevec[13];
			cevec[14] = 1.0;
			cevec[15] = 0.0;
			/*INVERSE*/
			double deno = 0.5/(2.*a + a2); double auxe1 = e1*deno; double auxe2 = e2*deno; double auxd = d*0.25/(a-a2);
			cinv_evec[1] = cinv_evec[9] = -0.25 - auxd;
			cinv_evec[5] = cinv_evec[13] = - cinv_evec[1];
			cinv_evec[2] =  -auxe1;
			cinv_evec[10] = - cinv_evec[2];
			cinv_evec[6] = -0.5 - auxe2;
			cinv_evec[14] = -cinv_evec[6];
			cinv_evec[3] =  -0.5 -auxe1;
			cinv_evec[11] = - cinv_evec[3];
			cinv_evec[7] = -auxe2;
			cinv_evec[15] = -cinv_evec[7];
			cinv_evec[0] = 0.25 + auxe1 + auxd;
			cinv_evec[4] = 0.25 + auxe2 -auxd;
			cinv_evec[8] = 0.25 - auxe1 + auxd;
			cinv_evec[12] = 0.25 - auxe2 - auxd;
        }
    }

     else if (name.find("5.7a") != string::npos) {
		a = -(rate_matrix[0] + rate_matrix[10])/6.;
		e1 = rate_matrix[0] + 3.*a;
		e2 = rate_matrix[5] + 3.*a;
		a2 = (rate_matrix[2] - a + e1)/2.;
		b = rate_matrix[1] - a + a2 - e2;
		double deno1 = 2.0*a + a2 - b + e1 + e2;
		double deno2 = 2.0* a + a2 + b - e1 + e2;
		double deno3 = 2.*a + a2 - b;
		double deno4 = 2.*a + a2 + b;
		if (abs(deno1) < 1.0e-7 || abs(deno2) < 1.0e-7|| abs(deno3) < 1.0e-7 || abs(deno4) < 1.0e-7) {
			nondiagonalizable = true;
		} else {
			nondiagonalizable = false;
			/*cout <<"Los parametros son a = " << a << " e1 =  " << e1 << "e2 = " << e2 << " a2 =  " << a2 << " b = " << b << endl;*/
			/******** eigenvalues *********/
			//Eigenvalues = {0, -4 (a - a2), -2 (2 a + a2 - b), -2 (2 a + a2 + b)}
			ceval[0] = 0.0; ceval[1] = -4.0*(a - a2); ceval[2] = -2.0*(2.0*a + a2 - b); ceval[3] = ceval[2]-4.0*b;

			/******** right eigenvectors *********/
			// {{1, 1, 1, 1},
			// {-1, 1, -1, 1},
			// {-((2 a + a2 - b - e1 - e2)/(2 a + a2 - b + e1 + e2)), -((2 a + a2 - b - e1 - e2)/(2 a + a2 - b + e1 + e2)), 1, 1},
			//{1, -((2 a + a2 + b + e1 - e2)/(2 a + a2 + b - e1 + e2)), -((2 a + a2 + b + e1 - e2)/(2 a + a2 + b - e1 + e2)), 1}}

			//v0
			cevec[0] = cevec[1] = cevec[2] = cevec[3] = 1.0;
			//v1
			cevec[4] = -1.0;
			cevec[5] = 1.0;
			cevec[6] = -1.0;
			cevec[7] = 1.0;
			//v2
			cevec[8] = -((2.0*a + a2 - b - e1 - e2)/deno1);
			cevec[9] = cevec[8];
			cevec[10] = 1.0;
			cevec[11] = 1.0;
			//v3
			cevec[12] = 1.0;
			cevec[13] = -((2.0*a + a2 + b + e1 - e2)/deno2);
			cevec[14] = cevec[13];
			cevec[15] = 1.0;
			/******** INVERSE ****** There is a sum zero property ***/
			double auxs = 0.25*(e1 + e2)/deno3; double auxr = 0.25*(-e1 + e2)/deno4;
			cinv_evec[2] = cinv_evec[6] =-0.25 - auxs;
			cinv_evec[10] = cinv_evec[14] = - cinv_evec[2];
			cinv_evec[3] = cinv_evec[15] = 0.25 + auxr;
			cinv_evec[7] = cinv_evec[11] =-cinv_evec[3];
			cinv_evec[0] = 0.25 + auxs - auxr;
			cinv_evec[4] = 0.25 + auxs + auxr;
			cinv_evec[8] = 0.5 - cinv_evec[0];
			cinv_evec[12] = 0.5 - cinv_evec[4];
			cinv_evec[1] = cinv_evec[9] = -0.25;
			cinv_evec[5] = cinv_evec[13] = 0.25;
		}
     } else if (name.find("5.7b") != string::npos) {
		f1 = (rate_matrix[0] - rate_matrix[10])*0.5;
		f2 = (rate_matrix[15] - rate_matrix[5])*0.5;
		a = (-rate_matrix[0] + f1)/3.;
		a2 = (rate_matrix[2] - a + f1)*0.5;
		b = rate_matrix[1] - a + a2 -f2;
		double deno57br1 = (3.0*a2 - b - f1 - f2);
		double deno57br2 = (3.0*a2 + b + f1 - f2);
		double deno3 = 3.*a2 - b;
		double deno4 = 3.*a2 + b;
		if (abs(deno57br1) < 1.0e-7 || abs(deno57br2) < 1.0e-7 || abs(deno3) < 1.0e-7 || abs(deno4) < 1.0e-7) {
			nondiagonalizable = true;
		} else {
			nondiagonalizable = false;
			deno57br1 = 1/deno57br1;
			deno57br2 = 1/deno57br2;
			 /*cout << "These are the parameters  a  = " << a << " a2 =  " << a2 << "   b = " << b << "  f1 =  " << f1 << "  f2 = " << f2 << endl;*/
			/******** eigenvalues *********/

			//Eigenvalues = {0, -4 (a - a2), -2 (2 a + a2 - b), -2 (2 a + a2 + b)}
			ceval[0] = 0.0; ceval[1] = -4.0*(a - a2); ceval[2] = -2.0*(2.0*a + a2 - b); ceval[3] = ceval[2] - 4.0*b;

			/******** right eigenvectors *********/
			// {{1, 1, 1, 1},
			//{-1, 1, -1, 1},
			//{-1, -((3 a2 - b + f1 + f2)/(3 a2 - b - f1 - f2)), -((-3 a2 + b - f1 - f2)/(3 a2 - b - f1 - f2)), 1},
			// {-((-3 a2 - b + f1 - f2)/(3 a2 + b + f1 - f2)), -((3 a2 + b - f1 + f2)/(3 a2 + b + f1 - f2)), -1, 1}}

			//v0
			cevec[0] = cevec[1] = cevec[2] = cevec[3] = 1.0;
			//v1
			cevec[4] = -1.0;
			cevec[5] = 1.0;
			cevec[6] = -1.0;
			cevec[7] = 1.0;
			//v4
			cevec[8] = -1.0;
			cevec[9] = -(3.0*a2 - b + f1 + f2)*deno57br1;
			cevec[10] = -cevec[9];
			cevec[11] = 1.0;
			//v3
			cevec[12] = -(-3.0*a2 - b + f1 - f2)*deno57br2;
			cevec[13] = -cevec[12];
			cevec[14] = -1.0;
			cevec[15] = 1.0;
			/*INVERSE : They satisfy a sumzero property*/

			cinv_evec[0] =cinv_evec[4] =cinv_evec[8] =cinv_evec[12] = 0.25;
			cinv_evec[2] =cinv_evec[6] = -0.25 + 0.25*(f1 + f2)/(3.*a2 - b);
			cinv_evec[10] =cinv_evec[14] = - cinv_evec[2];
			cinv_evec[3] =cinv_evec[15] = 0.25 + 0.25*(f1 - f2)/(3.*a2 + b);
			cinv_evec[11] =cinv_evec[7] = -cinv_evec[3];
			cinv_evec[1] = -cinv_evec[0]-cinv_evec[2]-cinv_evec[3];
			cinv_evec[5] = -cinv_evec[4]-cinv_evec[6] - cinv_evec[7];
			cinv_evec[9] = -cinv_evec[1] - 0.5;
			cinv_evec[13] = -cinv_evec[5] +0.5;
		}

	} else if (name.find("5.7c") != string::npos) {
		g1 = (rate_matrix[0] - rate_matrix[10])*0.5;
		g2 = (rate_matrix[15] - rate_matrix[5])*0.5;
		a = (-rate_matrix[0] + g1)/3.;
		a2 = (rate_matrix[2] - a - g1)*0.5;
		b = rate_matrix[1] - a + a2 + g1;
		double aux57b1 = 9.0*a2*a2 - b*b;
		double aux57b2 = 6.0*a2*g1 + 2.0 * b * g2;
		double aux57b3 = 6.0*a2*g2 + 2.0 * b * g1;
		double deno57bl = aux57b1 + aux57b3;
		double deno2 = 3.*a2 - b;
		double deno3 = 3.*a2 + b;

		if (abs(deno57bl) < 1.0e-7 || abs(deno2) < 1.0e-7 || abs(deno3) < 1.0e-7) {
			nondiagonalizable = true;
		} else {
			nondiagonalizable = false;
			 deno57bl = 1./deno57bl;
			//Eigenvalues = {0, -4 (a - a2), -2 (2 a + a2 - b), -2 (2 a + a2 + b)}
			ceval[0] = 0.0; ceval[1] = -4.0*(a - a2); ceval[2] = -2.0*(2.0*a + a2 - b); ceval[3] = ceval[2] - 4.0*b;

			/******** right eigenvectors *********/
			// {{1, 1, 1, 1},
			//{-((9 a2^2 - b^2 + 6 a2 g1 + 2 b g2)/(9 a2^2 - b^2 + 2 b g1 + 6 a2 g2)),
			// -((-9 a2^2 + b^2 + 2 b g1 + 6 a2 g2)/(9 a2^2 - b^2 + 2 b g1 + 6 a2 g2)),
			// -((9 a2^2 - b^2 - 6 a2 g1 - 2 b g2)/(9 a2^2 - b^2 + 2 b g1 + 6 a2 g2)),
			// 1},
			// {-1, -1, 1, 1},
			// {1, -1, -1, 1}}

			//v0
			cevec[0] = cevec[1] = cevec[2] = cevec[3] = 1.0;
			//v1
			cevec[4] = -(aux57b1 + aux57b2 )* deno57bl;
			cevec[5] = -(-aux57b1 + aux57b3)* deno57bl;
			cevec[6] = -(aux57b1 - aux57b2)* deno57bl;
			cevec[7] = 1.0;
			//v2
			cevec[8] = -1.0;
			cevec[9] = -1.0;
			cevec[10] = 1.0;
			cevec[11] = 1.0;
			//v3
			cevec[12] = 1.0;
			cevec[13] = -1.0;
			cevec[14] = -1.0;
			cevec[15] = 1.0;
			/*INVERSE*/
			cinv_evec[0] =cinv_evec[4] =cinv_evec[8] =cinv_evec[12] = 0.25;

			cinv_evec[2] = -0.25 + 0.25*(g1 + g2)/deno2;
			cinv_evec[6] =-0.5 - cinv_evec[2];
			cinv_evec[10] =-cinv_evec[6];
			cinv_evec[14] = -cinv_evec[2];

			cinv_evec[3] = 0.25 - 0.25*(g1 - g2)/deno3;
			cinv_evec[7] =-cinv_evec[3];
			cinv_evec[11] =cinv_evec[3] - 0.5;
			cinv_evec[15] = - cinv_evec[11];

			cinv_evec[1] = cinv_evec[9] = -cinv_evec[0]-cinv_evec[2]-cinv_evec[3];
			cinv_evec[5] = cinv_evec[13] = -cinv_evec[1];
		}

	} else if (name.find("5.11a") != string::npos) {

		e2 = (rate_matrix[1] - rate_matrix[3])*0.5;
		e1 = (rate_matrix[4] - rate_matrix[6])*0.5;
		d1 = (rate_matrix[5] - rate_matrix[0] - e2 + e1)*0.5;
		a = (-rate_matrix[0] -d1 + e1)/3.;
		a2 = -rate_matrix[4] + a + e1;
		double deno511a = 2.* a + a2 - d1 + e2;
		double deno2 = 2.*a + a2 +d1;
		double deno3 = 2.*a + a2 - d1;

		if (abs(deno511a) < 1.0e-7 || abs(deno2) < 1.0e-7 || abs(deno3) < 1.0e-7 || abs(e1) < 1.0e-7) {
			nondiagonalizable = true;
		} else {
			nondiagonalizable = false;
			/******** eigenvalues REPASAR *********/
			//Eigenvalues = {0, -4 (a - a2), -2 (2 a + a2 - d1), -2 (2 a + a2 + d1)}
			ceval[0] = 0.0; ceval[1] = -4.0*(a - a2); ceval[2] = -2.0*(2.0*a + a2 - d1); ceval[3] = ceval[2] - 4.*d1;

			/******** right eigenvectors *********/
			//{{1, 1, 1, 1},
			//{-1, 1, -1, 1},
			//{e2/(2 a + a2 - d1 + e2), -((2 a + a2 - d1 - e2)/(2 a + a2 - d1 + e2)), e2/(2 a + a2 - d1 + e2), 1},
			//{-((2 a + a2 + d1 - e1)/e1), 1, -((-2 a - a2 - d1 - e1)/e1), 1}}

			//v0
			cevec[0] = cevec[1] = cevec[2] = cevec[3] = 1.0;
			//v1
			cevec[4] = -1.0;
			cevec[5] = 1.0;
			cevec[6] = -1.0;
			cevec[7] = 1.0;
			//v2

			cevec[8] = e2/deno511a;
			cevec[9] =-1. + 2.*cevec[8];
			cevec[10] = cevec[8];
			cevec[11] = 1.0;
			//v3
			double aux511a = 2.* a + a2 + d1;
			cevec[12] = 1.- aux511a/e1;
			cevec[13] = 1.0;
			cevec[14] = -cevec[12] +2.;
			cevec[15] = 1.0;
			/*INVERSE*/
			cinv_evec[2] =cinv_evec[10] = 0.;
			cinv_evec[5] =cinv_evec[13] = 0.25;
			cinv_evec[1] =cinv_evec[9] = -0.25;
			cinv_evec[3] =-0.5*e1/(2.*a + a2 +d1);
			cinv_evec[0] = 0.25 -cinv_evec[3];
			cinv_evec[4] = 0.25 + 0.5*e2/(2.*a + a2 - d1);
			cinv_evec[6] = -cinv_evec[4] - 0.25;
			cinv_evec[8] = -cinv_evec[0] + 0.5;
			cinv_evec[12] = -cinv_evec[4] + 0.5;
			cinv_evec[11] = -cinv_evec[3];
			cinv_evec[14] = -cinv_evec[6];
		}
	} else if (name.find("5.11b") != string::npos) {
		f2 = (rate_matrix[1] - rate_matrix[3])*0.5;
		f1 = (rate_matrix[6] - rate_matrix[4])*0.5;
		a = -(rate_matrix[0] + rate_matrix[5] -f1 + f2)/6.;
		a2 = -rate_matrix[1] + f2 + a;
		d1 = -rate_matrix[0] - 3.*a + f1;
		double deno1 = 3.* a2 - d1 - f2;

		if (abs(deno1) < 1.0e-7 || abs(f1) < 1.0e-7|| abs(3.*a2 -d1) < 1.0e-7 || abs(3.*a2 + d1) < 1.0e-7) {
			nondiagonalizable = true;
		} else {
			nondiagonalizable = false;
			/******** eigenvalues *********/
			//Eigenvalues = {0, -4 (a - a2), -2 (2 a + a2 - d1), -2 (2 a + a2 + d1)}
			ceval[0] = 0.0; ceval[1] = -4.0*(a - a2); ceval[2] = -2.0*(2.0*a + a2 - d1); ceval[3] = ceval[2] - 4.0*d1;

			/******** right eigenvectors  *********/
			// {{1, 1, 1, 1},
			// {-1, 1, -1, 1},
			// {f2/(3 a2 - d1 - f2), -((3 a2 - d1 + f2)/(3 a2 - d1 - f2)), f2/(3 a2 - d1 - f2), 1},
			// {-((-3 a2 - d1 + f1)/f1), 1, -((3 a2 + d1 + f1)/f1), 1}}

			//v0
			cevec[0] = cevec[1] = cevec[2] = cevec[3] = 1.0;
			//v1
			cevec[4] = -1.0;
			cevec[5] = 1.0;
			cevec[6] = -1.0;
			cevec[7] = 1.0;
			//v2
			cevec[8] = f2/deno1;
			cevec[9] = -(1. + 2.* f2/deno1);
			cevec[10] = cevec[8];
			cevec[11] = 1.0;
			//v3
			cevec[12] = (3.0*a2 + d1)/f1 -1.;
			cevec[13] = 1.0;
			cevec[14] = -cevec[12]-2.;
			cevec[15] = 1.0;
			/*INVERSE*/
			cinv_evec[0] =cinv_evec[4] =cinv_evec[8] =cinv_evec[12] = 0.25;
			cinv_evec[2] =cinv_evec[7] =cinv_evec[10] =cinv_evec[15] = 0.0;

			cinv_evec[3] = 0.5*f1/(3.*a2 + d1);
			cinv_evec[1] = - cinv_evec[3] - 0.25;
			cinv_evec[11] =-cinv_evec[3];
			cinv_evec[9] =cinv_evec[3] - 0.25;

			double auxf2 = 0.5*f2/(3.*a2 -d1);
			cinv_evec[5] = 0.25 - auxf2;
			cinv_evec[6] = -0.5 + auxf2;
			cinv_evec[13] = 0.25 + auxf2;
			cinv_evec[14] = 0.5 - auxf2;
		}
	} else if (name.find("5.11c") != string::npos) {
		g1 = (rate_matrix[0] - rate_matrix[10])*0.5;
		g2 = (rate_matrix[15] - rate_matrix[5])*0.5;
		d1 =(rate_matrix[5] - rate_matrix[0] + g2 + g1)*0.5;
		a = (-rate_matrix[0] - d1 + g1)/3.;
		a2 = -rate_matrix[1] + a -g1;
		double auxdeno1 = 3.*a2 + d1;
		double auxnum1 = auxdeno1 - 2.*d1;
		double auxdeno2 = auxnum1 + 2.*g2;
		double deno2 = 3.*a2 - d1;
		double deno3 = 3.*a2 + d1;
		if (abs(auxdeno1) < 1.0e-7 || abs(auxdeno2) < 1.0e-7 || abs(deno2) < 1.0e-7 || abs(deno3) < 1.0e-7) {
			nondiagonalizable = true;
		} else {
			nondiagonalizable = false;
			double deno = 1/(auxdeno1 * auxdeno2);
			/*cout <<"The parameters are a = " << a << " a2 =  " << a2 << " d1 = " << d1 << " g1 = " << g1 << "  g2 = " << g2 << endl;*/
			/******** eigenvalues *********/
			//Eigenvalues = {0, -4 (a - a2), -2 (2 a + a2 - d1), -2 (2 a + a2 + d1)}
			ceval[0] = 0.0; ceval[1] = -4.0*(a - a2); ceval[2] = -2.0*(2.0*a + a2 - d1); ceval[3] = ceval[2] - 4.0*d1;

			/******** right eigenvectors *********/
			// {{1, 1, 1, 1},
			// {-(((3 a2 - d1) (3 a2 + d1 + 2 g1))/((3 a2 + d1) (3 a2 - d1 + 2 g2))),
			//-((-3 a2 + d1 + 2 g2)/(3 a2 - d1 + 2 g2)),
			//(-9 a2^2 + d1 (d1 - 2 g1) + 6 a2 g1)/((3 a2 + d1) (3 a2 - d1 + 2 g2)),
			//1},
			 //0, -1, 0, 1},
			 //{-1, 0, 1, 0}}
			//v0
			cevec[0] = cevec[1] = cevec[2] = cevec[3] = 1.0;
			//v1

			cevec[4] = -auxnum1*(auxdeno1 + 2.*g1)*deno;
			cevec[5] = -(-auxnum1 + 2.*g2)/auxdeno2;
			cevec[6] = (-auxdeno1 + 2.*g1)*auxnum1*deno;
			cevec[7] = 1.0;
			//v2
			cevec[8] = 0.0;
			cevec[9] = -1.0;
			cevec[10] = 0.0;
			cevec[11] = 1.0;
			//v3
			cevec[12] = -1.0;
			cevec[13] = 0.0;
			cevec[14] = 1.0;
			cevec[15] = 0.0;
			/*INVERSE*/
			cinv_evec[0] =cinv_evec[4] =cinv_evec[8] =cinv_evec[12] = 0.25;
			cinv_evec[2] =cinv_evec[10] = 0.5*g2 /deno2;
			cinv_evec[7] =cinv_evec[15] = -0.5*g1 /deno3;
			cinv_evec[1] = cinv_evec[9] = -0.25 - cinv_evec[2];
			cinv_evec[3] = -0.5 -cinv_evec[7];
			cinv_evec[11] = 0.5 - cinv_evec[7];
			cinv_evec[5] = cinv_evec[13] = 0.25 + cinv_evec[2];
			cinv_evec[6] = -0.5 -cinv_evec[2];
			cinv_evec[14] = 0.5 - cinv_evec[2];
		}

	} else if (name.find("5.16") != string::npos) {
		g1 = (rate_matrix[0] - rate_matrix[10])*0.5;
		g2 = (rate_matrix[15] - rate_matrix[5])*0.5;
		d = (-rate_matrix[5] + rate_matrix[0] - g2 - g1)*0.5;
		a = (-rate_matrix[0] + d + g1)/3.;
		a2 = -rate_matrix[1] + a - d - g1;
		double aux516r1 = 3.*a2*(a -a2);
		double aux516r2 = 2.*g1*(a -a2);
		double aux516r3 = 2.*g2*(a -a2);
		double aux513r4 = 3.*a2*d;
		double deno516r = aux516r1 + aux513r4 + aux516r3;
		if (abs(deno516r) < 1.0e-7 || abs(a2) < 1.0e-7 || abs(a-a2) < 1.0e-7) {
			nondiagonalizable = true;
		} else {
			/*cout <<"The parameters are a = " << a << " a2 =  " << a2 << " d = " << d << " g1 = " << g1 << "  g2 = " << g2 << endl;*/
			/******** eigenvalues *********/
			//Eigenvalues = {0, -4 (a - a2), -2 (2 a + a2), -2 (2 a + a2)}
			ceval[0] = 0.0; ceval[1] = -4.0*(a - a2); ceval[2] = -2.0*(2.0*a + a2); ceval[3] = ceval[2];

			/******** right eigenvectors *********/
			// {{1, 1, 1, 1},
			// {-((3 a a2 - 3 a2^2 - 3 a2 d + 2 a g1 - 2 a2 g1)/(3 a a2 - 3 a2^2 + 3 a2 d + 2 a g2 - 2 a2 g2)),
			// (3 a a2 - 3 a2^2 + 3 a2 d - 2 a g2 + 2 a2 g2)/(3 a a2 - 3 a2^2 + 3 a2 d + 2 a g2 - 2 a2 g2),
			// -((3 a a2 - 3 a2^2 - 3 a2 d - 2 a g1 + 2 a2 g1)/(3 a a2 - 3 a2^2 + 3 a2 d + 2 a g2 - 2 a2 g2)),
			// 1},
			// {0, -1, 0, 1},
			// {-1, 0, 1, 0}}

			//v0
			cevec[0] = cevec[1] = cevec[2] = cevec[3] = 1.0;
			//v1

			cevec[4] = -(aux516r1 - aux513r4 + aux516r2)/deno516r;
			cevec[5] = (aux516r1 + aux513r4  - aux516r3)/deno516r;
			cevec[6] = -(aux516r1 - aux513r4 - aux516r2)/deno516r;
			cevec[7] = 1.0;
			//v2
			cevec[8] = 0.0;
			cevec[9] = -1.0;
			cevec[10] = 0.0;
			cevec[11] = 1.0;
			//v3
			cevec[12] = -1.0;
			cevec[13] = 0.0;
			cevec[14] = 1.0;
			cevec[15] = 0.0;

			/***INVERSE***/
			double auxg1 = g1/(6.*a2); double auxg2 = g2/(6.*a2);
			cinv_evec[0] =cinv_evec[8] = 0.25 + 0.25*d/(a-a2);
			cinv_evec[4] =cinv_evec[12] = 0.5 - cinv_evec[0];
			cinv_evec[2] =cinv_evec[10] =auxg2;
			cinv_evec[6] = -0.5 - auxg2;
			cinv_evec[14] = 0.5 - auxg2;
			cinv_evec[3] = -0.5 + auxg1;
			cinv_evec[11] = 0.5 + auxg1;
			cinv_evec[7] =cinv_evec[15] = - auxg1;
			cinv_evec[1] =cinv_evec[9] = -cinv_evec[0] -cinv_evec[2];
			cinv_evec[13] =cinv_evec[5] =-cinv_evec[1];
		}

	    }
	else {
    	cout << "This line should not appear." << endl;
    }
   if (nondiagonalizable == false) {
		/* check eigenvalue equation*/
		std::complex<double> zero;
		int error = 0;
		for (j = 0; j < num_states; j++) {
			for (i = 0, zero = 0.0; i < num_states; i++) {
				for (int k = 0; k < num_states; k++)
					zero += rate_matrix[i*num_states+k] * cevec[j*num_states+k];
				zero -= ceval[j] * cevec[j*num_states+i];
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
					zero += cinv_evec[i*num_states+k] * cevec[k*num_states+j];
				double deltaij = 0;
				if (i == j) deltaij = 1;
				zero -= deltaij;
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
			  cout << "Here we start: complex eigenvalues:";
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
		}
    }
}

void ModelLieMarkov::computeTransMatrix(double time, double *trans_matrix, int mixture, int selected_row) {
    return ModelMarkov::computeTransMatrix(time, trans_matrix, mixture, selected_row);
    /*
  MatrixExpTechnique technique = phylo_tree->params->matrix_exp_technique;
  if (technique == MET_SCALING_SQUARING || nondiagonalizable ) {
        Matrix4d A = Map<Matrix4d>(rate_matrix);
        A = (A.transpose() * time).exp();
        Map<Matrix4d> P(trans_matrix);
        P = A.transpose();

        int i, j;
        for (i = 0; i < 4; i++) {
            double sum = 0.0;
            for (j = 0; j < 4; j++)
                sum += (trans_matrix[i*4+j]);
            ASSERT(fabs(sum-1.0) < 1e-4);
        }
    } else if (technique == MET_EIGEN3LIB_DECOMPOSITION || technique == MET_LIE_MARKOV_DECOMPOSITION) {
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
        if (technique == MET_EIGEN3LIB_DECOMPOSITION) {
			for (i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					trans_matrix[i*4+j] = res(j, i).real();
					ASSERT(fabs(res(j,i).imag()) < 1e-6);
					ASSERT(trans_matrix[i*4+j] >= -0.000001);
					ASSERT(trans_matrix[i*4+j] <=  1.000001);
					if (trans_matrix[i*4+j] < 0)
						trans_matrix[i*4+j] = 0.0;
					if (trans_matrix[i*4+j] > 1)
						trans_matrix[i*4+j] = 1.0;
				}

				ASSERT(fabs(trans_matrix[i*4]+trans_matrix[i*4+1]+trans_matrix[i*4+2]+trans_matrix[i*4+3]-1.0) < 1e-4);
			}
        } else {
        	for (i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					trans_matrix[i*4+j] = res(i,j).real();
					ASSERT(fabs(res(j,i).imag()) < 1e-6);
					ASSERT(trans_matrix[i*4+j] >= -0.000001);
					ASSERT(trans_matrix[i*4+j] <=  1.000001);
					if (trans_matrix[i*4+j] < 0)
						trans_matrix[i*4+j] = 0.0;
					if (trans_matrix[i*4+j] > 1)
						trans_matrix[i*4+j] = 1.0;
				}
				ASSERT(fabs(trans_matrix[i*4]+trans_matrix[i*4+1]+trans_matrix[i*4+2]+trans_matrix[i*4+3]-1.0) < 1e-4);
			}
        }

    } else
        ModelMarkov::computeTransMatrix(time, trans_matrix);
     */
}

