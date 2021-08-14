//
//  alisimulatorheterogeneity.h
//  iqtree
//
//  Created by Nhan Ly-Trong on 23/03/2021.
//

#ifndef alisimulatorheterogeneity_h
#define alisimulatorheterogeneity_h

#include "alisimulator.h"

class AliSimulatorHeterogeneity : public AliSimulator
{
protected:
    /**
        get site-specific rates based on Continuous Gamma Distribution
    */
    virtual void getSiteSpecificRatesContinuousGamma(vector<double> &site_specific_rates, int sequence_length);
    
    /**
        get site-specific rates based on Discrete Distribution (Gamma/FreeRate)
    */
    void getSiteSpecificRatesDiscrete(vector<short int> &new_site_specific_rate_index, vector<double> &site_specific_rates, int sequence_length);
    
    /**
      estimate the state from accumulated trans_matrices
    */
    virtual int estimateStateFromAccumulatedTransMatrices(double *cache_trans_matrix, double site_specific_rate, int site_index, int num_rate_categories, int max_num_states, int dad_state);
    
    /**
      estimate the state from an original trans_matrix
    */
    virtual int estimateStateFromOriginalTransMatrix(ModelSubst *model, int model_component_index, double rate, double *trans_matrix, int max_num_states, double branch_length, int dad_state);
    
    /**
        initialize site specific model index based on its weights in the mixture model
    */
    void intializeSiteSpecificModelIndex(int length, vector<short int> &new_site_specific_model_index);
    
    /**
        initialize caching accumulated_trans_matrix
    */
    void intializeCachingAccumulatedTransMatrices(double *cache_trans_matrix, int num_models, int num_rate_categories, int max_num_states, DoubleVector branch_lengths, double *trans_matrix, ModelSubst* model);
    
    /**
        regenerate sequence based on mixture model component base fequencies
    */
    vector<short int> regenerateSequenceMixtureModel(int length, vector<short int> new_site_specific_model_index);
    
    /**
        simulate a sequence for a node from a specific branch after all variables has been initializing
    */
    virtual void simulateASequenceFromBranchAfterInitVariables(ModelSubst *model, int sequence_length, vector<double> site_specific_rates, double *trans_matrix, int max_num_states, Node *node, NeighborVec::iterator it, string lengths = "");
    
    /**
        initialize variables (e.g., site-specific rate)
    */
    virtual void initVariables(int sequence_length, vector<double> &site_specific_rates);
    
    /**
    *  insert a new sequence into the current sequence
    *
    */
    virtual void insertNewSequenceForInsertionEvent(Node *node, int position, vector<short int> &new_sequence, vector<double> &site_specific_rates);
    
public:
    
    RateHeterogeneity *rate_heterogeneity;
    vector<short int> site_specific_model_index;
    vector<short int> site_specific_rate_index;
    const int RATE_ZERO_INDEX = -1;
    const int RATE_ONE_INDEX = 0;
    
    /**
        constructor
    */
    AliSimulatorHeterogeneity(Params *params);
    
    /**
        constructor
    */
    AliSimulatorHeterogeneity(AliSimulator *alisimulator);
    
    /**
        get site-specific rates
    */
    void getSiteSpecificRates(vector<short int> &new_site_specific_rate_index, vector<double> &site_specific_rates, int sequence_length);

    /**
    *  simulate sequences for all nodes in the tree
    */
    virtual void simulateSeqsForTree(map<string,string> input_msa, string output_filepath);
};

#endif /* alisimulatorheterogeneity_h */
