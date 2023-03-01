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
        get site-specific on Posterior Mean Rates (Discrete Gamma/FreeRate)
    */
    void getSiteSpecificPosteriorRateHeterogeneity(vector<short int> &new_site_specific_rate_index, vector<double> &site_specific_rates, int sequence_length, IntVector &site_to_patternID);
    
    /**
      estimate the state from accumulated trans_matrices
    */
    virtual int estimateStateFromAccumulatedTransMatrices(double *cache_trans_matrix, double site_specific_rate, int site_index, int num_rate_categories, int dad_state, int* rstream);
    
    /**
      estimate the state from an original trans_matrix
    */
    virtual int estimateStateFromOriginalTransMatrix(ModelSubst *model, int model_component_index, double rate, double *trans_matrix, double branch_length, int dad_state, int site_index, int* rstream);
    
    /**
        initialize site specific model index based on its weights in the mixture model
    */
    void intializeSiteSpecificModelIndex(int length, vector<short int> &new_site_specific_model_index, IntVector &site_to_patternID);
    
    /**
        initialize site specific model index based on posterior model probability
    */
    void intSiteSpecificModelIndexPosteriorProb(int length, vector<short int> &new_site_specific_model_index, IntVector &site_to_patternID);
    
    /**
        initialize caching accumulated_trans_matrix
    */
    void intializeCachingAccumulatedTransMatrices(double *cache_trans_matrix, int num_models, int num_rate_categories, DoubleVector &branch_lengths, double *trans_matrix, ModelSubst* model);
    
    /**
        regenerate sequence based on mixture model component base fequencies
    */
    vector<short int> regenerateSequenceMixtureModel(int length, vector<short int> &new_site_specific_model_index);
    
    /**
        regenerate sequence based on posterior mean state frequencies (for mixture models)
    */
    vector<short int> regenerateSequenceMixtureModelPosteriorMean(int length, IntVector &site_to_patternID);
    
    /**
        simulate a sequence for a node from a specific branch after all variables has been initializing
    */
    virtual void simulateASequenceFromBranchAfterInitVariables(int segment_start, ModelSubst *model, double *trans_matrix, vector<short int> &dad_seq_chunk, vector<short int> &node_seq_chunk, Node *node, NeighborVec::iterator it, int* rstream, string lengths = "");
    
    /**
        initialize variables (e.g., site-specific rate)
    */
    virtual void initVariablesRateHeterogeneity(int sequence_length, bool regenerate_root_sequence = false);
    
    /**
    *  insert a new sequence into the current sequence
    *
    */
    virtual void insertNewSequenceForInsertionEvent(vector<short int> &indel_sequence, int position, vector<short int> &new_sequence);
    
    /**
        initialize variables for Rate_matrix approach: total_sub_rate, accumulated_rates, num_gaps
    */
    virtual void initVariables4RateMatrix(int segment_start, double &total_sub_rate, int &num_gaps, vector<double> &sub_rate_by_site, vector<short int> &sequence);
    
    /**
        extract pattern- posterior mean state frequencies and posterior model probability
    */
    void extractPatternPosteriorFreqsAndModelProb();
    
public:
    
    RateHeterogeneity *rate_heterogeneity;
    
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
    void getSiteSpecificRates(vector<short int> &new_site_specific_rate_index, vector<double> &new_site_specific_rates, vector<short int> &new_site_specific_model_index, int sequence_length, IntVector &site_to_patternID);
};

#endif /* alisimulatorheterogeneity_h */
