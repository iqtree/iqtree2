//
//  alisimulatorheterogeneity.cpp
//  model
//
//  Created by Nhan Ly-Trong on 23/03/2021.
//

#include "alisimulatorheterogeneity.h"

AliSimulatorHeterogeneity::AliSimulatorHeterogeneity(Params *params) :
AliSimulator(params) {
    rate_heterogeneity = tree->getRate();
}

AliSimulatorHeterogeneity::AliSimulatorHeterogeneity(AliSimulator *alisimulator){
    tree = alisimulator->tree;
    params = alisimulator->params;
    num_sites_per_state = alisimulator->num_sites_per_state;
    length_ratio = alisimulator->length_ratio;
    inverse_length_ratio = alisimulator->inverse_length_ratio;
    expected_num_sites = alisimulator->expected_num_sites;
    partition_rate = alisimulator->partition_rate;
    rate_heterogeneity = tree->getRate();
    max_length_taxa_name = alisimulator->max_length_taxa_name;
    fundi_items = alisimulator->fundi_items;
    STATE_UNKNOWN = alisimulator->STATE_UNKNOWN;
    max_num_states = alisimulator->max_num_states;
    seq_length_indels = alisimulator->seq_length_indels;
    map_seqname_node = alisimulator->map_seqname_node;
    latest_insertion = alisimulator->latest_insertion;
    first_insertion = alisimulator->first_insertion;
    starting_pos = alisimulator->starting_pos;
    output_line_length = alisimulator->output_line_length;
    num_threads = alisimulator->num_threads;
    force_output_PHYLIP = alisimulator->force_output_PHYLIP;
}

/**
    initialize site specific model index based on its weights in the mixture model
*/
void AliSimulatorHeterogeneity::intializeSiteSpecificModelIndex(int sequence_length, vector<short int> &new_site_specific_model_index, IntVector &site_to_patternID)
{
    new_site_specific_model_index.resize(sequence_length);
    
    // if a mixture model is used -> randomly select a model for each site
    if (tree->getModel()->isMixture())
    {
        // if inference_mode is used -> randomly select a model for each site based on the posterior model probability
        if (tree->params->alisim_inference_mode)
            intSiteSpecificModelIndexPosteriorProb(sequence_length, new_site_specific_model_index, site_to_patternID);
        // otherwise, randomly select a model for each site based on the weights of model components
        else
        {
            // get/init variables
            ModelSubst* model = tree->getModel();
            int num_models = model->getNMixtures();
            mixture_accumulated_weight = new double[num_models];
            
            // get the weights of model components
            bool isFused = model->isFused();
            mixture_max_weight_pos = 0;
            
            // fused model, take the weight from site_rate
            if (isFused)
            {
                double fused_denominator = 1.0 / (1.0 - tree->getRate()->getPInvar());
                for (int i = 0; i < num_models; i++)
                {
                    mixture_accumulated_weight[i] = tree->getRate()->getProp(i) * fused_denominator;
                    
                    // finding the max probability position
                    if (mixture_accumulated_weight[i] > mixture_accumulated_weight[mixture_max_weight_pos])
                        mixture_max_weight_pos = i;
                }
            }
            // otherwise, non-fused models -> take the mixture weight
            else
            {
                for (int i = 0; i < num_models; i++)
                {
                    mixture_accumulated_weight[i] = model->getMixtureWeight(i);
                    
                    // finding the max probability position
                    if (mixture_accumulated_weight[i] > mixture_accumulated_weight[mixture_max_weight_pos])
                        mixture_max_weight_pos = i;
                }
            }
                
            // convert the model_prop into an accumulated model_prop
            convertProMatrixIntoAccumulatedProMatrix(mixture_accumulated_weight, 1, num_models);
            
            for (int i = 0; i < sequence_length; i++)
            {
                // randomly select a model from the set of model components, considering its probability array.
                new_site_specific_model_index[i] = getRandomItemWithAccumulatedProbMatrixMaxProbFirst(mixture_accumulated_weight, 0, num_models, mixture_max_weight_pos, NULL);
            }
            
            // delete the mixture_accumulated_weight if mixture model at substitution level is not used
            if (!params->alisim_mixture_at_sub_level)
            {
                delete[] mixture_accumulated_weight;
                mixture_accumulated_weight = NULL;
            }
        }
    }
    // otherwise, if it's not a mixture model -> set model index = 0 for all sites
    else
    {
        // set model index = 0 for all sites
        for (int i = 0; i < sequence_length; i++)
        {
            new_site_specific_model_index[i] = 0;
        }
    }
}

/**
    initialize site specific model index based on posterior model probability
*/
void AliSimulatorHeterogeneity::intSiteSpecificModelIndexPosteriorProb(int sequence_length, vector<short int> &new_site_specific_model_index, IntVector &site_to_patternID)
{
    // dummy variables
    int nmixture = tree->getModel()->getNMixtures();
    
    // extract pattern- posterior mean state frequencies and posterior model probability
    extractPatternPosteriorFreqsAndModelProb();
    
    ASSERT(site_to_patternID.size() >= sequence_length);
    for (int i = 0; i < sequence_length; i++)
    {
        double rand_num = random_double();
        // extract pattern id from site id
        int site_pattern_id = site_to_patternID[i];
        
        int starting_index = site_pattern_id * nmixture;
        new_site_specific_model_index[i] = binarysearchItemWithAccumulatedProbabilityMatrix(ptn_model_dis, rand_num, starting_index, starting_index + nmixture - 1, starting_index) - starting_index;
    }
    
    // delete ptn_model_dis if we don't need to use it for handling insertions (in Indels)
    if (tree->params->alisim_insertion_ratio + tree->params->alisim_deletion_ratio == 0)
    {
        delete [] ptn_model_dis;
        ptn_model_dis = NULL;
    }
}

/**
    regenerate ancestral sequence based on mixture model component base fequencies
*/
vector<short int> AliSimulatorHeterogeneity::regenerateSequenceMixtureModel(int length, vector<short int> &new_site_specific_model_index){
    // dummy variables
    ModelSubst* model = tree->getModel();
    int num_models = model->getNMixtures();
    int num_states = tree->aln->getMaxNumStates();
    
    // initialize base frequencies maxtrix
    double * base_freqs_all_components = new double[num_models * num_states];
    double * base_freqs_one_component = new double[num_states];
    
    // retrieve base frequencies of each model component
    double* base_freqs_all_components_pointer = base_freqs_all_components;
    for (int i = 0; i < num_models; i++, base_freqs_all_components_pointer += num_states)
    {
        model->getStateFrequency(base_freqs_one_component, i);
        
        // copy base_freqs_one_component to base_freqs_all_components
        for (int j = 0; j < num_states; j++)
            base_freqs_all_components_pointer[j] = base_freqs_one_component[j];
    }
    
    // delete base_freqs_one_component
    delete [] base_freqs_one_component;
    
    // convert base_freqs_all_components to accummulated matrix
    convertProMatrixIntoAccumulatedProMatrix(base_freqs_all_components, num_models, num_states);
    
    // re-generate the sequence
    vector <short int> new_sequence(length, num_states);
    int num_states_minus_one = num_states - 1;
    for (int i = 0; i < length; i++)
    {
        double rand_num = random_double();
        // NHANLT: potential improvement
        // cache new_site_specific_model_index[i]*num_states
        int starting_index = new_site_specific_model_index[i] * num_states;
        new_sequence[i] = binarysearchItemWithAccumulatedProbabilityMatrix(base_freqs_all_components, rand_num, starting_index, starting_index + num_states_minus_one, starting_index) - starting_index;
    }
    
    // delete base_freqs_one_component
    delete [] base_freqs_all_components;
    
    return new_sequence;
}

/**
    extract pattern- posterior mean state frequencies and posterior model probability
*/
void AliSimulatorHeterogeneity::extractPatternPosteriorFreqsAndModelProb()
{
    // get pattern-specific state frequencies (ptn_state_freq)
    int nptn = tree->aln->getNPattern();
    int nmixture = tree->getModel()->getNMixtures();
    if (!ptn_state_freq)
    {
        ptn_state_freq = new double[nptn * max_num_states];
        
        SiteFreqType tmp_site_freq_type = tree->params->print_site_state_freq;
        tree->params->print_site_state_freq = WSF_POSTERIOR_MEAN;
        tree->computePatternStateFreq(ptn_state_freq);
        // get pattern-specific posterior model probability
        int nptn_times_nmixture = nptn * nmixture;
        ptn_model_dis = new double[nptn_times_nmixture];
        memcpy(ptn_model_dis, tree->getPatternLhCatPointer(), nptn_times_nmixture * sizeof(double));
        tree->params->print_site_state_freq = tmp_site_freq_type;
        
        // convert ptn_model_dis to accummulated matrix
        convertProMatrixIntoAccumulatedProMatrix(ptn_model_dis, nptn, nmixture);
    }
}

/**
    regenerate sequence based on posterior mean state frequencies (for mixture models)
*/
vector<short int> AliSimulatorHeterogeneity::regenerateSequenceMixtureModelPosteriorMean(int length, IntVector &site_to_patternID)
{
    ASSERT(tree->params->alisim_stationarity_heterogeneity == POSTERIOR_MEAN);
    
    // extract pattern- posterior mean state frequencies and posterior model probability
    extractPatternPosteriorFreqsAndModelProb();
    
    // init ptn_accumulated_state_freq
    if (!ptn_accumulated_state_freq)
    {
        int nptn = tree->aln->getNPattern();
        int nptn_times_max_num_states = nptn * max_num_states;
        ptn_accumulated_state_freq = new double[nptn_times_max_num_states];
        memcpy(ptn_accumulated_state_freq, ptn_state_freq, nptn_times_max_num_states * sizeof(double));
        
        // convert ptn_state_freq to ptn_accumulated_state_freq
        convertProMatrixIntoAccumulatedProMatrix(ptn_accumulated_state_freq, nptn, max_num_states);
    }
    
    // re-generate the sequence
    vector <short int> new_sequence(length, max_num_states);
    int max_num_states_minus_one = max_num_states - 1;
    for (int i = 0; i < length; i++)
    {
        double rand_num = random_double();
        // extract pattern id from site id
        int site_pattern_id = site_to_patternID[i];
        
        int starting_index = site_pattern_id * max_num_states;
        new_sequence[i] = binarysearchItemWithAccumulatedProbabilityMatrix(ptn_accumulated_state_freq, rand_num, starting_index, starting_index + max_num_states_minus_one, starting_index) - starting_index;
    }
    
    // delete ptn_accumulated_state_freq if we don't need to use it for handling insertions (in Indels)
    if (tree->params->alisim_insertion_ratio + tree->params->alisim_deletion_ratio == 0)
    {
        delete [] ptn_accumulated_state_freq;
        ptn_accumulated_state_freq = NULL;
    }
    
    return new_sequence;
}

/**
    initialize caching accumulated_trans_matrix
*/
void AliSimulatorHeterogeneity::intializeCachingAccumulatedTransMatrices(double *cache_trans_matrix, int num_models, int num_rate_categories, DoubleVector &branch_lengths, double *trans_matrix, ModelSubst* model)
{
    bool fuse_mixture_model = (model->isMixture() && model->isFused());
    
    // initialize the cache_trans_matrix
    double combine_rate = partition_rate * params->alisim_branch_scale;
    double* cache_trans_matrix_pointer = cache_trans_matrix;
    int num_state_square = max_num_states * max_num_states;
    for (int model_index = 0; model_index < num_models; model_index++)
    {
        // Bug fixed
        // in mixture model where each mixture component has a specific rate then we need to use that rate to compute the transition matrix
        // we need to use total_num_subst * total_num_subst (instead of total_num_subst) to cancel "/total_num_subst" in the computeTrans function
        if (model->isMixture())
        {
            double total_num_subst = ((ModelMarkov*) model->getMixtureClass(model_index))->total_num_subst;
            if (fabs(total_num_subst - 1.0) > 1e-6)
                combine_rate *=  total_num_subst * total_num_subst;
        }
        
        for (int category_index = 0; category_index < num_rate_categories; category_index++, cache_trans_matrix_pointer += num_state_square)
        {
            // skip computing unused trans_matrices if a mixture with fused site rate is used
            if (fuse_mixture_model && model_index != category_index)
                continue;
            
            double rate = rate_heterogeneity->getNRate() == 1?1:rate_heterogeneity->getRate(category_index);
            double branch_length_by_category = rate_heterogeneity->isHeterotachy()?branch_lengths[category_index]:branch_lengths[0];
            
            // compute the transition matrix
            model->computeTransMatrix(combine_rate * branch_length_by_category * rate, trans_matrix, model_index);
            
            // copy the transition matrix to the cache_trans_matrix
            for (int trans_index = 0; trans_index < num_state_square; trans_index++)
                cache_trans_matrix_pointer[trans_index] = trans_matrix[trans_index];
        }
    }
    
    // convert cache_trans_matrix into an accumulated cache_trans_matrix
    convertProMatrixIntoAccumulatedProMatrix(cache_trans_matrix, num_models * num_rate_categories * max_num_states, max_num_states);
}

/**
  estimate the state from accumulated trans_matrices
*/
int AliSimulatorHeterogeneity::estimateStateFromAccumulatedTransMatrices(double *cache_trans_matrix, double site_specific_rate, int site_index, int num_rate_categories, int dad_state, int* rstream)
{
    // randomly select the state, considering it's dad states, and the accumulated trans_matrices
    int model_index_times_num_rate_categories = site_specific_model_index[site_index];
    if (model_index_times_num_rate_categories > 0)
        model_index_times_num_rate_categories *= num_rate_categories;
    
    int starting_index = site_specific_rate_index[site_index];
    ASSERT(starting_index > RATE_ZERO_INDEX);
    starting_index += model_index_times_num_rate_categories;
    if (starting_index > 0)
        starting_index *= max_num_states;
    
    starting_index = (starting_index + dad_state) * max_num_states;
  
    return getRandomItemWithAccumulatedProbMatrixMaxProbFirst(cache_trans_matrix, starting_index, max_num_states, dad_state, rstream);
}

/**
  estimate the state from an original trans_matrix
*/
int AliSimulatorHeterogeneity::estimateStateFromOriginalTransMatrix(ModelSubst *model, int model_component_index, double rate, double *trans_matrix, double branch_length, int dad_state, int site_index, int* rstream)
{
    double combine_rate = partition_rate * params->alisim_branch_scale;
    // Bug fixed
    // in mixture model where each mixture component has a specific rate then we need to use that rate to compute the transition matrix
    // we need to use total_num_subst * total_num_subst (instead of total_num_subst) to cancel "/total_num_subst" in the computeTrans function
    if (model->isMixture())
    {
        double total_num_subst = ((ModelMarkov*) model->getMixtureClass(model_component_index))->total_num_subst;
        if (fabs(total_num_subst - 1.0) > 1e-6)
            combine_rate *=  total_num_subst * total_num_subst;
    }
    
    // update state freqs of the model component based on posterior mean site_freqs if needed
    if (model->isMixture() && model->isMixtureSameQ() && params->alisim_stationarity_heterogeneity == POSTERIOR_MEAN)
    {
        ASSERT(site_to_patternID.size() > site_index && ptn_state_freq);
        double *tmp_state_freqs = ptn_state_freq + site_to_patternID[site_index] * max_num_states;
        #ifdef _OPENMP
        #pragma omp critical
        #endif
        {
            model->setStateFrequency(tmp_state_freqs);
            // compute the transition matrix
            model->computeTransMatrix(combine_rate * branch_length * rate, trans_matrix, model_component_index, dad_state);
        }
    }
    // otherwise, only need to compute the transition matrix
    else
        model->computeTransMatrix(combine_rate * branch_length * rate, trans_matrix, model_component_index, dad_state);
    
    // NHANLT: potential improvement
    // cache dad_state * max_num_states
    // iteratively select the state, considering it's dad states, and the transition_probability_matrix
    int starting_index = dad_state * max_num_states;
    return getRandomItemWithProbabilityMatrix(trans_matrix, starting_index, max_num_states, rstream);
}

/**
    get site-specific rates based on Continuous Gamma Distribution
*/
void AliSimulatorHeterogeneity::getSiteSpecificRatesContinuousGamma(vector<double> &site_specific_rates, int sequence_length)
{
    RateContinuousGamma *rate_continuous_gamma = new RateContinuousGamma(rate_heterogeneity->getGammaShape());
    
    rate_continuous_gamma->getSiteSpecificRates(site_specific_rates, sequence_length);
    
    // delete rate_continuous_gamma
    delete rate_continuous_gamma;
}

/**
    get site-specific rates based on Discrete Distribution (Gamma/FreeRate)
*/
void AliSimulatorHeterogeneity::getSiteSpecificRatesDiscrete(vector<short int> &new_site_specific_rate_index, vector<double> &site_specific_rates, int sequence_length)
{
    int num_rate_categories = rate_heterogeneity->getNDiscreteRate();
    
    // initialize the probability array of rate categories
    double *category_probability_matrix = new double[num_rate_categories];
    int max_prob_pos = 0;
    for (int i = 0; i < num_rate_categories; i++)
    {
        category_probability_matrix[i] = rate_heterogeneity->getProp(i);
        
        // finding the max probability position
        if (category_probability_matrix[i] > category_probability_matrix[max_prob_pos])
            max_prob_pos = i;
    }
    
    // convert the probability matrix of rate categories into an accumulated probability matrix of rate categories
    convertProMatrixIntoAccumulatedProMatrix(category_probability_matrix, 1, num_rate_categories, false);
    
    // BUG FIXED
    // normallize the accumulated probability matrix if sum of the weights of all categories is less than 1
    if (category_probability_matrix[num_rate_categories - 1] + tree->getRate()->getPInvar() < 1)
    {
        outWarning("Normalizing weights of rate categories so that sum of them is 1!");
        double inverse_sum = 1.0 / (category_probability_matrix[num_rate_categories - 1] + tree->getRate()->getPInvar());
        double prev_accumulated_weight = 0;
        for (int i = 0; i < num_rate_categories - 1; i++)
        {
            category_probability_matrix[i] *= inverse_sum;
            
            // update the weight of the rate category
            tree->getRate()->setProp(i, category_probability_matrix[i] - prev_accumulated_weight);
            prev_accumulated_weight = category_probability_matrix[i];
        }
        
        // update the weight for the last rate category
        category_probability_matrix[num_rate_categories - 1] = 1 - tree->getRate()->getPInvar();
        tree->getRate()->setProp(num_rate_categories - 1, category_probability_matrix[num_rate_categories - 1] - prev_accumulated_weight);
        
    }
    
    // initialize the site-specific rates
    for (int i = 0; i < sequence_length; i++)
    {
        // randomly select a rate from the set of rate categories, considering its probability array.
        int rate_category = getRandomItemWithAccumulatedProbMatrixMaxProbFirst(category_probability_matrix, 0, num_rate_categories, max_prob_pos, NULL);
        
        // if rate_category == -1 <=> this site is invariant -> return dad's state
        if (rate_category == -1)
        {
            site_specific_rates[i] = 0;
            new_site_specific_rate_index[i] = RATE_ZERO_INDEX;
        }
        else // otherwise, get the rate of that rate_category
        {
            site_specific_rates[i] = rate_heterogeneity->getRate(rate_category);
            new_site_specific_rate_index[i] = rate_category;
        }
    }
    
    // delete the probability array of rate categories
    delete[] category_probability_matrix;
}

/**
    get site-specific on Posterior Mean Rates (Discrete Gamma/FreeRate)
*/
void AliSimulatorHeterogeneity::getSiteSpecificPosteriorRateHeterogeneity(vector<short int> &new_site_specific_rate_index, vector<double> &site_specific_rates, int sequence_length, IntVector &site_to_patternID)
{
    int num_rates = rate_heterogeneity->getNDiscreteRate();
    
    // get pattern-specific rate (ptn_rate)
    if (pattern_rates.size() == 0)
    {
        IntVector pattern_cat;
        tree->getRate()->computePatternRates(pattern_rates, pattern_cat);
        
        // extract pattern rate distribution if the user wants to sample a rate for each site from posterior distribution
        if (tree->params->alisim_rate_heterogeneity == POSTERIOR_DIS)
        {
            // init variables
            int num_ptns = pattern_rates.size();
            int num_ptns_times_num_rates = num_ptns * num_rates;
            ptn_accumulated_rate_dis = new double[num_ptns_times_num_rates];
            
            // clone ptn_rate_dis
            memcpy(ptn_accumulated_rate_dis, tree->getPatternLhCatPointer(), num_ptns_times_num_rates * sizeof(double));
        
            // convert ptn_rate_dis into ptn_accumulated_rate_dis
            convertProMatrixIntoAccumulatedProMatrix(ptn_accumulated_rate_dis, num_ptns, num_rates);
            
            // normalize ptn_accumulated_rate_dis
            int i_times_num_rates = 0;
            int num_rates_minus_one = num_rates - 1;
            for (int i = 0; i < num_ptns; i++, i_times_num_rates += num_rates)
            {
                double inverse_row_sum = 1.0 / ptn_accumulated_rate_dis[i_times_num_rates + num_rates_minus_one];
                for (int j = 0; j < num_rates; j++)
                    ptn_accumulated_rate_dis[i_times_num_rates + j] *= inverse_row_sum;
            }
        }
    }
    
    // init site-specific posterior rate
    // extract posterior mean rate from pattern
    if (tree->params->alisim_rate_heterogeneity == POSTERIOR_MEAN)
        for (int i = 0; i < sequence_length; i++)
        {
            // extract pattern id from site id
            int site_pattern_id = site_to_patternID[i];

            ASSERT(site_pattern_id < pattern_rates.size());
            site_specific_rates[i] = pattern_rates[site_pattern_id];
        }
    // otherwise, sample site rate from posterior distribution
    else if (tree->params->alisim_rate_heterogeneity == POSTERIOR_DIS)
        for (int i = 0; i < sequence_length; i++)
        {
            // extract pattern id from site id
            int site_pattern_id = site_to_patternID[i];
            int starting_index = site_pattern_id * num_rates;
            double rand_num = random_double();
            int rate_cat = binarysearchItemWithAccumulatedProbabilityMatrix(ptn_accumulated_rate_dis, rand_num, starting_index, starting_index + num_rates - 1, starting_index) - starting_index;
            site_specific_rates[i] = rate_heterogeneity->getRate(rate_cat);
            new_site_specific_rate_index[i] = rate_cat;
        }
    
    // delete ptn_accumulated_rate_dis if we don't need to use it for handling insertions (in Indels)
    if (ptn_accumulated_rate_dis && tree->params->alisim_insertion_ratio + tree->params->alisim_deletion_ratio == 0)
    {
        delete [] ptn_accumulated_rate_dis;
        ptn_accumulated_rate_dis = NULL;
    }
}

/**
    get site-specific rates
*/
void AliSimulatorHeterogeneity::getSiteSpecificRates(vector<short int> &new_site_specific_rate_index, vector<double> &site_specific_rates, vector<short int> &new_site_specific_model_index, int sequence_length, IntVector &site_to_patternID)
{
    new_site_specific_rate_index.resize(sequence_length);
    site_specific_rates.resize(sequence_length, 1);
    
    // if a mixture model is supplied and it's fused with site rates -> set site_specific_rate_index equals to site_specific_model_index
    if (tree->getModel()->isMixture() && tree->getModel()->isFused())
    {
        // get invariant_probability
        double invariant_prop = tree->getRate()->getPInvar();

        for (int i = 0; i < sequence_length; i++)
        {
            // handle invariant sites
            if (random_double() <= invariant_prop)
            {
                new_site_specific_rate_index[i] = RATE_ZERO_INDEX;
                site_specific_rates[i] = 0;
            }
            // or set the rate index equal to the model index
            else
            {
                new_site_specific_rate_index[i] = new_site_specific_model_index[i];
                site_specific_rates[i] = rate_heterogeneity->getRate(new_site_specific_rate_index[i]);
            }
        }
        return;
    }
    
    string rate_name = tree->getRateName();
    
    // mixture model without site rate heterogeneity
    if (rate_name.empty())
    {
        // initialize all site's rate equally at 1
        for (int i = 0; i < sequence_length; i++)
        {
            site_specific_rates[i] = 1;
            new_site_specific_rate_index[i] = RATE_ONE_INDEX;
        }
    }
    // otherwise, it's the case with site rate heterogeneity
    else
    {
        // initalize rates based on continuous gamma distribution
        if ((rate_name.find("+G") != std::string::npos) && tree->getModelFactory()->is_continuous_gamma)
        {
            getSiteSpecificRatesContinuousGamma(site_specific_rates, sequence_length);
        }
        // initalize rates based on discrete distribution (gamma/freerate)
        else
        {
            if (applyPosRateHeterogeneity)
                getSiteSpecificPosteriorRateHeterogeneity(new_site_specific_rate_index, site_specific_rates, sequence_length, site_to_patternID);
            else
                getSiteSpecificRatesDiscrete(new_site_specific_rate_index, site_specific_rates, sequence_length);
        }
    }
}

/**
    simulate a sequence for a node from a specific branch after all variables has been initializing
*/
void AliSimulatorHeterogeneity::simulateASequenceFromBranchAfterInitVariables(int segment_start, ModelSubst *model, double *trans_matrix, vector<short int> &dad_seq_chunk, vector<short int> &node_seq_chunk, Node *node, NeighborVec::iterator it, int* rstream, string lengths){
    
    // estimate the sequence for the current neighbor
    // check if trans_matrix could be caching (without rate_heterogeneity or the num of rate_categories is lowr than the threshold (5)) or not
    if ((tree->getRateName().empty()
        || (!tree->getModelFactory()->is_continuous_gamma && !(applyPosRateHeterogeneity && params->alisim_rate_heterogeneity == POSTERIOR_MEAN) && rate_heterogeneity && rate_heterogeneity->getNDiscreteRate() <= params->alisim_max_rate_categories_for_applying_caching))
        && !(model->isMixture() && model->isMixtureSameQ() && params->alisim_stationarity_heterogeneity == POSTERIOR_MEAN))
    {
        int num_models = tree->getModel()->isMixture()?tree->getModel()->getNMixtures():1;
        int num_rate_categories  = tree->getRateName().empty()?1:rate_heterogeneity->getNDiscreteRate();
        double *cache_trans_matrix = new double[num_models * num_rate_categories * max_num_states * max_num_states];
        
        // initialize a set of branch_lengths
        DoubleVector branch_lengths;
        // if heterotachy model is used in branch-specific model -> parse multiple lengths from the input string 'lengths'
        if (rate_heterogeneity->isHeterotachy() && lengths.length() > 0)
        {
            // parse lengths
            convert_double_vec_with_distributions(lengths.c_str(), branch_lengths, '/');
            
            if (num_rate_categories != branch_lengths.size())
                outError("The number of lengths ("+convertIntToString(branch_lengths.size())+") is different from the number of caterogies ("+convertIntToString(num_rate_categories)+"). Please check and try again!");
        }
        // otherwise, get branch-length from the tree
        else
        {
            branch_lengths.resize(num_rate_categories);
            for (int i = 0; i < num_rate_categories; i++)
                branch_lengths[i] = (*it)->getLength(i);
        }
        
        // initialize caching accumulated trans_matrices
        intializeCachingAccumulatedTransMatrices(cache_trans_matrix, num_models, num_rate_categories, branch_lengths, trans_matrix, model);

        // estimate the sequence
        for (int i = 0 ; i < node_seq_chunk.size(); i++)
        {
            // if the parent's state is a gap -> the children's state should also be a gap
            if (dad_seq_chunk[i] == STATE_UNKNOWN)
                node_seq_chunk[i] = STATE_UNKNOWN;
            else
            {
                node_seq_chunk[i] = estimateStateFromAccumulatedTransMatrices(cache_trans_matrix, site_specific_rates[segment_start + i] , segment_start + i, num_rate_categories, dad_seq_chunk[i], rstream);
            }
        }
        
        // delete cache_trans_matrix
        delete [] cache_trans_matrix;
    }
    // otherwise, estimating the sequence without trans_matrix caching
    else
    {
        for (int i = 0 ; i < node_seq_chunk.size(); i++)
        {
            // if the parent's state is a gap -> the children's state should also be a gap
            if (dad_seq_chunk[i] == STATE_UNKNOWN)
                node_seq_chunk[i] = STATE_UNKNOWN;
            else
            {
                // randomly select the state, considering it's dad states, and the transition_probability_matrix
                node_seq_chunk[i] = estimateStateFromOriginalTransMatrix(model, site_specific_model_index[segment_start + i], site_specific_rates[segment_start + i], trans_matrix, (*it)->length, dad_seq_chunk[i], segment_start + i, rstream);
            }
        }
    }
}

/**
    initialize variables (e.g., site-specific rate)
*/
void AliSimulatorHeterogeneity::initVariablesRateHeterogeneity(int sequence_length, bool regenerate_root_sequence)
{    
    // initialize site specific model index based on its weights (in the mixture model)
    intializeSiteSpecificModelIndex(sequence_length, site_specific_model_index, site_to_patternID);
    
    // only regenerate the ancestral sequence if mixture model is used and the ancestral sequence is not specified by the user.
    if (regenerate_root_sequence && tree->getModel()->isMixture() && !tree->params->alisim_ancestral_sequence_aln_filepath)
    {
        // re-generate sequence based on posterior mean/distribution state frequencies if users want to do so
        if (tree->getModel()->isMixtureSameQ() && tree->params->alisim_stationarity_heterogeneity == POSTERIOR_MEAN)
            tree->root->sequence->sequence_chunks[0] = regenerateSequenceMixtureModelPosteriorMean(expected_num_sites, site_to_patternID);
        // otherwise re-generate sequence based on the state frequencies the model component for each site
        else
            tree->root->sequence->sequence_chunks[0] = regenerateSequenceMixtureModel(expected_num_sites, site_specific_model_index);
    
        // separate root sequence into chunks
        separateSeqIntoChunks(tree->root);
    }

    
    // initialize site-specific rates
    getSiteSpecificRates(site_specific_rate_index, site_specific_rates, site_specific_model_index, sequence_length, site_to_patternID);
}

/**
*  insert a new sequence into the current sequence
*
*/
void AliSimulatorHeterogeneity::insertNewSequenceForInsertionEvent(vector<short int> &indel_sequence, int position, vector<short int> &new_sequence)
{
    // init new_site_to_patternID
    IntVector new_site_to_patternID;
    if (tree->params->alisim_inference_mode && (tree->params->alisim_rate_heterogeneity!=UNSPECIFIED || tree->params->alisim_stationarity_heterogeneity!=UNSPECIFIED))
    {
        new_site_to_patternID.resize(new_sequence.size());
        
        // randomly pick a pattern for each site
        int site_id;
        for (int i = 0; i < new_sequence.size(); i++)
        {
            site_id = random_int(site_to_patternID.size());
            new_site_to_patternID[i] = site_to_patternID[site_id];
        }
        
        // insert new_site_to_patternID into site_to_patternID
        site_to_patternID.insert(site_to_patternID.begin()+position, new_site_to_patternID.begin(), new_site_to_patternID.end());
    }
    
    // initialize new_site_specific_model_index
    vector<short int> new_site_specific_model_index;
    intializeSiteSpecificModelIndex(new_sequence.size(), new_site_specific_model_index, new_site_to_patternID);
    
    // insert new_site_specific_model_index into site_specific_model_index
    site_specific_model_index.insert(site_specific_model_index.begin()+position, new_site_specific_model_index.begin(), new_site_specific_model_index.end());
    
    // initialize new_site_specific_rates, and new_site_specific_rate_index for new sequence
    vector<double> new_site_specific_rates;
    vector<short int> new_site_specific_rate_index;
    getSiteSpecificRates(new_site_specific_rate_index, new_site_specific_rates, new_site_specific_model_index, new_sequence.size(), new_site_to_patternID);
    
    // insert new_site_specific_rates into site_specific_rates
    site_specific_rates.insert(site_specific_rates.begin()+position, new_site_specific_rates.begin(), new_site_specific_rates.end());
    
    // insert new_site_specific_rate_index into site_specific_rate_index
    site_specific_rate_index.insert(site_specific_rate_index.begin()+position, new_site_specific_rate_index.begin(), new_site_specific_rate_index.end());
    
    // regenerate new_sequence if mixture model is used
    if (tree->getModel()->isMixture())
    {
        // re-generate sequence based on posterior mean/distribution state frequencies if users want to do so
        if (tree->getModel()->isMixtureSameQ() && tree->params->alisim_stationarity_heterogeneity == POSTERIOR_MEAN)
            new_sequence = regenerateSequenceMixtureModelPosteriorMean(new_site_specific_model_index.size(), new_site_to_patternID);
        // otherwise re-generate sequence based on the state frequencies the model component for each site
        else
            new_sequence = regenerateSequenceMixtureModel(new_site_specific_model_index.size(), new_site_specific_model_index);
    }
    
    // insert new_sequence into the current sequence
    AliSimulator::insertNewSequenceForInsertionEvent(indel_sequence, position, new_sequence);
}

/**
    initialize variables for Rate_matrix approach: total_sub_rate, accumulated_rates, num_gaps
*/
void AliSimulatorHeterogeneity::initVariables4RateMatrix(int segment_start, double &total_sub_rate, int &num_gaps, vector<double> &sub_rate_by_site, vector<short int> &sequence)
{
    // initialize variables
    total_sub_rate = 0;
    num_gaps = 0;
    sub_rate_by_site.resize(sequence.size(), 0);
    
    // check if sub_rates could be caching (without continuous gamma and not use Posterior Mean Rates) -> compute sub_rate_by_site efficiently using cache_sub_rates
    if (!tree->getModelFactory()->is_continuous_gamma && !applyPosRateHeterogeneity)
    {
        int num_models = tree->getModel()->isMixture()?tree->getModel()->getNMixtures():1;
        int num_rate_categories  = tree->getRateName().empty()?1:rate_heterogeneity->getNDiscreteRate();
        int num_categories_times_num_states = num_rate_categories * max_num_states;
        int total_elements = num_models * num_categories_times_num_states;
        double *cache_sub_rates = new double[total_elements];
        vector<int> sub_rate_count(total_elements, 0);
        
        // initialize cache_sub_rates
        bool fuse_mixture_model = (tree->getModel()->isMixture() && tree->getModel()->isFused());
        int overall_index = 0;
        int model_index_times_num_states = 0;

        for (int model_index = 0; model_index < num_models; model_index++, model_index_times_num_states += max_num_states)
        {
            for (int rate_category_index = 0; rate_category_index < num_rate_categories; rate_category_index++, overall_index += max_num_states)
            {
                // skip computing unused cache_sub_rates if a mixture with fused site rate is used
                if (fuse_mixture_model && model_index != rate_category_index)
                    continue;
                
                // extract site's rate
                double rate = rate_heterogeneity->getNRate() == 1 ? 1 : rate_heterogeneity->getRate(rate_category_index);
                
                // compute cache_sub_rates
                if (rate == 1)
                {
                    for (int state = 0; state < max_num_states; state++)
                        cache_sub_rates[overall_index + state] = sub_rates[model_index_times_num_states + state];
                }
                else
                {
                    for (int state = 0; state < max_num_states; state++)
                        cache_sub_rates[overall_index + state] = sub_rates[model_index_times_num_states + state] * rate;
                }
            }
        }
        
        // NHANLT: potential improvement
        // cache rate_index * max_num_states
        // compute sub_rate_by_site
        int segment_start_plus_i = segment_start;
        for (int i = 0; i < sequence.size(); i++, ++segment_start_plus_i)
        {
            // not compute the substitution rate for gaps/deleted sites or constant sites
            if (sequence[i] != STATE_UNKNOWN && site_specific_rates[segment_start_plus_i] != 0)
            {
                // get the mixture model index and site_specific_rate_index
                int model_index = site_specific_model_index[segment_start_plus_i];
                int rate_index = site_specific_rate_index[segment_start_plus_i];
                
                // update sub_rate_by_site for the current site
                int index = (model_index == 0 ? 0 : (model_index * num_categories_times_num_states)) + (rate_index == 0 ? 0 : (rate_index * max_num_states)) + sequence[i];
                sub_rate_count[index]++;
                sub_rate_by_site[i] = cache_sub_rates[index];
            }
            else
            {
                sub_rate_by_site[i] = 0;
                
                if (sequence[i] == STATE_UNKNOWN)
                    num_gaps++;
            }
        }
        
        // update total_sub_rate
        for (int i = 0; i < total_elements; i++)
            total_sub_rate += sub_rate_count[i] * cache_sub_rates[i];
        
        // delete cache_sub_rates
        delete[] cache_sub_rates;
    }
    // otherwise, sub_rate_by_site for all sites one by one
    else
    {
        int segment_start_plus_i = segment_start;
        for (int i = 0; i < sequence.size(); i++, ++segment_start_plus_i)
        {
            // not compute the substitution rate for gaps/deleted sites or constant sites
            if (sequence[i] != STATE_UNKNOWN && site_specific_rates[segment_start_plus_i] != 0)
            {
                // get the mixture model index
                int model_index_times_num_states = site_specific_model_index[segment_start_plus_i];
                model_index_times_num_states = model_index_times_num_states == 0 ? 0 : (model_index_times_num_states * max_num_states);

                
                // update sub_rate_by_site for the current site
                sub_rate_by_site[i] = sub_rates[model_index_times_num_states + sequence[i]] * site_specific_rates[segment_start_plus_i];
            }
            else
            {
                sub_rate_by_site[i] = 0;
                
                if (sequence[i] == STATE_UNKNOWN)
                    num_gaps++;
            }
            
            // update total_sub_rate
            total_sub_rate += sub_rate_by_site[i];
        }
    }
}
