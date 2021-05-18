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
    expected_num_sites = alisimulator->expected_num_sites;
    partition_rate = alisimulator->partition_rate;
    rate_heterogeneity = tree->getRate();
}

/**
    initialize site specific model index based on its weights in the mixture model
*/
void AliSimulatorHeterogeneity::intializeSiteSpecificModelIndex()
{
    int sequence_length = expected_num_sites;
    site_specific_model_index.resize(sequence_length);
    
    // if a mixture model is used -> randomly select a model for each site based on the weights of model components
    if (tree->getModel()->isMixture())
    {
        // get/init variables
        ModelSubst* model = tree->getModel();
        int num_models = model->getNMixtures();
        double *model_prop = new double[num_models];
        
        // get the weights of model components
        bool isFused = model->isFused();
        for (int i = 0; i < num_models; i++)
        {
            // fused model, take the weight from site_rate
            if (isFused)
                model_prop[i] = tree->getRate()->getProp(i) / (1.0 - tree->getRate()->getPInvar());
            else
                model_prop[i] = model->getMixtureWeight(i);
        }
            
        // convert the model_prop into an accumulated model_prop
        convertProMatrixIntoAccumulatedProMatrix(model_prop, 1, num_models);
        
        for (int i = 0; i < sequence_length; i++)
        {
            // randomly select a model from the set of model components, considering its probability array.
            site_specific_model_index[i] = getRandomItemWithAccumulatedProbabilityMatrix(model_prop, 0, num_models);
        }
        
        // delete the probability array of rate categories
        delete[] model_prop;
    }
    // otherwise, if it's not a mixture model -> set model index = 0 for all sites
    else
    {
        // set model index = 0 for all sites
        for (int i = 0; i < sequence_length; i++)
        {
            site_specific_model_index[i] = 0;
        }
    }
}

/**
    initialize state freqs for all model components (of a mixture model)
*/
void AliSimulatorHeterogeneity::intializeStateFreqsMixtureModel()
{
    // get/init variables
    ModelSubst* model = tree->getModel();
    
    // only initialize state freqs if it's a mixture model && the state freqs have not been estimated by an inference process yet
    if (model->isMixture() && tree->aln->aln_file.length() == 0 && model->getFreqType() == FREQ_EMPIRICAL)
    {
        // get max_num_bases
        int max_num_states = tree->aln->getMaxNumStates();
        
        // initialize state freqs
        double *state_freq = new double[max_num_states];
        
        // get the weights of model components
        for (int i = 0; i < model->getNMixtures(); i++)
            if (model->getMixtureClass(i)->getFreqType() == FREQ_EMPIRICAL)
            {
                generateRandomBaseFrequencies(state_freq, max_num_states);
                model->getMixtureClass(i)->setStateFrequency(state_freq);
            }
        
        // delete state_freq
        delete [] state_freq;
    }
}

/**
    initialize caching accumulated_trans_matrix
*/
void AliSimulatorHeterogeneity::intializeCachingAccumulatedTransMatrices(double *cache_trans_matrix, int num_models, int num_rate_categories, int max_num_states, DoubleVector branch_lengths, double *trans_matrix, ModelSubst* model)
{
    bool fuse_mixture_model = (model->isMixture() && model->isFused());
    
    // initialize the cache_trans_matrix
    for (int model_index = 0; model_index < num_models; model_index++)
        for (int category_index = 0; category_index < num_rate_categories; category_index++)
        {
            // skip computing unused trans_matrices if a mixture with fused site rate is used
            if (fuse_mixture_model && model_index != category_index)
                continue;
            
            double rate = tree->getRateName().empty()?1:rate_heterogeneity->getRate(category_index);
            double branch_length_by_category = rate_heterogeneity->isHeterotachy()?branch_lengths[category_index]:branch_lengths[0];
            
            // compute the transition matrix
            model->computeTransMatrix(partition_rate*branch_length_by_category*rate, trans_matrix, model_index);
            
            // copy the transition matrix to the cache_trans_matrix
            for (int trans_index = 0; trans_index < max_num_states*max_num_states; trans_index++)
            {
                int cache_index = model_index*num_rate_categories*max_num_states*max_num_states + category_index*max_num_states*max_num_states + trans_index;
                cache_trans_matrix[cache_index] = trans_matrix[trans_index];
            }
        }
    
    // convert cache_trans_matrix into an accumulated cache_trans_matrix
    convertProMatrixIntoAccumulatedProMatrix(cache_trans_matrix, num_models*num_rate_categories*max_num_states, max_num_states);
}

/**
*  simulate sequences for all nodes in the tree by DFS
*
*/
void AliSimulatorHeterogeneity::simulateSeqs(int sequence_length, double *site_specific_rates, ModelSubst *model, double *trans_matrix, int max_num_states, Node *node, Node *dad, int thread_id, int num_threads)
{
    // process its neighbors/children
    NeighborVec::iterator it;
    FOR_NEIGHBOR(node, dad, it) {
        
        // estimate the sequence for the current neighbor
        if ((*it)->node->sequence.size() != sequence_length)
        {
#ifdef _OPENMP
#pragma omp critical
#endif
            if ((*it)->node->sequence.size() != sequence_length)
                (*it)->node->sequence.resize(sequence_length);
        }
        
        // check if trans_matrix could be caching (without rate_heterogeneity or the num of rate_categories is lowr than the threshold (5)) or not
        if (tree->getRateName().empty()
            || (!params->alisim_continuous_gamma && rate_heterogeneity && rate_heterogeneity->getNDiscreteRate() <= params->alisim_max_rate_categories_for_applying_caching))
        {
            int num_models = tree->getModel()->isMixture()?tree->getModel()->getNMixtures():1;
            int num_rate_categories  = tree->getRateName().empty()?1:rate_heterogeneity->getNDiscreteRate();
            double *cache_trans_matrix = new double[num_models*num_rate_categories*max_num_states*max_num_states];
            
            // initialize a set of branch_lengths
            DoubleVector branch_lengths;
            branch_lengths.resize(num_rate_categories);
            for (int i = 0; i < num_rate_categories; i++)
                branch_lengths[i] = (*it)->getLength(i);
            
            // initialize caching accumulated trans_matrices
            intializeCachingAccumulatedTransMatrices(cache_trans_matrix, num_models, num_rate_categories, max_num_states, branch_lengths, trans_matrix, model);

            // estimate the sequence
            for (int i = thread_id; i < sequence_length; i += num_threads)
                (*it)->node->sequence[i] = estimateStateFromAccumulatedTransMatrices(cache_trans_matrix, site_specific_rates[i] , i, num_rate_categories, max_num_states, node->sequence[i]);
            
            // delete cache_trans_matrix
            delete [] cache_trans_matrix;
        }
        // otherwise, estimating the sequence without trans_matrix caching
        else
        {
            for (int i = thread_id; i < sequence_length; i += num_threads)
               // randomly select the state, considering it's dad states, and the transition_probability_matrix
                (*it)->node->sequence[i] = estimateStateFromOriginalTransMatrix(model, site_specific_model_index[i], site_specific_rates[i], trans_matrix, max_num_states, (*it)->length, node->sequence[i]);
        }
        
        // browse 1-step deeper to the neighbor node
        simulateSeqs(sequence_length, site_specific_rates, model, trans_matrix, max_num_states, (*it)->node, node, thread_id, num_threads);
    }
}

/**
  estimate the state from accumulated trans_matrices
*/
int AliSimulatorHeterogeneity::estimateStateFromAccumulatedTransMatrices(double *cache_trans_matrix, double site_specific_rate, int site_index, int num_rate_categories, int max_num_states, int dad_state)
{
    // randomly select the state, considering it's dad states, and the accumulated trans_matrices
    int model_index = site_specific_model_index[site_index];
    int category_index = site_specific_rate_index[site_index];
    int starting_index = model_index*num_rate_categories*max_num_states*max_num_states + category_index*max_num_states*max_num_states + max_num_states*dad_state;
    
    ASSERT(category_index > RATE_ZERO_INDEX);
    
    return getRandomItemWithAccumulatedProbabilityMatrix(cache_trans_matrix, starting_index, max_num_states);
}

/**
  estimate the state from an original trans_matrix
*/
int AliSimulatorHeterogeneity::estimateStateFromOriginalTransMatrix(ModelSubst *model, int model_component_index, double rate, double *trans_matrix, int max_num_states, double branch_length, int dad_state)
{
    // compute the transition matrix
    model->computeTransMatrix(partition_rate*branch_length*rate, trans_matrix, model_component_index);
    
    // iteratively select the state, considering it's dad states, and the transition_probability_matrix
    int starting_index = dad_state*max_num_states;
    return getRandomItemWithProbabilityMatrix(trans_matrix, starting_index, max_num_states);
}

/**
    get site-specific rates based on Continuous Gamma Distribution
*/
void AliSimulatorHeterogeneity::getSiteSpecificRatesContinuousGamma(double *site_specific_rates, int sequence_length)
{
    RateContinuousGamma *rate_continuous_gamma = new RateContinuousGamma(rate_heterogeneity->getGammaShape(), params->ran_seed);
    
    rate_continuous_gamma->getSiteSpecificRates(site_specific_rates, sequence_length);
    
    // delete rate_continuous_gamma
    delete rate_continuous_gamma;
}

/**
    get site-specific rates based on Discrete Distribution (Gamma/FreeRate)
*/
void AliSimulatorHeterogeneity::getSiteSpecificRatesDiscrete(double *site_specific_rates, int sequence_length)
{
    int num_rate_categories = rate_heterogeneity->getNDiscreteRate();
    
    // initialize the probability array of rate categories
    double *category_probability_matrix = new double[num_rate_categories];
    for (int i = 0; i < num_rate_categories; i++)
        category_probability_matrix[i] = rate_heterogeneity->getProp(i);
    
    // convert the probability matrix of rate categories into an accumulated probability matrix of rate categories
    convertProMatrixIntoAccumulatedProMatrix(category_probability_matrix, 1, num_rate_categories);
    
    // initialize the site-specific rates
    for (int i = 0; i < sequence_length; i++)
    {
        // randomly select a rate from the set of rate categories, considering its probability array.
        int rate_category = getRandomItemWithAccumulatedProbabilityMatrix(category_probability_matrix, 0, num_rate_categories);
        
        // if rate_category == -1 <=> this site is invariant -> return dad's state
        if (rate_category == -1)
        {
            site_specific_rates[i] = 0;
            site_specific_rate_index[i] = RATE_ZERO_INDEX;
        }
        else // otherwise, get the rate of that rate_category
        {
            site_specific_rates[i] = rate_heterogeneity->getRate(rate_category);
            site_specific_rate_index[i] = rate_category;
        }
    }
    
    // delete the probability array of rate categories
    delete[] category_probability_matrix;
}

/**
    get site-specific rates
*/
void AliSimulatorHeterogeneity::getSiteSpecificRates(double *site_specific_rates, int sequence_length)
{
    site_specific_rate_index.resize(sequence_length);
    
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
                site_specific_rate_index[i] = RATE_ZERO_INDEX;
                site_specific_rates[i] = 0;
            }
            // or set the rate index equal to the model index
            else
            {
                site_specific_rate_index[i] = site_specific_model_index[i];
                site_specific_rates[i] = rate_heterogeneity->getRate(site_specific_rate_index[i]);
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
            site_specific_rate_index[i] = RATE_ONE_INDEX;
        }
    }
    // otherwise, it's the case with site rate heterogeneity
    else
    {
        // initalize rates based on continuous gamma distribution
        if ((rate_name.find("+G") != std::string::npos) && params->alisim_continuous_gamma)
        {
            getSiteSpecificRatesContinuousGamma(site_specific_rates, sequence_length);
        }
        // initalize rates based on discrete distribution (gamma/freerate)
        else
        {
            getSiteSpecificRatesDiscrete(site_specific_rates, sequence_length);
        }
    }
}

/**
*  simulate sequences for all nodes in the tree
*/
void AliSimulatorHeterogeneity::simulateSeqsForTree(){
    // get variables
    int sequence_length = expected_num_sites;
    ModelSubst *model = tree->getModel();
    int max_num_states = tree->aln->getMaxNumStates();
    
    // initialize site specific model index based on its weights (in the mixture model)
    intializeSiteSpecificModelIndex();
    
    // initialize site-specific rates
    double *site_specific_rates = new double[sequence_length];
    getSiteSpecificRates(site_specific_rates, sequence_length);
    
    // initialize trans_matrix
    double *trans_matrix;
    
    // simulate sequences
    int num_threads = 1;
    int thread_id = 0;
#ifdef _OPENMP
#pragma omp parallel private(thread_id, trans_matrix)
#endif
    {
#ifdef _OPENMP
#pragma omp single
        num_threads = omp_get_num_threads();
        
        thread_id = omp_get_thread_num();
#endif

        trans_matrix = new double[max_num_states*max_num_states];
        
        simulateSeqs(sequence_length, site_specific_rates, model, trans_matrix, max_num_states, tree->MTree::root, tree->MTree::root, thread_id, num_threads);
    
        // delete trans_matrix array
        delete[] trans_matrix;
    }
    
    // delete the site-specific rates
    delete[] site_specific_rates;
    
    // removing constant states if it's necessary
    if (params->alisim_length_ratio > 1)
        removeConstantSites();
}
