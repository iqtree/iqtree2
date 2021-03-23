//
//  ratecontinuousgamma.cpp
//  model
//
//  Created by Nhan Ly-Trong on 22/03/2021.
//

#include "ratecontinuousgamma.h"

RateContinuousGamma::RateContinuousGamma(): RateHeterogeneity()
{
    
}

RateContinuousGamma::RateContinuousGamma(double shape, int rand_seed): RateHeterogeneity()
{
    gamma_shape = shape;
    seed = rand_seed;
}

void RateContinuousGamma::startCheckpoint() {
    checkpoint->startStruct("RateContinuousGamma");
}

void RateContinuousGamma::saveCheckpoint() {
    startCheckpoint();
    CKP_SAVE(gamma_shape);
    checkpoint->endStruct();
    RateContinuousGamma::saveCheckpoint();
}

void RateContinuousGamma::restoreCheckpoint() {
    RateContinuousGamma::restoreCheckpoint();
    startCheckpoint();
    CKP_RESTORE(gamma_shape);
    checkpoint->endStruct();
}

void RateContinuousGamma::setGammaShape(double gs) {
    gamma_shape = gs;
}

void RateContinuousGamma::getSiteSpecificRates(double * site_specific_rates, int sequence_length)
{
    // initialize gamma distribution
    default_random_engine generator;
    generator.seed(seed);
    gamma_distribution<double> distribution(gamma_shape, gamma_shape);
    
    double sum_rate = 0;
    
    for (int i = 0; i < sequence_length; i++)
    {
        site_specific_rates[i] = distribution(generator);
        
        // update sum_rate
        sum_rate += site_specific_rates[i];
    }
    
    // compute mean_rate
    double mean_rate = sum_rate/sequence_length;
    
    // normalize the rates
    for (int i = 0; i < sequence_length; i++)
    {
        site_specific_rates[i] /= mean_rate;
    }
}

void RateContinuousGamma::writeInfo(ostream &out) {
    out << "Gamma shape alpha: " << gamma_shape << endl;
}

void RateContinuousGamma::writeParameters(ostream &out) {
    out << "\t" << gamma_shape;
}

