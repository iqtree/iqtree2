//
//  ratecontinuousgammainvar.cpp
//  model
//
//  Created by Nhan Ly-Trong on 22/03/2021.
//

#include "ratecontinuousgammainvar.h"

RateContinuousGammaInvar::RateContinuousGammaInvar(double shape, int seed, double p_invar_sites) :
RateInvar(p_invar_sites), RateContinuousGamma(shape,seed) {
}

void RateContinuousGammaInvar::setPInvar(double pInvar) {
    p_invar = pInvar;
}

void RateContinuousGammaInvar::startCheckpoint() {
    checkpoint->startStruct("RateGammaInvar");
}

void RateContinuousGammaInvar::saveCheckpoint() {
    RateInvar::saveCheckpoint();
    RateContinuousGamma::saveCheckpoint();
}

void RateContinuousGammaInvar::restoreCheckpoint() {
    // should restore p_invar first before gamma, because RateGamma will call computeRates()
    RateInvar::restoreCheckpoint();
    RateContinuousGamma::restoreCheckpoint();
}

string RateContinuousGammaInvar::getNameParams() {
    return RateInvar::getNameParams() + RateContinuousGamma::getNameParams();
}

void RateContinuousGammaInvar::getSiteSpecificRates(double * site_specific_rates, int sequence_length)
{
    // initialize gamma distribution
    default_random_engine generator;
    generator.seed(seed);
    gamma_distribution<double> distribution(gamma_shape, gamma_shape);
    
    double sum_rate = 0;
    
    for (int i = 0; i < sequence_length; i++)
    {
        // if this site is invariant -> its rate is zero
        if (random_double() <= p_invar)
            site_specific_rates[i] = 0;
        else
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

void RateContinuousGammaInvar::writeInfo(ostream &out) {
    RateInvar::writeInfo(out);
    RateContinuousGamma::writeInfo(out);
}

void RateContinuousGammaInvar::writeParameters(ostream &out) {
    RateInvar::writeParameters(out);
    RateContinuousGamma::writeParameters(out);
}
