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

void RateContinuousGammaInvar::getSiteSpecificRates(vector<double> &site_specific_rates, int sequence_length)
{
    // initialize gamma distribution
    gamma_distribution<double> distribution(gamma_shape, 1/gamma_shape);
    
    for (int i = 0; i < sequence_length; i++)
    {
        // if this site is invariant -> its rate is zero
        if (random_double() <= p_invar)
            site_specific_rates[i] = 0;
        else
            site_specific_rates[i] = distribution(generator);
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
