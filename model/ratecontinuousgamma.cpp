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

RateContinuousGamma::RateContinuousGamma(double shape): RateHeterogeneity()
{
    gamma_shape = shape;
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

void RateContinuousGamma::getSiteSpecificRates(vector<double> &site_specific_rates, int sequence_length)
{
    // initialize gamma distribution
    gamma_distribution<double> distribution(gamma_shape, 1/gamma_shape);
    default_random_engine generator = Params::getInstance().generator;
    
    for (int i = 0; i < sequence_length; i++)
    {
        site_specific_rates[i] = distribution(generator);
    }
}

void RateContinuousGamma::writeInfo(ostream &out) {
    out << "Gamma shape alpha: " << gamma_shape << endl;
}

void RateContinuousGamma::writeParameters(ostream &out) {
    out << "\t" << gamma_shape;
}

