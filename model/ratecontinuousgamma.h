//
//  ratecontinuousgamma.h
//  iqtree
//
//  Created by Nhan Ly-Trong on 22/03/2021.
//

#ifndef ratecontinuousgamma_h
#define ratecontinuousgamma_h

#include "rateheterogeneity.h"
#include <random>

class RateContinuousGamma: virtual public RateHeterogeneity
{

    friend class RateContinuousGammaInvar;
    
public:
    /**
        constructor
    */
    RateContinuousGamma();
    /**
        constructor
        @param shape Gamma shape parameter
    */
    RateContinuousGamma(double shape);

    /**
        start structure for checkpointing
    */
    virtual void startCheckpoint();

    /**
        save object into the checkpoint
    */
    virtual void saveCheckpoint();

    /**
        restore object from the checkpoint
    */
    virtual void restoreCheckpoint();

    virtual double getGammaShape() { return gamma_shape; }

    virtual void setGammaShape(double gs);

    /**
        @return site-specific rates
    */

    virtual void getSiteSpecificRates(vector<double> &site_specific_rates, int sequence_length);

    /**
        write information
        @param out output stream
    */
    virtual void writeInfo(ostream &out);

    /**
        write parameters, used with modeltest
        @param out output stream
    */
    virtual void writeParameters(ostream &out);

protected:
    /**
        the gamma shape parameter 'alpha'
    */
    double gamma_shape;

};

#endif /* ratecontinuousgamma_h */
