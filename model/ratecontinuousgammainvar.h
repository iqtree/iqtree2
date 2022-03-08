//
//  ratecontinuousgammainvar.h
//  iqtree
//
//  Created by Nhan Ly-Trong on 22/03/2021.
//

#ifndef ratecontinuousgammainvar_h
#define ratecontinuousgammainvar_h

#include "rateinvar.h"
#include "ratecontinuousgamma.h"

class RateContinuousGammaInvar : public RateInvar, public RateContinuousGamma
{
public:
     /**
        constructor
        @param alpha and p_invar
    */
    RateContinuousGammaInvar(double shape, double p_invar_sites);

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

    /**
        set the proportion of invariable sites. Default: do nothing
        @param pinv the proportion of invariable sites
    */
    virtual void setPInvar(double pInvar);

    /**
     * @return model name with parameters in form of e.g. GTR{a,b,c,d,e,f}
     */
    virtual string getNameParams();
    
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

};

#endif /* ratecontinuousgammainvar_h */
