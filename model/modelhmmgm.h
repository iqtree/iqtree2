//
//  modelhmmgm.h
//  model
//
//  Created by Thomas Wong on 03/02/23.
//

#ifndef modelhmmgm_h
#define modelhmmgm_h

#include "modelhmm.h"

#define MIN_VALUE 1e-5
#define MAX_VALUE 10000

using namespace std;

class PhyloHmm;

class ModelHmmGm: public ModelHmm
{
public:
    /**
     constructor
     @param ncat : number of categories
     */
    ModelHmmGm(int numcat);
    
    /**
     destructor
     */
    ~ModelHmmGm();
    
    /**
     @return the number of dimensions
     */
    virtual int getNDim() { return ndim; }
    
    /**
     * @return HMM model name
     */
    virtual string getName() { return "GRM"; }
    
    /**
     * @return HMM model full name
     */
    virtual string getFullName() { return "General Reversible Model"; }

    // initialize transitLog array
    virtual void initialize_transitLog();

    /**
     Optimize the model parameters
     @param gradient_epsilon: epsilon for optimization
     @return log-likelihood value
     */
    virtual double optimizeParameters(double gradient_epsilon);

    /**
     Show parameters
     */
    virtual void showParameters(ostream& out);

    virtual int getNParameters();

protected:
    
    int ndim; // number of dimensions
    
    // read the parameters and write into "variables"
    virtual void setVariables(double *variables);
    
    // read the "variables" and write into the parameters
    virtual void getVariables(double *variables);
    
    // set the bounds
    virtual void setBounds(double *lower_bound, double *upper_bound, bool* bound_check);
    
private:
    
    // Parameter: transition matrix
    // dimension: ncat x ncat
    // transit[i*ncat+j] : transition probability for a site with cat i to the next site with cat j
    double* transit;
    double* transit_normalize;
    // sum_j(transit_normalize[i*ncat+j]) = 1.0

    // compute the normalized values of transition matrix
    virtual void computeNormalizedTransits();

    // compute the log values of transition matrix
    virtual void computeLogTransits();

    /**
     Optimize the transition matrix by EM algorithm
     @return log-likelihood value
     */
    virtual double optimizeParametersByEM();

};


#endif
