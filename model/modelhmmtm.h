//
//  modelhmmtm.h
//  model
//
//  Created by Thomas Wong on 13/10/23.
//

#ifndef modelhmmtm_h
#define modelhmmtm_h

#include "modelhmm.h"

#define MIN_VALUE 1e-5
#define MAX_VALUE 10000

using namespace std;

class PhyloHmm;

// HMM type-dependent model
class ModelHmmTm: public ModelHmm
{
public:
    /**
     constructor
     @param ncat : number of categories
     @param numtype : number of types
     */
    ModelHmmTm(int numcat, int numtype, int* siteTypes, int nsite, string* typeDescs);
    
    /**
     destructor
     */
    ~ModelHmmTm();
    
    /**
     @return the number of dimensions
     */
    virtual int getNDim() { return ntypepair; }
    
    /**
     * @return HMM model name
     */
    virtual string getName() { return "TM"; }
    
    /**
     * @return HMM model full name
     */
    virtual string getFullName() { return "Site-Type dependent Model"; }
    
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
    
    /**
     * @return log values of transition matrix
     */
    virtual double* getTransitLog(int site_i);
    
    // probability of transition between the same category for each typepair
    double* tranSameCats;
    
protected:
    
    int ndim; // number of dimensions
    int ntype; // number of types
    int ntypepair; // number of type pairs
    int nsite; // number of sites
    
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
    
    // which transit matrix should be used for each site
    int* transit_fwd; // site i-1 -> site i
    
    // representation of the IDs of the transition matrix
    vector<string> transit_id_descs;
    
    // update the transition matrix according to the array tranSameCats
    void updateTransits();
    
    // compute the normalized values of transition matrix
    virtual void computeNormalizedTransits();
    
    // compute the log values of transition matrix
    virtual void computeLogTransits();
    
    /**
     Optimize the transition matrix by EM algorithm
     @return log-likelihood value
     */
    virtual double optimizeParametersByEM();
    
    // create the mapping between site ID and which matrix should be used for the transition between the previous site and the current site
    void initialize_transit_id(int* siteTypes, int nsite, string* typeDescs);
};


#endif
