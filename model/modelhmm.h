//
//  modelhmm.h
//  model
//
//  Created by Thomas Wong on 02/02/23.
//

#ifndef modelhmm_h
#define modelhmm_h

#include <string>
#include "utils/tools.h"
#include "utils/optimization.h"
#include "tree/phylotree.h"
#include "tree/phylohmm.h"

#define MIN_TRAN_PROB 1e-10
#define INITIAL_PROB_SAME_CAT 0.9999

using namespace std;

class PhyloHmm;

/**
HMM transition Simple model
*/
class ModelHmm: public Optimization
{
    
public:
    /**
     constructor
     @param ncat : number of categories
     */
    ModelHmm(int numcat);
    
    /**
     destructor
     */
    ~ModelHmm();
    
    /**
     @return the number of dimensions
     */
    virtual int getNDim() { return 1; }
    
    /**
     * @return HMM model name
     */
    virtual string getName() { return "SM"; }

    /**
     * @return HMM model full name
     */
    virtual string getFullName() { return "Simple Model"; }

    // initialize transitLog array
    virtual void initialize_transitLog();
    /**
        set the associated PhyloHmm
        @param phyloHmm the associated PhyloHmm
    */
    void setPhyloHmm(PhyloHmm *phyloHmm);

    /**
     Optimize the model parameters
     @param gradient_epsilon: epsilon for optimization
     @return log-likelihood value
     */
    virtual double optimizeParameters(double gradient_epsilon);
    
    /**
     * @return log values of transition matrix
     */
    virtual double* getTransitLog(int site_i) { return transitLog; }
    
    /**
     Show parameters
     */
    virtual void showParameters(ostream& out);
    
    virtual int getNParameters();

    // probability of transition between the same category
    double tranSameCat;

protected:
    
    int ncat; // number of categories
    int sq_ncat; // ncat * ncat

    /**
        PhyloHmm associated
    */
    PhyloHmm *phylo_hmm;

    // compute the log values of transition matrix
    virtual void computeLogTransits();
    
    // Parameter: log values of transition matrix
    // dimension: ncat x ncat
    // transitLog[i*ncat+j] : log of transition probability for a site with cat i to the next site with cat j
    double* transitLog;

private:
    
    // for optimization
    double computeFunction(double tran_same_cat);

    /**
     Optimize the transition matrix by EM algorithm
     @return log-likelihood value
     */
    virtual double optimizeParametersByEM();
};

#endif
