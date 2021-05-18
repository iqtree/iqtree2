//
//  alisimulatorinvar.h
//  iqtree
//
//  Created by Nhan Ly-Trong on 23/03/2021.
//

#ifndef alisimulatorinvar_h
#define alisimulatorinvar_h

#include "alisimulator.h"

class AliSimulatorInvar : public AliSimulator
{
protected:
    double invariant_proportion;

    /**
    *  simulate sequences for all nodes in the tree by DFS
    *
    */
    virtual void simulateSeqs(int sequence_length, double *site_specific_rates, ModelSubst *model, double *trans_matrix, int max_num_states, Node *node, Node *dad, int thread_id = 0, int num_threads = 1);
    
public:
    
    /**
        constructor
    */
    AliSimulatorInvar(Params *params, double invar_prop);
    
    /**
        constructor
    */
    AliSimulatorInvar(AliSimulator *alisimulator, double invar_prop);

    /**
    *  simulate sequences for all nodes in the tree
    */
    virtual void simulateSeqsForTree();
};

#endif /* alisimulatorinvar_h */
