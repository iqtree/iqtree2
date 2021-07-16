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
    virtual void simulateSeqs(int sequence_length, double *site_specific_rates, ModelSubst *model, double *trans_matrix, int max_num_states, Node *node, Node *dad, ostream &out, vector<string> state_mapping);
    
    /**
        simulate a sequence for a node from a specific branch
    */
    virtual void simulateASequenceFromBranch(ModelSubst *model, int sequence_length, double *trans_matrix, int max_num_states, Node *node, NeighborVec::iterator it);
    
    /**
        simulate a sequence for a node from a specific branch after all variables has been initializing
    */
    void simulateASequenceFromBranchAfterInitVariables(ModelSubst *model, int sequence_length, double *site_specific_rates, double *trans_matrix, int max_num_states, Node *node, NeighborVec::iterator it);
    
    /**
        initialize variables (e.g., site-specific rate)
    */
    void initVariables(int sequence_length, double *site_specific_rates);
    
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
    virtual void simulateSeqsForTree(string output_filepath);
};

#endif /* alisimulatorinvar_h */
