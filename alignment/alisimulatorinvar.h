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
        simulate a sequence for a node from a specific branch after all variables has been initializing
    */
    virtual void simulateASequenceFromBranchAfterInitVariables(ModelSubst *model, int sequence_length, double *trans_matrix, Node *node, NeighborVec::iterator it, string lengths = "");
    
    /**
        initialize variables (e.g., site-specific rate)
    */
    virtual void initVariables(int sequence_length, bool regenerate_root_sequence = false);
    
    /**
    *  insert a new sequence into the current sequence
    *
    */
    virtual void insertNewSequenceForInsertionEvent(vector<short int> &indel_sequence, int position, vector<short int> &new_sequence);
    
    /**
      initialize site_specific_rates
    */
    void initSiteSpecificRates(vector<double> &site_specific_rates, int sequence_length);
    
public:
    
    /**
        constructor
    */
    AliSimulatorInvar(Params *params, double invar_prop);
    
    /**
        constructor
    */
    AliSimulatorInvar(AliSimulator *alisimulator, double invar_prop);
};

#endif /* alisimulatorinvar_h */
