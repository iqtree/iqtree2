/*
 *  alisim.h
 *  implemetation of AliSim (Alignment Simulator)
 *  Created on: Mar 13, 2021
 *      Author: Nhan Ly-Trong
 */

#include "alisim.h"

void runAliSim(Params params)
{
    cout << "[Alignment Simulator] Executing" <<"\n";
    
    // case 1 (default): without rate heterogeneity
    AliSimulator *alisimulator = new AliSimulator(&params);
    
    // get variables
    string rate_name = alisimulator->tree->getRateName();
    double invariant_proportion = alisimulator->tree->getRate()->getPInvar();
    
    // case 2: with rate heterogeneity
    if (!rate_name.empty())
    {
        if((rate_name.find("+G") != std::string::npos) || (rate_name.find("+R") != std::string::npos))
        {
            // case 2.1: with rate heterogeneity (gamma/freerate model with invariant sites)
            if (invariant_proportion > 0)
            {
                alisimulator = new AliSimulatorHeterogeneityInvar(alisimulator, invariant_proportion);
            }
            // case 2.2: with rate heterogeneity (gamma/freerate model without invariant sites)
            else
            {
                alisimulator = new AliSimulatorHeterogeneity(alisimulator);
                
            }
        }
        // case 2.3: without gamma/freerate model with only invariant sites
        else if (rate_name.find("+I") != std::string::npos)
        {
            alisimulator = new AliSimulatorInvar(alisimulator, invariant_proportion);
        }
    }
    
    // show parameters
    alisimulator->showParameters();
    
    // iteratively generate multiple/a single  alignment(s) for each tree
    alisimulator->generateMultipleAlignmentsFromSingleTree();
    
    cout << "[Alignment Simulator] Done"<<"\n";
    
    // delete alisimulator
    delete alisimulator;
}
