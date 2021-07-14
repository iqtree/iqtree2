//
//  randomdistribution.h
//  iqtree
//
//  Created by Nhan Ly-Trong on 14/07/2021.
//


#ifndef randomdistribution_h
#define randomdistribution_h

#include "tools.h"

class RandomDistribution
{
private:
    map<string, string> distributions;
    int num_rand_numbers;
    
    /**
            randomly select a base frequency from the pool of freqs
            @param distribution_name storing name of distribution
     */
    double random_a_base_freq(string distribution_name);

public:
    /**
        constructor
    */
    RandomDistribution();
    
    /**
            read distributions from built-in string or user-specified file
     */

    void read_distributions(char* filepath);
    
    /**
            randomly select the frequency for each nucleotide base from the pools of freqs
            @param freqs storing the output base frequencies
            @param list_distribution_names specifying a list of distributions corresponding for each state
            @param num_states the number of states
     */
    void random_base_frequencies(double *freqs, int num_states = 4, string list_distribution_names = "");

    
};

#endif /* randomdistribution_h */
