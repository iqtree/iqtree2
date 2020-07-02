//
//  alignmentsummary.hpp
//  alignment
//
//  Created by James Barbetti on 1/7/20.
//

#ifndef alignmentsummary_hpp
#define alignmentsummary_hpp

#include <vector>
#include <map>

/**
Summary (for an Alignment) of sites where there are variations
        @author James Barbetti
 */

class Alignment;

struct AlignmentSummary
{
public:
    AlignmentSummary(const Alignment* a, bool keepConstSites);
    ~AlignmentSummary();
    const Alignment*   alignment;
    std::vector<int>   siteNumbers;      //of sites with variation
    std::vector<int>   siteFrequencies;  //ditto
    std::map<int, int> stateToSumOfConstantSiteFrequencies;
    int                totalFrequency;   //sum of frequencies (*including* constant sites!)
    int                totalFrequencyOfNonConstSites; //ditto (*excluding* constant sites!)
    StateType          minState; //found on any site where there is variation
    StateType          maxState; //ditto
    char*              sequenceMatrix;
    size_t             sequenceLength;  //Sequence length
    size_t             sequenceCount;   //The number of sequences
    int                getSumOfConstantSiteFrequenciesForState(int state);
    bool constructSequenceMatrix(bool treatAllAmbiguousStatesAsUnknown);
};

#endif /* alignmentsummary_hpp */
