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
#include <utils/progress.h> //for progress_display

/**
Summary (for an Alignment) of sites where there are variations
        @author James Barbetti
 */

class Alignment;

struct AlignmentSummary
{
public:
    AlignmentSummary(const Alignment* a, bool keepConstSites, bool keepBoringSites);
    ~AlignmentSummary();
    const Alignment*   alignment;
    std::vector<int>   siteNumbers;      //of sites with variation
    std::vector<int>   siteFrequencies;  //ditto
    std::vector<int>   nonConstSiteFrequencies; //ditto, but zeroed if site
                                                //isConst according to alignment
    std::map<int, int> stateToSumOfConstantSiteFrequencies;
    size_t             totalFrequency;    //sum of frequencies (*including* constant sites!)
    size_t             totalFrequencyOfNonConstSites; //ditto (*excluding* constant sites!)
    StateType          minState; //found on any site where there is variation
    StateType          maxState; //ditto
    char*              sequenceMatrix;
    size_t             sequenceLength;  //Sequence length
    size_t             sequenceCount;   //The number of sequences
    size_t             getSumOfConstantSiteFrequenciesForState(int state);
    bool constructSequenceMatrix ( bool treatAllAmbiguousStatesAsUnknown
                                 , progress_display *progress = nullptr);
    bool constructSequenceMatrixNoisily ( bool treatAllAmbiguousStatesAsUnknown, 
        const char* taskName, const char* verb);
};

#endif /* alignmentsummary_hpp */
