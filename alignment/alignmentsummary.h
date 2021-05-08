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
protected:
    const Alignment*      alignment;
    std::vector<intptr_t> siteNumbers;      //of sites with variation
    IntVector             siteFrequencies;  //ditto
    IntVector             nonConstSiteFrequencies; //ditto, but zeroed if site
                                                //isConst according to alignment
    std::map<int, int>    stateToSumOfConstantSiteFrequencies;
    size_t                totalFrequency;    //sum of frequencies (*including* constant sites!)
    size_t                totalFrequencyOfNonConstSites; //ditto (*excluding* constant sites!)
    StateType             minState;        //found on any site where there is variation
    StateType             maxState;        //ditto
    char*                 sequenceMatrix;
    intptr_t              sequenceLength;  //Sequence length (or: count of sites per sequence)
    intptr_t              sequenceCount;   //The number of sequences

public:
    AlignmentSummary(const Alignment* a, bool keepConstSites, bool keepBoringSites);
    ~AlignmentSummary();
    bool                hasSequenceMatrix() const;
    intptr_t            getSequenceCount() const;
    intptr_t            getSumOfConstantSiteFrequenciesForState(int state) const;
    const    IntVector& getSiteFrequencies() const;
    intptr_t            getTotalFrequency() const;
    const    IntVector& getNonConstSiteFrequencies() const;
    size_t              getTotalFrequencyOfNonConstSites() const;
    
    const char* getSequenceMatrix() const;
    const char* getSequence(int sequence_id) const;
    intptr_t    getSequenceLength() const;
    intptr_t    getStateCount() const;
    bool        constructSequenceMatrix ( bool treatAllAmbiguousStatesAsUnknown,
                                          progress_display_ptr progress = nullptr);
    bool        constructSequenceMatrixNoisily ( bool treatAllAmbiguousStatesAsUnknown,
                                                 const char* taskName, const char* verb);
};

#endif /* alignmentsummary_hpp */
