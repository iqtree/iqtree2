//
//  alignmentsummary.cpp
//  alignment
//
//  Created by James Barbetti on 1/7/20.
//

#include "alignment.h"
#include "alignmentsummary.h"

AlignmentSummary::AlignmentSummary(const Alignment* a
                                   , bool keepConstSites
                                   , bool keepBoringSites) {
    alignment      = a;
    sequenceMatrix = nullptr;
    sequenceCount  = a->getNSeq();
    totalFrequency = 0;
    totalFrequencyOfNonConstSites = 0;
    if (sequenceCount==0) {
        minState = a->STATE_UNKNOWN;
        maxState = a->STATE_UNKNOWN;
        return;
    }
    
    struct SiteSummary
    {
    public:
        bool      isConst;
        int       frequency;
        StateType minState;
        StateType maxState;
        SiteSummary(): isConst(false), frequency(0), minState(0), maxState(0) {}
    };
    
    size_t siteCount = alignment->size();
    std::vector<SiteSummary> sites;
    sites.resize(siteCount);
    #ifdef _OPENMP
        #pragma omp parallel for
    #endif
    for (size_t site=0; site<siteCount; ++site) { //per site
        auto        itSite = alignment->begin() + site;
        SiteSummary &s     = sites[site];
        StateType minStateForSite = (*itSite)[0];
        StateType maxStateForSite = minStateForSite;
        s.isConst   = itSite->isConst();
        s.frequency = itSite->frequency;
        for (size_t seq=1; seq<sequenceCount; ++seq) {
            auto state = (*itSite)[seq] ;
            if (state < minStateForSite ) {
                minStateForSite = state;
            }
            else if (maxStateForSite < state ) {
                maxStateForSite = state;
            }
        }
        s.minState  = minStateForSite;
        s.maxState  = maxStateForSite;
    }
    sequenceLength = 0; //Number sites where there's some variability
    std::map<int, int> & map = stateToSumOfConstantSiteFrequencies;
    for (size_t site=0; site<siteCount; ++site) {
        SiteSummary &s      = sites[site];
        bool alreadyCounted = false;
        totalFrequency     += s.frequency;
        totalFrequencyOfNonConstSites += s.isConst ? 0 : s.frequency;
        if ( keepConstSites || !s.isConst) {
            if ( keepBoringSites || ( 0 < s.frequency && s.minState < s.maxState) ) {
                alreadyCounted = true;
                if (sequenceLength==0) {
                    minState = s.minState;
                    maxState = s.maxState;
                } else {
                    if (s.minState < minState) {
                        minState = s.minState;
                    }
                    if (maxState < s.maxState) {
                        maxState = s.maxState;
                    }
                }
                siteNumbers.emplace_back(site);
                siteFrequencies.emplace_back(s.frequency);
                nonConstSiteFrequencies.emplace_back(s.isConst ? 0 : s.frequency);
                ++sequenceLength;
            }
        }
        if (s.minState == s.maxState && !alreadyCounted) {
            if (map.find(s.minState)==map.end()) {
                map[s.minState] = s.frequency;
            } else {
                map[s.minState] += s.frequency;
            }
        }
    }
}

AlignmentSummary::~AlignmentSummary() {
    delete [] sequenceMatrix;
    sequenceMatrix = nullptr;
    sequenceLength = 0;
    sequenceCount  = 0;
}

size_t AlignmentSummary::getSumOfConstantSiteFrequenciesForState(int state) {
    auto it = stateToSumOfConstantSiteFrequencies.find(state);
    return (it == stateToSumOfConstantSiteFrequencies.end()) ? 0 : it->second;
}

bool AlignmentSummary::constructSequenceMatrixNoisily(bool treatAllAmbiguousStatesAsUnknown
    , const char* taskName, const char* verb) {
    progress_display progress(sequenceCount, taskName, verb, "sequence");
    return constructSequenceMatrix(treatAllAmbiguousStatesAsUnknown, &progress);
}

bool AlignmentSummary::constructSequenceMatrix ( bool treatAllAmbiguousStatesAsUnknown
                                                , progress_display* progress) {
    delete [] sequenceMatrix;
    sequenceMatrix = nullptr;
    if ( minState<0 || 127 < maxState  ) {
        return false;
    }
    sequenceMatrix =  new char[ sequenceCount * sequenceLength ];
    const int* posToSite = siteNumbers.data();
    if (treatAllAmbiguousStatesAsUnknown)
    {
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (size_t seq=0; seq<sequenceCount; ++seq) { //sequence
            char *sequence = sequenceMatrix + seq * sequenceLength;
            for (size_t seqPos = 0; seqPos < sequenceLength; ++seqPos) { //position in reduced sequence
                auto state = alignment->at(posToSite[seqPos])[seq] ;
                if ( this->alignment->num_states <= state ) {
                    state = this->alignment->STATE_UNKNOWN;
                }
                sequence[seqPos] = static_cast<char> ( state );
                //the state at the (seqPos)th non-constant site, in the (seq)th sequence
            }
            if (progress!=nullptr && (seq % 100) == 0) {
                (*progress) += 100;
            }
        }
    }
    else
    {
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (size_t seq=0; seq<sequenceCount; ++seq) { //sequence
            char *sequence = sequenceMatrix + seq * sequenceLength;
            for (size_t seqPos = 0; seqPos < sequenceLength; ++seqPos) { //position in reduced sequence
                sequence[seqPos] = static_cast<char> ( alignment->at(posToSite[seqPos])[seq] );
                //the state at the (seqPos)th non-constant site, in the (seq)th sequence
            }
            if (progress != nullptr && (seq % 100) == 0) {
                (*progress) += 100;
            }
        }
    }
    return true;
}
