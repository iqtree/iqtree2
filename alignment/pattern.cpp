//
// C++ Implementation: pattern
//
// Description:
//
//
// Author: BUI Quang Minh, Steffen Klaere, Arndt von Haeseler <minh.bui@univie.ac.at>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "pattern.h"
#include "alignment/alignment.h"
#include <vectorclass/vectorclass.h>

Pattern::Pattern()
        : vector<StateType>()
{
    frequency = 0;
//    is_const = false;
//    is_informative = false;
    flag = 0;
    const_char = -1;
    num_chars = 0;
}

Pattern::Pattern(int nseq, int freq)
: vector<StateType>(nseq)
{
    frequency = freq;
    //    is_const = false;
    //    is_informative = false;
    flag = 0;
    const_char = -1;
    num_chars = 0;
}

Pattern::Pattern(const Pattern &pat)
        : vector<StateType>(pat)
{
    frequency = pat.frequency;
//    is_const = pat.is_const;
//    is_informative = pat.is_informative;
    flag = pat.flag;
    const_char = pat.const_char;
    num_chars = pat.num_chars;
}

Pattern::~Pattern()
{
}

int Pattern::computeAmbiguousChar(int num_states) const {
    int num = 0;
    for (auto it = begin(); it != end(); ++it) {
        if (static_cast<int>(*it) >= num_states) {
            num++;
        }
    }
    return num;
}

#define VECTORIZE_GAPCHAR_COUNT 1
int Pattern::computeGapChar(int num_states, int STATE_UNKNOWN) const {
    int num = 0;
#if VECTORIZE_GAPCHAR_COUNT
    //This won't compile unless value_type is based on uint32_t
    //(nor should it! You'd need to use different vector types!)
    const uint32_t* dataStart = data();
    size_type  count   = size();
    size_type  vecSize = Vec8ui::size();
    Vec8ui     unknown = STATE_UNKNOWN;
    const uint32_t* dataStop   = dataStart + count;
    const uint32_t* blockStop  = dataStop - (count & (vecSize-1));
    for (const uint32_t* block=dataStart; block<blockStop; block+=vecSize) {
        Vec8ui a;
        a.load(block);
        num -= horizontal_add( Vec8ui(a == unknown) );
    }
    for (const uint32_t* single=blockStop; single<dataStop; ++single) {
        if (*single == STATE_UNKNOWN) {
            ++num;
        }
    }
#else
    for (iterator i = begin(); i != end(); i++)
        if (*i == STATE_UNKNOWN) num++;
#endif
    return num;
}

bool Pattern::isAllGaps(int STATE_UNKNOWN) const {
    //Probably not worth vectorizing this function, since it will usually
    //return early (it can return false as soon as a known state  is
    //seen, which will probably happen very early in the sequence).
    for (auto it = begin(); it != end(); ++it) {
        if (*it != STATE_UNKNOWN) {
            return false;
        }
    }
    return true;
}

void Pattern::countAppearances(Alignment* aln) {
    is_const = true;
    is_invariant = false;
    is_informative = false;
    // critical fix: const_char was set wrongly to num_states
    // in some data type (binary, codon), causing wrong
    // log-likelihood computation for +I or +I+G model

    const_char = aln->STATE_UNKNOWN+1;
//    if (STATE_UNKNOWN == num_states)
//    	pat.const_char = STATE_UNKNOWN+1;
//    else
//    	pat.const_char = STATE_UNKNOWN;

    state_app.reset();
    for (int j = 0; j < aln->num_states; j++) {
        state_app[j] = 1;
    }

    num_app.resize (aln->num_states, 0);
    last_app.resize(aln->num_states, 0);

    auto pat_data = data();
    int  pat_len  = static_cast<int>(size());
    for (int i = 0; i < pat_len; ++i) {
        Pattern::value_type j = pat_data[i];
    	StateBitset this_app;
    	aln->getAppearance(j, this_app);
    	state_app &= this_app;
        if (static_cast<int>(j) < aln->num_states) {
            auto state = (int)j;
            num_app [state]++;
            last_app[state] = i;
        }
//        else if (*i != STATE_UNKNOWN) {
//            // ambiguous characters
//            is_const = false;
//        }
    }
    count_multi        = 0; // number of states with >= 2 appearances
    count_singleton    = 0; // number of singleton states (with 1 appearance only)
    num_chars          = 0; // number of states with >= 1 appearance
    for (int j = 0; j < aln->num_states; j++) { 
        if (num_app[j]) {
            num_chars++;
            if (num_app[j] >= 2) {
                ++count_multi;
            }
            else if (num_app[j] == 1) {
                ++count_singleton;                           
            }
        }
    }
}

void Pattern::setConstCharForAlignment(Alignment* aln) {
    if (state_app.count() == aln->num_states) {
        const_char = aln->STATE_UNKNOWN;
    } else if (state_app.count() == 1) {
        for (int j = 0; j < aln->num_states; j++) {
            if (state_app[j]) {
                const_char = j;
                break;
            }
        }
    } else if (aln->seq_type == SeqType::SEQ_DNA) {
        const_char = aln->num_states-1;
        for (int j = 0; j < aln->num_states; j++)
            if (state_app[j])
                const_char += (1<<j);
    } else if (aln->seq_type == SeqType::SEQ_PROTEIN) {
        if (state_app[2] && state_app[3]) //4+8, // B = N or D
            const_char = aln->num_states;
        else if (state_app[5] && state_app[6]) //32+64, // Z = Q or E
            const_char = aln->num_states+1;
        else if (state_app[9] && state_app[10]) // 512+1024 // U = I or L
            const_char = aln->num_states+2;
        else ASSERT(0);
    } else {
        ASSERT(0);
    }    
}

void Pattern::setInformativeFlags(Alignment* aln) {
    is_informative = (count_multi > 1);
    is_const       = (state_app.count() >= 1);
    if (is_const) {
        setConstCharForAlignment(aln);
    }

    // compute is_invariant
    is_invariant = (state_app.count() >= 1);
    ASSERT(is_invariant >= is_const);

    // Wed Jun 28 16:01:30 BST 2017. The calculation of these properties seems
    // to be OKish. They are only used for reports and to calculate the
    // parsimony tree in the beginning anyways.

    // if (seq_type == SeqType::SEQ_POMO) {
    //     // For PoMo most sites are informative (ambiguous map from data to state space)
    //     is_informative = true;
    //     // For PoMo there are hardly any constant sites
    //     is_const = false;
    //     is_invariant = false;
    // }

    flag  = (is_const       ? PAT_CONST       : 0);
    flag |= (is_invariant   ? PAT_INVARIANT   : 0);
    flag |= (is_informative ? PAT_INFORMATIVE : 0);
}

void Pattern::countTowardSingletonParsimonyStates(std::vector<UINT>& singleton_parsimony_states) const {
    int num_states = static_cast<int>(num_app.size());
    if (!is_informative && 0<count_singleton) {
        for (int j = 0; j < num_states; j++) {
            if (num_app[j]==1) {
                size_t last_taxon_id = last_app[j];
                singleton_parsimony_states[last_taxon_id] += frequency;
            }
        }
    }

}