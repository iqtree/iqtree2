#ifndef SUBSTITUTION_H
#define SUBSTITUTION_H

#include "alignment.h"

/** class storing a substition (evolving from one state to another state) at a site */
class Substitution {
private:
    /**
     The old state
     */
    short int old_state;
    
    /**
     The new state
     */
    short int new_state;
    
    /**
     Position where the substitution occurs
     */
    int position;
    
    /**
     Parse a state from readable character(s)
     */
    int parseState(const std::string& old_state_str, Alignment* const aln) const;
    
public:
    /**
     Default constructor
     */
    // Substitution();
    
    /**
     Custom constructor from a string
     */
    Substitution(const std::string& sub_str, Alignment* const aln, const int& seq_length);
    
    /**
     Get the old state
     */
    short int getOldState() const;
    
    /**
     Get the new state
     */
    short int getNewState() const;
    
    /**
     Get the position
     */
    int getPosition() const;
};

/** class storing a set of substitions */
class Substitutions: public std::vector<Substitution> {
public:
    /**
     Custom constructor from a string {<Sub>/.../<Sub>}
     */
    Substitutions(const std::string& sub_str, Alignment* const aln, const int& seq_length);
};

#endif
