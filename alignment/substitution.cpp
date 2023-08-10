#include "substitution.h"

using namespace std;

Substitution::Substitution(const std::string& sub_str, Alignment* const aln)
{
    // validate the input
    if (!aln)
        outError("Null alignment found when parsing the predefined mutation: " + sub_str);
    const int num_sites_per_state = aln->seq_type == SEQ_CODON ? 3 : 1;
    const int input_length = sub_str.length();
    if (input_length < num_sites_per_state + num_sites_per_state + 1)
        outError("Failed to parse the predefined mutation: " + sub_str);
    
    // Parse the old state
    old_state = parseState(sub_str.substr(0, num_sites_per_state), aln);
    if (old_state >= aln->num_states)
        outError("Failed to parse the predefined mutation: " + sub_str + ". The old state is invalid.");
    
    // Parse the new state
    new_state = parseState(sub_str.substr(input_length - num_sites_per_state, num_sites_per_state), aln);
    if (new_state >= aln->num_states)
        outError("Failed to parse the predefined mutation: " + sub_str + ". The new state is invalid.");
    
    // Parse the position
    position = convert_int(sub_str.substr(num_sites_per_state, input_length - (num_sites_per_state + num_sites_per_state)).c_str()) - 1;
    if (aln->seq_type == SEQ_CODON)
        position = position * ONE_THIRD;
    if (position < 0)
        outError("Failed to parse the predefined mutation: " + sub_str + ". Position must be positive!");
}

int Substitution::parseState(const std::string& old_state_str, Alignment* const aln) const
{
    // init a state
    int state = 0;
    
    // parse the state regarding their sequence type
    // CODON
    if (aln->seq_type == SEQ_CODON)
    {
        // validate input
        ASSERT(old_state_str.length() == 3);
        
        // dummy variables
        std:string sequence_name;
        ostringstream err_str;
        int num_error = 0;
        
        // parse the codon
        state = aln->getCodonStateTypeFromSites(aln->convertState(old_state_str[0], SEQ_DNA), aln->convertState(old_state_str[1], SEQ_DNA), aln->convertState(old_state_str[2], SEQ_DNA), sequence_name, 0, err_str, num_error);
    }
    // Other data types
    else
        state = aln->convertState(old_state_str[0], aln->seq_type);
    
    // return the state
    return state;
}

short int Substitution::getOldState() const
{
    return old_state;
}

short int Substitution::getNewState() const
{
    return new_state;
}


int Substitution::getPosition() const
{
    return position;
}

Substitutions::Substitutions(const std::string& sub_str, Alignment* const aln)
{
    const int length = sub_str.length();
    // Validate the input
    if (length < 2)
        outError("Failed to parse a list of predefined mutations: '" + sub_str + "'. It should start with { and end with }.");
    if (sub_str[0] != '{')
        outError("List of predefined mutations must start with {");
    if (sub_str[length - 1] != '}')
        outError("List of predefined mutations must end with }");
    
    // Remove {}
    std::string list_mut_str = sub_str.substr(1, length - 2);
    
    // parse mutations one by one
    const std::string delimiter = "/";
    const int max_subs = std::count(list_mut_str.begin(), list_mut_str.end(), delimiter[0]) + 1;
    reserve(max_subs);
    size_t pos = 0;
    while ((pos = list_mut_str.find(delimiter)) != std::string::npos) {
        emplace_back(list_mut_str.substr(0, pos), aln);
        list_mut_str.erase(0, pos + delimiter.length());
    }
    if (list_mut_str.length())
        emplace_back(list_mut_str, aln);
}
