//
// C++ Implementation: alignment
//
// Description: 
//
//
// Author: BUI Quang Minh, Steffen Klaere, Arndt von Haeseler <minh.bui@univie.ac.at>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "alignment.h"
#include "myreader.h"

char symbols_protein[] = "ARNDCQEGHILKMFPSTWYVX"; // X for unknown AA
char symbols_dna[]     = "ACGT";
char symbols_rna[]     = "ACGU";
char symbols_binary[]  = "01";

Alignment::Alignment()
 : vector<Pattern>()
{
	num_states = 0;
	frac_const_sites = 0.0;
}

string &Alignment::getSeqName(int i) {
	assert(i >= 0 && i < seq_names.size());
	return seq_names[i];
}

int Alignment::getSeqID(string &seq_name) {
	for (int i = 0; i < getNSeq(); i++)
		if (seq_name == getSeqName(i)) return i;
	return -1;
}

void checkSeqName(StrVector &seq_names) {
	ostringstream warn_str;
	for (StrVector::iterator it = seq_names.begin(); it != seq_names.end(); it++) {
		string orig_name = (*it);
		for (string::iterator i = it->begin(); i != it->end(); i++) {
			if (!isalnum(*i) && (*i) != '_' && (*i) != '-' && (*i) != '.') {
				(*i) = '_';
			}
		}
		if (orig_name != (*it)) 
			warn_str << orig_name << " -> " << (*it) << endl;
	}
	if (warn_str.str() != "") {
		string str = "Some sequence names are changed as follows:\n";
		outWarning(str + warn_str.str());
	}
}

Alignment::Alignment(char *filename, char *sequence_type, InputType &intype) : vector<Pattern>() {

	cout << "Reading alignment file " << filename << "..." << endl;
	intype = detectInputFile(filename);

	try {
	
		if (intype == IN_NEXUS) {
			readNexus(filename);
		} else {
			readPhylip(filename, sequence_type);
			//outError("in progress");
			//outError("Alignment format not supported, use NEXUS format.");
		}
	} catch (ios::failure) {
		outError(ERR_READ_INPUT);
	} catch (const char *str) {
		outError(str);
	} catch (string str) {
		outError(str);
	}

	if (getNSeq() < 3) 
		outError("Alignment must have at least 3 sequences");
		
	checkSeqName(seq_names);
	cout << "Alignment contains " << getNSeq() << " sequences with " << getNSite() << " characters" << endl;
	cout << "Number of character states is " << num_states << endl;
	cout << "Number of patterns = " << size() << endl;
	countConstSite();
	cout << "Fraction of constant sites: " << frac_const_sites << endl;

}

int Alignment::readNexus(char *filename) {
	NxsTaxaBlock *taxa_block;
	NxsAssumptionsBlock *assumptions_block;
	NxsDataBlock *data_block = NULL;
	NxsTreesBlock *trees_block = NULL;

	taxa_block = new NxsTaxaBlock();
	assumptions_block = new NxsAssumptionsBlock(taxa_block);
	data_block = new NxsDataBlock(taxa_block, assumptions_block);
	trees_block = new TreesBlock(taxa_block);

	MyReader nexus(filename);

	nexus.Add(taxa_block);
	nexus.Add(assumptions_block);
	nexus.Add(data_block);
	nexus.Add(trees_block);

	MyToken token(nexus.inf);
	nexus.Execute(token);

	if (data_block->GetNTax() == 0) {
		outError("No data is given in the input file");	
		return 0;
	}
	if (verbose_mode >= VB_DEBUG)
		data_block->Report(cout);
	
	extractDataBlock(data_block);

	return 1;
}

void Alignment::extractDataBlock(NxsCharactersBlock *data_block) {
	int nseq = data_block->GetNTax();
	int nsite = data_block->GetNCharTotal();
	char *symbols;
	//num_states = strlen(symbols);
	char char_to_state[NUM_CHAR];
	char state_to_char[NUM_CHAR];

	NxsCharactersBlock::DataTypesEnum data_type = (NxsCharactersBlock::DataTypesEnum)data_block->GetDataType();
	if (data_type == NxsCharactersBlock::continuous) {
		outError("Continuous characters not supported");
	} else if (data_type == NxsCharactersBlock::dna || data_type == NxsCharactersBlock::rna || 
		data_type == NxsCharactersBlock::nucleotide) 
	{
		num_states = 4;
		if (data_type == NxsCharactersBlock::rna) 
			symbols = symbols_rna;
		else
			symbols = symbols_dna;
	} else if (data_type == NxsCharactersBlock::protein) {
		num_states = 20;
		symbols = symbols_protein;
	} else {
		num_states = 2;
		symbols = symbols_binary;
	}

	memset(char_to_state, STATE_UNKNOWN, NUM_CHAR);
	memset(state_to_char, '?', NUM_CHAR);
	for (int i = 0; i < strlen(symbols); i++) {
		char_to_state[(int)symbols[i]] = i;	
		state_to_char[i] = symbols[i];
	}
	state_to_char[(int)STATE_UNKNOWN] = '-';


	int seq, site;

	for (seq = 0; seq < nseq; seq++) {
		seq_names.push_back(data_block->GetTaxonLabel(seq));
	}

	site_pattern.resize(nsite, -1);

	int num_gaps_only = 0;

	for (site = 0; site < nsite; site++) {
 		Pattern pat;
		for (seq = 0; seq < nseq; seq++) {
			int nstate = data_block->GetNumStates(seq, site);
			if (nstate == 0) 
				pat += STATE_UNKNOWN;
			else if (nstate == 1) {
				pat += char_to_state[(int)data_block->GetState(seq, site, 0)];
			} else {
				assert(data_type != NxsCharactersBlock::dna || data_type != NxsCharactersBlock::rna || data_type != NxsCharactersBlock::nucleotide);
				char pat_ch = 0;
				for (int state = 0; state < nstate; state++) {
					pat_ch |= (1 << char_to_state[(int)data_block->GetState(seq, site, state)]);
				}
				pat_ch += 3;
				pat += pat_ch;
			} 
		}
		num_gaps_only += addPattern(pat, site);
	}
	if (num_gaps_only)
		cout << "WARNING: " << num_gaps_only << " sites contain only gaps or unknown chars." << endl;
	if (verbose_mode >= VB_MAX)
		for (site = 0; site < size(); site++) {
			for (seq = 0; seq < nseq; seq++)
				cout << state_to_char[(int)(*this)[site][seq]];
 			cout << "  " << (*this)[site].frequency << endl;
		}
}

bool Alignment::addPattern(Pattern &pat, int site, int freq) {
	// check if pattern contains only gaps
	bool gaps_only = true;
	for (Pattern::iterator it = pat.begin(); it != pat.end(); it++)
		if ((*it) != STATE_UNKNOWN) { 
			gaps_only = false; 
			break;
		}
	if (gaps_only) {
		if (verbose_mode >= VB_DEBUG)
			cout << "Site " << site << " contains only gaps or unknown characters" << endl;
		//return true;
	}
	PatternIntMap::iterator pat_it = pattern_index.find(pat);
	if (pat_it == pattern_index.end()) { // not found
		pat.frequency = freq;
		pat.computeConst();
		push_back(pat);
		pattern_index[pat] = size()-1;
		site_pattern[site] = size()-1;
	} else {
		int index = pat_it->second;
		at(index).frequency += freq;
		site_pattern[site] = index;
	} 
	return gaps_only;
}

enum SeqType {SEQ_DNA, SEQ_PROTEIN, SEQ_BINARY, SEQ_UNKNOWN};

/**
	detect the data type of the input sequences
	@param sequences vector of strings
	@return the data type of the input sequences
*/
SeqType detectSequenceType(StrVector &sequences) {
	int num_nuc = 0;
	int num_ungap = 0;
	int num_bin = 0;
	int num_alphabet = 0;

	for (StrVector::iterator it = sequences.begin(); it != sequences.end(); it++)
		for (string::iterator i = it->begin(); i != it->end(); i++) {
			if ((*i) != '?' && (*i) != '-') num_ungap++;
			if ((*i) == 'A' || (*i) == 'C' || (*i) == 'G' || (*i) == 'T' || (*i) == 'U')
				num_nuc++;
			if ((*i) == '0' || (*i) == '1')
				num_bin++;
			if (isalpha(*i)) num_alphabet++;
		}
	if (((double)num_nuc) / num_ungap > 0.9)
		return SEQ_DNA;
	if (((double)num_bin) / num_ungap > 0.9)
		return SEQ_BINARY;
	if (((double)num_alphabet) / num_ungap < 0.5)
		return SEQ_UNKNOWN;
	return SEQ_PROTEIN;
}

/**
	convert a raw characer state into ID, indexed from 0
	@param state input raw state
	@param seq_type data type (SEQ_DNA, etc.)
	@return state ID
*/
char convertState(char state, SeqType seq_type) {
	if (state == '?' || state == '-')
		return STATE_UNKNOWN;

	char *loc;

	switch (seq_type) {
	case SEQ_BINARY:
		switch (state) {
		case '0': return 0;
		case '1': return 1;
		default: return STATE_INVALID;
		}
	case SEQ_DNA: // DNA
		switch (state) {
			case 'A': return 0;
			case 'C': return 1;
			case 'G': return 2;
			case 'T': return 3;
			case 'U': return 3;
			case 'R': return 1+4+3; // A or G, Purine
			case 'Y': return 2+8+3; // C or T, Pyrimidine
			case 'N': return STATE_UNKNOWN;
			case 'W': return 1+8+3; // A or T, Weak
			case 'S': return 2+4+3; // G or C, Strong
			case 'M': return 1+2+3; // A or C, Amino
			case 'K': return 4+8+3; // G or T, Keto
			case 'B': return 2+4+8+3; // C or G or T
			case 'H': return 1+2+8+3; // A or C or T
			case 'D': return 1+4+8+3; // A or G or T
			case 'V': return 1+2+4+3; // A or G or C
			default: return STATE_INVALID; // unrecognize character
		}
		return state;
	case SEQ_PROTEIN: // Protein
		loc = strchr(symbols_protein, state);
		
		if (!loc) return STATE_INVALID; // unrecognize character
		state = loc - symbols_protein;
		if (state < 20) 
			return state;
		else 
			return STATE_UNKNOWN;
	default:
		return STATE_INVALID;
	}
}

char Alignment::convertStateBack(char state) {
	if (state == STATE_UNKNOWN) return '-';
	if (state == STATE_INVALID) return '?';

	switch (num_states) {
	case 2:
		switch (state) {
		case 0: return '0';
		case 1: return '1';
		default: return STATE_INVALID;
		}
	case 4: // DNA
		switch (state) {
			case 0: return 'A';
			case 1: return 'C';
			case 2: return 'G';
			case 3: return 'T';
			case 1+4+3: return 'R'; // A or G, Purine
			case 2+8+3: return 'Y'; // C or T, Pyrimidine
			case 1+8+3: return 'W'; // A or T, Weak
			case 2+4+3: return 'S'; // G or C, Strong
			case 1+2+3: return 'M'; // A or C, Amino
			case 4+8+3: return 'K'; // G or T, Keto
			case 2+4+8+3: return 'B'; // C or G or T
			case 1+2+8+3: return 'H'; // A or C or T
			case 1+4+8+3: return 'D'; // A or G or T
			case 1+2+4+3: return 'V'; // A or G or C
			default: return '?'; // unrecognize character
		}
		return state;
	case 20: // Protein
		if (state < 20) 
			return symbols_protein[state];
		else 
			return '-';
	default:
		return '?';
	}
}

int Alignment::readPhylip(char *filename, char *sequence_type) {
	
	StrVector sequences;
	ostringstream err_str;
	ifstream in;
	int line_num = 1;
	// set the failbit and badbit
	in.exceptions(ios::failbit | ios::badbit);
	in.open(filename);
	int nseq = 0, nsite = 0;
	int seq_id = 0;
	string line;
	// remove the failbit
	in.exceptions(ios::badbit);
	
	for (; !in.eof(); line_num++) {
		getline(in, line);
		if (line == "") continue;

		//cout << line << endl;		
		if (nseq == 0) { // read number of sequences and sites
			istringstream line_in(line);
			if (!(line_in >> nseq >> nsite))
				throw "Invalid PHYLIP format. First line must contain number of sequences and sites";
			//cout << "nseq: " << nseq << "  nsite: " << nsite << endl;
			if (nseq < 3)
				throw "There must be at least 3 sequences";
			if (nsite < 1)
				throw "No alignment columns";

			seq_names.resize(nseq, "");
			site_pattern.resize(nsite, -1);
			sequences.resize(nseq, "");
			clear();
			pattern_index.clear();

		} else { // read sequence contents
			if (seq_names[seq_id] == "") { // cut out the sequence name
				string::size_type pos = line.find(' ');
				if (pos == string::npos) pos = 10; //  assume standard phylip
				seq_names[seq_id] = line.substr(0, pos);
				line.erase(0, pos);
			}
			for (string::iterator it = line.begin(); it != line.end(); it++) {
				if ((*it) <= ' ') continue;
				if (isalpha(*it) || (*it) == '-' || (*it) == '?')
					sequences[seq_id].append(1, toupper(*it));
				else {
					err_str << "Unrecognized character " << *it << " on line " << line_num;
					throw err_str.str();
				}
			}
			seq_id++;
			if (seq_id == nseq) seq_id = 0;
		} 
		//sequences.	
	}
	in.clear();
	// set the failbit again
	in.exceptions(ios::failbit | ios::badbit);
	in.close();

	/* now check that all sequence names are correct */
	for (seq_id = 0; seq_id < nseq; seq_id ++) {
		ostringstream err_str;
		if (seq_names[seq_id] == "") 
			err_str << "Sequence number " << seq_id+1 << " has no names\n";
		// check that all the names are different
		for (int i = 0; i < seq_id; i++)
			if (seq_names[i] == seq_names[seq_id])
				err_str << "The sequence name " << seq_names[seq_id] << " is dupplicated\n";
	}
	if (err_str.str() != "")
		throw err_str.str();


	/* now check that all sequences have the same length */
	for (seq_id = 0; seq_id < nseq; seq_id ++) {
		if (sequences[seq_id].length() != nsite) {
			err_str << "Sequence " << seq_names[seq_id] << " contains ";
			if (sequences[seq_id].length() < nsite)
				err_str << "not enough";
			else
				err_str << "too many";
				
			err_str << " characters (" << sequences[seq_id].length() << ")\n";
		}
	}

	if (err_str.str() != "")
		throw err_str.str();

	/* now check data type */		
	SeqType seq_type = SEQ_UNKNOWN;
	seq_type = detectSequenceType(sequences);
	switch (seq_type) {
	case SEQ_BINARY: 
		num_states = 2;
		cout << "Alignment most likely contains binary sequences" << endl;
		break;
	case SEQ_DNA: 
		num_states = 4;
		cout << "Alignment most likely contains DNA/RNA sequences" << endl;
		break;
	case SEQ_PROTEIN:
		num_states = 20;
		cout << "Alignment most likely contains protein sequences" << endl;
		break;
	default: 
		if (!sequence_type)
			throw "Unknown sequence type.";
	}
	SeqType user_seq_type;
	if (sequence_type) {
		if (strcmp(sequence_type, "B") == 0) {
			num_states = 2;
			user_seq_type = SEQ_BINARY;
		} else if (strcmp(sequence_type, "D") == 0) {
			num_states = 4;
			user_seq_type = SEQ_DNA;
		} else if (strcmp(sequence_type, "P") == 0) {
			num_states = 20;
			user_seq_type = SEQ_PROTEIN;
		} else
			throw "Invalid sequence type.";
		if (user_seq_type != seq_type && seq_type != SEQ_UNKNOWN)
			outWarning("Your specified sequence type is different from the detected one");
		seq_type = user_seq_type;
	}

	// now convert to patterns
	int site, seq, num_gaps_only = 0;

	for (site = 0; site < nsite; site++) {
 		Pattern pat;
		for (seq = 0; seq < nseq; seq++) {
			char state = convertState(sequences[seq][site], seq_type);
			if (state == STATE_INVALID) 
				err_str << "Sequence " << seq_names[seq] << " has invalid character " << 
					sequences[seq][site] << " at site " << site+1 << "\n";
			pat.push_back(state);
		}
		num_gaps_only += addPattern(pat, site);
	}
	if (num_gaps_only)
		cout << "WARNING: " << num_gaps_only << " sites contain only gaps or unknown chars." << endl;
	if (err_str.str() != "")
		throw err_str.str();

}


void Alignment::printPhylip(char *file_name) {
	try {
		ofstream out;
		out.exceptions(ios::failbit | ios::badbit);
		out.open(file_name);
		out << seq_names.size() << " " << site_pattern.size() << endl;
		StrVector::iterator it;
		int max_seq_len = 10;
		for (it = seq_names.begin(); it != seq_names.end(); it++)
			if ((*it).length() > max_seq_len) max_seq_len = (*it).length();
		int seq_id = 0;
		for (it = seq_names.begin(); it != seq_names.end(); it++, seq_id++) {
			out.width(max_seq_len);
			out << left << (*it) << " ";
			for (IntVector::iterator i = site_pattern.begin();  i != site_pattern.end(); i++)
				out << convertStateBack(at(*i)[seq_id]);
			out << endl;
		}
		out.close();
		cout << "Alignment was printed to " << file_name << endl;
	} catch (ios::failure) {
		outError(ERR_WRITE_OUTPUT, file_name);
	}	
}


void Alignment::extractSubAlignment(Alignment *aln, IntVector &seq_id) {
	IntVector::iterator it;
	for (it = seq_id.begin(); it != seq_id.end(); it++) {
		assert(*it >= 0 && *it < aln->getNSeq());
		seq_names.push_back(aln->getSeqName(*it));
	}
	num_states = aln->num_states;
	site_pattern.resize(aln->getNPattern(), -1);
	clear();
	pattern_index.clear();
	int site = 0;
	VerboseMode save_mode = verbose_mode; 
	verbose_mode = VB_MIN; // to avoid printing gappy sites in addPattern
	for (iterator pit = aln->begin(); pit != aln->end(); pit++, site++) {
		Pattern pat;
		for (it = seq_id.begin(); it != seq_id.end(); it++)
			pat.push_back((*pit)[*it]);
		addPattern(pat, site, (*pit).frequency);
	}
	verbose_mode = save_mode;
	countConstSite();
}

void Alignment::countConstSite() {
	int num_const_sites = 0;
	for (iterator it = begin(); it != end(); it++)
		if ((*it).is_const) num_const_sites += (*it).frequency;
	frac_const_sites = ((double)num_const_sites) / getNSite();
}

Alignment::~Alignment()
{
}

double Alignment::computeObsDist(int seq1, int seq2) {
	int diff_pos = 0, total_pos = 0;
	for (iterator it = begin(); it != end(); it++) 
		if  ((*it)[seq1] < num_states && (*it)[seq2] < num_states) {
		//if ((*it)[seq1] != STATE_UNKNOWN && (*it)[seq2] != STATE_UNKNOWN) {
			total_pos += (*it).frequency;
			if ((*it)[seq1] != (*it)[seq2] )
				diff_pos += (*it).frequency;
		}
	if (!total_pos) total_pos = 1;
	return ((double)diff_pos) / total_pos;
}

double Alignment::computeJCDist(int seq1, int seq2) {
	double obs_dist = computeObsDist(seq1, seq2);
	double x = 1.0 - ((double)num_states * obs_dist / (num_states-1));
	if (x <= 0) {
		string str = "Too long distance between two sequences ";
		str += getSeqName(seq1);
		str += " and ";
		str += getSeqName(seq2);
		outWarning(str);
		return MAX_GENETIC_DIST;
	}
	return -log(x) * (num_states - 1) / num_states;
}

void Alignment::printDist(ostream &out, double *dist_mat) {
	int nseqs = getNSeq();
	out << nseqs << endl;
	int pos = 0;
	out.precision(6);
	out << fixed;
	for (int seq1 = 0; seq1 < nseqs; seq1 ++)  {
		out.width(10);
		out << left << getSeqName(seq1) << " ";
		out.width(8);
		for (int seq2 = 0; seq2 < nseqs; seq2 ++) {
			out << dist_mat[pos++] << " ";
		}	
		out << endl;
	}
}

void Alignment::printDist(const char *file_name, double *dist_mat) {
	try {
		ofstream out;
		out.exceptions(ios::failbit | ios::badbit);
		out.open(file_name);
		printDist(out, dist_mat);
		out.close();
		cout << "Distance matrix was printed to " << file_name << endl;
	} catch (ios::failure) {
		outError(ERR_WRITE_OUTPUT, file_name);
	}	
}

void Alignment::readDist(istream &in, double *dist_mat) {
	int nseqs;
	in >> nseqs;
	if (nseqs != getNSeq())
		throw "Distance file has different number of taxa";
	int pos = 0, seq1, seq2;
	for (seq1 = 0; seq1 < nseqs; seq1 ++)  {
		string seq_name;
		in >> seq_name;
		if (seq_name != getSeqName(seq1))
			throw "Sequence name " + seq_name + " is different from " + getSeqName(seq1);
		for (seq2 = 0; seq2 < nseqs; seq2 ++) {
			in >> dist_mat[pos++];
		}	
	}
	// check for symmetric matrix
	for (seq1 = 0; seq1 < nseqs-1; seq1++) {
		if (dist_mat[seq1*nseqs+seq1] != 0.0)
			throw "Diagonal elements of distance matrix is not ZERO";
		for (seq2 = seq1+1; seq2 < nseqs; seq2++)
			if (dist_mat[seq1*nseqs+seq2] != dist_mat[seq2*nseqs+seq1])
				throw "Distance between " + getSeqName(seq1) + " and " + getSeqName(seq2) + " is not symmetric";
	}
}

void Alignment::readDist(const char *file_name, double *dist_mat) {
	try {
		ifstream in;
		in.exceptions(ios::failbit | ios::badbit);
		in.open(file_name);
		readDist(in, dist_mat);
		in.close();
		cout << "Distance matrix was read from " << file_name << endl;
	} catch (const char *str) {
		outError(str);
	} catch (string str) {
		outError(str);
	} catch (ios::failure) {
		outError(ERR_READ_INPUT, file_name);
	}
	
}


void Alignment::computeStateFreq (double *stateFrqArr) {
	int stateNo_;
	int nState_ = num_states;
	int nseqs = getNSeq();
	double timeAppArr_[num_states];
	double siteAppArr_[num_states]; //App = appearance
	double newSiteAppArr_[num_states];

	for (stateNo_ = 0; stateNo_ < nState_; stateNo_ ++)
		stateFrqArr [ stateNo_ ] = 1.0 / nState_;

	int NUM_TIME = 8;
	//app = appeareance
	for (int time_ = 0; time_ < NUM_TIME; time_ ++) 
	{
		for (stateNo_ = 0; stateNo_ < nState_; stateNo_ ++)
			timeAppArr_[stateNo_] = 0.0;

		for (iterator it = begin(); it != end(); it++) 
			for (int i = 0; i < (*it).frequency; i++)	
			{
			for (int seq = 0; seq < nseqs; seq++) {
				int stateNo_ = (*it)[seq];

				getAppearance (stateNo_, siteAppArr_);

				double totalSiteApp_ = 0.0;
				for (stateNo_ = 0; stateNo_ < nState_; stateNo_ ++) {
					newSiteAppArr_[stateNo_] = stateFrqArr[stateNo_] * siteAppArr_[stateNo_];
					totalSiteApp_ += newSiteAppArr_[stateNo_];
				}

				for (stateNo_ = 0; stateNo_ < nState_; stateNo_ ++)
					timeAppArr_[stateNo_] += newSiteAppArr_[stateNo_] / totalSiteApp_;
			}
		}

		double totalTimeApp_ = 0.0;
		int stateNo_;
		for (stateNo_ = 0; stateNo_ < nState_; stateNo_ ++)
			totalTimeApp_ += timeAppArr_[stateNo_];


		for (stateNo_ = 0; stateNo_ < nState_; stateNo_ ++)
			stateFrqArr[stateNo_] = timeAppArr_[stateNo_] / totalTimeApp_;

	} //end of for time_

	//  std::cout << "state frequency ..." << endl;
	// for (stateNo_ = 0; stateNo_ < nState_; stateNo_ ++)
	// std::cout << stateFrqArr[stateNo_] << endl;

	if (verbose_mode >= VB_DEBUG) {
		cout << "Empirical state frequencies: ";
		for (stateNo_ = 0; stateNo_ < nState_; stateNo_ ++)
			cout << stateFrqArr[stateNo_] << " ";
		cout << endl;
	}

}

void Alignment::getAppearance(char state, double *state_app) {
	int i;
	if (state == STATE_UNKNOWN) {
		for (i = 0; i < num_states; i++)
			state_app[i] = 1.0;
		return;
	}

	memset(state_app, 0, num_states * sizeof(double));
	if (state < num_states) {
		state_app[state] = 1.0;
		return;
	}
	state -= (num_states-1);
	for (i = 0; i < num_states; i++) 
	if (state & (1 << i)) {
		state_app[i] = 1.0;
	}
}


void Alignment::computeEmpiricalRate (double *rates) {
	int i, j, k;
	assert(rates);
	int nseqs = getNSeq();
	double **pair_rates = (double**) new double[num_states];
	for (i = 0; i < num_states; i++) {
		pair_rates[i] = new double[num_states];
		memset(pair_rates[i], 0, sizeof(double)*num_states);
	}

	for (iterator it = begin(); it != end(); it++) {
		for (i = 0; i < nseqs-1; i++) {
			char state1 = (*it)[i];
			if (state1 >= num_states) continue;
			for (j = i+1; j < nseqs; j++) {
				char state2 = (*it)[j];
				if (state2 < num_states) pair_rates[state1][state2] += (*it).frequency;
			}
		}
	}

	k = 0;
	double last_rate = pair_rates[num_states-2][num_states-1] + pair_rates[num_states-1][num_states-2];
	for (i = 0; i < num_states-1; i++)
		for (j = i+1; j < num_states; j++)
			rates[k++] = (pair_rates[i][j] + pair_rates[j][i]) / last_rate;
	if (verbose_mode >= VB_DEBUG) {
		cout << "Empirical rates: ";
		for (k = 0; k < num_states*(num_states-1)/2; k++)
			cout << rates[k] << " ";
		cout << endl;
	}

	for (i = num_states-1; i >= 0; i--) {
		delete [] pair_rates[i];
	}
	delete [] pair_rates;
}

double Alignment::computeUnconstrainedLogL() {
	int nptn = size();
	double logl = 0.0;
	int nsite = getNSite(), i;
	double lognsite = log(nsite);
	for (i = 0; i < nptn; i++)
		logl += (log(at(i).frequency) - lognsite) * at(i).frequency;
	return logl;
}

