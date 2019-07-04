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
#include "utils/tools.h"
#include "alignment.h"
#include "nclextra/myreader.h"
#include <numeric>
#include <sstream>
#include "model/rategamma.h"
#include "gsl/mygsl.h"
#include "utils/gzstream.h"

using namespace std;

char symbols_protein[] = "ARNDCQEGHILKMFPSTWYVX"; // X for unknown AA
char symbols_dna[]     = "ACGT";
char symbols_rna[]     = "ACGU";
//char symbols_binary[]  = "01";
char symbols_morph[] = "0123456789ABCDEFGHIJKLMNOPQRSTUV";
// genetic code from tri-nucleotides (AAA, AAC, AAG, AAT, ..., TTT) to amino-acids
// Source: http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
// Base1:                AAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTT
// Base2:                AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT
// Base3:                ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
char genetic_code1[]  = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"; // Standard
char genetic_code2[]  = "KNKNTTTT*S*SMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"; // Vertebrate Mitochondrial
char genetic_code3[]  = "KNKNTTTTRSRSMIMIQHQHPPPPRRRRTTTTEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"; // Yeast Mitochondrial
char genetic_code4[]  = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"; // Mold, Protozoan, etc.
char genetic_code5[]  = "KNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"; // Invertebrate Mitochondrial
char genetic_code6[]  = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSS*CWCLFLF"; // Ciliate, Dasycladacean and Hexamita Nuclear
// note: tables 7 and 8 are not available in NCBI
char genetic_code9[]  = "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"; // Echinoderm and Flatworm Mitochondrial
char genetic_code10[] = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSCCWCLFLF"; // Euplotid Nuclear
char genetic_code11[] = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"; // Bacterial, Archaeal and Plant Plastid
char genetic_code12[] = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLSLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"; // Alternative Yeast Nuclear
char genetic_code13[] = "KNKNTTTTGSGSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"; // Ascidian Mitochondrial
char genetic_code14[] = "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYY*YSSSSWCWCLFLF"; // Alternative Flatworm Mitochondrial
char genetic_code15[] = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YQYSSSS*CWCLFLF"; // Blepharisma Nuclear
char genetic_code16[] = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLYSSSS*CWCLFLF"; // Chlorophycean Mitochondrial
// note: tables 17-20 are not available in NCBI
char genetic_code21[] = "NNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"; // Trematode Mitochondrial
char genetic_code22[] = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLY*SSS*CWCLFLF"; // Scenedesmus obliquus mitochondrial
char genetic_code23[] = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWC*FLF"; // Thraustochytrium Mitochondrial
char genetic_code24[] = "KNKNTTTTSSKSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"; // Pterobranchia mitochondrial
char genetic_code25[] = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSGCWCLFLF"; // Candidate Division SR1 and Gracilibacteria

Alignment::Alignment()
        : vector<Pattern>()
{
    num_states = 0;
    frac_const_sites = 0.0;
    frac_invariant_sites = 0.0;
    codon_table = NULL;
    genetic_code = NULL;
    non_stop_codon = NULL;
    seq_type = SEQ_UNKNOWN;
    STATE_UNKNOWN = 126;
    pars_lower_bound = NULL;
}

string &Alignment::getSeqName(int i) {
    ASSERT(i >= 0 && i < (int)seq_names.size());
    return seq_names[i];
}

vector<string>& Alignment::getSeqNames() {
	return seq_names;
}

int Alignment::getSeqID(string &seq_name) {
    for (int i = 0; i < getNSeq(); i++)
        if (seq_name == getSeqName(i)) return i;
    return -1;
}

int Alignment::getMaxSeqNameLength() {
    int len = 0;
    for (int i = 0; i < getNSeq(); i++)
        if (getSeqName(i).length() > len)
            len = getSeqName(i).length();
    return len;
}

/**
   probability that the observed chi-square exceeds chi2 even if model is correct
   @param deg degree of freedom
   @param chi2 chi-square value
   @return p-value
   */
double chi2prob (int deg, double chi2)
{
    double a = 0.5*deg;
    double x = 0.5*chi2;
    return 1.0-RateGamma::cmpIncompleteGamma (x, a, RateGamma::cmpLnGamma(a));
//	return IncompleteGammaQ (0.5*deg, 0.5*chi2);
} /* chi2prob */


int Alignment::checkAbsentStates(string msg) {
    double *state_freq = new double[num_states];
    computeStateFreq(state_freq);
    string absent_states, rare_states;
    int count = 0;
    // Skip check for PoMo.
    if (seq_type == SEQ_POMO)
      return 0;
    for (int i = 0; i < num_states; i++)
        if (state_freq[i] == 0.0) {
            if (!absent_states.empty())
                absent_states += ", ";
            absent_states += convertStateBackStr(i);
            count++;
        } else if (state_freq[i] <= MIN_FREQUENCY) {
            if (!rare_states.empty())
                rare_states += ", ";
            rare_states += convertStateBackStr(i);
        }
    if (count >= num_states-1 && Params::getInstance().fixed_branch_length != BRLEN_FIX)
        outError("Only one state is observed in " + msg);
    if (!absent_states.empty())
        cout << "NOTE: State(s) " << absent_states << " not present in " << msg << " and thus removed from Markov process to prevent numerical problems" << endl;
    if (!rare_states.empty())
        cout << "WARNING: States(s) " << rare_states << " rarely appear in " << msg << " and may cause numerical problems" << endl;
    delete[] state_freq;
    return count;
}

void Alignment::checkSeqName() {
    ostringstream warn_str;
    StrVector::iterator it;
    for (it = seq_names.begin(); it != seq_names.end(); it++) {
        string orig_name = (*it);
        if (renameString(*it))
            warn_str << orig_name << " -> " << (*it) << endl;
    }
    if (warn_str.str() != "") {
        string str = "Some sequence names are changed as follows:\n";
        outWarning(str + warn_str.str());
    }
    // now check that sequence names are different
    StrVector names;
    names.insert(names.begin(), seq_names.begin(), seq_names.end());
    sort(names.begin(), names.end());
    bool ok = true;
    for (it = names.begin(); it != names.end(); it++) {
        if (it+1==names.end()) break;
        if (*it == *(it+1)) {
            cout << "ERROR: Duplicated sequence name " << *it << endl;
            ok = false;
        }
    }
    if (!ok) outError("Please rename sequences listed above!");

    double *state_freq = new double[num_states];
//    double *freq_per_sequence = new double[num_states*getNSeq()];
    double *freq_per_sequence = new double[num_states];
    unsigned *count_per_seq = new unsigned[num_states*getNSeq()];
    computeStateFreq(state_freq);
//    computeStateFreqPerSequence(freq_per_sequence);
    countStatePerSequence(count_per_seq);

    int i, df = -1;
    for (i = 0; i < num_states; i++)
        if (state_freq[i] > 0.0) df++;
    
    if (seq_type == SEQ_POMO)
        cout << "NOTE: The composition test for PoMo only tests the proportion of fixed states!" << endl;

    int max_len = getMaxSeqNameLength()+1;
    cout.width(max_len+14);
    cout << right << "Gap/Ambiguity" << "  Composition  p-value"<< endl;
    int num_problem_seq = 0;
    int total_gaps = 0;
    cout.precision(2);
    int num_failed = 0;
    for (i = 0; i < seq_names.size(); i++) {
        int j;
        int num_gaps = getNSite() - countProperChar(i);
        total_gaps += num_gaps;
        double percent_gaps = ((double)num_gaps / getNSite())*100.0;
        cout.width(4);
        cout << right << i+1 << "  ";
        cout.width(max_len);
        cout << left << seq_names[i] << " ";
        cout.width(6);
        cout << right << percent_gaps << "%";
        if (percent_gaps > 50) {
            num_problem_seq++;
        }

        double chi2 = 0.0;
        unsigned sum_count = 0;

        if (seq_type == SEQ_POMO) {
            // FIXME: Number of nucleotides hardcoded here.
            int nnuc = 4;
            df = nnuc-1;
            // Have to normalize allele frequencies.
            double state_freq_norm[nnuc];
            double sum_freq = 0.0;
            for (j = 0; j < nnuc; j++) {
                sum_freq += state_freq[j];
                state_freq_norm[j] = state_freq[j];
            }
            for (j = 0; j < nnuc; j++) {
                state_freq_norm[j] /= sum_freq;
            }

            for (j = 0; j < nnuc; j++)
                sum_count += count_per_seq[i*num_states+j];
            double sum_inv = 1.0/sum_count;
            for (j = 0; j < nnuc; j++)
                freq_per_sequence[j] = count_per_seq[i*num_states+j]*sum_inv;
            for (j = 0; j < nnuc; j++)
                chi2 += (state_freq_norm[j] - freq_per_sequence[j]) * (state_freq_norm[j] - freq_per_sequence[j]) / state_freq_norm[j];

            // chi2 *= getNSite();
            chi2 *= sum_count;
            double pvalue = chi2prob(nnuc-1, chi2);
            if (pvalue < 0.05) {
                cout << "    failed ";
                num_failed++;
            } else
                cout << "    passed ";
            cout.width(9);
            cout << right << pvalue*100 << "%";
        } else {
            for (j = 0; j < num_states; j++)
                sum_count += count_per_seq[i*num_states+j];
            double sum_inv = 1.0/sum_count;
            for (j = 0; j < num_states; j++)
                freq_per_sequence[j] = count_per_seq[i*num_states+j]*sum_inv;
            for (j = 0; j < num_states; j++)
                if (state_freq[j] > 0.0)
                    chi2 += (state_freq[j] - freq_per_sequence[j]) * (state_freq[j] - freq_per_sequence[j]) / state_freq[j];
            
            chi2 *= sum_count;
            double pvalue = chi2prob(df, chi2);
            if (pvalue < 0.05) {
                cout << "    failed ";
                num_failed++;
            } else
                cout << "    passed ";
            cout.width(9);
            cout << right << pvalue*100 << "%";
        }
        cout << endl;
    }
    if (num_problem_seq) cout << "WARNING: " << num_problem_seq << " sequences contain more than 50% gaps/ambiguity" << endl;
    cout << "**** ";
    cout.width(max_len+2);
    cout << left << " TOTAL  ";
    cout.width(6);
    cout << right << ((double)total_gaps/getNSite())/getNSeq()*100 << "% ";
    cout << " " << num_failed << " sequences failed composition chi2 test (p-value<5%; df=" << df << ")" << endl;
    cout.precision(3);

    delete [] count_per_seq;
    delete [] freq_per_sequence;
    delete [] state_freq;
}

int Alignment::checkIdenticalSeq()
{
	int num_identical = 0;
    IntVector checked;
    checked.resize(getNSeq(), 0);
	for (int seq1 = 0; seq1 < getNSeq(); seq1++) {
        if (checked[seq1]) continue;
		bool first = true;
		for (int seq2 = seq1+1; seq2 < getNSeq(); seq2++) {
			bool equal_seq = true;
			for (iterator it = begin(); it != end(); it++)
				if  ((*it)[seq1] != (*it)[seq2]) {
					equal_seq = false;
					break;
				}
			if (equal_seq) {
				if (first)
					cout << "WARNING: Identical sequences " << getSeqName(seq1);
				cout << ", " << getSeqName(seq2);
				num_identical++;
				checked[seq2] = 1;
				first = false;
			}
		}
		checked[seq1] = 1;
		if (!first) cout << endl;
	}
	if (num_identical)
		outWarning("Some identical sequences found that should be discarded before the analysis");
	return num_identical;
}

Alignment *Alignment::removeIdenticalSeq(string not_remove, bool keep_two, StrVector &removed_seqs, StrVector &target_seqs)
{
    IntVector checked;
    vector<bool> removed;
    checked.resize(getNSeq(), 0);
    removed.resize(getNSeq(), false);
    int seq1;

	for (seq1 = 0; seq1 < getNSeq(); seq1++) {
        if (checked[seq1]) continue;
        bool first_ident_seq = true;
		for (int seq2 = seq1+1; seq2 < getNSeq(); seq2++) {
			if (getSeqName(seq2) == not_remove || removed[seq2]) continue;
			bool equal_seq = true;
			for (iterator it = begin(); it != end(); it++)
				if  ((*it)[seq1] != (*it)[seq2]) {
					equal_seq = false;
					break;
				}
			if (equal_seq) {
				if (removed_seqs.size() < getNSeq()-3 && (!keep_two || !first_ident_seq)) {
					removed_seqs.push_back(getSeqName(seq2));
					target_seqs.push_back(getSeqName(seq1));
					removed[seq2] = true;
				} else {
                    cout << "NOTE: " << getSeqName(seq2) << " is identical to " << getSeqName(seq1) << " but kept for subsequent analysis" << endl;
                }
				checked[seq2] = 1;
				first_ident_seq = false;
			}
		}
		checked[seq1] = 1;
	}

	if (removed_seqs.size() > 0) {
		if (removed_seqs.size() >= getNSeq()-3)
			outWarning("Your alignment contains too many identical sequences!");
		IntVector keep_seqs;
		for (seq1 = 0; seq1 < getNSeq(); seq1++)
			if (!removed[seq1]) keep_seqs.push_back(seq1);
		Alignment *aln = new Alignment;
		aln->extractSubAlignment(this, keep_seqs, 0);
		return aln;
	} else return this;
}


bool Alignment::isGapOnlySeq(int seq_id) {
    ASSERT(seq_id < getNSeq());
    for (iterator it = begin(); it != end(); it++)
        if ((*it)[seq_id] != STATE_UNKNOWN) {
            return false;
        }
    return true;
}

Alignment *Alignment::removeGappySeq() {
	IntVector keep_seqs;
	int i, nseq = getNSeq();
	for (i = 0; i < nseq; i++)
		if (! isGapOnlySeq(i)) {
			keep_seqs.push_back(i);
		}
	if (keep_seqs.size() == nseq)
		return this;
    // 2015-12-03: if resulting alignment has too few seqs, try to add some back
    if (keep_seqs.size() < 3 && getNSeq() >= 3) {
        for (i = 0; i < nseq && keep_seqs.size() < 3; i++)
            if (isGapOnlySeq(i))
                keep_seqs.push_back(i);
    }
	Alignment *aln = new Alignment;
	aln->extractSubAlignment(this, keep_seqs, 0);
	return aln;
}

void Alignment::checkGappySeq(bool force_error) {
    int nseq = getNSeq(), i;
    int wrong_seq = 0;
    for (i = 0; i < nseq; i++)
        if (isGapOnlySeq(i)) {
            outWarning("Sequence " + getSeqName(i) + " contains only gaps or missing data");
            wrong_seq++;
        }
    if (wrong_seq && force_error) {
        outError("Some sequences (see above) are problematic, please check your alignment again");
    }
}

Alignment::Alignment(char *filename, char *sequence_type, InputType &intype, string model) : vector<Pattern>() {
    name = "Noname";
    this->model_name = model;
    if (sequence_type)
        this->sequence_type = sequence_type;
    aln_file = filename;
    num_states = 0;
    frac_const_sites = 0.0;
    frac_invariant_sites = 0.0;
    codon_table = NULL;
    genetic_code = NULL;
    non_stop_codon = NULL;
    seq_type = SEQ_UNKNOWN;
    STATE_UNKNOWN = 126;
    pars_lower_bound = NULL;
    cout << "Reading alignment file " << filename << " ... ";
    intype = detectInputFile(filename);

    try {

        if (intype == IN_NEXUS) {
            cout << "Nexus format detected" << endl;
            readNexus(filename);
        } else if (intype == IN_FASTA) {
            cout << "Fasta format detected" << endl;
            readFasta(filename, sequence_type);
        } else if (intype == IN_PHYLIP) {
            cout << "Phylip format detected" << endl;
            if (Params::getInstance().phylip_sequential_format)
                readPhylipSequential(filename, sequence_type);
            else
                readPhylip(filename, sequence_type);
        } else if (intype == IN_COUNTS) {
            cout << "Counts format (PoMo) detected" << endl;
            readCountsFormat(filename, sequence_type);
        } else if (intype == IN_CLUSTAL) {
            cout << "Clustal format detected" << endl;
            readClustal(filename, sequence_type);
        } else if (intype == IN_MSF) {
            cout << "MSF format detected" << endl;
            readMSF(filename, sequence_type);
        } else {
            outError("Unknown sequence format, please use PHYLIP, FASTA, CLUSTAL, MSF, or NEXUS format");
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

    countConstSite();

    cout << "Alignment has " << getNSeq() << " sequences with " << getNSite()
         << " columns, " << getNPattern() << " distinct patterns" << endl
         << num_informative_sites << " parsimony-informative, "
         << num_variant_sites-num_informative_sites << " singleton sites, "
         << (int)(frac_const_sites*getNSite()) << " constant sites" << endl;
    buildSeqStates();
    checkSeqName();
    // OBSOLETE: identical sequences are handled later
//	checkIdenticalSeq();
    //cout << "Number of character states is " << num_states << endl;
    //cout << "Number of patterns = " << size() << endl;
    //cout << "Fraction of constant sites: " << frac_const_sites << endl;

}

bool Alignment::isStopCodon(int state) {
    // 2017-05-27: all stop codon removed from Markov process
    return false;

	if (seq_type != SEQ_CODON || state >= num_states) return false;
	ASSERT(genetic_code);
	return (genetic_code[state] == '*');
}

int Alignment::getNumNonstopCodons() {
    if (seq_type != SEQ_CODON) return num_states;
	ASSERT(genetic_code);
	int c = 0;
	for (char *ch = genetic_code; *ch != 0; ch++)
		if (*ch != '*') c++;
	return c;
}

bool Alignment::isStandardGeneticCode() {
    if (seq_type != SEQ_CODON) return false;
	return (genetic_code == genetic_code1 || genetic_code == genetic_code11);
}

void Alignment::buildSeqStates(bool add_unobs_const) {
	string unobs_const;
	if (add_unobs_const) unobs_const = getUnobservedConstPatterns();
	seq_states.clear();
	seq_states.resize(getNSeq());
	for (int seq = 0; seq < getNSeq(); seq++) {
		vector<bool> has_state;
		has_state.resize(STATE_UNKNOWN+1, false);
		for (int site = 0; site < getNPattern(); site++)
			has_state[at(site)[seq]] = true;
		for (string::iterator it = unobs_const.begin(); it != unobs_const.end(); it++)
			has_state[*it] = true;
        seq_states[seq].clear();
		for (int state = 0; state < STATE_UNKNOWN; state++)
			if (has_state[state])
				seq_states[seq].push_back(state);
	}
}

int Alignment::readNexus(char *filename) {
    NxsTaxaBlock *taxa_block;
    NxsAssumptionsBlock *assumptions_block;
    NxsDataBlock *data_block = NULL;
    NxsTreesBlock *trees_block = NULL;
    NxsCharactersBlock *char_block = NULL;

    taxa_block = new NxsTaxaBlock();
    assumptions_block = new NxsAssumptionsBlock(taxa_block);
    data_block = new NxsDataBlock(taxa_block, assumptions_block);
    char_block = new NxsCharactersBlock(taxa_block, assumptions_block);
    trees_block = new TreesBlock(taxa_block);

    MyReader nexus(filename);

    nexus.Add(taxa_block);
    nexus.Add(assumptions_block);
    nexus.Add(data_block);
	nexus.Add(char_block);
    nexus.Add(trees_block);

    MyToken token(nexus.inf);
    nexus.Execute(token);

	if (data_block->GetNTax() && char_block->GetNTax()) {
		outError("I am confused since both DATA and CHARACTERS blocks were specified");
		return 0;
	}

    if (data_block->GetNTax() == 0 && char_block->GetNTax() == 0) {
        outError("No DATA or CHARACTERS blocks found");
        return 0;
    }

    if (char_block->GetNTax() > 0) {
        extractDataBlock(char_block);
        if (verbose_mode >= VB_DEBUG)
            char_block->Report(cout);
    } else {
        extractDataBlock(data_block);
        if (verbose_mode >= VB_DEBUG)
            data_block->Report(cout);
    }

    delete trees_block;
    delete char_block;
    delete data_block;
    delete assumptions_block;
    delete taxa_block;
    return 1;
}

void Alignment::computeUnknownState() {
    switch (seq_type) {
    case SEQ_DNA: STATE_UNKNOWN = 18; break;
    case SEQ_PROTEIN: STATE_UNKNOWN = 23; break;
    case SEQ_POMO: {
        if (pomo_sampling_method == SAMPLING_SAMPLED) STATE_UNKNOWN = num_states;
        else STATE_UNKNOWN = 0xffffffff;
        break;
    }
    default: STATE_UNKNOWN = num_states; break;
    }
}

int getDataBlockMorphStates(NxsCharactersBlock *data_block) {
    int nseq = data_block->GetNTax();
    int nsite = data_block->GetNCharTotal();
    int seq, site;
    char ch;
    int nstates = 0;

    for (seq = 0; seq < nseq; seq++)
        for (site = 0; site < nsite; site++) {
            int nstate = data_block->GetNumStates(seq, site);
            if (nstate == 0)
                continue;
            if (nstate == 1) {
                ch = data_block->GetState(seq, site, 0);
                if (!isalnum(ch)) continue;
                if (ch >= '0' && ch <= '9')
                    ch = ch - '0' + 1;
                else if (ch >= 'A' && ch <= 'Z')
                    ch = ch - 'A' + 11;
                else
                    outError(data_block->GetTaxonLabel(seq) + " has invalid single state " + ch + " at site " + convertIntToString(site+1));
                if (ch > nstates) nstates = ch;
                continue;
            }
            //cout << "NOTE: " << data_block->GetTaxonLabel(seq) << " has ambiguous state at site " << site+1 << " which is treated as unknown" << endl;
        }
    return nstates;
}

void Alignment::extractDataBlock(NxsCharactersBlock *data_block) {
    int nseq = data_block->GetNTax();
    int nsite = data_block->GetNCharTotal();
    char *symbols = NULL;
    //num_states = strlen(symbols);
    char char_to_state[NUM_CHAR];
    char state_to_char[NUM_CHAR];
    
    if (!data_block->GetMatrix())
        outError("MATRIX command undeclared or invalid");
    
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
        seq_type = SEQ_DNA;
    } else if (data_type == NxsCharactersBlock::protein) {
        num_states = 20;
        symbols = symbols_protein;
        seq_type = SEQ_PROTEIN;
    } else {
    	// standard morphological character
//        num_states = data_block->GetMaxObsNumStates();
        num_states = getDataBlockMorphStates(data_block);
        if (num_states > 32)
        	outError("Number of states can not exceed 32");
        if (num_states < 2)
        	outError("Number of states can not be below 2");
        if (num_states == 2)
        	seq_type = SEQ_BINARY;
        else
    		seq_type = SEQ_MORPH;
        symbols = symbols_morph;
    }

    computeUnknownState();
    memset(char_to_state, STATE_UNKNOWN, NUM_CHAR);
    memset(state_to_char, '?', NUM_CHAR);
    for (int i = 0; i < strlen(symbols); i++) {
        char_to_state[(int)symbols[i]] = i;
        state_to_char[i] = symbols[i];
    }
    state_to_char[(int)STATE_UNKNOWN] = '-';


    int seq, site;

    if (data_block->taxa->GetNumTaxonLabels() == 0)
        outError("MATRIX not found, make sure nexus command before MATRIX ends with semi-colon (;)");

    if (data_block->taxa->GetNumTaxonLabels() != nseq)
        outError("ntax is different from number of matrix rows");

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
                pat.push_back(STATE_UNKNOWN);
            else if (nstate == 1) {
                pat.push_back(char_to_state[(int)data_block->GetState(seq, site, 0)]);
            } else if (data_type == NxsCharactersBlock::dna || data_type == NxsCharactersBlock::rna || data_type == NxsCharactersBlock::nucleotide) {
                // 2018-06-07: correctly interpret ambiguous nucleotide
                char pat_ch = 0;
                for (int state = 0; state < nstate; state++) {
                    pat_ch |= (1 << char_to_state[(int)data_block->GetState(seq, site, state)]);
                }
                pat_ch += 3;
                pat.push_back(pat_ch);
            } else {
                // other ambiguous characters are treated as unknown
                stringstream str;
                str << "Sequence " << seq_names[seq] << " site " << site+1 << ": {";
                for (int state = 0; state < nstate; state++) {
                    str << data_block->GetState(seq, site, state);
                }
                str << "} treated as unknown character";
                outWarning(str.str());
                pat.push_back(STATE_UNKNOWN);
            }
        }
        num_gaps_only += addPattern(pat, site);
    }
    if (num_gaps_only)
        cout << "WARNING: " << num_gaps_only << " sites contain only gaps or ambiguous characters." << endl;
    if (verbose_mode >= VB_MAX)
        for (site = 0; site < size(); site++) {
            for (seq = 0; seq < nseq; seq++)
                cout << state_to_char[(int)(*this)[site][seq]];
            cout << "  " << (*this)[site].frequency << endl;
        }
}

/**
	determine if the pattern is constant. update the is_const variable.
*/
void Alignment::computeConst(Pattern &pat) {
    bool is_const = true;
    bool is_invariant = false;
    bool is_informative = false;
    // critical fix: const_char was set wrongly to num_states in some data type (binary, codon),
    // causing wrong log-likelihood computation for +I or +I+G model
    pat.const_char = STATE_UNKNOWN+1;
//    if (STATE_UNKNOWN == num_states)
//    	pat.const_char = STATE_UNKNOWN+1;
//    else
//    	pat.const_char = STATE_UNKNOWN;
    StateBitset state_app;
    state_app.reset();
    int j;
    for (j = 0; j < num_states; j++)
    	state_app[j] = 1;

    // number of appearance for each state, to compute is_informative
    size_t num_app[num_states];
    memset(num_app, 0, num_states*sizeof(size_t));

    for (Pattern::iterator i = pat.begin(); i != pat.end(); i++) {
    	StateBitset this_app;
    	getAppearance(*i, this_app);
    	state_app &= this_app;
        if (*i < num_states) {
            num_app[(int)(*i)]++;
        }
//        else if (*i != STATE_UNKNOWN) {
//            // ambiguous characters
//            is_const = false;
//        }
    }
    int count = 0; // number of states with >= 2 appearances
    pat.num_chars = 0; // number of states with >= 1 appearance
    for (j = 0; j < num_states; j++) if (num_app[j]) {
        pat.num_chars++;
        if (num_app[j] >= 2) {
            count++;
        }
    }

    // at least 2 states, each appearing at least twice
    is_informative = (count >= 2);

    // compute is_const
    /*
    is_const = is_const && (pat.num_chars <= 1);
    if (is_const) {
        if (pat.num_chars == 0) // all-gap pattern
            pat.const_char = num_states;
        else {
            // pat.num_chars is 1
            for (j = 0; j < num_states; j++)
                if (num_app[j]) {
                    pat.const_char = j;
                    break;
                }
        }
    }
    */
    is_const = (state_app.count() >= 1);
    if (is_const) {
        if (state_app.count() == num_states) {
            pat.const_char = STATE_UNKNOWN;
        } else if (state_app.count() == 1) {
            for (j = 0; j < num_states; j++)
                if (state_app[j]) {
                    pat.const_char = j;
                    break;
                }
        } else if (seq_type == SEQ_DNA) {
            pat.const_char = num_states-1;
            for (j = 0; j < num_states; j++)
                if (state_app[j])
                    pat.const_char += (1<<j);
        } else if (seq_type == SEQ_PROTEIN) {
            if (state_app[2] && state_app[3]) //4+8, // B = N or D
                pat.const_char = num_states;
            else if (state_app[5] && state_app[6]) //32+64, // Z = Q or E
                pat.const_char = num_states+1;
            else if (state_app[9] && state_app[10]) // 512+1024 // U = I or L
                pat.const_char = num_states+2;
            else ASSERT(0);
        } else {
            ASSERT(0);
        }
    }

//    delete [] num_app;

    // compute is_invariant
    is_invariant = (state_app.count() >= 1);
    ASSERT(is_invariant >= is_const);


    // Wed Jun 28 16:01:30 BST 2017. The calculation of these properties seems
    // to be OKish. They are only used for reports and to calculate the
    // parsimony tree in the beginning anyways.

    // if (seq_type == SEQ_POMO) {
    //     // For PoMo most sites are informative (ambiguous map from data to state space)
    //     is_informative = true;
    //     // For PoMo there are hardly any constant sites
    //     is_const = false;
    //     is_invariant = false;
    // }

    pat.flag = 0;
    if (is_const) pat.flag |= PAT_CONST;
    if (is_invariant) pat.flag |= PAT_INVARIANT;
    if (is_informative) pat.flag |= PAT_INFORMATIVE;
}


void Alignment::printSiteInfo(ostream &out, int part_id) {
    int nsite = getNSite();
    for (int site = 0; site != nsite; site++) {
        Pattern ptn = getPattern(site);
        if (part_id >= 0)
            out << part_id << "\t";
        out << site+1 << "\t";
        if (ptn.isInformative())
            out << "I";
        else if (ptn.isConst()) {
            if (ptn.const_char == STATE_UNKNOWN)
                out << "-";
            else if (ptn.const_char < num_states)
                out << "C";
            else
                out << "c";
        } else
            out << "U";
        out << endl;
    }
}

void Alignment::printSiteInfoHeader(ostream &out, const char* filename, bool partition) {
    out << "# Alignment site statistics" << endl
        << "# This file can be read in MS Excel or in R with command:" << endl
        << "#   tab=read.table('" <<  filename << "',header=TRUE)" << endl
        << "# Columns are tab-separated with following meaning:" << endl;
    if (partition)
        out << "#   Part:   Partition ID" << endl
            << "#   Site:   Site ID within partition (starting from 1 for each partition)" << endl;
    else
        out << "#   Site:   Site ID" << endl;

    out << "#   Stat:   Statistic, I=informative, C=constant, c=constant+ambiguous," << endl
        << "#           U=Uninformative but not constant, -=all-gaps" << endl;
    if (partition)
        out << "Part\t";
    out << "Site\tStat" << endl;
}

void Alignment::printSiteInfo(const char* filename) {
    try {
        ofstream out(filename);
        printSiteInfoHeader(out, filename);
        printSiteInfo(out, -1);
        out.close();
    } catch (...) {
        outError(ERR_WRITE_OUTPUT, filename);
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
            cout << "Site " << site << " contains only gaps or ambiguous characters" << endl;
        //return true;
    }
    PatternIntMap::iterator pat_it = pattern_index.find(pat);
    if (pat_it == pattern_index.end()) { // not found
        pat.frequency = freq;
        computeConst(pat);
        push_back(pat);
        pattern_index[back()] = size()-1;
        site_pattern[site] = size()-1;
    } else {
        int index = pat_it->second;
        at(index).frequency += freq;
        site_pattern[site] = index;
    }
    return gaps_only;
}

void Alignment::addConstPatterns(char *freq_const_patterns) {
	IntVector vec;
	convert_int_vec(freq_const_patterns, vec);
	if (vec.size() != num_states)
		outError("Const pattern frequency vector has different number of states: ", freq_const_patterns);

	int nsite = getNSite(), orig_nsite = getNSite();
	int i;
	for (i = 0; i < vec.size(); i++) {
		nsite += vec[i];
		if (vec[i] < 0)
			outError("Const pattern frequency must be non-negative");
	}
    site_pattern.resize(nsite, -1);
	int nseq = getNSeq();
	nsite = orig_nsite;
	for (i = 0; i < vec.size(); i++) if (vec[i] > 0) {
		Pattern pat;
		pat.resize(nseq, i);
//		if (pattern_index.find(pat) != pattern_index.end()) {
//			outWarning("Constant pattern of all " + convertStateBackStr(i) + " already exists");
//		}
		for (int j = 0; j < vec[i]; j++)
			addPattern(pat, nsite++, 1);
	}
    countConstSite();
    buildSeqStates();
}

void Alignment::orderPatternByNumChars(int pat_type) {
    int nptn = getNPattern();
    int ptn, site, i = 0;
    int *num_chars = new int[nptn];
    int *ptn_order = new int[nptn];
    const int UINT_BITS = sizeof(UINT)*8;
    if (pat_type == PAT_INFORMATIVE)
        num_parsimony_sites = num_informative_sites;
    else
        num_parsimony_sites = num_variant_sites;

    int maxi = (num_parsimony_sites+UINT_BITS-1)/UINT_BITS;
    pars_lower_bound = new UINT[maxi+1];
    UINT sum = 0;
    memset(pars_lower_bound, 0, (maxi+1)*sizeof(UINT));
    for (ptn = 0; ptn < nptn; ptn++) {
        num_chars[ptn] =  -at(ptn).num_chars + (at(ptn).isInvariant())*1024;
        ptn_order[ptn] = ptn;
    }
    quicksort(num_chars, 0, nptn-1, ptn_order);
    ordered_pattern.clear();
    for (ptn = 0, site = 0, i = 0; ptn < nptn; ptn++) {
        if (pat_type == PAT_INFORMATIVE) {
            if (!at(ptn_order[ptn]).isInformative())
                break;
        } else {
            if (at(ptn_order[ptn]).isInvariant())
                break;
        }
        ordered_pattern.push_back(at(ptn_order[ptn]));
        int freq = ordered_pattern.back().frequency;
        UINT num = ordered_pattern.back().num_chars - 1;
        for (int j = 0; j < freq; j++, site++) {
            if (site == UINT_BITS) {
                sum += pars_lower_bound[i];
                i++;
                site = 0;
            }
            pars_lower_bound[i] += num;
        }
    }
    sum += pars_lower_bound[i];
    // now transform lower_bound
//    assert(i == maxi-1);

    for (int j = 0; j <= i; j++) {
        UINT newsum = sum - pars_lower_bound[j];
        pars_lower_bound[j] = sum;
        sum = newsum;
    }

    if (verbose_mode >= VB_MAX) {
//        for (ptn = 0; ptn < nptn; ptn++)
//            cout << at(ptn_order[ptn]).num_chars << " ";
        for (int j = 0; j <= i; j++) {
            cout << pars_lower_bound[j] << " ";
        }
        cout << endl << sum << endl;
    }
    delete [] ptn_order;
    delete [] num_chars;
//    cout << ordered_pattern.size() << " ordered_pattern" << endl;
}

void Alignment::ungroupSitePattern()
{
	vector<Pattern> stored_pat = (*this);
	clear();
	for (int i = 0; i < getNSite(); i++) {
		Pattern pat = stored_pat[getPatternID(i)];
		pat.frequency = 1;
		push_back(pat);
		site_pattern[i] = i;
	}
	pattern_index.clear();
}

void Alignment::regroupSitePattern(int groups, IntVector& site_group)
{
	vector<Pattern> stored_pat = (*this);
	IntVector stored_site_pattern = site_pattern;
	clear();
	site_pattern.clear();
	site_pattern.resize(stored_site_pattern.size(), -1);
	int count = 0;
	for (int g = 0; g < groups; g++) {
		pattern_index.clear();
		for (int i = 0; i < site_group.size(); i++)
		if (site_group[i] == g) {
			count++;
			Pattern pat = stored_pat[stored_site_pattern[i]];
			addPattern(pat, i);
		}
	}
	ASSERT(count == stored_site_pattern.size());
	count = 0;
	for (iterator it = begin(); it != end(); it++)
		count += it->frequency;
	ASSERT(count == getNSite());
	pattern_index.clear();
	//printPhylip("/dev/stdout");
}


/**
	detect the data type of the input sequences
	@param sequences vector of strings
	@return the data type of the input sequences
*/
SeqType Alignment::detectSequenceType(StrVector &sequences) {
    int num_nuc = 0;
    int num_ungap = 0;
    int num_bin = 0;
    int num_alpha = 0;
    int num_digit = 0;

    for (StrVector::iterator it = sequences.begin(); it != sequences.end(); it++)
        for (string::iterator i = it->begin(); i != it->end(); i++) {
            if ((*i) != '?' && (*i) != '-' && (*i) != '.' && *i != 'N' && *i != 'X' &&  (*i) != '~') num_ungap++;
            if ((*i) == 'A' || (*i) == 'C' || (*i) == 'G' || (*i) == 'T' || (*i) == 'U')
                num_nuc++;
            if ((*i) == '0' || (*i) == '1')
                num_bin++;
            if (isalpha(*i)) num_alpha++;
            if (isdigit(*i)) num_digit++;
        }
    if (((double)num_nuc) / num_ungap > 0.9)
        return SEQ_DNA;
    if (((double)num_bin) / num_ungap > 0.9)
        return SEQ_BINARY;
    if (((double)num_alpha) / num_ungap > 0.9)
        return SEQ_PROTEIN;
    if (((double)(num_alpha+num_digit)) / num_ungap > 0.9)
        return SEQ_MORPH;
    return SEQ_UNKNOWN;
}

void Alignment::buildStateMap(char *map, SeqType seq_type) {
    memset(map, STATE_INVALID, NUM_CHAR);
    ASSERT(STATE_UNKNOWN < 126);
    map[(unsigned char)'?'] = STATE_UNKNOWN;
    map[(unsigned char)'-'] = STATE_UNKNOWN;
    map[(unsigned char)'~'] = STATE_UNKNOWN;
    map[(unsigned char)'.'] = STATE_UNKNOWN;
    int len;
    switch (seq_type) {
    case SEQ_BINARY:
        map[(unsigned char)'0'] = 0;
        map[(unsigned char)'1'] = 1;
        return;
    case SEQ_DNA: // DNA
	case SEQ_CODON:
        map[(unsigned char)'A'] = 0;
        map[(unsigned char)'C'] = 1;
        map[(unsigned char)'G'] = 2;
        map[(unsigned char)'T'] = 3;
        map[(unsigned char)'U'] = 3;
        map[(unsigned char)'R'] = 1+4+3; // A or G, Purine
        map[(unsigned char)'Y'] = 2+8+3; // C or T, Pyrimidine
        map[(unsigned char)'N'] = STATE_UNKNOWN;
        map[(unsigned char)'X'] = STATE_UNKNOWN;
        map[(unsigned char)'W'] = 1+8+3; // A or T, Weak
        map[(unsigned char)'S'] = 2+4+3; // G or C, Strong
        map[(unsigned char)'M'] = 1+2+3; // A or C, Amino
        map[(unsigned char)'K'] = 4+8+3; // G or T, Keto
        map[(unsigned char)'B'] = 2+4+8+3; // C or G or T
        map[(unsigned char)'H'] = 1+2+8+3; // A or C or T
        map[(unsigned char)'D'] = 1+4+8+3; // A or G or T
        map[(unsigned char)'V'] = 1+2+4+3; // A or G or C
        return;
    case SEQ_PROTEIN: // Protein
        for (int i = 0; i < 20; i++)
            map[(int)symbols_protein[i]] = i;
        map[(int)symbols_protein[20]] = STATE_UNKNOWN;
//		map[(unsigned char)'B'] = 4+8+19; // N or D
//		map[(unsigned char)'Z'] = 32+64+19; // Q or E
        map[(unsigned char)'B'] = 20; // N or D
        map[(unsigned char)'Z'] = 21; // Q or E
        map[(unsigned char)'J'] = 22; // I or L
        map[(unsigned char)'*'] = STATE_UNKNOWN; // stop codon
        map[(unsigned char)'U'] = STATE_UNKNOWN; // 21st amino acid
        map[(unsigned char)'O'] = STATE_UNKNOWN; // 22nd amino acid

        return;
    case SEQ_MULTISTATE:
        for (int i = 0; i <= STATE_UNKNOWN; i++)
            map[i] = i;
        return;
    case SEQ_MORPH: // Protein
    	len = strlen(symbols_morph);
        for (int i = 0; i < len; i++)
            map[(int)symbols_morph[i]] = i;
        return;
    default:
        return;
    }
}


/**
	convert a raw characer state into ID, indexed from 0
	@param state input raw state
	@param seq_type data type (SEQ_DNA, etc.)
	@return state ID
*/
char Alignment::convertState(char state, SeqType seq_type) {
    if (state == '?' || state == '-' || state == '.' || state == '~')
        return STATE_UNKNOWN;

    char *loc;

    switch (seq_type) {
    case SEQ_BINARY:
        switch (state) {
        case '0':
            return 0;
        case '1':
            return 1;
        default:
            return STATE_INVALID;
        		}
		break;
    case SEQ_DNA: // DNA
        switch (state) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        case 'U':
            return 3;
        case 'R':
            return 1+4+3; // A or G, Purine
        case 'Y':
            return 2+8+3; // C or T, Pyrimidine
        case 'O':
        case 'N':
        case 'X':
            return STATE_UNKNOWN;
        case 'W':
            return 1+8+3; // A or T, Weak
        case 'S':
            return 2+4+3; // G or C, Strong
        case 'M':
            return 1+2+3; // A or C, Amino
        case 'K':
            return 4+8+3; // G or T, Keto
        case 'B':
            return 2+4+8+3; // C or G or T
        case 'H':
            return 1+2+8+3; // A or C or T
        case 'D':
            return 1+4+8+3; // A or G or T
        case 'V':
            return 1+2+4+3; // A or G or C
        default:
            return STATE_INVALID; // unrecognize character
        }
        return state;
    case SEQ_PROTEIN: // Protein
//		if (state == 'B') return 4+8+19;
//		if (state == 'Z') return 32+64+19;
		if (state == 'B') return 20;
		if (state == 'Z') return 21;
		if (state == 'J') return 22;
        if (state == '*') return STATE_UNKNOWN; // stop codon
        if (state == 'U') return STATE_UNKNOWN; // 21st amino-acid
        if (state == 'O') return STATE_UNKNOWN; // 22nd amino-acid
        loc = strchr(symbols_protein, state);

        if (!loc) return STATE_INVALID; // unrecognize character
        state = loc - symbols_protein;
        if (state < 20)
            return state;
        else
            return STATE_UNKNOWN;
    case SEQ_MORPH: // Standard morphological character
        loc = strchr(symbols_morph, state);

        if (!loc) return STATE_INVALID; // unrecognize character
        state = loc - symbols_morph;
	    return state;
    default:
        return STATE_INVALID;
    }
}

// TODO: state should int
char Alignment::convertState(char state) {
	return convertState(state, seq_type);
}


// TODO: state should int
char Alignment::convertStateBack(char state) {
    if (state == STATE_UNKNOWN) return '-';
    if (state == STATE_INVALID) return '?';

    switch (seq_type) {
    case SEQ_BINARY:
        switch (state) {
        case 0:
            return '0';
        case 1:
            return '1';
        default:
            return STATE_INVALID;
        }
    case SEQ_DNA: // DNA
        switch (state) {
        case 0:
            return 'A';
        case 1:
            return 'C';
        case 2:
            return 'G';
        case 3:
            return 'T';
        case 1+4+3:
            return 'R'; // A or G, Purine
        case 2+8+3:
            return 'Y'; // C or T, Pyrimidine
        case 1+8+3:
            return 'W'; // A or T, Weak
        case 2+4+3:
            return 'S'; // G or C, Strong
        case 1+2+3:
            return 'M'; // A or C, Amino
        case 4+8+3:
            return 'K'; // G or T, Keto
        case 2+4+8+3:
            return 'B'; // C or G or T
        case 1+2+8+3:
            return 'H'; // A or C or T
        case 1+4+8+3:
            return 'D'; // A or G or T
        case 1+2+4+3:
            return 'V'; // A or G or C
        default:
            return '?'; // unrecognize character
        }
        return state;
    case SEQ_PROTEIN: // Protein
        if (state < 20)
            return symbols_protein[(int)state];
		else if (state == 20) return 'B';
		else if (state == 21) return 'Z';
		else if (state == 22) return 'J';
//		else if (state == 4+8+19) return 'B';
//		else if (state == 32+64+19) return 'Z';
        else
            return '-';
    case SEQ_MORPH:
    	// morphological state
        if (state < strlen(symbols_morph))
            return symbols_morph[(int)state];
        else
            return '-';
    default:
    	// unknown
    	return '*';
    }
}

string Alignment::convertStateBackStr(StateType state) {
	string str;
    if (seq_type == SEQ_POMO)
        return string("POMO")+convertIntToString(state);
    if (seq_type == SEQ_MULTISTATE)
        return " " + convertIntToString(state);
	if (seq_type == SEQ_CODON) {
        // codon data
        if (state >= num_states) return "???";
        assert(codon_table);
        state = codon_table[(int)state];
        str = symbols_dna[state/16];
        str += symbols_dna[(state%16)/4];
        str += symbols_dna[state%4];
        return str;
	}
    // all other data types
    str = convertStateBack(state);
	return str;
}

void Alignment::convertStateStr(string &str, SeqType seq_type) {
    for (string::iterator it = str.begin(); it != str.end(); it++)
        (*it) = convertState(*it, seq_type);
}

void Alignment::initCodon(char *gene_code_id) {
    // build index from 64 codons to non-stop codons
	int transl_table = 1;
	if (strlen(gene_code_id) > 0) {
		try {
			transl_table = convert_int(gene_code_id);
		} catch (string &str) {
			outError("Wrong genetic code ", gene_code_id);
		}
		switch (transl_table) {
		case 1: genetic_code = genetic_code1; break;
		case 2: genetic_code = genetic_code2; break;
		case 3: genetic_code = genetic_code3; break;
		case 4: genetic_code = genetic_code4; break;
		case 5: genetic_code = genetic_code5; break;
		case 6: genetic_code = genetic_code6; break;
		case 9: genetic_code = genetic_code9; break;
		case 10: genetic_code = genetic_code10; break;
		case 11: genetic_code = genetic_code11; break;
		case 12: genetic_code = genetic_code12; break;
		case 13: genetic_code = genetic_code13; break;
		case 14: genetic_code = genetic_code14; break;
		case 15: genetic_code = genetic_code15; break;
		case 16: genetic_code = genetic_code16; break;
		case 21: genetic_code = genetic_code21; break;
		case 22: genetic_code = genetic_code22; break;
		case 23: genetic_code = genetic_code23; break;
		case 24: genetic_code = genetic_code24; break;
		case 25: genetic_code = genetic_code25; break;
		default:
			outError("Wrong genetic code ", gene_code_id);
			break;
		}
	} else {
		genetic_code = genetic_code1;
	}
	ASSERT(strlen(genetic_code) == 64);

	int codon;
	num_states = 0;
	for (codon = 0; codon < strlen(genetic_code); codon++)
		if (genetic_code[codon] != '*')
			num_states++; // only count non-stop codons
	codon_table = new char[num_states];
	non_stop_codon = new char[strlen(genetic_code)];
	int state = 0;
	for (int codon = 0; codon < strlen(genetic_code); codon++) {
		if (genetic_code[codon] != '*') {
			non_stop_codon[codon] = state++;
			codon_table[(int)non_stop_codon[codon]] = codon;
		} else {
			non_stop_codon[codon] = STATE_INVALID;
		}
	}
//	num_states = strlen(genetic_code);
//	codon_table = new char[num_states];
//	non_stop_codon = new char[strlen(genetic_code)];
//	int state = 0;
//	for (int codon = 0; codon < strlen(genetic_code); codon++) {
//		non_stop_codon[codon] = state++;
//		codon_table[(int)non_stop_codon[codon]] = codon;
//	}
//	cout << "num_states = " << num_states << endl;
}

int getMorphStates(StrVector &sequences) {
	char maxstate = 0;
	for (StrVector::iterator it = sequences.begin(); it != sequences.end(); it++)
		for (string::iterator pos = it->begin(); pos != it->end(); pos++)
			if ((*pos) > maxstate && isalnum(*pos)) maxstate = *pos;
	if (maxstate >= '0' && maxstate <= '9') return (maxstate - '0' + 1);
	if (maxstate >= 'A' && maxstate <= 'V') return (maxstate - 'A' + 11);
	return 0;
}

SeqType Alignment::getSeqType(const char *sequence_type) {
    SeqType user_seq_type = SEQ_UNKNOWN;
    if (strcmp(sequence_type, "BIN") == 0) {
        user_seq_type = SEQ_BINARY;
    } else if (strcmp(sequence_type, "NT") == 0 || strcmp(sequence_type, "DNA") == 0) {
        user_seq_type = SEQ_DNA;
    } else if (strcmp(sequence_type, "AA") == 0 || strcmp(sequence_type, "PROT") == 0) {
        user_seq_type = SEQ_PROTEIN;
    } else if (strncmp(sequence_type, "NT2AA", 5) == 0) {
        user_seq_type = SEQ_PROTEIN;
    } else if (strcmp(sequence_type, "NUM") == 0 || strcmp(sequence_type, "MORPH") == 0) {
        user_seq_type = SEQ_MORPH;
    } else if (strcmp(sequence_type, "TINA") == 0 || strcmp(sequence_type, "MULTI") == 0) {
        user_seq_type = SEQ_MULTISTATE;
    } else if (strncmp(sequence_type, "CODON", 5) == 0) {
        user_seq_type = SEQ_CODON;
    }
    return user_seq_type;
}

int Alignment::buildPattern(StrVector &sequences, char *sequence_type, int nseq, int nsite) {
    int seq_id;
    ostringstream err_str;
    codon_table = NULL;
    genetic_code = NULL;
    non_stop_codon = NULL;


    if (nseq != seq_names.size()) throw "Different number of sequences than specified";

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
    case SEQ_MORPH:
        num_states = getMorphStates(sequences);
        if (num_states < 2 || num_states > 32) throw "Invalid number of states.";
        cout << "Alignment most likely contains " << num_states << "-state morphological data" << endl;
        break;
    case SEQ_POMO:
        throw "Counts Format pattern is built in Alignment::readCountsFormat().";
        break;
    default:
        if (!sequence_type)
            throw "Unknown sequence type.";
    }
    bool nt2aa = false;
    if (sequence_type && strcmp(sequence_type,"") != 0) {
        SeqType user_seq_type;
        if (strcmp(sequence_type, "BIN") == 0) {
            num_states = 2;
            user_seq_type = SEQ_BINARY;
        } else if (strcmp(sequence_type, "NT") == 0 || strcmp(sequence_type, "DNA") == 0) {
            num_states = 4;
            user_seq_type = SEQ_DNA;
        } else if (strcmp(sequence_type, "AA") == 0 || strcmp(sequence_type, "PROT") == 0) {
            num_states = 20;
            user_seq_type = SEQ_PROTEIN;
        } else if (strncmp(sequence_type, "NT2AA", 5) == 0) {
            if (seq_type != SEQ_DNA)
                outWarning("Sequence type detected as non DNA!");
            initCodon(&sequence_type[5]);
            seq_type = user_seq_type = SEQ_PROTEIN;
            num_states = 20;
            nt2aa = true;
            cout << "Translating to amino-acid sequences with genetic code " << &sequence_type[5] << " ..." << endl;
        } else if (strcmp(sequence_type, "NUM") == 0 || strcmp(sequence_type, "MORPH") == 0) {
            num_states = getMorphStates(sequences);
            if (num_states < 2 || num_states > 32) throw "Invalid number of states";
            user_seq_type = SEQ_MORPH;
        } else if (strcmp(sequence_type, "TINA") == 0 || strcmp(sequence_type, "MULTI") == 0) {
            cout << "Multi-state data with " << num_states << " alphabets" << endl;
            user_seq_type = SEQ_MULTISTATE;
        } else if (strncmp(sequence_type, "CODON", 5) == 0) {
            if (seq_type != SEQ_DNA)
				outWarning("You want to use codon models but the sequences were not detected as DNA");
            seq_type = user_seq_type = SEQ_CODON;
        	initCodon(&sequence_type[5]);
            cout << "Converting to codon sequences with genetic code " << &sequence_type[5] << " ..." << endl;
        } else
            throw "Invalid sequence type.";
        if (user_seq_type != seq_type && seq_type != SEQ_UNKNOWN)
            outWarning("Your specified sequence type is different from the detected one");
        seq_type = user_seq_type;
    }

    // now convert to patterns
    int site, seq, num_gaps_only = 0;

    char char_to_state[NUM_CHAR];
    char AA_to_state[NUM_CHAR];
    computeUnknownState();
    if (nt2aa) {
        buildStateMap(char_to_state, SEQ_DNA);
        buildStateMap(AA_to_state, SEQ_PROTEIN);
    } else
        buildStateMap(char_to_state, seq_type);

    Pattern pat;
    pat.resize(nseq);
    int step = ((seq_type == SEQ_CODON || nt2aa) ? 3 : 1);
    if (nsite % step != 0)
    	outError("Number of sites is not multiple of 3");
    site_pattern.resize(nsite/step, -1);
    clear();
    pattern_index.clear();
    int num_error = 0;
    for (site = 0; site < nsite; site+=step) {
        for (seq = 0; seq < nseq; seq++) {
            //char state = convertState(sequences[seq][site], seq_type);
            char state = char_to_state[(int)(sequences[seq][site])];
            if (seq_type == SEQ_CODON || nt2aa) {
            	// special treatment for codon
            	char state2 = char_to_state[(int)(sequences[seq][site+1])];
            	char state3 = char_to_state[(int)(sequences[seq][site+2])];
            	if (state < 4 && state2 < 4 && state3 < 4) {
//            		state = non_stop_codon[state*16 + state2*4 + state3];
            		state = state*16 + state2*4 + state3;
            		if (genetic_code[(int)state] == '*') {
                        err_str << "Sequence " << seq_names[seq] << " has stop codon " <<
                        		sequences[seq][site] << sequences[seq][site+1] << sequences[seq][site+2] <<
                        		" at site " << site+1 << endl;
                        num_error++;
                        state = STATE_UNKNOWN;
            		} else if (nt2aa) {
                        state = AA_to_state[(int)genetic_code[(int)state]];
                    } else {
                        state = non_stop_codon[(int)state];
                    }
            	} else if (state == STATE_INVALID || state2 == STATE_INVALID || state3 == STATE_INVALID) {
            		state = STATE_INVALID;
            	} else {
            		if (state != STATE_UNKNOWN || state2 != STATE_UNKNOWN || state3 != STATE_UNKNOWN) {
            			ostringstream warn_str;
                        warn_str << "Sequence " << seq_names[seq] << " has ambiguous character " <<
                        		sequences[seq][site] << sequences[seq][site+1] << sequences[seq][site+2] <<
                        		" at site " << site+1;
                        outWarning(warn_str.str());
            		}
            		state = STATE_UNKNOWN;
            	}
            }
            if (state == STATE_INVALID) {
                if (num_error < 100) {
                    err_str << "Sequence " << seq_names[seq] << " has invalid character " << sequences[seq][site];
                    if (seq_type == SEQ_CODON)
                        err_str << sequences[seq][site+1] << sequences[seq][site+2];
                    err_str << " at site " << site+1 << endl;
                } else if (num_error == 100)
                    err_str << "...many more..." << endl;
                num_error++;
            }
            pat[seq] = state;
        }
        if (!num_error)
            num_gaps_only += addPattern(pat, site/step);
    }
    if (num_gaps_only)
        cout << "WARNING: " << num_gaps_only << " sites contain only gaps or ambiguous characters." << endl;
    if (err_str.str() != "")
        throw err_str.str();
    return 1;
}

void processSeq(string &sequence, string &line, int line_num) {
    for (string::iterator it = line.begin(); it != line.end(); it++) {
        if ((*it) <= ' ') continue;
        if (isalnum(*it) || (*it) == '-' || (*it) == '?'|| (*it) == '.' || (*it) == '*' || (*it) == '~')
            sequence.append(1, toupper(*it));
        else if (*it == '(' || *it == '{') {
            auto start_it = it;
            while (*it != ')' && *it != '}' && it != line.end())
                it++;
            if (it == line.end())
                throw "Line " + convertIntToString(line_num) + ": No matching close-bracket ) or } found";
            sequence.append(1, '?');
            cout << "NOTE: Line " << line_num << ": " << line.substr(start_it-line.begin(), (it-start_it)+1) << " is treated as unknown character" << endl;
        } else {
            throw "Line " + convertIntToString(line_num) + ": Unrecognized character "  + *it;
        }
    }
}

int Alignment::readPhylip(char *filename, char *sequence_type) {

    StrVector sequences;
    ostringstream err_str;
    igzstream in;
    int line_num = 1;
    // set the failbit and badbit
    in.exceptions(ios::failbit | ios::badbit);
    in.open(filename);
    int nseq = 0, nsite = 0;
    int seq_id = 0;
    string line;
    // remove the failbit
    in.exceptions(ios::badbit);
    bool tina_state = (sequence_type && (strcmp(sequence_type,"TINA") == 0 || strcmp(sequence_type,"MULTI") == 0));
    num_states = 0;

    for (; !in.eof(); line_num++) {
        safeGetline(in, line);
        line = line.substr(0, line.find_first_of("\n\r"));
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
            sequences.resize(nseq, "");

        } else { // read sequence contents
            if (seq_names[seq_id] == "") { // cut out the sequence name
                string::size_type pos = line.find_first_of(" \t");
                if (pos == string::npos) pos = 10; //  assume standard phylip
                seq_names[seq_id] = line.substr(0, pos);
                line.erase(0, pos);
            }
            int old_len = sequences[seq_id].length();
            if (tina_state) {
                stringstream linestr(line);
                int state;
                while (!linestr.eof() ) {
                    state = -1;
                    linestr >> state;
                    if (state < 0) break;
                    sequences[seq_id].append(1, state);
                    if (num_states < state+1) num_states = state+1;
                }
            } else processSeq(sequences[seq_id], line, line_num);
            if (sequences[seq_id].length() != sequences[0].length()) {
                err_str << "Line " << line_num << ": Sequence " << seq_names[seq_id] << " has wrong sequence length " << sequences[seq_id].length() << endl;
                throw err_str.str();
            }
            if (sequences[seq_id].length() > old_len)
                seq_id++;
            if (seq_id == nseq) {
                seq_id = 0;
                // make sure that all sequences have the same length at this moment
            }
        }
        //sequences.
    }
    in.clear();
    // set the failbit again
    in.exceptions(ios::failbit | ios::badbit);
    in.close();

    return buildPattern(sequences, sequence_type, nseq, nsite);
}

int Alignment::readPhylipSequential(char *filename, char *sequence_type) {

    StrVector sequences;
    ostringstream err_str;
    igzstream in;
    int line_num = 1;
    // set the failbit and badbit
    in.exceptions(ios::failbit | ios::badbit);
    in.open(filename);
    int nseq = 0, nsite = 0;
    int seq_id = 0;
    string line;
    // remove the failbit
    in.exceptions(ios::badbit);
    num_states = 0;

    for (; !in.eof(); line_num++) {
        safeGetline(in, line);
        line = line.substr(0, line.find_first_of("\n\r"));
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
            sequences.resize(nseq, "");

        } else { // read sequence contents
            if (seq_id >= nseq)
                throw "Line " + convertIntToString(line_num) + ": Too many sequences detected";

            if (seq_names[seq_id] == "") { // cut out the sequence name
                string::size_type pos = line.find_first_of(" \t");
                if (pos == string::npos) pos = 10; //  assume standard phylip
                seq_names[seq_id] = line.substr(0, pos);
                line.erase(0, pos);
            }
            processSeq(sequences[seq_id], line, line_num);
            if (sequences[seq_id].length() > nsite)
                throw ("Line " + convertIntToString(line_num) + ": Sequence " + seq_names[seq_id] + " is too long (" + convertIntToString(sequences[seq_id].length()) + ")");
            if (sequences[seq_id].length() == nsite) {
                seq_id++;
            }
        }
        //sequences.
    }
    in.clear();
    // set the failbit again
    in.exceptions(ios::failbit | ios::badbit);
    in.close();

    return buildPattern(sequences, sequence_type, nseq, nsite);
}

int Alignment::readFasta(char *filename, char *sequence_type) {

    StrVector sequences;
    ostringstream err_str;
    igzstream in;
    // ifstream in;
    int line_num = 1;
    string line;

    // PoMo with Fasta files is not supported yet.
    // if (sequence_type) {
    //     string st (sequence_type);
    //     if (st.substr(0,2) != "CF")
    //         throw "PoMo does not support reading fasta files yet, please use a Counts File.";
    // }

    // set the failbit and badbit
    in.exceptions(ios::failbit | ios::badbit);
    in.open(filename);
    // remove the failbit
    in.exceptions(ios::badbit);

    for (; !in.eof(); line_num++) {
        safeGetline(in, line);
        if (line == "") continue;

        //cout << line << endl;
        if (line[0] == '>') { // next sequence
            string::size_type pos = line.find_first_of("\n\r");
            seq_names.push_back(line.substr(1, pos-1));
            trimString(seq_names.back());
            sequences.push_back("");
            continue;
        }
        // read sequence contents
        if (sequences.empty()) throw "First line must begin with '>' to define sequence name";
        processSeq(sequences.back(), line, line_num);
    }
    in.clear();
    // set the failbit again
    in.exceptions(ios::failbit | ios::badbit);
    in.close();

    // now try to cut down sequence name if possible
    int i, j, step = 0;
    StrVector new_seq_names, remain_seq_names;
    new_seq_names.resize(seq_names.size());
    remain_seq_names = seq_names;

    for (step = 0; step < 4; step++) {
        bool duplicated = false;
        for (i = 0; i < seq_names.size(); i++) {
            if (remain_seq_names[i].empty()) continue;
            size_t pos = remain_seq_names[i].find_first_of(" \t");
            if (pos == string::npos) {
                new_seq_names[i] += remain_seq_names[i];
                remain_seq_names[i] = "";
            } else {
                new_seq_names[i] += remain_seq_names[i].substr(0, pos);
                remain_seq_names[i] = "_" + remain_seq_names[i].substr(pos+1);
            }
            // now check for duplication
            if (!duplicated)
            for (j = 0; j < i-1; j++)
                if (new_seq_names[j] == new_seq_names[i]) {
                    duplicated = true;
                    break;
                }
        }
        if (!duplicated) break;
    }

    if (step > 0) {
        for (i = 0; i < seq_names.size(); i++)
            if (seq_names[i] != new_seq_names[i]) {
                cout << "NOTE: Change sequence name '" << seq_names[i] << "' -> " << new_seq_names[i] << endl;
            }
    }

    seq_names = new_seq_names;

    return buildPattern(sequences, sequence_type, seq_names.size(), sequences.front().length());
}

int Alignment::readClustal(char *filename, char *sequence_type) {


    StrVector sequences;
    igzstream in;
    int line_num = 1;
    string line;
    num_states = 0;


    // set the failbit and badbit
    in.exceptions(ios::failbit | ios::badbit);
    in.open(filename);
    // remove the failbit
    in.exceptions(ios::badbit);
    safeGetline(in, line);
    if (line.substr(0, 7) != "CLUSTAL") {
        throw "ClustalW file does not start with 'CLUSTAL'";
    }

    int seq_count = 0;
    for (line_num = 2; !in.eof(); line_num++) {
        safeGetline(in, line);
        trimString(line);
        if (line == "") {
            seq_count = 0;
            continue;
        }
        if (line[0] == '*' || line[0] == ':' || line[0] == '.') continue; // ignore conservation line

        size_t pos = line.find_first_of(" \t");
        if (pos == string::npos) {
            throw "Line " + convertIntToString(line_num) + ": whitespace not found between sequence name and content";
        }
        string seq_name = line.substr(0, pos);
        if (seq_count == seq_names.size()) {
            seq_names.push_back(seq_name);
            sequences.push_back("");
        } else if (seq_count > seq_names.size()){
            throw "Line " + convertIntToString(line_num) + ": New sequence name is not allowed here";
        } else if (seq_name != seq_names[seq_count]) {
            throw "Line " + convertIntToString(line_num) + ": Sequence name " + seq_name + " does not match previously declared " +seq_names[seq_count];
        }

        line = line.substr(pos+1);
        trimString(line);
        pos = line.find_first_of(" \t");
        line = line.substr(0, pos);
        // read sequence contents
        processSeq(sequences[seq_count], line, line_num);
        seq_count++;
    }
    in.clear();
    // set the failbit again
    in.exceptions(ios::failbit | ios::badbit);
    in.close();

    if (sequences.empty())
        throw "No sequences found. Please check input (e.g. newline character)";

    return buildPattern(sequences, sequence_type, seq_names.size(), sequences.front().length());


}


int Alignment::readMSF(char *filename, char *sequence_type) {


    StrVector sequences;
    igzstream in;
    int line_num = 1;
    string line;
    num_states = 0;


    // set the failbit and badbit
    in.exceptions(ios::failbit | ios::badbit);
    in.open(filename);
    // remove the failbit
    in.exceptions(ios::badbit);
    safeGetline(in, line);
    if (line.find("MULTIPLE_ALIGNMENT") == string::npos) {
        throw "MSF file must start with header line MULTIPLE_ALIGNMENT";
    }

    int seq_len = 0, seq_count = 0;
    bool seq_started = false;

    for (line_num = 2; !in.eof(); line_num++) {
        safeGetline(in, line);
        trimString(line);
        if (line == "") {
            continue;
        }
        size_t pos;

        if (line.substr(0, 2) == "//") {
            seq_started = true;
            continue;
        }

        if (line.substr(0,5) == "Name:") {
            if (seq_started)
                throw "Line " + convertIntToString(line_num) + ": Cannot declare sequence name here";
            line = line.substr(5);
            trimString(line);
            pos = line.find_first_of(" \t");
            if (pos == string::npos)
                throw "Line " + convertIntToString(line_num) + ": No whitespace found after sequence name";
            string seq_name = line.substr(0,pos);
            seq_names.push_back(seq_name);
            sequences.push_back("");
            pos = line.find("Len:");
            if (pos == string::npos)
                throw "Line " + convertIntToString(line_num) + ": Sequence description does not contain 'Len:'";
            line = line.substr(pos+4);
            trimString(line);
            pos = line.find_first_of(" \t");
            if (pos == string::npos)
                throw "Line " + convertIntToString(line_num) + ": No whitespace found after sequence length";

            int len;
            line = line.substr(0, pos);
            try {
                len = convert_int(line.c_str());
            } catch (string &str) {
                throw "Line " + convertIntToString(line_num) + ": " + str;
            }
            if (len <= 0)
                throw "Line " + convertIntToString(line_num) + ": Non-positive sequence length not allowed";
            if (seq_len == 0)
                seq_len = len;
            else if (seq_len != len)
                throw "Line " + convertIntToString(line_num) + ": Sequence length " + convertIntToString(len) + " is different from previously defined " + convertIntToString(seq_len);
            continue;
        }

        if (!seq_started) continue;

        if (seq_names.empty())
            throw "No sequence name declared in header";

        if (isdigit(line[0])) continue;
        pos = line.find_first_of(" \t");
        if (pos == string::npos)
            throw "Line " + convertIntToString(line_num) + ": whitespace not found between sequence name and content - " + line;

        string seq_name = line.substr(0, pos);
        if (seq_name != seq_names[seq_count])
            throw "Line " + convertIntToString(line_num) + ": Sequence name " + seq_name + " does not match previously declared " +seq_names[seq_count];

        line = line.substr(pos+1);
        // read sequence contents
        processSeq(sequences[seq_count], line, line_num);
        seq_count++;
        if (seq_count == seq_names.size())
            seq_count = 0;
    }
    in.clear();
    // set the failbit again
    in.exceptions(ios::failbit | ios::badbit);
    in.close();
    return buildPattern(sequences, sequence_type, seq_names.size(), sequences.front().length());


}

// TODO: Use outWarning to print warnings.
int Alignment::readCountsFormat(char* filename, char* sequence_type) {
    int npop = 0;                // Number of populations.
    int nsites = 0;              // Number of sites.
    int N = 9;                   // Virtual population size; defaults
                                 // to 9.  If `-st CFXX` is given, it
                                 // will be set to XX below.
    int nnuc = 4;                // Number of nucleotides (base states).
    ostringstream err_str;
    igzstream in;

    // Variables to stream the data.
    string line;
    string field;
    string val_str;             // String of counts value.
    int value;                  // Actual int value.
    int line_num = 0;           // Line number counter.
    int field_num;              // Field number counter.
    int site_count = 0;         // Site / base counter.
    // Delimiters
//    char const field_delim = '\t';
    char const value_delim = ',';

    // Vector of nucleotide base counts in order A, C, G, T.
    IntVector values;
    // Sampled vector of nucleotide base counts (N individuals are
    // sampled out of =values=).
    IntVector sampled_values;
    // Iterator to loop over bases.
    IntVector::iterator i;

    // Variables to convert sampled_values to a state in the pattern.
    int sum;
    int count;
    int id1, id2;
    int r_int;
    // Index of polymorphism type; ranges from 0 to 5: [AC], [AG],
    // [AT], [CG], [CT], [GT].
    int j;
    // String with states.  Each character represents an integer state
    // value ranging from 0 to 4+(4 choose 2)*(N-1)-1.  E.g., 0 to 57
    // if N is 10.
    Pattern pattern;
    // The state a population is in at a specific site.
    // 0 ... 3 = fixed A,C,G,T
    // 4 + j*(N-2)+j ... 4 + (j+1)*(N-2)+j = polymorphism of type j
    // E.g., 4 = [1A,9C]; 5 = [2A,8C]; 12 = [9A,1C]; 13 = [1A,9G]
    int state;

    // Strings to check counts-file identification line.
    string ftype, npop_str, nsites_str;

    bool everything_ok = true;
    int fails = 0;

    // Access model_name in global parameters; needed to get N and
    // sampling method.
    Params params = Params::getInstance();
    // TODO DS: Do not temper with params; use another way to set PoMo
    // flag.
    params.pomo = true;

    // Initialize sampling method.
    pomo_sampling_method = SAMPLING_WEIGHTED_BINOM;

    // Check for custom virtual population size or sampling method.
    size_t n_pos_start = params.model_name.find("+N");
    size_t n_pos_end   = params.model_name.find_first_of("+", n_pos_start+1);
    if (n_pos_start != string::npos) {
        int length;
        if (n_pos_end != string::npos)
            length = n_pos_end - n_pos_start - 2;
        else
            length = params.model_name.length() - n_pos_start - 2;
        try {
            N = convert_int(params.model_name.substr(n_pos_start+2,length).c_str());
        }
        catch (string str) {
            cout << "The model string is faulty." << endl;
            cout << "The virtual population size N is not clear when reading in data." << endl;
            cout << "Use, e.g., \"+N7\"." << endl;
            cout << "For each run, N can only be set once." << endl;
            outError(str);
        }

        if (((N != 10) && (N != 2) && (N % 2 == 0)) || (N < 2) || (N > 19))
            outError("Custom virtual population size of PoMo not 2, 10 or any other odd number between 3 and 19.");
    }
    // TODO: probably remove virtual_pop_size and use N only.
    params.pomo_pop_size = N;
    virtual_pop_size = N;

    size_t w_pos = params.model_name.find("+WB");
    size_t h_pos = params.model_name.find("+WH");
    size_t s_pos = params.model_name.find("+S");
    int count_sampling_methods = 0;
    if (w_pos != string::npos) {
      pomo_sampling_method = SAMPLING_WEIGHTED_BINOM;
      count_sampling_methods += 1;
    }
    if (h_pos != string::npos) {
      pomo_sampling_method = SAMPLING_WEIGHTED_HYPER;
      count_sampling_methods += 1;
    }
    if (s_pos != string::npos) {
      pomo_sampling_method = SAMPLING_SAMPLED;
      count_sampling_methods += 1;
    }

    if (count_sampling_methods > 1)
      outError("Multiple sampling methods specified.");

    // Print error if sequence type is given (not supported anymore).
    if (sequence_type) {
        cout << "Counts files are auto detected." << endl;
        cout << "PoMo does not support -st flag." << endl;
        cout << "Please use model string to specifcy virtual population size and sampling method." << endl;
        outError("Abort.");
    }

    // if (sequence_type) {
    //     string st (sequence_type);
    //     if (st.substr(0,2) == "CR")
    //         pomo_random_sampling = true;
    //     else if (st.substr(0,2) == "CF")
    //         pomo_random_sampling = false;
    //     else
    //         throw "Counts File detected but sequence type (-st) is neither 'CF' nor 'CR'.";
    //     string virt_pop_size_str = st.substr(2);
    //     if (virt_pop_size_str != "") {
    //         int virt_pop_size = atoi(virt_pop_size_str.c_str());
    //         N = virt_pop_size;
    //     }
    // }

    // Set the number of states.  If nnuc=4:
    // 4 + (4 choose 2)*(N-1) = 58.
    num_states = nnuc + nnuc*(nnuc-1)/2*(N-1);
    seq_type = SEQ_POMO;

    // Set UNKNOWN_STATE.  This state is set if no information is in
    // the alignment.  If we use partial likelihood we do not know the
    // number of different patterns in the alignment yet and hence,
    // cannot set the variable STATE_UNKNOWN yet (see
    // `state_unknown_buffer`).
    computeUnknownState();

    // Use a buffer for STATE_UNKNOWN.  I.e., if an unknown state is
    // encountered, the pattern is added to this buffer.  Only after
    // all sites have been read in, the patterns from this temporal
    // buffer are added to the normal alignment because then, the
    // value of STATE_UNKNOWN is known.
    vector<Pattern> su_buffer;
    // The site numbers of the patterns that include unknown states.
    IntVector su_site_counts;
    int su_number = 0;

    // BQM: not neccessary, su_site_count will be equal to su_site_counts.size()
    //    int su_site_count = 0;
    bool includes_state_unknown = false;

    // Variables to calculate mean number of samples per population.
    // If N is way above the average number of samples, PoMo has been
    // ovserved to be unstable and a big warning is printed.
    int n_samples_sum = 0;
    int n_sites_sum = 0;
    // Average number of samples.
    double n_samples_bar = 0;

    // Open counts file.
    // Set the failbit and badbit.
    in.exceptions(ios::failbit | ios::badbit);
    in.open(filename);
    // Remove the failbit.
    in.exceptions(ios::badbit);

    // Skip comments.
    do {
        getline(in, line);
        line_num++;
    }
    while (line[0] == '#');

    // Read in npop and nsites;
    istringstream ss1(line);
    // Read and check counts file headerline.
    if (!(ss1 >> ftype >> npop_str >> npop >> nsites_str >> nsites)) {
        err_str << "Counts-File identification line could not be read.";
        throw err_str.str();
    }
    if ((ftype.compare("COUNTSFILE") ||
         npop_str.compare("NPOP") ||
         nsites_str.compare("NSITES")) != 0) {
        err_str << "Counts-File identification line could not be read.";
        throw err_str.str();
    }
    cout << endl;
    cout << "----------------------------------------------------------------------" << endl;
    cout << "Number of populations:     " << npop << endl;
    cout << "Number of sites:           " << nsites << endl;

    if (nsites > 0)
        site_pattern.resize(nsites);
    else {
        err_str << "Number of sites is 0.";
        throw err_str.str();
    }

    // Skip comments.
    do {
        getline(in, line);
        line_num++;
    }
    while (line[0] == '#');

    // Headerline.
    istringstream ss2(line);

    for (field_num = 0; (ss2 >> field); field_num++) {
        if (field_num == 0) {
            if ((field.compare("Chrom") != 0) && (field.compare("CHROM") != 0)) {
                err_str << "Unrecognized header field " << field << ".";
                throw err_str.str();
            }
        }
        else if (field_num == 1) {
            if ((field.compare("Pos") != 0) && (field.compare("POS") != 0)) {
                err_str << "Unrecognized header field " << field << ".";
                throw err_str.str();
            }
        }
        else {
            //Read in sequence names.
            seq_names.push_back(field);
        }
    }
    if ((int) seq_names.size() != npop) {
                err_str << "Number of populations in headerline doesn't match NPOP.";
                throw err_str.str();
    }

    // Data.
    // Loop over sites.
    for ( ; getline(in, line); ) {
        line_num++;
    	field_num = 0;
        pattern.clear();
        everything_ok = true;
        includes_state_unknown = false;
        istringstream fieldstream(line);
        // Loop over populations / individuals.
        for ( ; (fieldstream >> field); ) {
            // Skip Chrom and Pos columns.
            if ( (field_num == 0) || (field_num == 1)) {
                field_num++;
                continue;
            }
            // Clear value vectors.
            values.clear();
            sampled_values.clear();
            sampled_values.resize(nnuc,0);
            istringstream valuestream(field);
            // Loop over bases within one population.
            for (; getline(valuestream, val_str, value_delim);) {
            	try {
            		value = convert_int(val_str.c_str());
            	} catch(string &str) {
            		err_str << "Could not read value " << val_str << " on line " << line_num << ".";
            		throw err_str.str();
            	}
            	values.push_back(value);
            }
            if (values.size() != nnuc) {
                err_str << "Number of bases does not match on line " << line_num << ".";
                throw err_str.str();
            }

            // Read in the data.
            sum = 0;
            count = 0;
            id1 = -1;
            id2 = -1;
            // Sum over elements and count non-zero elements.
            for(i = values.begin(); i != values.end(); ++i) {
                // `i` is an iterator object that points to some
                // element of `value`.
            	if (*i != 0) {
                    // `i - values.begin()` ranges from 0 to 3 and
                    // determines the nucleotide or allele type.
            		if (id1 == -1) id1 = i - values.begin();
            		else id2 = i - values.begin();
            		count++;
                	sum += *i;
            	}
            }
            // Determine state (cf. above).
            if (count == 1) {
                n_samples_sum += values[id1];
                n_sites_sum++;
                if (pomo_sampling_method == SAMPLING_SAMPLED) {
                    // Fixed state, state ID is just id1.
                    state = id1;
                } else {
                    if (values[id1] >= 16384) {
                        cout << "WARNING: Pattern on line " <<
                            line_num << " exceeds count limit of 16384." << endl;
                        everything_ok = false;
                    }
                    uint32_t pomo_state = (id1 | (values[id1]) << 2);
                    IntIntMap::iterator pit = pomo_sampled_states_index.find(pomo_state);
                    if (pit == pomo_sampled_states_index.end()) { // not found
                        state = pomo_sampled_states_index[pomo_state] = pomo_sampled_states.size();
                        pomo_sampled_states.push_back(pomo_state);
                    } else {
                        state = pit->second;
                    }
                    state += num_states; // make the state larger than num_states
                }
            }
            else if (count == 0) {
                state = STATE_UNKNOWN;
                su_number++;
                includes_state_unknown = true;
            }
            else if (count > 2) {
                if (verbose_mode >= VB_MAX) {
                    std::cout << "WARNING: More than two bases are present on line ";
                    std::cout << line_num << "." << std::endl;
                }
                everything_ok = false;
            	// err_str << "More than 2 bases are present on line " << line_num << ".";
            	// throw err_str.str();
            }
            // Now we deal with the important polymorphic states with two alleles.
            else if (count == 2) {
                n_samples_sum += values[id1];
                n_samples_sum += values[id2];
                n_sites_sum++;
                if (pomo_sampling_method == SAMPLING_SAMPLED) {
                     // Binomial sampling.  2 bases are present.
                    for(int k = 0; k < N; k++) {
                        r_int = random_int(sum);
                        if (r_int < values[id1]) sampled_values[id1]++;
                        else sampled_values[id2]++;
                    }
                    if (sampled_values[id1] == 0) state = id2;
                    else if (sampled_values[id2] == 0) state = id1;
                    else {
                        if (id1 == 0) j = id2 - 1;
                        else j = id1 + id2;
                        state = nnuc + j*(N-2) + j + sampled_values[id1] - 1;
                    }
                } else {
                    /* BQM 2015-07: store both states now */
                    if (values[id1] >= 16384 || values[id2] >= 16384)
                        // Cannot add sites where more than 16384
                        // individuals have the same base within one
                        // population.
                        everything_ok = false;
                    uint32_t pomo_state = (id1 | (values[id1]) << 2) | ((id2 | (values[id2]<<2))<<16);
                    IntIntMap::iterator pit = pomo_sampled_states_index.find(pomo_state);
                    if (pit == pomo_sampled_states_index.end()) { // not found
                        state = pomo_sampled_states_index[pomo_state] = pomo_sampled_states.size();
                        pomo_sampled_states.push_back(pomo_state);
                    } else {
                        state = pit->second;
                    }
                    state += num_states; // make the state larger than num_states
                }
            }
            else {
                err_str << "Unexpected error on line number " << line_num << ".";
                throw err_str.str();
            }

            // Now we have the state to build a pattern ;-).
            pattern.push_back(state);
        }
        if ((int) pattern.size() != npop) {
            err_str << "Number of species does not match on line " << line_num << ".";
            throw err_str.str();
        }
        // Pattern has been built and is now added to the vector of
        // patterns.
        if (everything_ok == true) {
            if (includes_state_unknown) {
//                su_site_count++;
                if (pomo_sampling_method == SAMPLING_WEIGHTED_BINOM ||
                    pomo_sampling_method == SAMPLING_WEIGHTED_HYPER) {
                    su_buffer.push_back(pattern);
                    su_site_counts.push_back(site_count);
                }
                // Add pattern if we use random sampling because then,
                // STATE_UNKNOWN = num_states is well defined already at
                // this stage.
                else
                    addPattern(pattern, site_count);

                // BQM: it is neccessary to always increase site_count
                site_count++;
            }
            else {
                addPattern(pattern, site_count);
                site_count++;
            }
        }
        else {
            fails++;
            if (verbose_mode >= VB_MAX) {
                cout << "WARNING: Pattern on line " <<
                    line_num << " was not added." << endl;
            }
        }
    }

    if (site_count + fails != nsites) {
        err_str << "Number of sites does not match NSITES.";
        throw err_str.str();
    }

    if (pomo_sampling_method == SAMPLING_WEIGHTED_BINOM ||
        pomo_sampling_method == SAMPLING_WEIGHTED_HYPER) {
        // Now we can correctly set STATE_UNKNOWN.
        STATE_UNKNOWN = pomo_sampled_states.size() + num_states;

        // Process sites that include an unknown state.
        for (vector<Pattern>::iterator pat_it = su_buffer.begin();
             pat_it != su_buffer.end(); pat_it++) {
            for (Pattern::iterator sp_it = pat_it->begin(); sp_it != pat_it->end(); sp_it++)
                if (*sp_it == 0xffffffff) *sp_it = STATE_UNKNOWN;
        }

        for (unsigned int i = 0; i < su_buffer.size(); i++)
                addPattern(su_buffer[i], su_site_counts[i]);
    }

    cout << "---" << endl;
    cout << "Normal sites:              " << site_count - su_site_counts.size() << endl;
    cout << "Sites with unknown states: " << su_site_counts.size() << endl;
    cout << "Total sites read:          " << site_count << endl;
    cout << "Fails:                     " << fails << endl;
    if (pomo_sampling_method == SAMPLING_WEIGHTED_BINOM ||
        pomo_sampling_method == SAMPLING_WEIGHTED_HYPER) {
        cout << "---" << endl;
        cout << "Compound states:           " << pomo_sampled_states.size() << endl;
    }
    cout << "----------------------------------------------------------------------" << endl << endl;

    // Check if N is not too large.
    n_samples_bar = n_samples_sum / (double) n_sites_sum;
    cout << "The average number of samples is " << n_samples_bar << endl;
    if ((pomo_sampling_method == SAMPLING_WEIGHTED_BINOM) &&
        (n_samples_bar * 3.0 <= N)) {
        cout << "----------------------------------------------------------------------" << endl;
        cout << "WARNING: The virtual population size N is much larger ";
        cout << "than the average number of samples." << endl;
        cout << "WARNING: This setting together with /weighted binomial/ sampling ";
        cout << "may be numerically unstable." << endl << endl;
        cout << "----------------------------------------------------------------------" << endl;
    }

    site_pattern.resize(site_count);

    in.clear();
    // set the failbit again
    in.exceptions(ios::failbit | ios::badbit);
    in.close();

    return 1;
}

bool Alignment::getSiteFromResidue(int seq_id, int &residue_left, int &residue_right) {
    int i, j;
    int site_left = -1, site_right = -1;
    for (i = 0, j = -1; i < getNSite(); i++) {
        if (at(site_pattern[i])[seq_id] != STATE_UNKNOWN) j++;
        if (j == residue_left) site_left = i;
        if (j == residue_right-1) site_right = i+1;
    }
    if (site_left < 0 || site_right < 0)
        cout << "Out of range: Maxmimal residue number is " << j+1 << endl;
    if (site_left == -1) outError("Left residue range is too high");
    if (site_right == -1) {
        outWarning("Right residue range is set to alignment length");
        site_right = getNSite();
    }
    residue_left = site_left;
    residue_right = site_right;
    return true;
}

int Alignment::buildRetainingSites(const char *aln_site_list, IntVector &kept_sites,
		bool exclude_gaps, bool exclude_const_sites, const char *ref_seq_name)
{
    if (aln_site_list) {
        int seq_id = -1;
        if (ref_seq_name) {
            string ref_seq = ref_seq_name;
            seq_id = getSeqID(ref_seq);
            if (seq_id < 0) outError("Reference sequence name not found: ", ref_seq_name);
        }
        cout << "Reading site position list " << aln_site_list << " ..." << endl;
        kept_sites.resize(getNSite(), 0);
        try {
            ifstream in;
            in.exceptions(ios::failbit | ios::badbit);
            in.open(aln_site_list);
            in.exceptions(ios::badbit);

            while (!in.eof()) {
                int left, right;
                left = right = 0;
                in >> left;
                if (in.eof()) break;
                in >> right;
                cout << left << "-" << right << endl;
                if (left <= 0 || right <= 0) throw "Range must be positive";
                if (left > right) throw "Left range is bigger than right range";
                left--;
                if (right > getNSite()) throw "Right range is bigger than alignment size";
                if (seq_id >= 0) getSiteFromResidue(seq_id, left, right);
                for (int i = left; i < right; i++)
                    kept_sites[i] = 1;
            }
            in.close();
        } catch (ios::failure) {
            outError(ERR_READ_INPUT, aln_site_list);
        } catch (const char* str) {
            outError(str);
        }
    } else {
        kept_sites.resize(getNSite(), 1);
    }

    int j;
    if (exclude_gaps) {
        for (j = 0; j < kept_sites.size(); j++)
            if (kept_sites[j] && at(site_pattern[j]).computeAmbiguousChar(num_states) > 0) {
                kept_sites[j] = 0;
            }
    }
    if (exclude_const_sites) {
        for (j = 0; j < kept_sites.size(); j++)
        	if (at(site_pattern[j]).isInvariant())
        		kept_sites[j] = 0;

    }

    int final_length = 0;
    for (j = 0; j < kept_sites.size(); j++)
        if (kept_sites[j]) final_length++;
    return final_length;
}

void Alignment::printPhylip(ostream &out, bool append, const char *aln_site_list,
                            bool exclude_gaps, bool exclude_const_sites, const char *ref_seq_name, bool print_taxid) {
    IntVector kept_sites;
    int final_length = buildRetainingSites(aln_site_list, kept_sites, exclude_gaps, exclude_const_sites, ref_seq_name);
    if (seq_type == SEQ_CODON)
        final_length *= 3;

	out << getNSeq() << " " << final_length << endl;
	int max_len = getMaxSeqNameLength();
    if (print_taxid) max_len = 10;
	if (max_len < 10) max_len = 10;
	int seq_id;
	for (seq_id = 0; seq_id < seq_names.size(); seq_id++) {
		out.width(max_len);
        if (print_taxid)
            out << left << seq_id << " ";
        else
            out << left << seq_names[seq_id] << " ";
		int j = 0;
		for (IntVector::iterator i = site_pattern.begin();  i != site_pattern.end(); i++, j++)
			if (kept_sites[j])
				out << convertStateBackStr(at(*i)[seq_id]);
		out << endl;
	}
}

void Alignment::printPhylip(const char *file_name, bool append, const char *aln_site_list,
                            bool exclude_gaps, bool exclude_const_sites, const char *ref_seq_name) {
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);

        if (append)
            out.open(file_name, ios_base::out | ios_base::app);
        else
            out.open(file_name);

        printPhylip(out, append, aln_site_list, exclude_gaps, exclude_const_sites, ref_seq_name);

        out.close();
        if (verbose_mode >= VB_MED)
        	cout << "Alignment was printed to " << file_name << endl;
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, file_name);
    }
}

void Alignment::printFasta(const char *file_name, bool append, const char *aln_site_list
                           , bool exclude_gaps, bool exclude_const_sites, const char *ref_seq_name)
{
    IntVector kept_sites;
    buildRetainingSites(aln_site_list, kept_sites, exclude_gaps, exclude_const_sites, ref_seq_name);
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        if (append)
            out.open(file_name, ios_base::out | ios_base::app);
        else
            out.open(file_name);
        StrVector::iterator it;
        int seq_id = 0;
        for (it = seq_names.begin(); it != seq_names.end(); it++, seq_id++) {
            out << ">" << (*it) << endl;
            int j = 0;
            for (IntVector::iterator i = site_pattern.begin();  i != site_pattern.end(); i++, j++)
                if (kept_sites[j])
                    out << convertStateBackStr(at(*i)[seq_id]);
            out << endl;
        }
        out.close();
        cout << "Alignment was printed to " << file_name << endl;
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, file_name);
    }
}


void Alignment::extractSubAlignment(Alignment *aln, IntVector &seq_id, int min_true_char, int min_taxa, IntVector *kept_partitions) {
    IntVector::iterator it;
    for (it = seq_id.begin(); it != seq_id.end(); it++) {
        ASSERT(*it >= 0 && *it < aln->getNSeq());
        seq_names.push_back(aln->getSeqName(*it));
    }
    name = aln->name;
    model_name = aln->model_name;
    sequence_type = aln->sequence_type;
    position_spec = aln->position_spec;
    aln_file = aln->aln_file;
    num_states = aln->num_states;
    seq_type = aln->seq_type;
    STATE_UNKNOWN = aln->STATE_UNKNOWN;
	genetic_code = aln->genetic_code;
    if (seq_type == SEQ_CODON) {
    	codon_table = new char[num_states];
    	memcpy(codon_table, aln->codon_table, num_states);
    	non_stop_codon = new char[strlen(genetic_code)];
    	memcpy(non_stop_codon, aln->non_stop_codon, strlen(genetic_code));
    }
    site_pattern.resize(aln->getNSite(), -1);
    clear();
    pattern_index.clear();
    int site = 0, removed_sites = 0;
    VerboseMode save_mode = verbose_mode;
    verbose_mode = min(verbose_mode, VB_MIN); // to avoid printing gappy sites in addPattern
//    for (iterator pit = aln->begin(); pit != aln->end(); pit++) {
    for (site = 0; site < aln->getNSite(); site++) {
        iterator pit = aln->begin() + (aln->getPatternID(site));
        Pattern pat;
        int true_char = 0;
        for (it = seq_id.begin(); it != seq_id.end(); it++) {
            char ch = (*pit)[*it];
            if (ch != STATE_UNKNOWN) true_char++;
            pat.push_back(ch);
        }
        if (true_char < min_true_char)
            removed_sites++;
        else
            addPattern(pat, site-removed_sites);
//        for (int i = 0; i < (*pit).frequency; i++)
//            site_pattern[site++] = size()-1;
    }
    site_pattern.resize(aln->getNSite() - removed_sites);
    verbose_mode = save_mode;
    countConstSite();
    buildSeqStates();
    ASSERT(size() <= aln->size());
    if (kept_partitions)
        kept_partitions->push_back(0);
}


void Alignment::extractPatterns(Alignment *aln, IntVector &ptn_id) {
    int i;
    for (i = 0; i < aln->getNSeq(); i++) {
        seq_names.push_back(aln->getSeqName(i));
    }
    name = aln->name;
    model_name = aln->model_name;
    sequence_type = aln->sequence_type;
    position_spec = aln->position_spec;
    aln_file = aln->aln_file;
    num_states = aln->num_states;
    seq_type = aln->seq_type;
    STATE_UNKNOWN = aln->STATE_UNKNOWN;
    genetic_code = aln->genetic_code;
    if (seq_type == SEQ_CODON) {
    	codon_table = new char[num_states];
    	memcpy(codon_table, aln->codon_table, num_states);
    	non_stop_codon = new char[strlen(genetic_code)];
    	memcpy(non_stop_codon, aln->non_stop_codon, strlen(genetic_code));
    }
    site_pattern.resize(aln->getNSite(), -1);
    clear();
    pattern_index.clear();
    int site = 0;
    VerboseMode save_mode = verbose_mode;
    verbose_mode = min(verbose_mode, VB_MIN); // to avoid printing gappy sites in addPattern
    for (i = 0; i != ptn_id.size(); i++) {
        ASSERT(ptn_id[i] >= 0 && ptn_id[i] < aln->getNPattern());
        Pattern pat = aln->at(ptn_id[i]);
        addPattern(pat, site, aln->at(ptn_id[i]).frequency);
        for (int j = 0; j < aln->at(ptn_id[i]).frequency; j++)
            site_pattern[site++] = size()-1;
    }
    site_pattern.resize(site);
    verbose_mode = save_mode;
    countConstSite();
    buildSeqStates();
    ASSERT(size() <= aln->size());
}

void Alignment::extractPatternFreqs(Alignment *aln, IntVector &ptn_freq) {
    int i;
    ASSERT(ptn_freq.size() <= aln->getNPattern());
    for (i = 0; i < aln->getNSeq(); i++) {
        seq_names.push_back(aln->getSeqName(i));
    }
    name = aln->name;
    model_name = aln->model_name;
    sequence_type = aln->sequence_type;
    position_spec = aln->position_spec;
    aln_file = aln->aln_file;
    num_states = aln->num_states;
    seq_type = aln->seq_type;
    genetic_code = aln->genetic_code;
    if (seq_type == SEQ_CODON) {
    	codon_table = new char[num_states];
    	memcpy(codon_table, aln->codon_table, num_states);
    	non_stop_codon = new char[strlen(genetic_code)];
    	memcpy(non_stop_codon, aln->non_stop_codon, strlen(genetic_code));
    }
    STATE_UNKNOWN = aln->STATE_UNKNOWN;
    site_pattern.resize(accumulate(ptn_freq.begin(), ptn_freq.end(), 0), -1);
    clear();
    pattern_index.clear();
    int site = 0;
    VerboseMode save_mode = verbose_mode;
    verbose_mode = min(verbose_mode, VB_MIN); // to avoid printing gappy sites in addPattern
    for (i = 0; i != ptn_freq.size(); i++)
        if (ptn_freq[i]) {
            ASSERT(ptn_freq[i] > 0);
            Pattern pat = aln->at(i);
            addPattern(pat, site, ptn_freq[i]);
            for (int j = 0; j < ptn_freq[i]; j++)
                site_pattern[site++] = size()-1;
        }
    site_pattern.resize(site);
    verbose_mode = save_mode;
    countConstSite();
    buildSeqStates();
    ASSERT(size() <= aln->size());
}

void Alignment::extractSites(Alignment *aln, IntVector &site_id) {
    int i;
    for (i = 0; i < aln->getNSeq(); i++) {
        seq_names.push_back(aln->getSeqName(i));
    }
    name = aln->name;
    model_name = aln->model_name;
    sequence_type = aln->sequence_type;
    position_spec = aln->position_spec;
    aln_file = aln->aln_file;
    num_states = aln->num_states;
    seq_type = aln->seq_type;
    STATE_UNKNOWN = aln->STATE_UNKNOWN;
    genetic_code = aln->genetic_code;
    if (seq_type == SEQ_CODON) {
    	codon_table = new char[num_states];
    	memcpy(codon_table, aln->codon_table, num_states);
    	non_stop_codon = new char[strlen(genetic_code)];
    	memcpy(non_stop_codon, aln->non_stop_codon, strlen(genetic_code));
    }
    site_pattern.resize(site_id.size(), -1);
    clear();
    pattern_index.clear();
    VerboseMode save_mode = verbose_mode;
    verbose_mode = min(verbose_mode, VB_MIN); // to avoid printing gappy sites in addPattern
    for (i = 0; i != site_id.size(); i++) {
        Pattern pat = aln->getPattern(site_id[i]);
        addPattern(pat, i);
    }
    verbose_mode = save_mode;
    countConstSite();
    buildSeqStates();
    // sanity check
    for (iterator it = begin(); it != end(); it++)
    	if (it->at(0) == -1)
    		ASSERT(0);

    //cout << getNSite() << " positions were extracted" << endl;
    //cout << __func__ << " " << num_states << endl;
}



void Alignment::convertToCodonOrAA(Alignment *aln, char *gene_code_id, bool nt2aa) {
    if (aln->seq_type != SEQ_DNA)
        outError("Cannot convert non-DNA alignment into codon alignment");
    if (aln->getNSite() % 3 != 0)
        outError("Sequence length is not divisible by 3 when converting to codon sequences");
    int i, site;
    char AA_to_state[NUM_CHAR];
    for (i = 0; i < aln->getNSeq(); i++) {
        seq_names.push_back(aln->getSeqName(i));
    }
    name = aln->name;
    model_name = aln->model_name;
    sequence_type = aln->sequence_type;
    position_spec = aln->position_spec;
    aln_file = aln->aln_file;
//    num_states = aln->num_states;
    seq_type = SEQ_CODON;
    initCodon(gene_code_id);
    if (nt2aa) {
        seq_type = SEQ_PROTEIN;
        num_states = 20;
    }

    computeUnknownState();

    if (nt2aa) {
        buildStateMap(AA_to_state, SEQ_PROTEIN);
    }

    site_pattern.resize(aln->getNSite()/3, -1);
    clear();
    pattern_index.clear();
    int step = ((seq_type == SEQ_CODON || nt2aa) ? 3 : 1);

    VerboseMode save_mode = verbose_mode;
    verbose_mode = min(verbose_mode, VB_MIN); // to avoid printing gappy sites in addPattern
    int nsite = aln->getNSite();
    int nseq = aln->getNSeq();
    Pattern pat;
    pat.resize(nseq);
    int num_error = 0;
    ostringstream err_str;

    for (site = 0; site < nsite; site+=step) {
        for (int seq = 0; seq < nseq; seq++) {
            //char state = convertState(sequences[seq][site], seq_type);
            char state = aln->at(aln->getPatternID(site))[seq];
            // special treatment for codon
            char state2 = aln->at(aln->getPatternID(site+1))[seq];
            char state3 = aln->at(aln->getPatternID(site+2))[seq];
            if (state < 4 && state2 < 4 && state3 < 4) {
//            		state = non_stop_codon[state*16 + state2*4 + state3];
                state = state*16 + state2*4 + state3;
                if (genetic_code[(int)state] == '*') {
                    err_str << "Sequence " << seq_names[seq] << " has stop codon "
                            << " at site " << site+1 << endl;
                    num_error++;
                    state = STATE_UNKNOWN;
                } else if (nt2aa) {
                    state = AA_to_state[(int)genetic_code[(int)state]];
                } else {
                    state = non_stop_codon[(int)state];
                }
            } else if (state == STATE_INVALID || state2 == STATE_INVALID || state3 == STATE_INVALID) {
                state = STATE_INVALID;
            } else {
                if (state != STATE_UNKNOWN || state2 != STATE_UNKNOWN || state3 != STATE_UNKNOWN) {
                    ostringstream warn_str;
                    warn_str << "Sequence " << seq_names[seq] << " has ambiguous character " <<
                        " at site " << site+1;
                    outWarning(warn_str.str());
                }
                state = STATE_UNKNOWN;
            }
            if (state == STATE_INVALID) {
                if (num_error < 100) {
                    err_str << "Sequence " << seq_names[seq] << " has invalid character ";
                    err_str << " at site " << site+1 << endl;
                } else if (num_error == 100)
                    err_str << "...many more..." << endl;
                num_error++;
            }
            pat[seq] = state;
        }
        if (!num_error)
            addPattern(pat, site/step);
    }
    if (num_error)
        outError(err_str.str());
    verbose_mode = save_mode;
    countConstSite();
    buildSeqStates();
    // sanity check
    for (iterator it = begin(); it != end(); it++)
    	if (it->at(0) == -1)
    		ASSERT(0);

}

void convert_range(const char *str, int &lower, int &upper, int &step_size, char* &endptr) throw (string) {
    //char *endptr;
    char *beginptr = (char*) str;

    // parse the lower bound of the range
    int d = strtol(str, &endptr, 10);
    if ((d == 0 && endptr == str) || abs(d) == HUGE_VALL) {
        string err = "Expecting integer, but found \"";
        err += str;
        err += "\" instead";
        throw err;
    }
    lower = d;
    //int d_save = d;
    upper = d;
    step_size = 1;
    // skip blank chars
    for (; *endptr == ' '; endptr++) {}
    if (*endptr != '-') return;

    // parse the upper bound of the range
    endptr++;
    // skip blank chars
    for (; *endptr == ' '; endptr++) {}
    str = endptr;
    d = strtol(str, &endptr, 10);
    if ((d == 0 && endptr == str) || abs(d) == HUGE_VALL) {
        if (str[0] == '.') {
            // 2019-06-03: special character '.' for whatever ending position
            d = lower-1;
            endptr++;
        } else {
            string err = "Expecting integer, but found \"";
            err += str;
            err += "\" instead";
            throw err;
        }
    }

    //lower = d_save;
    upper = d;
    // skip blank chars
    for (; *endptr == ' '; endptr++) {}

    if (*endptr != '\\') return;

    // parse the step size of the range
    str = endptr+1;
    d = strtol(str, &endptr, 10);
    if ((d == 0 && endptr == str) || abs(d) == HUGE_VALL) {
        string err = "Expecting integer, but found \"";
        err += str;
        err += "\" instead";
        throw err;
    }

    step_size = d;
    str = beginptr;

}

void extractSiteID(Alignment *aln, const char* spec, IntVector &site_id) {
    int i;
    char *str = (char*)spec;
    int nchars = 0;
    try {
        for (; *str != 0; ) {
            int lower, upper, step;
            convert_range(str, lower, upper, step, str);
            // 2019-06-03: special '.' character
            if (upper == lower-1)
                upper = aln->getNSite();
            lower--;
            upper--;
            nchars += (upper-lower+1)/step;
            if (aln->seq_type == SEQ_CODON) {
                lower /= 3;
                upper /= 3;
            }
            if (upper >= aln->getNSite()) throw "Too large site ID";
            if (lower < 0) throw "Negative site ID";
            if (lower > upper) throw "Wrong range";
            if (step < 1) throw "Wrong step size";
            for (i = lower; i <= upper; i+=step)
                site_id.push_back(i);
            if (*str == ',' || *str == ' ') str++;
            //else break;
        }
        if (aln->seq_type == SEQ_CODON && nchars % 3 != 0)
            throw (string)"Range " + spec + " length is not multiple of 3 (necessary for codon data)";
    } catch (const char* err) {
        outError(err);
    } catch (string err) {
        outError(err);
    }
}

void Alignment::extractSites(Alignment *aln, const char* spec) {
    IntVector site_id;
    extractSiteID(aln, spec, site_id);
    extractSites(aln, site_id);
}

void Alignment::createBootstrapAlignment(Alignment *aln, IntVector* pattern_freq, const char *spec) {
    if (aln->isSuperAlignment()) outError("Internal error: ", __func__);
    name = aln->name;
    model_name = aln->model_name;
    sequence_type = aln->sequence_type;
    position_spec = aln->position_spec;
    aln_file = aln->aln_file;
    int site, nsite = aln->getNSite();
    seq_names.insert(seq_names.begin(), aln->seq_names.begin(), aln->seq_names.end());
    num_states = aln->num_states;
    seq_type = aln->seq_type;
    genetic_code = aln->genetic_code;
    if (seq_type == SEQ_CODON) {
    	codon_table = new char[num_states];
    	memcpy(codon_table, aln->codon_table, num_states);
    	non_stop_codon = new char[strlen(genetic_code)];
    	memcpy(non_stop_codon, aln->non_stop_codon, strlen(genetic_code));
    }
    STATE_UNKNOWN = aln->STATE_UNKNOWN;
    site_pattern.resize(nsite, -1);
    clear();
    pattern_index.clear();

    // 2016-07-05: copy variables for PoMo
    pomo_sampled_states = aln->pomo_sampled_states;
    pomo_sampled_states_index = aln->pomo_sampled_states_index;
    pomo_sampling_method = aln->pomo_sampling_method;
    virtual_pop_size = aln->virtual_pop_size;

    VerboseMode save_mode = verbose_mode;
    verbose_mode = min(verbose_mode, VB_MIN); // to avoid printing gappy sites in addPattern
    if (pattern_freq) {
        pattern_freq->resize(0);
        pattern_freq->resize(aln->getNPattern(), 0);
    }

    if (!aln->site_state_freq.empty()) {
        // resampling also the per-site state frequency vector
        if (aln->site_state_freq.size() != aln->getNPattern() || spec)
            outError("Unsupported bootstrap feature, pls contact the developers");
    }
    
    if (Params::getInstance().jackknife_prop > 0.0 && spec) {
        outError((string)"Unsupported jackknife with sampling " + spec);
    }

	IntVector site_vec;
    if (!spec) {
		// standard bootstrap
        int added_sites = 0;
        IntVector sample;
        random_resampling(nsite, sample);
		for (site = 0; site < nsite; site++)
        for (int rep = 0; rep < sample[site]; rep++) {
			int ptn_id = aln->getPatternID(site);
			Pattern pat = aln->at(ptn_id);
            int nptn = getNPattern();
			addPattern(pat, added_sites);
            if (!aln->site_state_freq.empty() && getNPattern() > nptn) {
                // a new pattern is added, copy state frequency vector
                double *state_freq = new double[num_states];
                memcpy(state_freq, aln->site_state_freq[ptn_id], num_states*sizeof(double));
                site_state_freq.push_back(state_freq);
            }
			if (pattern_freq) ((*pattern_freq)[ptn_id])++;
            added_sites++;
		}
        if (added_sites < nsite)
            site_pattern.resize(added_sites);
    } else if (strncmp(spec, "GENESITE,", 9) == 0) {
		// resampling genes, then resampling sites within resampled genes
		convert_int_vec(spec+9, site_vec);
		int i;
		IntVector begin_site;
		for (i = 0, site = 0; i < site_vec.size(); i++) {
			begin_site.push_back(site);
			site += site_vec[i];
			//cout << "site = " << site_vec[i] << endl;
		}
		if (site > getNSite())
			outError("Sum of lengths exceeded alignment length");

		for (i = 0; i < site_vec.size(); i++) {
			int part = random_int(site_vec.size());
			for (int j = 0; j < site_vec[part]; j++) {
				site = random_int(site_vec[part]) + begin_site[part];
				int ptn = aln->getPatternID(site);
				Pattern pat = aln->at(ptn);
				addPattern(pat, site);
				if (pattern_freq) ((*pattern_freq)[ptn])++;
			}
		}
    } else if (strncmp(spec, "GENE,", 5) == 0) {
		// resampling genes instead of sites
		convert_int_vec(spec+5, site_vec);
		int i;
		IntVector begin_site;
		for (i = 0, site = 0; i < site_vec.size(); i++) {
			begin_site.push_back(site);
			site += site_vec[i];
			//cout << "site = " << site_vec[i] << endl;
		}
		if (site > getNSite())
			outError("Sum of lengths exceeded alignment length");

		for (i = 0; i < site_vec.size(); i++) {
			int part = random_int(site_vec.size());
			for (site = begin_site[part]; site < begin_site[part] + site_vec[part]; site++) {
				int ptn = aln->getPatternID(site);
				Pattern pat = aln->at(ptn);
				addPattern(pat, site);
				if (pattern_freq) ((*pattern_freq)[ptn])++;
			}
		}
    } else {
    	// special bootstrap
    	convert_int_vec(spec, site_vec);
    	if (site_vec.size() % 2 != 0)
    		outError("Bootstrap specification length is not divisible by 2");
    	nsite = 0;
    	int part, begin_site = 0, out_site = 0;
    	for (part = 0; part < site_vec.size(); part+=2)
    		nsite += site_vec[part+1];
    	site_pattern.resize(nsite, -1);
    	for (part = 0; part < site_vec.size(); part += 2) {
    		if (begin_site + site_vec[part] > aln->getNSite())
    			outError("Sum of lengths exceeded alignment length");
    		for (site = 0; site < site_vec[part+1]; site++) {
    			int site_id = random_int(site_vec[part]) + begin_site;
    			int ptn_id = aln->getPatternID(site_id);
    			Pattern pat = aln->at(ptn_id);
    			addPattern(pat, site + out_site);
    			if (pattern_freq) ((*pattern_freq)[ptn_id])++;
    		}
    		begin_site += site_vec[part];
    		out_site += site_vec[part+1];
    	}
    }
    if (!aln->site_state_freq.empty()) {
        site_model = site_pattern;
        ASSERT(site_state_freq.size() == getNPattern());
    }
    verbose_mode = save_mode;
    countConstSite();
    buildSeqStates();
}

void Alignment::createBootstrapAlignment(IntVector &pattern_freq, const char *spec) {
	int nptn = getNPattern();
    pattern_freq.resize(nptn, 0);
    int *internal_freq = new int [nptn];
    createBootstrapAlignment(internal_freq, spec);
    for (int i = 0; i < nptn; i++)
    	pattern_freq[i] = internal_freq[i];
    delete [] internal_freq;
}

void Alignment::createBootstrapAlignment(int *pattern_freq, const char *spec, int *rstream) {
    int site, nsite = getNSite();
    memset(pattern_freq, 0, getNPattern()*sizeof(int));
	IntVector site_vec;
    if (Params::getInstance().jackknife_prop > 0.0 && spec)
        outError((string)"Unsupported jackknife with " + spec);

    if (spec && strncmp(spec, "SCALE=", 6) == 0) {
        // multi-scale bootstrapping called by AU test
        int orig_nsite = nsite;
        double scale = convert_double(spec+6);
        nsite = (int)round(scale * nsite);
        for (site = 0; site < nsite; site++) {
            int site_id = random_int(orig_nsite, rstream);
            int ptn_id = getPatternID(site_id);
            pattern_freq[ptn_id]++;
        }
    } else if (!spec) {

        int nptn = getNPattern();

        if (nsite/8 < nptn || Params::getInstance().jackknife_prop > 0.0) {
            IntVector sample;
            random_resampling(nsite, sample, rstream);
            for (site = 0; site < nsite; site++)
            for (int rep = 0; rep < sample[site]; rep++) {
                int ptn_id = getPatternID(site);
                pattern_freq[ptn_id]++;
            }
        } else {
            // BQM 2015-12-27: use multinomial sampling for faster generation if #sites is much larger than #patterns
            int ptn;
            double *prob = new double[nptn];
            for (ptn = 0; ptn < nptn; ptn++)
                prob[ptn] = at(ptn).frequency;
            gsl_ran_multinomial(nptn, nsite, prob, (unsigned int*)pattern_freq, rstream);
            int sum = 0;
            for (ptn = 0; ptn < nptn; ptn++)
                sum += pattern_freq[ptn];
            ASSERT(sum == nsite);
            delete [] prob;
        }
    } else if (strncmp(spec, "GENESITE,", 9) == 0) {
		// resampling genes, then resampling sites within resampled genes
		convert_int_vec(spec+9, site_vec);
		int i;
		IntVector begin_site;
		for (i = 0, site = 0; i < site_vec.size(); i++) {
			begin_site.push_back(site);
			site += site_vec[i];
			//cout << "site = " << site_vec[i] << endl;
		}
		if (site > getNSite())
			outError("Sum of lengths exceeded alignment length");

		for (i = 0; i < site_vec.size(); i++) {
			int part = random_int(site_vec.size(), rstream);
			for (int j = 0; j < site_vec[part]; j++) {
				site = random_int(site_vec[part], rstream) + begin_site[part];
				int ptn = getPatternID(site);
				pattern_freq[ptn]++;
			}
		}
	} else if (strncmp(spec, "GENE,", 5) == 0) {
		// resampling genes instead of sites
		convert_int_vec(spec+5, site_vec);
		int i;
		IntVector begin_site;
		for (i = 0, site = 0; i < site_vec.size(); i++) {
			begin_site.push_back(site);
			site += site_vec[i];
			//cout << "site = " << site_vec[i] << endl;
		}
		if (site > getNSite())
			outError("Sum of lengths exceeded alignment length");

		for (i = 0; i < site_vec.size(); i++) {
			int part = random_int(site_vec.size(), rstream);
			for (site = begin_site[part]; site < begin_site[part] + site_vec[part]; site++) {
				int ptn = getPatternID(site);
				pattern_freq[ptn]++;
			}
		}
	} else {
		// resampling sites within genes
        try {
            convert_int_vec(spec, site_vec);
        } catch (...) {
            outError("-bsam not allowed for non-partition model");
        }
		if (site_vec.size() % 2 != 0)
			outError("Bootstrap specification length is not divisible by 2");
		int part, begin_site = 0, out_site = 0;
		for (part = 0; part < site_vec.size(); part += 2) {
			if (begin_site + site_vec[part] > getNSite())
				outError("Sum of lengths exceeded alignment length");
			for (site = 0; site < site_vec[part+1]; site++) {
				int site_id = random_int(site_vec[part], rstream) + begin_site;
				int ptn_id = getPatternID(site_id);
				pattern_freq[ptn_id]++;
			}
			begin_site += site_vec[part];
			out_site += site_vec[part+1];
		}
    }
}


void Alignment::buildFromPatternFreq(Alignment & aln, IntVector new_pattern_freqs){
	int nsite = aln.getNSite();
    seq_names.insert(seq_names.begin(), aln.seq_names.begin(), aln.seq_names.end());
    name = aln.name;
    model_name = aln.model_name;
    sequence_type = aln.sequence_type;
    position_spec = aln.position_spec;
    aln_file = aln.aln_file;
    num_states = aln.num_states;
    seq_type = aln.seq_type;

    genetic_code = aln.genetic_code;
    STATE_UNKNOWN = aln.STATE_UNKNOWN;
    site_pattern.resize(nsite, -1);

    clear();
    pattern_index.clear();

    int site = 0;
    std::vector<Pattern>::iterator it;
    int p;

    for(it = aln.begin(), p = 0; it != aln.end(); ++it, ++p) {
    	if(new_pattern_freqs[p] > 0){
	    	Pattern pat = *it;
			addPattern(pat, site, new_pattern_freqs[p]);
			for (int j = 0; j < new_pattern_freqs[p]; j++)
				site_pattern[site++] = size()-1;
    	}
    }
    if (!aln.site_state_freq.empty()) {
        site_model = site_pattern;
        ASSERT(site_state_freq.size() == getNPattern());
    }

    countConstSite();
    buildSeqStates();
//    checkSeqName();
}


void Alignment::createGapMaskedAlignment(Alignment *masked_aln, Alignment *aln) {
    if (masked_aln->getNSeq() != aln->getNSeq()) outError("Different number of sequences in masked alignment");
    if (masked_aln->getNSite() != aln->getNSite()) outError("Different number of sites in masked alignment");

    int site, nsite = aln->getNSite(), nseq = aln->getNSeq();
    seq_names.insert(seq_names.begin(), aln->seq_names.begin(), aln->seq_names.end());
    name = aln->name;
    model_name = aln->model_name;
    sequence_type = aln->sequence_type;
    position_spec = aln->position_spec;
    aln_file = aln->aln_file;
    num_states = aln->num_states;
    seq_type = aln->seq_type;
    genetic_code = aln->genetic_code;
    if (seq_type == SEQ_CODON) {
    	codon_table = new char[num_states];
    	memcpy(codon_table, aln->codon_table, num_states);
    	non_stop_codon = new char[strlen(genetic_code)];
    	memcpy(non_stop_codon, aln->non_stop_codon, strlen(genetic_code));
    }
    STATE_UNKNOWN = aln->STATE_UNKNOWN;
    site_pattern.resize(nsite, -1);
    clear();
    pattern_index.clear();
    IntVector name_map;
    for (StrVector::iterator it = seq_names.begin(); it != seq_names.end(); it++) {
        int seq_id = masked_aln->getSeqID(*it);
        if (seq_id < 0) outError("Masked alignment does not contain taxon ", *it);
        name_map.push_back(seq_id);
    }
    VerboseMode save_mode = verbose_mode;
    verbose_mode = min(verbose_mode, VB_MIN); // to avoid printing gappy sites in addPattern
    for (site = 0; site < nsite; site++) {
        int ptn_id = aln->getPatternID(site);
        Pattern pat = aln->at(ptn_id);
        Pattern masked_pat = masked_aln->at(masked_aln->getPatternID(site));
        for (int seq = 0; seq < nseq; seq++)
            if (masked_pat[name_map[seq]] == STATE_UNKNOWN) pat[seq] = STATE_UNKNOWN;
        addPattern(pat, site);
    }
    verbose_mode = save_mode;
    countConstSite();
    buildSeqStates();
}

void Alignment::shuffleAlignment() {
    if (isSuperAlignment()) outError("Internal error: ", __func__);
    my_random_shuffle(site_pattern.begin(), site_pattern.end());
}


void Alignment::concatenateAlignment(Alignment *aln) {
    if (getNSeq() != aln->getNSeq()) outError("Different number of sequences in two alignments");
    if (num_states != aln->num_states) outError("Different number of states in two alignments");
    if (seq_type != aln->seq_type) outError("Different data type in two alignments");
    int site, nsite = aln->getNSite();
    int cur_sites = getNSite();
    site_pattern.resize(cur_sites + nsite , -1);
    IntVector name_map;
    for (StrVector::iterator it = seq_names.begin(); it != seq_names.end(); it++) {
        int seq_id = aln->getSeqID(*it);
        if (seq_id < 0) outError("The other alignment does not contain taxon ", *it);
        name_map.push_back(seq_id);
    }
    VerboseMode save_mode = verbose_mode;
    verbose_mode = min(verbose_mode, VB_MIN); // to avoid printing gappy sites in addPattern
    for (site = 0; site < nsite; site++) {
        Pattern pat = aln->at(aln->getPatternID(site));
        Pattern new_pat = pat;
        for (int i = 0; i < name_map.size(); i++) new_pat[i] = pat[name_map[i]];
        addPattern(new_pat, site + cur_sites);
    }
    verbose_mode = save_mode;
    countConstSite();
    buildSeqStates();
}

void Alignment::copyAlignment(Alignment *aln) {
    int site, nsite = aln->getNSite();
    seq_names.insert(seq_names.begin(), aln->seq_names.begin(), aln->seq_names.end());
    name = aln->name;
    model_name = aln->model_name;
    sequence_type = aln->sequence_type;
    position_spec = aln->position_spec;
    aln_file = aln->aln_file;
    num_states = aln->num_states;
    seq_type = aln->seq_type;
    genetic_code = aln->genetic_code;
    if (seq_type == SEQ_CODON) {
    	codon_table = new char[num_states];
    	memcpy(codon_table, aln->codon_table, num_states);
    	non_stop_codon = new char[strlen(genetic_code)];
    	memcpy(non_stop_codon, aln->non_stop_codon, strlen(genetic_code));
    }
    STATE_UNKNOWN = aln->STATE_UNKNOWN;
    site_pattern.resize(nsite, -1);
    clear();
    pattern_index.clear();
    VerboseMode save_mode = verbose_mode;
    verbose_mode = min(verbose_mode, VB_MIN); // to avoid printing gappy sites in addPattern
    for (site = 0; site < nsite; site++) {
        int site_id = site;
        int ptn_id = aln->getPatternID(site_id);
        Pattern pat = aln->at(ptn_id);
        addPattern(pat, site);
    }
    verbose_mode = save_mode;
    countConstSite();
    buildSeqStates();
}

void Alignment::countConstSite() {
    int num_const_sites = 0;
    num_informative_sites = 0;
    num_variant_sites = 0;
    int num_invariant_sites = 0;
    num_parsimony_sites = 0;
    for (iterator it = begin(); it != end(); it++) {
        if ((*it).isConst())
            num_const_sites += (*it).frequency;
        if (it->isInformative())
            num_informative_sites += it->frequency;
        if (it->isInvariant())
            num_invariant_sites += it->frequency;
        else
            num_variant_sites += it->frequency;
    }
    frac_const_sites = ((double)num_const_sites) / getNSite();
    frac_invariant_sites = ((double)num_invariant_sites) / getNSite();
}

string Alignment::getUnobservedConstPatterns() {
	string ret = "";
	for (char state = 0; state < num_states; state++)
    if (!isStopCodon(state))
    {
		Pattern pat;
		pat.resize(getNSeq(), state);
		if (pattern_index.find(pat) == pattern_index.end()) {
			// constant pattern is unobserved
			ret.push_back(state);
		}
	}
	return ret;
}

int Alignment::countProperChar(int seq_id) {
    int num_proper_chars = 0;
    for (iterator it = begin(); it != end(); it++) {
        if ((*it)[seq_id] < num_states + pomo_sampled_states.size()) num_proper_chars+=(*it).frequency;
    }
    return num_proper_chars;
}

Alignment::~Alignment()
{
	if (codon_table) {
		delete [] codon_table;
		codon_table = NULL;
	}
	if (non_stop_codon) {
		delete [] non_stop_codon;
		non_stop_codon = NULL;
	}
    if (pars_lower_bound) {
        delete [] pars_lower_bound;
        pars_lower_bound = NULL;
    }
    for (vector<double*>::reverse_iterator it = site_state_freq.rbegin(); it != site_state_freq.rend(); it++)
        if (*it) delete [] (*it);
    site_state_freq.clear();
    site_model.clear();
}

double Alignment::computeObsDist(int seq1, int seq2) {
    int diff_pos = 0, total_pos = 0;
    for (iterator it = begin(); it != end(); it++) {
        int state1 = convertPomoState((*it)[seq1]);
        int state2 = convertPomoState((*it)[seq2]);
        if  (state1 < num_states && state2 < num_states) {
            //if ((*it)[seq1] != STATE_UNKNOWN && (*it)[seq2] != STATE_UNKNOWN) {
            total_pos += (*it).frequency;
            if ((*it)[seq1] != (*it)[seq2] )
                diff_pos += (*it).frequency;
        }
    }
    if (!total_pos) {
        if (verbose_mode >= VB_MED)
            outWarning("No overlapping characters between " + getSeqName(seq1) + " and " + getSeqName(seq2));
        return MAX_GENETIC_DIST; // return +INF if no overlap between two sequences
    }
    return ((double)diff_pos) / total_pos;
}

double Alignment::computeJCDist(int seq1, int seq2) {
    double obs_dist = computeObsDist(seq1, seq2);
    double z = (double)num_states / (num_states-1);
    double x = 1.0 - (z * obs_dist);

    if (x <= 0) {
        // string str = "Too long distance between two sequences ";
        // str += getSeqName(seq1);
        // str += " and ";
        // str += getSeqName(seq2);
        // outWarning(str);
        return MAX_GENETIC_DIST;
    }

    return -log(x) / z;
}

void Alignment::printDist(ostream &out, double *dist_mat) {
    int nseqs = getNSeq();
    int max_len = getMaxSeqNameLength();
    if (max_len < 10) max_len = 10;
    out << nseqs << endl;
    int pos = 0;
    out.precision(max((int)ceil(-log10(Params::getInstance().min_branch_length))+1, 6));
    out << fixed;
    for (int seq1 = 0; seq1 < nseqs; seq1 ++)  {
        out.width(max_len);
        out << left << getSeqName(seq1) << " ";
        for (int seq2 = 0; seq2 < nseqs; seq2 ++) {
            out << dist_mat[pos++];
            /*if (seq2 % 7 == 6) {
            	out << endl;
            	out.width(max_len+1);
            } */
            out << " ";
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
        //cout << "Distance matrix was printed to " << file_name << endl;
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, file_name);
    }
}

double Alignment::readDist(istream &in, double *dist_mat) {
    double longest_dist = 0.0;
    int nseqs;
    in >> nseqs;
    if (nseqs != getNSeq())
        throw "Distance file has different number of taxa";
    double *tmp_dist_mat = new double[nseqs * nseqs];
    std::map< string, int > map_seqName_ID;
    int pos = 0, seq1, seq2, id = 0;
    // read in distances to a temporary array
    for (seq1 = 0; seq1 < nseqs; seq1++)  {
        string seq_name;
        in >> seq_name;
        // assign taxa name to integer id
        map_seqName_ID[seq_name] = id++;
        /*
        if (seq_name != getSeqName(seq1))
            throw "Sequence name " + seq_name + " is different from " + getSeqName(seq1);
        for (seq2 = 0; seq2 < nseqs; seq2++) {
            in >> dist_mat[pos++];
            if (dist_mat[pos-1] > longest_dist)
                longest_dist = dist_mat[pos-1];
        }
         */
        for (seq2 = 0; seq2 < nseqs; seq2++) {
            in >> tmp_dist_mat[pos++];
            //cout << tmp_dist_mat[pos - 1] << "  ";
            if (tmp_dist_mat[pos - 1] > longest_dist)
                longest_dist = tmp_dist_mat[pos - 1];
        }
        //cout << endl;
    }
    //cout << "Internal distance matrix: " << endl;
    // Now initialize the internal distance matrix, in which the sequence order is the same
    // as in the alignment
    for (seq1 = 0; seq1 < nseqs; seq1++) {
        for (seq2 = 0; seq2 < nseqs; seq2++) {
            string seq1Name = getSeqName(seq1);
            string seq2Name = getSeqName(seq2);
            if (map_seqName_ID.count(seq1Name) == 0) {
                throw "Could not find taxa name " + seq1Name;
            }
            if (map_seqName_ID.count(seq2Name) == 0) {
                throw "Could not find taxa name " + seq2Name;
            }
            int seq1_tmp_id = map_seqName_ID[seq1Name];
            int seq2_tmp_id = map_seqName_ID[seq2Name];
            dist_mat[seq1 * nseqs + seq2] = tmp_dist_mat[seq1_tmp_id * nseqs + seq2_tmp_id];
            //cout << dist_mat[seq1 * nseqs + seq2] << "  ";
        }
        //cout << endl;
    }

    // check for symmetric matrix
    for (seq1 = 0; seq1 < nseqs-1; seq1++) {
        if (dist_mat[seq1*nseqs+seq1] != 0.0)
            throw "Diagonal elements of distance matrix is not ZERO";
        for (seq2 = seq1+1; seq2 < nseqs; seq2++)
            if (dist_mat[seq1*nseqs+seq2] != dist_mat[seq2*nseqs+seq1])
                throw "Distance between " + getSeqName(seq1) + " and " + getSeqName(seq2) + " is not symmetric";
    }

    /*
    string dist_file = params.out_prefix;
    dist_file += ".userdist";
    printDist(dist_file.c_str(), dist_mat);*/

    delete [] tmp_dist_mat;

    return longest_dist;
}

double Alignment::readDist(const char *file_name, double *dist_mat) {
    double longest_dist = 0.0;

    try {
        ifstream in;
        in.exceptions(ios::failbit | ios::badbit);
        in.open(file_name);
        longest_dist = readDist(in, dist_mat);
        in.close();
        cout << "Distance matrix was read from " << file_name << endl;
    } catch (const char *str) {
        outError(str);
    } catch (string str) {
        outError(str);
    } catch (ios::failure) {
        outError(ERR_READ_INPUT, file_name);
    }
    return longest_dist;
}

// TODO DS: This only works when the sampling method is SAMPLING_SAMPLED or when
// the virtual population size is also the sample size (for every species and
// every site).
void Alignment::computeStateFreq (double *state_freq, size_t num_unknown_states) {
    int i, j;
    double *states_app = new double[num_states*(STATE_UNKNOWN+1)];
    double *new_freq = new double[num_states];
    unsigned *state_count = new unsigned[STATE_UNKNOWN+1];
    double *new_state_freq = new double[num_states];


    memset(state_count, 0, sizeof(unsigned)*(STATE_UNKNOWN+1));
    state_count[(int)STATE_UNKNOWN] = num_unknown_states;

    for (i = 0; i <= STATE_UNKNOWN; i++)
        getAppearance(i, &states_app[i*num_states]);

    size_t aln_len = 0;

    for (iterator it = begin(); it != end(); it++) {
        aln_len += it->frequency;
        for (Pattern::iterator it2 = it->begin(); it2 != it->end(); it2++)
            state_count[convertPomoState((int)*it2)] += it->frequency;
    }

    for (i = 0; i < num_states; i++)
        state_freq[i] = 1.0/num_states;

    const int NUM_TIME = 8;
    if (aln_len > 0)
    for (int k = 0; k < NUM_TIME; k++) {
        memset(new_state_freq, 0, sizeof(double)*num_states);

        for (i = 0; i <= STATE_UNKNOWN; i++) {
            if (state_count[i] == 0) continue;
            double sum_freq = 0.0;
            for (j = 0; j < num_states; j++) {
                new_freq[j] = state_freq[j] * states_app[i*num_states+j];
                sum_freq += new_freq[j];
            }
            sum_freq = 1.0/sum_freq;
            for (j = 0; j < num_states; j++) {
                new_state_freq[j] += new_freq[j]*sum_freq*state_count[i];
            }
        }

        double sum_freq = 0.0;
        for (j = 0; j < num_states; j++)
            sum_freq += new_state_freq[j];
        sum_freq = 1.0/sum_freq;
        for (j = 0; j < num_states; j++)
            state_freq[j] = new_state_freq[j]*sum_freq;
    }

	convfreq(state_freq);

    if (verbose_mode >= VB_MED) {
        cout << "Empirical state frequencies: ";
        cout << setprecision(10);
        for (i = 0; i < num_states; i++)
            cout << state_freq[i] << " ";
        cout << endl;
    }

    delete [] new_state_freq;
    delete [] state_count;
    delete [] new_freq;
    delete [] states_app;
}

int Alignment::convertPomoState(int state) {
  // This map from an observed state to a PoMo state influences parsimony
  // construction and the +I likelihood computation. It should not make too much
  // of a difference though.

    if (seq_type != SEQ_POMO) return state;
    if (state < num_states) return state;
    if (state == STATE_UNKNOWN) return state;
    state -= num_states;
    if (pomo_sampled_states.size() <= 0)
        outError("Alignment file is too short.");
    if (state >= pomo_sampled_states.size()) {
        cout << "state:              " << state << endl;
        cout << "pomo_sampled_states.size(): " << pomo_sampled_states.size() << endl;
    }
    ASSERT(state < pomo_sampled_states.size());
    int id1 = pomo_sampled_states[state] & 3;
    int id2 = (pomo_sampled_states[state] >> 16) & 3;
    int value1 = (pomo_sampled_states[state] >> 2) & 16383;
    int value2 = pomo_sampled_states[state] >> 18;
    int N = virtual_pop_size;
    int M = value1 + value2;

    // Mon Jun 13 13:24:55 CEST 2016. This is a stochastic way to assign PoMo
    // states. This is important if the sample size is small.

    // double stoch = (double) rand() / RAND_MAX - 0.5;
    // stoch /= 2.0;
    // int pick = (int)round(((double) value1*N/M) + stoch);

    // Fri Aug 18 15:37:22 BST 2017 However, Minh prefers a deterministic way
    // that, necessarily, introduces some systematic error.

    // BQM: Prefer the state with highest likelihood.

    // TODO: How to break tie? E.g., 4A4C but N=9? This way always prefers the
    // first allele (which is equivalent to a bias towards A and C, kind of).

    int pick = (int)round(((double) value1*N/M));

    int real_state;
    if (pick <= 0)
        real_state = id2;
    else if (pick >= N)
        real_state = id1;
    else {
        int j;
        if (id1 == 0) j = id2 - 1;
        else j = id1 + id2;
        real_state = 3 + j*(N-1) + pick;
    }
    state = real_state;
    ASSERT(state < num_states);
    return state;
}

void Alignment::computeAbsoluteStateFreq(unsigned int *abs_state_freq) {
    memset(abs_state_freq, 0, num_states * sizeof(unsigned int));

    if (seq_type == SEQ_POMO) {
        for (iterator it = begin(); it != end(); it++)
            for (Pattern::iterator it2 = it->begin(); it2 != it->end(); it2++)
                abs_state_freq[convertPomoState((int)*it2)] += it->frequency;
    } else {
        for (iterator it = begin(); it != end(); it++)
            for (Pattern::iterator it2 = it->begin(); it2 != it->end(); it2++)
                abs_state_freq[(int)*it2] += it->frequency;
    }
}


void Alignment::countStatePerSequence (unsigned *count_per_sequence) {
    int i;
    int nseqs = getNSeq();
    memset(count_per_sequence, 0, sizeof(unsigned)*num_states*nseqs);
    for (iterator it = begin(); it != end(); it++)
        for (i = 0; i != nseqs; i++) {
            int state = convertPomoState(it->at(i));
            if (state < num_states) {
                count_per_sequence[i*num_states + state] += it->frequency;
            }
        }
}

void Alignment::computeStateFreqPerSequence (double *freq_per_sequence) {
    int i, j;
    int nseqs = getNSeq();
    double *states_app = new double[num_states*(STATE_UNKNOWN+1)];
    double *new_freq = new double[num_states];
    unsigned *state_count = new unsigned[(STATE_UNKNOWN+1)*nseqs];
    double *new_state_freq = new double[num_states];


    memset(state_count, 0, sizeof(unsigned)*(STATE_UNKNOWN+1)*nseqs);

    for (i = 0; i <= STATE_UNKNOWN; i++)
        getAppearance(i, &states_app[i*num_states]);

    for (iterator it = begin(); it != end(); it++)
        for (i = 0; i != nseqs; i++) {
            state_count[i*(STATE_UNKNOWN+1) + it->at(i)] += it->frequency;
        }
    double equal_freq = 1.0/num_states;
    for (i = 0; i < num_states*nseqs; i++)
        freq_per_sequence[i] = equal_freq;

    const int NUM_TIME = 8;
    for (int k = 0; k < NUM_TIME; k++) {
        for (int seq = 0; seq < nseqs; seq++) {
            double *state_freq = &freq_per_sequence[seq*num_states];
            memset(new_state_freq, 0, sizeof(double)*num_states);
            for (i = 0; i <= STATE_UNKNOWN; i++) {
                if (state_count[seq*(STATE_UNKNOWN+1)+i] == 0) continue;
                double sum_freq = 0.0;
                for (j = 0; j < num_states; j++) {
                    new_freq[j] = state_freq[j] * states_app[i*num_states+j];
                    sum_freq += new_freq[j];
                }
                sum_freq = 1.0/sum_freq;
                for (j = 0; j < num_states; j++) {
                    new_state_freq[j] += new_freq[j]*sum_freq*state_count[seq*(STATE_UNKNOWN+1)+i];
                }
            }

            double sum_freq = 0.0;
            for (j = 0; j < num_states; j++)
                sum_freq += new_state_freq[j];
            sum_freq = 1.0/sum_freq;
            for (j = 0; j < num_states; j++)
                state_freq[j] = new_state_freq[j]*sum_freq;
         }
    }

//	convfreq(state_freq);
//
//    if (verbose_mode >= VB_MED) {
//        cout << "Empirical state frequencies: ";
//        for (i = 0; i < num_states; i++)
//            cout << state_freq[i] << " ";
//        cout << endl;
//    }

    delete [] new_state_freq;
    delete [] state_count;
    delete [] new_freq;
    delete [] states_app;
}

//void Alignment::computeStateFreq (double *stateFrqArr) {
//    int stateNo_;
//    int nState_ = num_states;
//    int nseqs = getNSeq();
//    double *timeAppArr_ = new double[num_states];
//    double *siteAppArr_ = new double[num_states]; //App = appearance
//    double *newSiteAppArr_ = new double[num_states];
//
//    for (stateNo_ = 0; stateNo_ < nState_; stateNo_ ++)
//        stateFrqArr [ stateNo_ ] = 1.0 / nState_;
//
//    int NUM_TIME = 8;
//    //app = appeareance
//    if (verbose_mode >= VB_MED)
//        cout << "Computing state frequencies..." << endl;
//    for (int time_ = 0; time_ < NUM_TIME; time_ ++)
//    {
//        for (stateNo_ = 0; stateNo_ < nState_; stateNo_ ++)
//            timeAppArr_[stateNo_] = 0.0;
//
//        for (iterator it = begin(); it != end(); it++)
//            for (int i = 0; i < (*it).frequency; i++)
//            {
//                for (int seq = 0; seq < nseqs; seq++) {
//                    int stateNo_ = (*it)[seq];
//
//                    getAppearance (stateNo_, siteAppArr_);
//
//                    double totalSiteApp_ = 0.0;
//                    for (stateNo_ = 0; stateNo_ < nState_; stateNo_ ++) {
//                        newSiteAppArr_[stateNo_] = stateFrqArr[stateNo_] * siteAppArr_[stateNo_];
//                        totalSiteApp_ += newSiteAppArr_[stateNo_];
//                    }
//                    totalSiteApp_ = 1.0 / totalSiteApp_;
//
//                    for (stateNo_ = 0; stateNo_ < nState_; stateNo_ ++)
//                        timeAppArr_[stateNo_] += newSiteAppArr_[stateNo_] * totalSiteApp_;
//                }
//            }
//
//        double totalTimeApp_ = 0.0;
//        int stateNo_;
//        for (stateNo_ = 0; stateNo_ < nState_; stateNo_ ++)
//            totalTimeApp_ += timeAppArr_[stateNo_];
//
//
//        for (stateNo_ = 0; stateNo_ < nState_; stateNo_ ++)
//            stateFrqArr[stateNo_] = timeAppArr_[stateNo_] / totalTimeApp_;
//
//    } //end of for time_
//
//    //  std::cout << "state frequency ..." << endl;
//    // for (stateNo_ = 0; stateNo_ < nState_; stateNo_ ++)
//    // std::cout << stateFrqArr[stateNo_] << endl;
//
//	convfreq(stateFrqArr);
//
//    if (verbose_mode >= VB_MED) {
//        cout << "Empirical state frequencies: ";
//        for (stateNo_ = 0; stateNo_ < nState_; stateNo_ ++)
//            cout << stateFrqArr[stateNo_] << " ";
//        cout << endl;
//    }
//	delete [] newSiteAppArr_;
//	delete [] siteAppArr_;
//	delete [] timeAppArr_;
//
//}

void Alignment::getAppearance(StateType state, double *state_app) {
    int i;
    if (state == STATE_UNKNOWN) {
        for (i = 0; i < num_states; i++)
            state_app[i] = 1.0;
        return;
    }

    memset(state_app, 0, num_states * sizeof(double));
    if (state < num_states) {
        state_app[(int)state] = 1.0;
        return;
    }
	// ambiguous characters
	int ambi_aa[] = {4+8, 32+64, 512+1024};
	switch (seq_type) {
	case SEQ_DNA:
	    state -= (num_states-1);
		for (i = 0; i < num_states; i++)
			if (state & (1 << i)) {
				state_app[i] = 1.0;
			}
		break;
	case SEQ_PROTEIN:
		ASSERT(state<23);
		state -= 20;
		for (i = 0; i < 11; i++)
			if (ambi_aa[(int)state] & (1<<i)) {
				state_app[i] = 1.0;
			}
		break;
    case SEQ_POMO:
//        state -= num_states;
//        assert(state < pomo_sampled_states.size());
//        // count the number of nucleotides
//        state_app[pomo_sampled_states[state] & 3] = 1.0;
//        state_app[(pomo_sampled_states[state] >> 16) & 3] = 1.0;
        state_app[convertPomoState(state)] = 1.0;
        break;
	default: ASSERT(0); break;
	}
}

void Alignment::getAppearance(StateType state, StateBitset &state_app) {

	int i;
    if (state == STATE_UNKNOWN) {
    	state_app.set();
        return;
    }

    state_app.reset();
    if (state < num_states) {
        state_app[(int)state] = 1;
        return;
    }
	// ambiguous characters
	int ambi_aa[] = {4+8, 32+64, 512+1024};
	switch (seq_type) {
	case SEQ_DNA:
	    state -= (num_states-1);
		for (i = 0; i < num_states; i++)
			if (state & (1 << i)) {
				state_app[i] = 1;
			}
		break;
	case SEQ_PROTEIN:
		if (state >= 23) return;
		state -= 20;
		for (i = 0; i < 11; i++)
			if (ambi_aa[(int)state] & (1<<i)) {
				state_app[i] = 1;
			}
		break;
    case SEQ_POMO:
//        state -= num_states;
//        assert(state < pomo_sampled_states.size());
//        // count the number of nucleotides
//        state_app[pomo_sampled_states[state] & 3] = 1;
//        state_app[(pomo_sampled_states[state] >> 16) & 3] = 1;
        state_app[convertPomoState(state)] = 1;
        break;
	default: ASSERT(0); break;
	}
}

void Alignment::computeCodonFreq(StateFreqType freq, double *state_freq, double *ntfreq) {
	int nseqs = getNSeq();
	int i, j;

	if (freq == FREQ_CODON_1x4) {
		memset(ntfreq, 0, sizeof(double)*4);
		for (iterator it = begin(); it != end(); it++) {
			for (int seq = 0; seq < nseqs; seq++) if ((*it)[seq] != STATE_UNKNOWN) {
				int codon = codon_table[(int)(*it)[seq]];
//				int codon = (int)(*it)[seq];
				int nt1 = codon / 16;
				int nt2 = (codon % 16) / 4;
				int nt3 = codon % 4;
				ntfreq[nt1] += (*it).frequency;
				ntfreq[nt2] += (*it).frequency;
				ntfreq[nt3] += (*it).frequency;
			}
		}
		double sum = 0;
		for (i = 0; i < 4; i++)
			sum += ntfreq[i];
		for (i = 0; i < 4; i++)
			ntfreq[i] /= sum;
		if (verbose_mode >= VB_MED) {
			for (i = 0; i < 4; i++)
				cout << "  " << symbols_dna[i] << ": " << ntfreq[i];
			cout << endl;
		}
		memcpy(ntfreq+4, ntfreq, sizeof(double)*4);
		memcpy(ntfreq+8, ntfreq, sizeof(double)*4);
        sum = 0.0;
		for (i = 0; i < num_states; i++) {
            int codon = codon_table[i];
            state_freq[i] = ntfreq[codon/16] * ntfreq[(codon%16)/4] * ntfreq[codon%4];
			if (isStopCodon(i)) {
//                sum_stop += state_freq[i];
				state_freq[i] = MIN_FREQUENCY;
			} else {
                sum += state_freq[i];
            }
        }
//        sum = (1.0-sum)/(1.0-sum_stop);
        sum = 1.0/sum;
		for (i = 0; i < num_states; i++)
            if (!isStopCodon(i))
                state_freq[i] *= sum;
        sum = 0.0;
		for (i = 0; i < num_states; i++)
                sum += state_freq[i];
        ASSERT(fabs(sum-1.0)<1e-5);
	} else if (freq == FREQ_CODON_3x4) {
		// F3x4 frequency model
		memset(ntfreq, 0, sizeof(double)*12);
		for (iterator it = begin(); it != end(); it++) {
			for (int seq = 0; seq < nseqs; seq++) if ((*it)[seq] != STATE_UNKNOWN) {
				int codon = codon_table[(int)(*it)[seq]];
//				int codon = (int)(*it)[seq];
				int nt1 = codon / 16;
				int nt2 = (codon % 16) / 4;
				int nt3 = codon % 4;
				ntfreq[nt1] += (*it).frequency;
				ntfreq[4+nt2] += (*it).frequency;
				ntfreq[8+nt3] += (*it).frequency;
			}
		}
		for (j = 0; j < 12; j+=4) {
			double sum = 0;
			for (i = 0; i < 4; i++)
				sum += ntfreq[i+j];
			for (i = 0; i < 4; i++)
				ntfreq[i+j] /= sum;
			if (verbose_mode >= VB_MED) {
				for (i = 0; i < 4; i++)
					cout << "  " << symbols_dna[i] << ": " << ntfreq[i+j];
				cout << endl;
			}
		}

//        double sum_stop=0.0;
        double sum = 0.0;
		for (i = 0; i < num_states; i++) {
            int codon = codon_table[i];
            state_freq[i] = ntfreq[codon/16] * ntfreq[4+(codon%16)/4] * ntfreq[8+codon%4];
			if (isStopCodon(i)) {
//                sum_stop += state_freq[i];
				state_freq[i] = MIN_FREQUENCY;
			} else {
                sum += state_freq[i];
            }
        }
//        sum = (1.0-sum)/(1.0-sum_stop);
        sum = 1.0 / sum;
		for (i = 0; i < num_states; i++)
            if (!isStopCodon(i))
                state_freq[i] *= sum;
        sum = 0.0;
		for (i = 0; i < num_states; i++)
                sum += state_freq[i];
        ASSERT(fabs(sum-1.0)<1e-5);

//		double sum = 0;
//		for (i = 0; i < num_states; i++)
//			if (isStopCodon(i)) {
//				state_freq[i] = 0.0;
//			} else {
//				//int codon = codon_table[i];
//				int codon = i;
//				state_freq[i] = ntfreq[codon/16] * ntfreq[4+(codon%16)/4] * ntfreq[8+codon%4];
//				sum += state_freq[i];
//			}
//		for (i = 0; i < num_states; i++)
//			state_freq[i] /= sum;

        // now recompute ntfreq based on state_freq
//        memset(ntfreq, 0, 12*sizeof(double));
//        for (i = 0; i < num_states; i++)
//            if (!isStopCodon(i)) {
//				int nt1 = i / 16;
//				int nt2 = (i % 16) / 4;
//				int nt3 = i % 4;
//                ntfreq[nt1] += state_freq[i];
//                ntfreq[nt2+4] += state_freq[i];
//                ntfreq[nt3+8] += state_freq[i];
//            }
//		for (j = 0; j < 12; j+=4) {
//			double sum = 0;
//			for (i = 0; i < 4; i++)
//				sum += ntfreq[i+j];
//			for (i = 0; i < 4; i++)
//				ntfreq[i+j] /= sum;
//			if (verbose_mode >= VB_MED) {
//				for (i = 0; i < 4; i++)
//					cout << "  " << symbols_dna[i] << ": " << ntfreq[i+j];
//				cout << endl;
//			}
//		}
	} else if (freq == FREQ_CODON_3x4C) {
        outError("F3X4C not yet implemented. Contact authors if you really need it.");
	} else if (freq == FREQ_EMPIRICAL || freq == FREQ_ESTIMATE) {
		memset(state_freq, 0, num_states*sizeof(double));
        i = 0;
        for (iterator it = begin(); it != end(); it++, i++)
			for (int seq = 0; seq < nseqs; seq++) {
				int state = it->at(seq);
				if (state >= num_states) continue;
				state_freq[state] += it->frequency;
			}
        double sum = 0.0;
        for (i = 0; i < num_states; i++)
        	sum += state_freq[i];
        for (i = 0; i < num_states; i++)
        	state_freq[i] /= sum;
	} else {
        outError("Unsupported codon frequency");
    }
	convfreq(state_freq);
}

void Alignment::computeDivergenceMatrix(double *rates) {
    int i, j, k;
    ASSERT(rates);
    int nseqs = getNSeq();
    unsigned *pair_rates = new unsigned[num_states*num_states];
    memset(pair_rates, 0, sizeof(unsigned)*num_states*num_states);
//    for (i = 0; i < num_states; i++) {
//        pair_rates[i] = new double[num_states];
//        memset(pair_rates[i], 0, sizeof(double)*num_states);
//    }

    unsigned *state_freq = new unsigned[STATE_UNKNOWN+1];

    for (iterator it = begin(); it != end(); it++) {
        memset(state_freq, 0, sizeof(unsigned)*(STATE_UNKNOWN+1));
        for (i = 0; i < nseqs; i++) {
            state_freq[(int)it->at(i)]++;
        }
        for (i = 0; i < num_states; i++) {
            if (state_freq[i] == 0) continue;
            pair_rates[i*num_states+i] += (state_freq[i]*(state_freq[i]-1)/2)*it->frequency;
            for (j = i+1; j < num_states; j++)
                pair_rates[i*num_states+j] += state_freq[i]*state_freq[j]*it->frequency;
        }
//            int state1 = it->at(i);
//            if (state1 >= num_states) continue;
//            int *this_pair = pair_rates + state1*num_states;
//            for (j = i+1; j < nseqs; j++) {
//                int state2 = it->at(j);
//                if (state2 < num_states) this_pair[state2] += it->frequency;
//            }
//        }
    }

    k = 0;
    double last_rate = pair_rates[(num_states-2)*num_states+num_states-1] + pair_rates[(num_states-1)*num_states+num_states-2];
    if (last_rate == 0) last_rate = 1;
    for (i = 0; i < num_states-1; i++)
        for (j = i+1; j < num_states; j++) {
            rates[k++] = (pair_rates[i*num_states+j] + pair_rates[j*num_states+i]) / last_rate;
            // BIG WARNING: zero rates might cause numerical instability!
//            if (rates[k-1] <= 0.0001) rates[k-1] = 0.01;
//            if (rates[k-1] > 100.0) rates[k-1] = 50.0;
        }
    rates[k-1] = 1;
    if (verbose_mode >= VB_MAX) {
        cout << "Empirical rates: ";
        for (k = 0; k < num_states*(num_states-1)/2; k++)
            cout << rates[k] << " ";
        cout << endl;
    }

//    for (i = num_states-1; i >= 0; i--) {
//        delete [] pair_rates[i];
//    }
    delete [] state_freq;
    delete [] pair_rates;
}

void Alignment::computeDivergenceMatrixNonRev (double *rates) {
    double *rates_mat = new double[num_states*num_states];
    int i, j, k;

    computeDivergenceMatrix(rates);

    for (i = 0, k = 0; i < num_states-1; i++)
        for (j = i+1; j < num_states; j++)
            rates_mat[i*num_states+j] = rates_mat[j*num_states+i] = rates[k++];

    for (i = 0, k = 0; i < num_states; i++)
        for (j = 0; j < num_states; j++)
            if (j != i) rates[k++] = rates_mat[i*num_states+j];
	delete [] rates_mat;

}

void Alignment::convfreq(double *stateFrqArr) {

    if (Params::getInstance().keep_zero_freq)
        return;
	int i, maxi=0;
	double freq, maxfreq, sum;
	int zero_states = 0;

	sum = 0.0;
	maxfreq = 0.0;
	for (i = 0; i < num_states; i++)
	{
		freq = stateFrqArr[i];
        // Do not check for a minimum frequency with PoMo because very
        // low frequencies are expected for polymorphic sites.
		if ((freq < MIN_FREQUENCY) &&
            (seq_type != SEQ_POMO)) {
			stateFrqArr[i] = MIN_FREQUENCY;
		}
		if (freq > maxfreq) {
			maxfreq = freq;
			maxi = i;
		}

		sum += stateFrqArr[i];
	}
	stateFrqArr[maxi] += 1.0 - sum;

	// make state frequencies a bit different from each other
//	for (i = 0; i < num_states - 1; i++)
//		if (!isStopCodon(i))
//			for (j = i + 1; j < num_states; j++)
//				if (!isStopCodon(j))
//					if (stateFrqArr[i] == stateFrqArr[j]) {
//						stateFrqArr[i] += MIN_FREQUENCY_DIFF;
//						stateFrqArr[j] -= MIN_FREQUENCY_DIFF;
//					}
	if (zero_states) {
		cout << "WARNING: " << zero_states << " states not present in alignment that might cause numerical instability" << endl;
	}
} /* convfreq */

double Alignment::computeUnconstrainedLogL() {
    int nptn = size();
    double logl = 0.0;
    int nsite = getNSite(), i;
    double lognsite = log(nsite);
    for (i = 0; i < nptn; i++)
        logl += (log(at(i).frequency) - lognsite) * at(i).frequency;
    return logl;
}

void Alignment::printSiteGaps(const char *filename) {
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);

        out.open(filename);
        int nsite = getNSite();
        out << nsite << endl << "Site_Gap  ";
        for (int site = 0; site < getNSite(); site++) {
            out << " " << at(getPatternID(site)).computeGapChar(num_states, STATE_UNKNOWN);
        }
        out << endl << "Site_Ambi ";
        for (int site = 0; site < getNSite(); site++) {
            out << " " << at(getPatternID(site)).computeAmbiguousChar(num_states);
        }
        out << endl;
        out.close();
        cout << "Site gap-counts printed to " << filename << endl;
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, filename);
    }
}

void Alignment::getPatternFreq(IntVector &freq) {
	freq.resize(getNPattern());
	int cnt = 0;
	for (iterator it = begin(); it < end(); it++, cnt++)
		freq[cnt] = (*it).frequency;
}

void Alignment::getPatternFreq(int *freq) {
    int cnt = 0;
    for (iterator it = begin(); it < end(); it++, cnt++)
        freq[cnt] = (*it).frequency;
}

//added by MA
void Alignment::multinomialProb(Alignment refAlign, double &prob)
{
// 	cout << "Computing the probability of this alignment given the multinomial distribution determined by a reference alignment ..." << endl;
    //should we check for compatibility of sequence's names and sequence's order in THIS alignment and in the objectAlign??
    //check alignment length
    int nsite = getNSite();
    ASSERT(nsite == refAlign.getNSite());
    double sumFac = 0;
    double sumProb = 0;
    double fac = logFac(nsite);
    int index;
    for ( iterator it = begin(); it != end() ; it++)
    {
        PatternIntMap::iterator pat_it = refAlign.pattern_index.find((*it));
        if ( pat_it == refAlign.pattern_index.end() ) //not found ==> error
            outError("Pattern in the current alignment is not found in the reference alignment!");
        sumFac += logFac((*it).frequency);
        index = pat_it->second;
        sumProb += (double)(*it).frequency*log((double)refAlign.at(index).frequency/(double)nsite);
    }
    prob = fac - sumFac + sumProb;
}

void Alignment::multinomialProb (DoubleVector logLL, double &prob)
{
    //cout << "Function in Alignment: Compute probability of the expected alignment (determined by patterns log-likelihood under some tree and model) given THIS alignment." << endl;

    //The expected normalized requencies
    IntVector expectedNorFre;

    if ( logLL.empty())
        outError("Error: log likelihood of patterns are not given!");

    int patNum = getNPattern();

    ASSERT(logLL.size() == patNum);

    int alignLen = getNSite();
    //resize the expectedNorFre vector
    expectedNorFre.resize(patNum,-1);

    //Vector containing the 'relative' likelihood of the pattern p_i
    DoubleVector LL(patNum,-1.0);
    double sumLL = 0; //sum of the likelihood of the patterns in the alignment
    double max_logl = *max_element(logLL.begin(), logLL.end()); // to rescale the log-likelihood
    //Compute the `relative' (to the first pattern) likelihood from the logLL
    for ( int i = 0; i < patNum; i++ )
    {
        LL[i] = exp(logLL[i]-max_logl);
        //LL[i] = exp(logLL[i]);
        sumLL += LL[i];
    }

    //Vector containing l_i = p_i*ell/sum_i(p_i)
    DoubleVector ell(patNum, -1.0);
    //Compute l_i
    for ( int i = 0; i < patNum; i++ )
    {
        ell[i] = (double)alignLen * LL[i] / sumLL;
    }


    //Vector containing r_i where r_0 = ell_0; r_{i+1} = ell_{i+1} + r_i - ordinaryRounding(r_i)
    DoubleVector r(patNum, -1.0);
    //Compute r_i and the expected normalized frequencies
    r[0] = ell[0];
    expectedNorFre[0] = (int)floor(ell[0]+0.5); //note that floor(_number+0.5) returns the ordinary rounding of _number
    //int sum = expectedNorFre[0];
    for (int j = 1; j < patNum; j++ )
    {
        r[j] = ell[j] + r[j-1] - floor(r[j-1]+0.5);
        expectedNorFre[j] = (int)floor(r[j]+0.5);
        //sum += expectedNorFre[j];
    }

    //cout << "Number of patterns: " << patNum << ", sum of expected sites: " << sum << endl;
    //return expectedNorFre;
    //compute the probability of having expectedNorFre given the observed pattern frequencies of THIS alignment
    double sumFac = 0;
    double sumProb = 0;
    double fac = logFac(alignLen);
    for (int patID = 0; patID < patNum; patID++) {
        int patFre = expectedNorFre[patID];
        sumFac += logFac(patFre);
        sumProb += (double)patFre*log((double)at(patID).frequency/(double)alignLen);
    }
    prob = fac - sumFac + sumProb;
}

void Alignment::multinomialProb (double *logLL, double &prob)
{
    //cout << "Function in Alignment: Compute probability of the expected alignment (determined by patterns log-likelihood under some tree and model) given THIS alignment." << endl;

    //The expected normalized requencies
    IntVector expectedNorFre;

    /*	if ( logLL.empty())
    		outError("Error: log likelihood of patterns are not given!");*/

    int patNum = getNPattern();

    //assert(logLL.size() == patNum);

    int alignLen = getNSite();
    //resize the expectedNorFre vector
    expectedNorFre.resize(patNum,-1);

    //Vector containing the 'relative' likelihood of the pattern p_i
    DoubleVector LL(patNum,-1.0);
    double sumLL = 0; //sum of the likelihood of the patterns in the alignment
    double max_logl = *max_element(logLL, logLL + patNum); // to rescale the log-likelihood
    //Compute the `relative' (to the first pattern) likelihood from the logLL
    for ( int i = 0; i < patNum; i++ )
    {
        LL[i] = exp(logLL[i]-max_logl);
        //LL[i] = exp(logLL[i]);
        sumLL += LL[i];
    }

    //Vector containing l_i = p_i*ell/sum_i(p_i)
    DoubleVector ell(patNum, -1.0);
    //Compute l_i
    for ( int i = 0; i < patNum; i++ )
    {
        ell[i] = (double)alignLen * LL[i] / sumLL;
    }


    //Vector containing r_i where r_0 = ell_0; r_{i+1} = ell_{i+1} + r_i - ordinaryRounding(r_i)
    DoubleVector r(patNum, -1.0);
    //Compute r_i and the expected normalized frequencies
    r[0] = ell[0];
    expectedNorFre[0] = (int)floor(ell[0]+0.5); //note that floor(_number+0.5) returns the ordinary rounding of _number
    //int sum = expectedNorFre[0];
    for (int j = 1; j < patNum; j++ )
    {
        r[j] = ell[j] + r[j-1] - floor(r[j-1]+0.5);
        expectedNorFre[j] = (int)floor(r[j]+0.5);
        //sum += expectedNorFre[j];
    }

    //cout << "Number of patterns: " << patNum << ", sum of expected sites: " << sum << endl;
    //return expectedNorFre;
    //compute the probability of having expectedNorFre given the observed pattern frequencies of THIS alignment
    double sumFac = 0;
    double sumProb = 0;
    double fac = logFac(alignLen);
    for (int patID = 0; patID < patNum; patID++) {
        int patFre = expectedNorFre[patID];
        sumFac += logFac(patFre);
        sumProb += (double)patFre*log((double)at(patID).frequency/(double)alignLen);
    }
    prob = fac - sumFac + sumProb;
}

double Alignment::multinomialProb (IntVector &pattern_freq)
{
    //cout << "Function in Alignment: Compute probability of the expected alignment (determined by patterns log-likelihood under some tree and model) given THIS alignment." << endl;

    //The expected normalized requencies

    //cout << "Number of patterns: " << patNum << ", sum of expected sites: " << sum << endl;
    //return expectedNorFre;
    //compute the probability of having expectedNorFre given the observed pattern frequencies of THIS alignment
    ASSERT(size() == pattern_freq.size());
    int patNum = getNPattern();
    int alignLen = getNSite();
    double sumFac = 0;
    double sumProb = 0;
    double fac = logFac(alignLen);
    for (int patID = 0; patID < patNum; patID++) {
        int patFre = pattern_freq[patID];
        sumFac += logFac(patFre);
        sumProb += (double)patFre*log((double)at(patID).frequency/(double)alignLen);
    }
    return (fac - sumFac + sumProb);
}

bool Alignment::readSiteStateFreq(const char* site_freq_file)
{
	cout << endl << "Reading site-specific state frequency file " << site_freq_file << " ..." << endl;
	site_model.resize(getNSite(), -1);
    int i;
    IntVector pattern_to_site; // vector from pattern to the first site
    pattern_to_site.resize(getNPattern(), -1);
    for (i = 0; i < getNSite(); i++)
        if (pattern_to_site[getPatternID(i)] == -1)
            pattern_to_site[getPatternID(i)] = i;

    bool aln_changed = false;

	try {
		ifstream in;
		in.exceptions(ios::failbit | ios::badbit);
		in.open(site_freq_file);
		double freq;
		string site_spec;
		int specified_sites = 0;
		in.exceptions(ios::badbit);
		for (int model_id = 0; !in.eof(); model_id++) {
			// remove the failbit
			in >> site_spec;
			if (in.eof()) break;
			IntVector site_id;
			extractSiteID(this, site_spec.c_str(), site_id);
			specified_sites += site_id.size();
			if (site_id.size() == 0) throw "No site ID specified";
			for (IntVector::iterator it = site_id.begin(); it != site_id.end(); it++) {
				if (site_model[*it] != -1) throw "Duplicated site ID";
				site_model[*it] = site_state_freq.size();
			}
			double *site_freq_entry = new double[num_states];
			double sum = 0;
			for (i = 0; i < num_states; i++) {
				in >> freq;
				if (freq <= 0.0 || freq >= 1.0) throw "Frequencies must be strictly positive and smaller than 1";
				site_freq_entry[i] = freq;
				sum += freq;
			}
			if (fabs(sum-1.0) > 1e-4) {
                if (fabs(sum-1.0) > 1e-3)
                    outWarning("Frequencies of site " + site_spec + " do not sum up to 1 and will be normalized");
                sum = 1.0/sum;
                for (i = 0; i < num_states; i++)
                    site_freq_entry[i] *= sum;
            }
			convfreq(site_freq_entry); // regularize frequencies (eg if some freq = 0)

            // 2016-02-01: now check for equality of sites with same site-pattern and same freq
            int prev_site = pattern_to_site[getPatternID(site_id[0])];
            if (site_id.size() == 1 && prev_site < site_id[0] && site_model[prev_site] != -1) {
                // compare freq with prev_site
                bool matched_freq = true;
                double *prev_freq = site_state_freq[site_model[prev_site]];
                for (i = 0; i < num_states; i++) {
                    if (site_freq_entry[i] != prev_freq[i]) {
                        matched_freq = false;
                        break;
                    }
                }
                if (matched_freq) {
                    site_model[site_id[0]] = site_model[prev_site];
                } else
                    aln_changed = true;
            }

            if (site_model[site_id[0]] == site_state_freq.size())
                site_state_freq.push_back(site_freq_entry);
            else
                delete [] site_freq_entry;
		}
		if (specified_sites < site_model.size()) {
            aln_changed = true;
			// there are some unspecified sites
			cout << site_model.size() - specified_sites << " unspecified sites will get default frequencies" << endl;
			for (i = 0; i < site_model.size(); i++)
				if (site_model[i] == -1)
					site_model[i] = site_state_freq.size();
			site_state_freq.push_back(NULL);
		}
		in.clear();
		// set the failbit again
		in.exceptions(ios::failbit | ios::badbit);
		in.close();
	} catch (const char* str) {
		outError(str);
	} catch (string str) {
		outError(str);
	} catch(ios::failure) {
		outError(ERR_READ_INPUT);
	}

    if (aln_changed) {
        cout << "Regrouping alignment sites..." << endl;
        regroupSitePattern(site_state_freq.size(), site_model);
    }
    cout << site_state_freq.size() << " distinct per-site state frequency vectors detected" << endl;
    return aln_changed;
}
