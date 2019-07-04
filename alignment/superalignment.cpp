/***************************************************************************
 *   Copyright (C) 2009 by BUI Quang Minh   *
 *   minh.bui@univie.ac.at   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <stdarg.h>
#include "tree/phylotree.h"
#include "superalignment.h"
#include "tree/phylosupertree.h"
#include "nclextra/msetsblock.h"
#include "nclextra/myreader.h"
#include "main/phylotesting.h"

SuperAlignment::SuperAlignment() : Alignment() {
    max_num_states = 0;
}

SuperAlignment::SuperAlignment(Params &params) : Alignment()
{
    
    cout << "Reading partition model file " << params.partition_file << " ..." << endl;
    if (detectInputFile(params.partition_file) == IN_NEXUS) {
        readPartitionNexus(params);
        if (partitions.empty()) {
            outError("No partition found in SETS block. An example syntax looks like: \n#nexus\nbegin sets;\n  charset part1=1-100;\n  charset part2=101-300;\nend;");
        }
    } else
        readPartitionRaxml(params);
    if (partitions.empty())
        outError("No partition found");
    
    // check for duplicated partition names
    unordered_set<string> part_names;
    for (auto pit = partitions.begin(); pit != partitions.end(); pit++) {
        if (part_names.find((*pit)->name) != part_names.end())
            outError("Duplicated partition name ", (*pit)->name);
        part_names.insert((*pit)->name);
    }

    // Initialize the counter for evaluated NNIs on subtrees
    cout << "Subset\tType\tSeqs\tSites\tInfor\tInvar\tModel\tName" << endl;
    int part = 0;
    for (auto it = partitions.begin(); it != partitions.end(); it++, part++) {
        cout << part+1 << "\t" << (*it)->sequence_type << "\t" << (*it)->getNSeq()
        << "\t" << (*it)->getNSite() << "\t" << (*it)->num_informative_sites
        << "\t" << (*it)->getNSite()-(*it)->num_variant_sites << "\t"
        << (*it)->model_name << "\t" << (*it)->name << endl;
        if ((*it)->num_informative_sites == 0) {
            outWarning("No parsimony-informative sites in partition " + (*it)->name);
        }
    }
    init();

    if (params.print_conaln) {
        string str = params.out_prefix;
        str = params.out_prefix;
        str += ".conaln";
        printCombinedAlignment(str.c_str());
    }
    
    cout << "Degree of missing data: " << computeMissingData() << endl;
    
#ifdef _OPENMP
    if (params.num_threads > partitions.size()) {
        cout << "Info: multi-threading strategy over alignment sites" << endl;
    } else {
        cout << "Info: multi-threading strategy over partitions" << endl;
    }
#endif
    cout << endl;

}

void SuperAlignment::init(StrVector *sequence_names) {
    // start original code
    
    max_num_states = 0;
	// first build taxa_index and partitions
	int site, seq, nsite = partitions.size();

    // BUG FIX 2016-11-29: when merging partitions with -m TESTMERGE, sequence order is changed
    // get the taxa names from existing tree
    if (sequence_names && !sequence_names->empty()) {
        seq_names = *sequence_names;
        taxa_index.resize(seq_names.size());
        for (auto i = taxa_index.begin(); i != taxa_index.end(); i++)
            i->resize(nsite, -1);
    }

    site = 0;
	for (auto it = partitions.begin(); it != partitions.end(); it++, site++) {
//        partitions.push_back((*it)->aln);
		int nseq = (*it)->getNSeq();
		//cout << "nseq  = " << nseq << endl;
		for (seq = 0; seq < nseq; seq++) {
			int id = getSeqID((*it)->getSeqName(seq));
			if (id < 0) {
				seq_names.push_back((*it)->getSeqName(seq));
				id = seq_names.size()-1;
				IntVector vec(nsite, -1);
				vec[site] = seq;
				taxa_index.push_back(vec);
			} else
				taxa_index[id][site] = seq;
		}
	}
	// now the patterns of sequence-genes presence/absence
	buildPattern();
}

void SuperAlignment::buildPattern() {
	int site, seq, nsite = partitions.size();

	seq_type = SEQ_BINARY;
	num_states = 2; // binary type because the super alignment presents the presence/absence of taxa in the partitions
	STATE_UNKNOWN = 2;
	site_pattern.resize(nsite, -1);
	clear();
	pattern_index.clear();
	VerboseMode save_mode = verbose_mode; 
	verbose_mode = min(verbose_mode, VB_MIN); // to avoid printing gappy sites in addPattern
	int nseq = getNSeq();
	for (site = 0; site < nsite; site++) {
 		Pattern pat;
 		pat.resize(nseq, 0);
		for (seq = 0; seq < nseq; seq++)
			pat[seq] = (taxa_index[seq][site] >= 0)? 1 : 0;
		addPattern(pat, site);
	}
	verbose_mode = save_mode;
	countConstSite();
    buildSeqStates();
}

void SuperAlignment::readPartition(Params &params) {
    try {
        ifstream in;
        in.exceptions(ios::failbit | ios::badbit);
        in.open(params.partition_file);
        in.exceptions(ios::badbit);
        
        while (!in.eof()) {
            CharSet info;
            getline(in, info.name, ',');
            if (in.eof()) break;
            getline(in, info.model_name, ',');
            if (model_name == "") info.model_name = params.model_name;
            getline(in, info.aln_file, ',');
            if (info.aln_file == "" && params.aln_file) info.aln_file = params.aln_file;
            getline(in, info.sequence_type, ',');
            if (info.sequence_type=="" && params.sequence_type)
                info.sequence_type = params.sequence_type;
            safeGetline(in, info.position_spec);
            trimString(info.sequence_type);
            //            cout << endl << "Reading partition " << info.name << " (model=" << info.model_name << ", aln=" <<
            //                    info.aln_file << ", seq=" << info.sequence_type << ", pos=" << ((info.position_spec.length() >= 20) ? info.position_spec.substr(0,20)+"..." : info.position_spec) << ") ..." << endl;

            // TODO move this to supertree
//            info.nniMoves[0].ptnlh = NULL;
//            info.nniMoves[1].ptnlh = NULL;
//            info.cur_ptnlh = NULL;
//            part_info.push_back(info);
            Alignment *part_aln = new Alignment((char*)info.aln_file.c_str(), (char*)info.sequence_type.c_str(), params.intype, info.model_name);
            if (!info.position_spec.empty()) {
                Alignment *new_aln = new Alignment();
                new_aln->extractSites(part_aln, info.position_spec.c_str());
                delete part_aln;
                part_aln = new_aln;
            }
            part_aln->name = info.name;
            part_aln->model_name = info.model_name;
            part_aln->position_spec = info.position_spec;
            part_aln->aln_file = info.aln_file;
            part_aln->sequence_type = info.sequence_type;
            partitions.push_back(part_aln);
            // TODO move this to supertree
//            PhyloTree *tree = new PhyloTree(part_aln);
//            push_back(tree);
        }
        
        in.clear();
        // set the failbit again
        in.exceptions(ios::failbit | ios::badbit);
        in.close();
    } catch(ios::failure) {
        outError(ERR_READ_INPUT);
    } catch (string str) {
        outError(str);
    }
    
    
}

void SuperAlignment::readPartitionRaxml(Params &params) {
    try {
        ifstream in;
        in.exceptions(ios::failbit | ios::badbit);
        in.open(params.partition_file);
        in.exceptions(ios::badbit);
//        PartitionInfo info;
        Alignment *input_aln = NULL;
        if (!params.aln_file)
            outError("Please supply an alignment with -s option");
        
        input_aln = new Alignment(params.aln_file, params.sequence_type, params.intype, params.model_name);
        
        cout << endl << "Partition file is not in NEXUS format, assuming RAxML-style partition file..." << endl;
        
        size_t pos = params.model_name.find_first_of("+*");
        string rate_type = "";
        if (pos != string::npos) rate_type = params.model_name.substr(pos);
        
        while (!in.eof()) {
            CharSet info;
            getline(in, info.model_name, ',');
            if (in.eof()) break;
            trimString(info.model_name);
            //            std::transform(info.model_name.begin(), info.model_name.end(), info.model_name.begin(), ::toupper);
            
            bool is_ASC = info.model_name.substr(0,4) == "ASC_";
            if (is_ASC) info.model_name.erase(0, 4);
            StateFreqType freq = FREQ_UNKNOWN;
            if (info.model_name.find_first_of("*+{") == string::npos ) {
                if (*info.model_name.rbegin() == 'F' && info.model_name != "DAYHOFF") {
                    freq = FREQ_EMPIRICAL;
                    info.model_name.erase(info.model_name.length()-1);
                } else if (*info.model_name.rbegin() == 'X' && info.model_name != "LG4X") {
                    freq = FREQ_ESTIMATE;
                    info.model_name.erase(info.model_name.length()-1);
                }
            }
            
            if (info.model_name.empty())
                outError("Please give model names in partition file!");
            if (info.model_name == "BIN") {
                info.sequence_type = "BIN";
                info.model_name = "GTR2";
            } else if (info.model_name == "DNA") {
                info.sequence_type = "DNA";
                info.model_name = "GTR";
            } else if (info.model_name == "MULTI") {
                info.sequence_type = "MORPH";
                info.model_name = "MK";
            } else if (info.model_name.substr(0,5) == "CODON") {
                info.sequence_type = info.model_name;
                info.model_name = "GY";
            } else {
                info.sequence_type = "AA";
                if (*info.model_name.begin() == '[') {
                    if (*info.model_name.rbegin() != ']')
                        outError("User-defined protein model should be [myProtenSubstitutionModelFileName]");
                    info.model_name = info.model_name.substr(1, info.model_name.length()-2);
                }
            }
            
            if (freq == FREQ_EMPIRICAL)
                info.model_name += "+F";
            else if (freq == FREQ_ESTIMATE)
                info.model_name += "+FO";
            if (is_ASC)
                info.model_name += "+ASC";
            info.model_name += rate_type;
            
            getline(in, info.name, '=');
            trimString(info.name);
            if (info.name.empty())
                outError("Please give partition names in partition file!");
            
            safeGetline(in, info.position_spec);
            trimString(info.position_spec);
            if (info.position_spec.empty())
                outError("Please specify alignment positions for partition" + info.name);
            std::replace(info.position_spec.begin(), info.position_spec.end(), ',', ' ');
            
            //            cout << "Reading partition " << info.name << " (model=" << info.model_name << ", seq=" << info.sequence_type << ", pos=" << ((info.position_spec.length() >= 20) ? info.position_spec.substr(0,20)+"..." : info.position_spec) << ") ..." << endl;
            
            // TODO to supertree
//            info.nniMoves[0].ptnlh = NULL;
//            info.nniMoves[1].ptnlh = NULL;
//            info.cur_ptnlh = NULL;
//            part_info.push_back(info);
            Alignment *part_aln = new Alignment();
            part_aln->extractSites(input_aln, info.position_spec.c_str());
            
            Alignment *new_aln;
            if (params.remove_empty_seq)
                new_aln = part_aln->removeGappySeq();
            else
                new_aln = part_aln;
            // also rebuild states set of each sequence for likelihood computation
            new_aln->buildSeqStates();
            
            if (part_aln != new_aln) delete part_aln;

            new_aln->name = info.name;
            new_aln->model_name = info.model_name;
            new_aln->position_spec = info.position_spec;
            new_aln->aln_file = info.aln_file;
            new_aln->sequence_type = info.sequence_type;
            partitions.push_back(new_aln);
            // TODO move to supertree
//            PhyloTree *tree = new PhyloTree(new_aln);
//            push_back(tree);
            //            cout << new_aln->getNSeq() << " sequences and " << new_aln->getNSite() << " sites extracted" << endl;
            //            params = origin_params;
        }
        
        in.clear();
        // set the failbit again
        in.exceptions(ios::failbit | ios::badbit);
        in.close();
    } catch(ios::failure) {
        outError(ERR_READ_INPUT);
    } catch (string str) {
        outError(str);
    }
    
    
}

void SuperAlignment::readPartitionNexus(Params &params) {
//    Params origin_params = params;
    MSetsBlock *sets_block = new MSetsBlock();
    MyReader nexus(params.partition_file);
    nexus.Add(sets_block);
    MyToken token(nexus.inf);
    nexus.Execute(token);
    
    Alignment *input_aln = NULL;
    if (params.aln_file) {
        input_aln = new Alignment(params.aln_file, params.sequence_type, params.intype, params.model_name);
    }
    
    bool empty_partition = true;
    vector<CharSet*>::iterator it;
    for (it = sets_block->charsets.begin(); it != sets_block->charsets.end(); it++)
        if ((*it)->model_name != "") {
            empty_partition = false;
            break;
        }
    if (empty_partition) {
        cout << "NOTE: No CharPartition defined, use all CharSets" << endl;
    }
    
    cout << endl << "Loading " << sets_block->charsets.size() << " partitions..." << endl;
    
    for (it = sets_block->charsets.begin(); it != sets_block->charsets.end(); it++)
        if (empty_partition || (*it)->char_partition != "") {
            if ((*it)->model_name == "")
                (*it)->model_name = params.model_name;
            if ((*it)->aln_file == "" && !params.aln_file)
                outError("No input data for partition ", (*it)->name);
            if ((*it)->sequence_type=="" && params.sequence_type)
                (*it)->sequence_type = params.sequence_type;
            
            if ((*it)->sequence_type == "" && !(*it)->model_name.empty()) {
                // try to get sequence type from model
            //TODO: why compile error?
                (*it)->sequence_type = detectSeqTypeName((*it)->model_name.substr(0, (*it)->model_name.find_first_of("+*")));
            }
            if ((*it)->aln_file == "" && ((*it)->position_spec == "" || (*it)->position_spec == "*"))
                outError("Empty position range for partition ", (*it)->name);
            trimString((*it)->sequence_type);
            //            cout << endl << "Reading partition " << info.name << " (model=" << info.model_name << ", aln=" <<
            //                info.aln_file << ", seq=" << info.sequence_type << ", pos=" << ((info.position_spec.length() >= 20) ? info.position_spec.substr(0,20)+"..." : info.position_spec) << ") ..." << endl;
            if ((*it)->sequence_type != "" && Alignment::getSeqType((*it)->sequence_type.c_str()) == SEQ_UNKNOWN)
                outError("Unknown sequence type " + (*it)->sequence_type);

            // TODO move to supertree
//            info.nniMoves[0].ptnlh = NULL;
//            info.nniMoves[1].ptnlh = NULL;
//            info.cur_ptnlh = NULL;
//            part_info.push_back(info);
            Alignment *part_aln;
            if ((*it)->aln_file != "") {
                part_aln = new Alignment((char*)(*it)->aln_file.c_str(), (char*)(*it)->sequence_type.c_str(), params.intype, (*it)->model_name);
            } else {
                part_aln = input_aln;
            }
            if (!(*it)->position_spec.empty() && (*it)->position_spec != "*") {
                Alignment *new_aln = new Alignment();
                new_aln->extractSites(part_aln, (*it)->position_spec.c_str());
                if (part_aln != input_aln) delete part_aln;
                part_aln = new_aln;
            }
            if (part_aln->seq_type == SEQ_DNA && ((*it)->sequence_type.substr(0, 5) == "CODON" || (*it)->sequence_type.substr(0, 5) == "NT2AA")) {
                Alignment *new_aln = new Alignment();
                new_aln->convertToCodonOrAA(part_aln, &(*it)->sequence_type[5], (*it)->sequence_type.substr(0, 5) == "NT2AA");
                if (part_aln != input_aln) delete part_aln;
                part_aln = new_aln;
            }
            Alignment *new_aln;
            if (params.remove_empty_seq)
                new_aln = part_aln->removeGappySeq();
            else
                new_aln = part_aln;
            // also rebuild states set of each sequence for likelihood computation
            new_aln->buildSeqStates();
            
            if (part_aln != new_aln && part_aln != input_aln) delete part_aln;
            new_aln->name = (*it)->name;
            new_aln->model_name = (*it)->model_name;
            new_aln->aln_file = (*it)->aln_file;
            new_aln->position_spec = (*it)->position_spec;
            new_aln->sequence_type = (*it)->sequence_type;
            partitions.push_back(new_aln);
//            PhyloTree *tree = new PhyloTree(new_aln);
//            push_back(tree);
//            params = origin_params;
            //            cout << new_aln->getNSeq() << " sequences and " << new_aln->getNSite() << " sites extracted" << endl;
        }
    
    if (input_aln)
        delete input_aln;
    delete sets_block;
}

void SuperAlignment::printPartition(const char *filename) {
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(filename);
        out << "#nexus" << endl << "[ partition information for alignment written in .conaln file ]" << endl
        << "begin sets;" << endl;
        int part; int start_site;
        for (part = 0, start_site = 1; part < partitions.size(); part++) {
            string name = partitions[part]->name;
            replace(name.begin(), name.end(), '+', '_');
            int end_site = start_site + partitions[part]->getNSite();
            out << "  charset " << name << " = " << start_site << "-" << end_site-1 << ";" << endl;
            start_site = end_site;
        }
        out << "  charpartition mymodels =" << endl;
        for (part = 0; part < partitions.size(); part++) {
            string name = partitions[part]->name;
            replace(name.begin(), name.end(), '+', '_');
            if (part > 0) out << "," << endl;
//            out << "    " << at(part)->getModelNameParams() << ":" << name;
            out << "    " << partitions[part]->model_name << ":" << name;
        }
        out << ";" << endl;
        out << "end;" << endl;
        out.close();
        cout << "Partition information was printed to " << filename << endl;
    } catch (ios::failure &) {
        outError(ERR_WRITE_OUTPUT, filename);
    }
    
}

void SuperAlignment::printBestPartition(const char *filename) {
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(filename);
        out << "#nexus" << endl
        << "begin sets;" << endl;
        int part;
        for (part = 0; part < partitions.size(); part++) {
            string name = partitions[part]->name;
            replace(name.begin(), name.end(), '+', '_');
            out << "  charset " << name << " = ";
            if (!partitions[part]->aln_file.empty()) out << partitions[part]->aln_file << ": ";
            if (partitions[part]->seq_type == SEQ_CODON)
                out << "CODON, ";
            string pos = partitions[part]->position_spec;
            replace(pos.begin(), pos.end(), ',' , ' ');
            out << pos << ";" << endl;
        }
        out << "  charpartition mymodels =" << endl;
        for (part = 0; part < partitions.size(); part++) {
            string name = partitions[part]->name;
            replace(name.begin(), name.end(), '+', '_');
            if (part > 0) out << "," << endl;
            out << "    " << partitions[part]->model_name << ": " << name;
        }
        out << ";" << endl;
        out << "end;" << endl;
        out.close();
        cout << "Partition information was printed to " << filename << endl;
    } catch (ios::failure &) {
        outError(ERR_WRITE_OUTPUT, filename);
    }
    
}


void SuperAlignment::printPartitionRaxml(const char *filename) {
    int part;
//    for (part = 0; part < partitions.size(); part++) {
//        if (partitions[part]->aln_file != "") {
//            cout << "INFO: Printing partition in RAxML format is not possible" << endl;
//            return;
//        }
//    }
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(filename);
        int start_site;
        for (part = 0, start_site = 1; part < partitions.size(); part++) {
            string name = partitions[part]->name;
            replace(name.begin(), name.end(), '+', '_');
            int end_site = start_site + partitions[part]->getNSite();
            switch (partitions[part]->seq_type) {
                case SEQ_DNA: out << "DNA, "; break;
                case SEQ_BINARY: out << "BIN, "; break;
                case SEQ_MORPH: out << "MULTI, "; break;
                default: out << partitions[part]->model_name << ","; break;
            }
            out << name << " = " << start_site << "-" << end_site-1 << endl;
            start_site = end_site;
        }
        out.close();
        cout << "Partition information in Raxml format was printed to " << filename << endl;
    } catch (ios::failure &) {
        outError(ERR_WRITE_OUTPUT, filename);
    }
    
}

void SuperAlignment::printBestPartitionRaxml(const char *filename) {
    int part;
//    for (part = 0; part < partitions.size(); part++) {
//        if (partitions[part]->aln_file != "") {
//            cout << "INFO: Printing partition in RAxML format is not possible" << endl;
//            return;
//        }
//    }
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(filename);
        for (part = 0; part < partitions.size(); part++) {
            string name = partitions[part]->name;
            replace(name.begin(), name.end(), '+', '_');
            if (partitions[part]->model_name.find("+ASC") != string::npos)
                out << "ASC_";
            switch (partitions[part]->seq_type) {
                case SEQ_DNA: out << "DNA"; break;
                case SEQ_BINARY: out << "BIN"; break;
                case SEQ_MORPH: out << "MULTI"; break;
                case SEQ_PROTEIN:
                    out << partitions[part]->model_name.substr(0, partitions[part]->model_name.find_first_of("*{+"));
                    break;
                case SEQ_CODON:
                    out << "CODON_" << partitions[part]->model_name.substr(0, partitions[part]->model_name.find_first_of("*{+"));
                    break;
                default: out << partitions[part]->model_name; break;
            }
            if (partitions[part]->model_name.find("+FO") != string::npos)
                out << "X";
            else if (partitions[part]->model_name.find("+F") != string::npos)
                out << "F";
            
            out << ", " << name << " = " << partitions[part]->position_spec << endl;
        }
        out.close();
        cout << "Partition information in Raxml format was printed to " << filename << endl;
    } catch (ios::failure &) {
        outError(ERR_WRITE_OUTPUT, filename);
    }
    
}


void SuperAlignment::linkSubAlignment(int part) {
	ASSERT(taxa_index.size() == getNSeq());
	int nseq = getNSeq(), seq;
	vector<bool> checked;
	checked.resize(partitions[part]->getNSeq(), false);
	for (seq = 0; seq < nseq; seq++) {
		int id = partitions[part]->getSeqID(getSeqName(seq));
		if (id < 0)
			taxa_index[seq][part] = -1;
		else {
			taxa_index[seq][part] = id;
			checked[id] = true;
		}
	}
	if (verbose_mode >= VB_MED) {

	}
	// sanity check that all seqnames in partition must be present in superalignment
	for (seq = 0; seq < checked.size(); seq++) {
		ASSERT(checked[seq]);
	}
}

void SuperAlignment::extractSubAlignment(Alignment *aln, IntVector &seq_id, int min_true_char, int min_taxa, IntVector *kept_partitions) {
	ASSERT(aln->isSuperAlignment());
	SuperAlignment *saln = (SuperAlignment*)aln;
    name = aln->name;
    model_name = aln->model_name;
    sequence_type = aln->sequence_type;
    position_spec = aln->position_spec;
    aln_file = aln->aln_file;

    int i;
    IntVector::iterator it;
    for (it = seq_id.begin(); it != seq_id.end(); it++) {
        ASSERT(*it >= 0 && *it < aln->getNSeq());
        seq_names.push_back(aln->getSeqName(*it));
    }

	// BUG HERE!
	//Alignment::extractSubAlignment(aln, seq_id, 0);

	taxa_index.resize(getNSeq());
	for (i = 0; i < getNSeq(); i++)
		taxa_index[i].resize(saln->partitions.size(), -1);

	int part = 0;
//	partitions.resize(saln->partitions.size());
    partitions.resize(0);
	for (vector<Alignment*>::iterator ait = saln->partitions.begin(); ait != saln->partitions.end(); ait++, part++) {
		IntVector sub_seq_id;
		for (IntVector::iterator it = seq_id.begin(); it != seq_id.end(); it++)
			if (saln->taxa_index[*it][part] >= 0)
				sub_seq_id.push_back(saln->taxa_index[*it][part]);
        if (sub_seq_id.size() < min_taxa)
            continue;
		Alignment *subaln = new Alignment;
		subaln->extractSubAlignment(*ait, sub_seq_id, 0);
		partitions.push_back(subaln);
		linkSubAlignment(partitions.size()-1);
        if (kept_partitions) kept_partitions->push_back(part);
//		cout << subaln->getNSeq() << endl;
//		subaln->printPhylip(cout);
	}

    if (partitions.size() < saln->partitions.size()) {
        for (i = 0; i < getNSeq(); i++)
            taxa_index[i].resize(partitions.size());
    }

	// now build the patterns based on taxa_index
	buildPattern();
}

Alignment *SuperAlignment::removeIdenticalSeq(string not_remove, bool keep_two, StrVector &removed_seqs, StrVector &target_seqs) {
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
			int part = 0;
			// check if seq1 and seq2 are identical over all partitions
			for (vector<Alignment*>::iterator ait = partitions.begin(); ait != partitions.end(); ait++, part++) {
				int subseq1 = taxa_index[seq1][part];
				int subseq2 = taxa_index[seq2][part];
				if (subseq1 < 0 && subseq2 < 0) // continue if both seqs are absent in this partition
					continue;
				if (subseq1 < 0 && subseq2 > 0) {
					// if one sequence is present and the other is absent for a gene, we conclude that they are not identical
					equal_seq = false;
					break;
				}
				if (subseq1 > 0 && subseq2 < 0) {
					// if one sequence is present and the other is absent for a gene, we conclude that they are not identical
					equal_seq = false;
					break;
				}
				// now if both seqs are present, check sequence content
				for (iterator it = (*ait)->begin(); it != (*ait)->end(); it++)
					if  ((*it)[subseq1] != (*it)[subseq2]) {
						equal_seq = false;
						break;
					}
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

	if (removed_seqs.empty()) return this; // do nothing if the list is empty

    if (removed_seqs.size() >= getNSeq()-3)
        outWarning("Your alignment contains too many identical sequences!");

	// now remove identical sequences
	IntVector keep_seqs;
	for (seq1 = 0; seq1 < getNSeq(); seq1++)
		if (!removed[seq1]) keep_seqs.push_back(seq1);
	SuperAlignment *aln;
	aln = new SuperAlignment;
	aln->extractSubAlignment(this, keep_seqs, 0);
	return aln;
}

int SuperAlignment::checkAbsentStates(string msg) {
    int count = 0;
    for (auto it = partitions.begin(); it != partitions.end(); it++)
        count += (*it)->checkAbsentStates("partition " + convertIntToString((it-partitions.begin())+1));
    return count;
}

/*
void SuperAlignment::checkGappySeq() {
	int nseq = getNSeq(), part = 0, i;
	IntVector gap_only_seq;
	gap_only_seq.resize(nseq, 1);
	//cout << "Checking gaps..." << endl;
	for (vector<Alignment*>::iterator it = partitions.begin(); it != partitions.end(); it++, part++) {
		IntVector keep_seqs;
		for (i = 0; i < nseq; i++)
			if (taxa_index[i][part] >= 0)
			if (!(*it)->isGapOnlySeq(taxa_index[i][part])) {
				keep_seqs.push_back(taxa_index[i][part]);
				gap_only_seq[i] = 0;
			}
		if (keep_seqs.size() < (*it)->getNSeq()) {
			cout << "Discard " << (*it)->getNSeq() - keep_seqs.size() 
				 << " sequences from partition number " << part+1 << endl;
			Alignment *aln = new Alignment;
			aln->extractSubAlignment((*it), keep_seqs, 0);
			delete (*it);
			(*it) = aln;
			linkSubAlignment(part);
		}
		cout << __func__ << " num_states = " << (*it)->num_states << endl;
	}
	int wrong_seq = 0;
	for (i = 0; i < nseq; i++)
		if (gap_only_seq[i]) {
			cout << "ERROR: Sequence " << getSeqName(i) << " contains only gaps or missing data" << endl;
			wrong_seq++;
		}
	if (wrong_seq) {
		outError("Some sequences (see above) are problematic, please check your alignment again");
		}
}
*/
void SuperAlignment::getSitePatternIndex(IntVector &pattern_index) {
	int nptn = 0;
	for (vector<Alignment*>::iterator it = partitions.begin(); it != partitions.end(); it++) {
		int nsite = pattern_index.size();
		pattern_index.insert(pattern_index.end(), (*it)->site_pattern.begin(), (*it)->site_pattern.end());
		for (int i = nsite; i < pattern_index.size(); i++)
			pattern_index[i] += nptn;
		nptn += (*it)->getNPattern();
	}
}

void SuperAlignment::getPatternFreq(IntVector &pattern_freq) {
	ASSERT(isSuperAlignment());
	size_t offset = 0;
	if (!pattern_freq.empty()) pattern_freq.resize(0);
	for (vector<Alignment*>::iterator it = partitions.begin(); it != partitions.end(); it++) {
		IntVector freq;
		(*it)->getPatternFreq(freq);
		pattern_freq.insert(pattern_freq.end(), freq.begin(), freq.end());
		offset += freq.size();
	}
}

void SuperAlignment::getPatternFreq(int *pattern_freq) {
    ASSERT(isSuperAlignment());
    size_t offset = 0;
    for (vector<Alignment*>::iterator it = partitions.begin(); it != partitions.end(); it++) {
        (*it)->getPatternFreq(pattern_freq + offset);
        offset += (*it)->getNPattern();
    }
}

void SuperAlignment::printSiteInfo(const char* filename) {
    try {
        ofstream out(filename);
        printSiteInfoHeader(out, filename, true);
        int id = 1;
        for (auto it = partitions.begin(); it != partitions.end(); it++, id++)
            (*it)->printSiteInfo(out, id);
        out.close();
    } catch (...) {
        outError(ERR_WRITE_OUTPUT, filename);
    }
}

/*
void SuperAlignment::createBootstrapAlignment(Alignment *aln, IntVector* pattern_freq, const char *spec) {
	ASSERT(aln->isSuperAlignment());
	Alignment::copyAlignment(aln);
	SuperAlignment *super_aln = (SuperAlignment*) aln;
	ASSERT(partitions.empty());
    name = aln->name;
    model_name = aln->model_name;
    sequence_type = aln->sequence_type;
    position_spec = aln->position_spec;
    aln_file = aln->aln_file;

	if (spec && strncmp(spec, "GENE", 4) == 0) {
		// resampling whole genes
        partitions.resize(super_aln->partitions.size(), NULL);
        int i, ptn;
        for (i = 0; i < super_aln->partitions.size(); i++) {

            // get a random gene
			int part = random_int(super_aln->partitions.size());

            // ptn_freq stores pattern frequency of bootstrap aln

            IntVector ptn_freq;
            if (strncmp(spec,"GENESITE",8) == 0) {
                // resample sites of this gene
                super_aln->partitions[part]->createBootstrapAlignment(ptn_freq);
                ASSERT(ptn_freq.size() == super_aln->partitions[part]->size());
            } else {
                // copy ptn_freq from this gene
                for (ptn = 0; ptn < super_aln->partitions[part]->size(); ptn++)
                    ptn_freq.push_back(super_aln->partitions[part]->at(ptn).frequency);
            }

            if (!partitions[part]) {
                // allocate the partition
                partitions[part] = new Alignment;
                partitions[part]->copyAlignment(super_aln->partitions[part]);
                for (ptn = 0; ptn < super_aln->partitions[part]->size(); ptn++)
                    partitions[part]->at(ptn).frequency = ptn_freq[ptn];
            } else {
                // increase frequency if already existed
                for (ptn = 0; ptn < super_aln->partitions[part]->size(); ptn++)
                    partitions[part]->at(ptn).frequency += ptn_freq[ptn];
            }
        }

        // fulfill genes that are missing
        for (i = 0; i < partitions.size(); i++)
            if (!partitions[i]) {
                partitions[i] = new Alignment;
                partitions[i]->copyAlignment(super_aln->partitions[i]);
                // reset all frequency
                for (ptn = 0; ptn < partitions[i]->size(); ptn++)
                    partitions[i]->at(ptn).frequency = 0;
            }

        // fill up pattern_freq vector
        if (pattern_freq) {
            pattern_freq->resize(0);
            for (i = 0; i < partitions.size(); i++)
                for (ptn = 0; ptn < partitions[i]->size(); ptn++)
                    pattern_freq->push_back(partitions[i]->at(ptn).frequency);
        }
    } else if (!spec) {
		// resampling sites within genes
        for (vector<Alignment*>::iterator it = super_aln->partitions.begin(); it != super_aln->partitions.end(); it++) {
            Alignment *boot_aln = new Alignment;
            if (pattern_freq) {
                IntVector part_pattern_freq;
                boot_aln->createBootstrapAlignment(*it, &part_pattern_freq);
                pattern_freq->insert(pattern_freq->end(), part_pattern_freq.begin(), part_pattern_freq.end());
            } else {
                boot_aln->createBootstrapAlignment(*it);
            }
            partitions.push_back(boot_aln);
        }
    } else {
        outError("Wrong -bsam, either -bsam GENE or -bsam GENESITE");
    }
	taxa_index = super_aln->taxa_index;
    countConstSite();
}
*/

void SuperAlignment::createBootstrapAlignment(Alignment *aln, IntVector* pattern_freq, const char *spec) {
    ASSERT(aln->isSuperAlignment());
    SuperAlignment *super_aln = (SuperAlignment*) aln;
    ASSERT(partitions.empty());
    name = aln->name;
    model_name = aln->model_name;
    sequence_type = aln->sequence_type;
    position_spec = aln->position_spec;
    aln_file = aln->aln_file;
    
    if (!spec) {
        // resampling sites within genes
        Alignment::copyAlignment(aln);
        partitions.reserve(super_aln->partitions.size());
        for (vector<Alignment*>::iterator it = super_aln->partitions.begin(); it != super_aln->partitions.end(); it++) {
            Alignment *boot_aln = new Alignment;
            if (pattern_freq) {
                IntVector part_pattern_freq;
                boot_aln->createBootstrapAlignment(*it, &part_pattern_freq);
                pattern_freq->insert(pattern_freq->end(), part_pattern_freq.begin(), part_pattern_freq.end());
            } else {
                boot_aln->createBootstrapAlignment(*it);
            }
            partitions.push_back(boot_aln);
        }
        taxa_index = super_aln->taxa_index;
        countConstSite();
    } else if (strcmp(spec, "GENE") == 0) {
        ASSERT(!pattern_freq);
        // resampling whole genes
        IntVector gene_freq;
        random_resampling(super_aln->partitions.size(), gene_freq);
        for (int i = 0; i < gene_freq.size(); i++)
            if (gene_freq[i] > 0) {
                Alignment *boot_aln = new Alignment;
                boot_aln->copyAlignment(super_aln->partitions[i]);
                if (gene_freq[i] > 1) {
                    for (auto it = boot_aln->begin(); it != boot_aln->end(); it++)
                        it->frequency *= gene_freq[i];
                    auto site_pattern = boot_aln->site_pattern;
                    for (int j = 1; j < gene_freq[i]; j++)
                        boot_aln->site_pattern.insert(boot_aln->site_pattern.end(), site_pattern.begin(), site_pattern.end());
                    boot_aln->countConstSite();
                }
                partitions.push_back(boot_aln);
            }
        init();
    } else if (strcmp(spec, "GENESITE") == 0) {
        ASSERT(!pattern_freq);
        // resampling whole genes then sites within resampled genes
        IntVector gene_freq;
        random_resampling(super_aln->partitions.size(), gene_freq);
        for (int i = 0; i < gene_freq.size(); i++)
            for (int rep = 0; rep < gene_freq[i]; rep++) {
            Alignment *boot_aln = new Alignment;
            boot_aln->createBootstrapAlignment(super_aln->partitions[i]);
            boot_aln->name = boot_aln->name + "." + convertIntToString(rep);
            partitions.push_back(boot_aln);
        }
        init();
    } else {
        outError("Wrong -bsam, either -bsam GENE or -bsam GENESITE");
    }
}

void SuperAlignment::createBootstrapAlignment(IntVector &pattern_freq, const char *spec) {
	ASSERT(isSuperAlignment());
	int nptn = 0;
	for (vector<Alignment*>::iterator it = partitions.begin(); it != partitions.end(); it++) {
		nptn += (*it)->getNPattern();
	}
	pattern_freq.resize(0);
	int *internal_freq = new int[nptn];
	createBootstrapAlignment(internal_freq, spec);
	pattern_freq.insert(pattern_freq.end(), internal_freq, internal_freq + nptn);
	delete [] internal_freq;

}


void SuperAlignment::createBootstrapAlignment(int *pattern_freq, const char *spec, int *rstream) {
	ASSERT(isSuperAlignment());
//	if (spec && strncmp(spec, "GENE", 4) != 0) outError("Unsupported yet. ", __func__);

	if (spec && strncmp(spec, "GENE", 4) == 0) {
		// resampling whole genes
		int nptn = 0;
		IntVector part_pos;
		for (vector<Alignment*>::iterator it = partitions.begin(); it != partitions.end(); it++) {
			part_pos.push_back(nptn);
			nptn += (*it)->getNPattern();
		}
		memset(pattern_freq, 0, nptn * sizeof(int));
        IntVector gene_freq;
        random_resampling(partitions.size(), gene_freq, rstream);
		for (int part = 0; part < partitions.size(); part++)
        for (int rep = 0; rep < gene_freq[part]; rep++){
			Alignment *aln = partitions[part];
			if (strncmp(spec,"GENESITE",8) == 0) {
				// then resampling sites in resampled gene
                IntVector sample;
                random_resampling(aln->getNSite(), sample, rstream);
				for (int site = 0; site < sample.size(); site++)
                for (int rep2 = 0; rep2 < sample[site]; rep2++) {
					int ptn_id = aln->getPatternID(site);
					pattern_freq[ptn_id + part_pos[part]]++;
				}
			} else {
				for (int ptn = 0; ptn < aln->getNPattern(); ptn++)
					pattern_freq[ptn + part_pos[part]] += aln->at(ptn).frequency;
			}
		}
	} else {
		// resampling sites within genes
		int offset = 0;
		for (vector<Alignment*>::iterator it = partitions.begin(); it != partitions.end(); it++) {
            if (spec && strncmp(spec, "SCALE=", 6) == 0)
                (*it)->createBootstrapAlignment(pattern_freq + offset, spec, rstream);
            else
                (*it)->createBootstrapAlignment(pattern_freq + offset, NULL, rstream);
			offset += (*it)->getNPattern();
		}
	}
}

/**
 * shuffle alignment by randomizing the order of sites
 */
void SuperAlignment::shuffleAlignment() {
	ASSERT(isSuperAlignment());
	for (vector<Alignment*>::iterator it = partitions.begin(); it != partitions.end(); it++) {
		(*it)->shuffleAlignment();
	}
}


double SuperAlignment::computeObsDist(int seq1, int seq2) {
	int site;
	int diff_pos = 0, total_pos = 0;
	for (site = 0; site < getNSite(); site++) {
		int id1 = taxa_index[seq1][site];
		int id2 = taxa_index[seq2][site];
		if (id1 < 0 || id2 < 0) continue;
		int num_states = partitions[site]->num_states;
		for (Alignment::iterator it = partitions[site]->begin(); it != partitions[site]->end(); it++) 
			if  ((*it)[id1] < num_states && (*it)[id2] < num_states) {
				total_pos += (*it).frequency;
				if ((*it)[id1] != (*it)[id2] )
					diff_pos += (*it).frequency;
			}
	}
	if (!total_pos)
		return MAX_GENETIC_DIST; // return +INF if no overlap between two sequences
	return ((double)diff_pos) / total_pos;
}


double SuperAlignment::computeDist(int seq1, int seq2) {
	if (partitions.empty()) return 0.0;
	double obs_dist = computeObsDist(seq1, seq2);
    int num_states = partitions[0]->num_states;
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
    //return computeObsDist(seq1, seq2);
	//  AVERAGE DISTANCE

	double dist = 0;
	int part = 0, num = 0;
	for (vector<Alignment*>::iterator it = partitions.begin(); it != partitions.end(); it++, part++) {
		int id1 = taxa_index[seq1][part];
		int id2 = taxa_index[seq2][part];
		if (id1 < 0 || id2 < 0) continue;
		dist += (*it)->computeDist(id1, id2);
	}
	if (num == 0) // two sequences are not overlapping at all!
		return MAX_GENETIC_DIST;
	return dist / num;
}

SuperAlignment::~SuperAlignment()
{
	for (vector<Alignment*>::reverse_iterator it = partitions.rbegin(); it != partitions.rend(); it++)
		delete (*it);
	partitions.clear();
}

void SuperAlignment::printCombinedAlignment(ostream &out, bool print_taxid) {
	vector<Alignment*>::iterator pit;
	int final_length = 0;
	for (pit = partitions.begin(); pit != partitions.end(); pit++)
        if ((*pit)->seq_type == SEQ_CODON)
            final_length += 3*(*pit)->getNSite();
        else
            final_length += (*pit)->getNSite();

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
		int part = 0;
		for (pit = partitions.begin(); pit != partitions.end(); pit++, part++) {
			int part_seq_id = taxa_index[seq_id][part];
			int nsite = (*pit)->getNSite();
			if (part_seq_id >= 0) {
				for (int i = 0; i < nsite; i++)
					out << (*pit)->convertStateBackStr((*pit)->getPattern(i) [part_seq_id]);
			} else {
				string str(nsite, '?');
				out << str;
			}
		}
		out << endl;
	}
}

void SuperAlignment::printCombinedAlignment(const char *file_name, bool append) {
	try {
		ofstream out;
		out.exceptions(ios::failbit | ios::badbit);

		if (append)
			out.open(file_name, ios_base::out | ios_base::app);
		else
			out.open(file_name);
        printCombinedAlignment(out);
		out.close();
		cout << "Concatenated alignment was printed to " << file_name << endl;
	} catch (ios::failure) {
		outError(ERR_WRITE_OUTPUT, file_name);
	}	
}

void SuperAlignment::printSubAlignments(Params &params) {
	vector<Alignment*>::iterator pit;
	string filename;
	int part;
	for (pit = partitions.begin(), part = 0; pit != partitions.end(); pit++, part++) {
		if (params.aln_output)
			filename = params.aln_output;
		else
			filename = params.out_prefix;
		filename += "." + (*pit)->name;
		 if (params.aln_output_format == ALN_PHYLIP)
			(*pit)->printPhylip(filename.c_str(), false, NULL, params.aln_nogaps, false, NULL);
		else if (params.aln_output_format == ALN_FASTA)
			(*pit)->printFasta(filename.c_str(), false, NULL, params.aln_nogaps, false, NULL);
	}
}

double SuperAlignment::computeUnconstrainedLogL() {
	double logl = 0.0;
	vector<Alignment*>::iterator pit;
	for (pit = partitions.begin(); pit != partitions.end(); pit++)
		logl += (*pit)->computeUnconstrainedLogL();
	return logl;
}

double SuperAlignment::computeMissingData() {
	double ret = 0.0;
	int len = 0;
	vector<Alignment*>::iterator pit;
	for (pit = partitions.begin(); pit != partitions.end(); pit++) {
		ret += (*pit)->getNSeq() * (*pit)->getNSite();
		len += (*pit)->getNSite();
	}
	ret /= getNSeq() * len;
	return 1.0 - ret;

}

Alignment *SuperAlignment::concatenateAlignments(set<int> &ids) {
	string union_taxa;
	int nsites = 0, nstates = 0;
    set<int>::iterator it;
	SeqType sub_type = SEQ_UNKNOWN;
	for (it = ids.begin(); it != ids.end(); it++) {
		int id = *it;
		ASSERT(id >= 0 && id < partitions.size());
		if (nstates == 0) nstates = partitions[id]->num_states;
		if (sub_type == SEQ_UNKNOWN) sub_type = partitions[id]->seq_type;
		if (sub_type != partitions[id]->seq_type)
			outError("Cannot concatenate sub-alignments of different type");
		if (nstates != partitions[id]->num_states)
			outError("Cannot concatenate sub-alignments of different #states");

		string taxa_set;
        Pattern taxa_pat = getPattern(id);
        taxa_set.insert(taxa_set.begin(), taxa_pat.begin(), taxa_pat.end());
		nsites += partitions[id]->getNSite();
		if (it == ids.begin()) union_taxa = taxa_set; else {
			for (int j = 0; j < union_taxa.length(); j++)
				if (taxa_set[j] == 1) union_taxa[j] = 1;
		}
	}

	Alignment *aln = new Alignment;
	for (int i = 0; i < union_taxa.length(); i++)
		if (union_taxa[i] == 1) {
			aln->seq_names.push_back(getSeqName(i));
		}
	aln->num_states = nstates;
	aln->seq_type = sub_type;
	aln->site_pattern.resize(nsites, -1);
    aln->clear();
    aln->pattern_index.clear();
    aln->STATE_UNKNOWN = partitions[*ids.begin()]->STATE_UNKNOWN;
    aln->genetic_code = partitions[*ids.begin()]->genetic_code;
    if (aln->seq_type == SEQ_CODON) {
    	aln->codon_table = new char[aln->num_states];
    	memcpy(aln->codon_table, partitions[*ids.begin()]->codon_table, aln->num_states);
    	aln->non_stop_codon = new char[strlen(aln->genetic_code)];
    	memcpy(aln->non_stop_codon, partitions[*ids.begin()]->non_stop_codon, strlen(aln->genetic_code));
    }

    int site = 0;
    for (it = ids.begin(); it != ids.end(); it++) {
    	int id = *it;
        // 2018-08-23: important bugfix in v1.6: taxa_set has wrong correspondance
		//string taxa_set;
        //Pattern taxa_pat = getPattern(id);
        //taxa_set.insert(taxa_set.begin(), taxa_pat.begin(), taxa_pat.end());
    	for (Alignment::iterator it = partitions[id]->begin(); it != partitions[id]->end(); it++) {
    		Pattern pat;
    		//int part_seq = 0;
    		for (int seq = 0; seq < union_taxa.size(); seq++)
    			if (union_taxa[seq] == 1) {
    				char ch = aln->STATE_UNKNOWN;
                    int seq_part = taxa_index[seq][id];
                    if (seq_part >= 0)
                        ch = (*it)[seq_part];
                    //if (taxa_set[seq] == 1) {
                    //    ch = (*it)[part_seq++];
                    //}
    				pat.push_back(ch);
    			}
    		//ASSERT(part_seq == partitions[id]->getNSeq());
    		aln->addPattern(pat, site, (*it).frequency);
    		// IMPORTANT BUG FIX FOLLOW
    		int ptnindex = aln->pattern_index[pat];
            for (int j = 0; j < (*it).frequency; j++)
                aln->site_pattern[site++] = ptnindex;

    	}
    }
    aln->countConstSite();
    aln->buildSeqStates();

	return aln;
}

Alignment *SuperAlignment::concatenateAlignments() {
    vector<SeqType> seq_types;
    vector<char*> genetic_codes;
    vector<set<int> > ids;
    for (int i = 0; i < partitions.size(); i++) {
        bool found = false;
        for (int j = 0; j < seq_types.size(); j++)
            if (partitions[i]->seq_type == seq_types[j] && partitions[i]->genetic_code == genetic_codes[j]) {
                ids[j].insert(i);
                found = true;
                break;
            }
        if (found)
            continue;
        // create a new partition
        seq_types.push_back(partitions[i]->seq_type);
        genetic_codes.push_back(partitions[i]->genetic_code);
        ids.push_back(set<int>());
        ids.back().insert(i);
    }
    if (seq_types.size() == 1)
        return concatenateAlignments(ids[0]);

    // mixed data with >= 2 partitions
    SuperAlignment *saln = new SuperAlignment();
    saln->max_num_states = 0;
    // first build taxa_index and partitions
    int site, seq, nsite = ids.size();
    
    // BUG FIX 2016-11-29: when merging partitions with -m TESTMERGE, sequence order is changed
    // get the taxa names from existing tree
    
    saln->seq_names = seq_names;
    saln->taxa_index.resize(saln->seq_names.size());
    for (auto it = saln->taxa_index.begin(); it != saln->taxa_index.end(); it++)
        it->resize(nsite, -1);
    
    for (site = 0; site != nsite; site++) {
        Alignment *part_aln = concatenateAlignments(ids[site]);
        saln->partitions.push_back(part_aln);
        int nseq = part_aln->getNSeq();
        //cout << "nseq  = " << nseq << endl;
        for (seq = 0; seq < nseq; seq++) {
            int id = saln->getSeqID(part_aln->getSeqName(seq));
            ASSERT(id >= 0);
            saln->taxa_index[id][site] = seq;
        }
    }
    // now the patterns of sequence-genes presence/absence
    saln->buildPattern();
    return saln;
}

void SuperAlignment::countConstSite() {
    num_informative_sites = 0;
    num_variant_sites = 0;
    max_num_states = 0;
    frac_const_sites = 0;
    frac_invariant_sites = 0;
    num_parsimony_sites = 0;
    size_t nsites = 0;
    for (vector<Alignment*>::iterator it = partitions.begin(); it != partitions.end(); it++) {
        (*it)->countConstSite();
        num_informative_sites += (*it)->num_informative_sites;
        num_variant_sites += (*it)->num_variant_sites;
        if ((*it)->num_states > max_num_states)
            max_num_states = (*it)->num_states;
        nsites += (*it)->getNSite();
        frac_const_sites += (*it)->frac_const_sites * (*it)->getNSite();
        frac_invariant_sites += (*it)->frac_invariant_sites * (*it)->getNSite();
    }
    frac_const_sites /= nsites;
    frac_invariant_sites /= nsites;
}

void SuperAlignment::orderPatternByNumChars(int pat_type) {
    const int UINT_BITS = sizeof(UINT)*8;
    if (pat_type == PAT_INFORMATIVE)
        num_parsimony_sites = num_informative_sites;
    else
        num_parsimony_sites = num_variant_sites;

    int maxi = (num_parsimony_sites+UINT_BITS-1)/UINT_BITS;
    pars_lower_bound = new UINT[maxi+1];
    memset(pars_lower_bound, 0, (maxi+1)*sizeof(UINT));
    int part, nseq = getNSeq();
    
    // compute ordered_pattern
    ordered_pattern.clear();
//    UINT sum_scores[npart];
    for (part  = 0; part != partitions.size(); part++) {
        partitions[part]->orderPatternByNumChars(pat_type);
        // partial_partition
        for (vector<Pattern>::iterator pit = partitions[part]->ordered_pattern.begin(); pit != partitions[part]->ordered_pattern.end(); pit++) {
            Pattern pattern(*pit);
            pattern.resize(nseq); // maximal unknown states
            for (int j = 0; j < nseq; j++)
                if (taxa_index[j][part] >= 0)
                    pattern[j] = (*pit)[taxa_index[j][part]];
                else
                    pattern[j] = partitions[part]->STATE_UNKNOWN;
            ordered_pattern.push_back(pattern);
        }
//        sum_scores[part] = partitions[part]->pars_lower_bound[0];
    }
    // TODO compute pars_lower_bound (lower bound of pars score for remaining patterns)
}
