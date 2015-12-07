/*
 * checkpoint.cpp
 *
 *  Created on: Jun 12, 2014
 *      Author: minh
 */

#include "checkpoint.h"
#include "tools.h"
#include "timeutil.h"

/*
 * The following parameters have been saved for checkpoint in IQPNNI
 *
Number iterations: 200
Maximum number iterations: 2000
Current number iterations: 139
Probability of deleting a sequence: 0.5
Number representatives: 4
Stopping rule (0: YES, 1: YES_MIN_ITER, 2: YES_MAX_ITER, 3: NO): 3
Type of data (0:NUCLEOTIDE, 1:AMINO_ACID): 0
Substitution model (0:HKY85, 1: TN93, 2:GTR, 3:WAG, 4:JTT, 5:VT, 6:MtREV24, 7:Blosum62, 8:Dayhoff, 9:rtREV, 10: User-defined): 0
Frequency of Base A: 0.248672
Frequency of Base C: 0.261687
Frequency of Base G: 0.250996
Frequency of Base T: 0.238645
Type of parameters (0:ESTIMATE,  1:USER_DEFINED, 2: EQUAL): 0
Transition/transversion ratito: 0.766912
Type of parameters (0:ESTIMATE,  1:USER_DEFINED): 0
Pyridimine/purine ratito: 1
Type of parameters (0:ESTIMATE,  1:USER_DEFINED): 0
Transition rate from A to G: -1
Transition rate from C to T: -1
Transversion rate from A to C: -1
Transversion rate from A to T: -1
Transversion rate from C to G: -1
Transversion rate from G to T: -1
Type of parameters (0:ESTIMATE,  1:USER_DEFINED): 0
Type of rate heterogeneity (0:UNIFORM, 1:SITE_SPECIFIC, 2:GAMMA): 0
Number rates: 1
Gamma distribution parameter alpha: 1
Type of parameters (0:ESTIMATE,  1:USER_DEFINED): 0
Invariant type (0: NONE, 1:ESTIMATE, 2: USER_DEFINED): 0
Proportion of invariable sites: 0
Out group sequence: 0
Bootstrap sample: 0
Current bootstrap sample: 0
Build consensus: 0
Current best log-likelihood: -11833.35062
Elapsed time: 23
Finished: 0
 */

Checkpoint::Checkpoint() {
	filename = "";
    prev_dump_time = 0;
    dump_interval = 30; // dumping at most once per 30 seconds
    struct_name = "";
    list_element = -1;
}


Checkpoint::~Checkpoint() {
}


const char* CKP_HEADER = "--- # IQ-TREE Checkpoint";

void Checkpoint::setFileName(string filename) {
	this->filename = filename;
}
void Checkpoint::load() {
	assert(filename != "");
    if (!fileExists(filename)) return;
    try {
        ifstream in;
        // set the failbit and badbit
        in.exceptions(ios::failbit | ios::badbit);
        in.open(filename.c_str());
        cout << "Loading checkpoint file " << filename << "..." << endl;
        string line;
        getline(in, line);
        if (line != CKP_HEADER)
        	throw ("Invalid checkpoint file");
        // remove the failbit
        in.exceptions(ios::badbit);
        string struct_name;
        size_t pos;
        while (!in.eof()) {
        	getline(in, line);
            pos = line.find('#');
            if (pos != string::npos)
                line.erase(pos);
            line.erase(line.find_last_not_of("\n\r\t")+1);
//            trimString(line);
            if (line.empty()) continue;
            if (line[0] != ' ') {
                struct_name = "";
            }
//            trimString(line);
            line.erase(0, line.find_first_not_of(" \n\r\t"));
            if (line.empty()) continue;
        	pos = line.find(": ");
        	if (pos != string::npos) {
                (*this)[struct_name + line.substr(0, pos)] = line.substr(pos+2);
            } else if (line[line.length()-1] == ':') {
                line.erase(line.length()-1);
                trimString(line);
                struct_name = line + '.';
                continue;
            } else {
        		throw "':' is expected between key and value";
            }
        }
        in.clear();
        // set the failbit again
        in.exceptions(ios::failbit | ios::badbit);
        in.close();
    } catch (ios::failure &) {
        outError(ERR_READ_INPUT);
    } catch (const char *str) {
        outError(str);
    } catch (string &str) {
        outError(str);
    }
}

void Checkpoint::setDumpInterval(double interval) {
    dump_interval = interval;
}


void Checkpoint::dump(bool force) {
	assert(filename != "");
    if (!force && getRealTime() < prev_dump_time + dump_interval) {
        return;
    }
    prev_dump_time = getRealTime();
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(filename.c_str());
        out << CKP_HEADER << endl;
        string struct_name;
        size_t pos;
        for (iterator i = begin(); i != end(); i++) {
            if ((pos = i->first.find('.')) != string::npos) {
                if (struct_name != i->first.substr(0, pos)) {
                    struct_name = i->first.substr(0, pos);
                    out << struct_name << ":" << endl;
                }
                out << "  " << i->first.substr(pos+1) << ": " << i->second << endl;
            } else
                out << i->first << ": " << i->second << endl;
        }
        out.close();
        cout << "Checkpoint dumped" << endl;
    } catch (ios::failure &) {
        outError(ERR_WRITE_OUTPUT, filename.c_str());
    }
}

bool Checkpoint::containsKey(string key) {
	return (find(key) != end());
}

/*-------------------------------------------------------------
 * series of get function to get value of a key
 *-------------------------------------------------------------*/

bool Checkpoint::getBool(string key) {
    string value;
    if (!get(key, value)) return false;
	if (value == "true") return true;
    if (value == "false") return false;
    outError("Invalid boolean value " + value + " for key " + key);
    return false;
}

double Checkpoint::getDouble(string key) {
    string value;
    if (!get(key, value)) return -DBL_MAX;
	return convert_double(value.c_str());
}

int Checkpoint::getInt(string key) {
    string value;
    if (!get(key, value)) return -INT_MAX;
	return convert_int(value.c_str());
}


/*-------------------------------------------------------------
 * series of put function to put pair of (key,value)
 *-------------------------------------------------------------*/

void Checkpoint::putBool(string key, bool value) {
    if (value)
        put(key, "true");
    else
        put(key, "false");
}


/*-------------------------------------------------------------
 * nested structures
 *-------------------------------------------------------------*/
void Checkpoint::startStruct(string name) {
    struct_name = struct_name + name + '.';
    list_element = -1;
}

/**
    end the current struct
*/
void Checkpoint::endStruct() {
    size_t pos = struct_name.find_last_of('.', struct_name.length()-2);
    if (pos == string::npos)
        struct_name = "";
    else
        struct_name.erase(pos+1);
    list_element = -1;
}

void Checkpoint::startListElement() {
    list_element++;    
    struct_name += convertIntToString(list_element) + ".";
}

void Checkpoint::endListElement() {
    size_t pos = struct_name.find_last_of('.', struct_name.length()-2);
    assert(pos != string::npos);
    struct_name.erase(pos+1);
}

/*-------------------------------------------------------------
 * CheckpointFactory
 *-------------------------------------------------------------*/

CheckpointFactory::CheckpointFactory() {
    checkpoint = NULL;
}

void CheckpointFactory::setCheckpoint(Checkpoint *checkpoint) {
    this->checkpoint = checkpoint;
}

Checkpoint *CheckpointFactory::getCheckpoint() {
    return checkpoint;
}

void CheckpointFactory::saveCheckpoint() {
    if (!checkpoint) return;
    checkpoint->dump();
}

void CheckpointFactory::restoreCheckpoint() {
    // do nothing
}

