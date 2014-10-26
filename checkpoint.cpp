/*
 * checkpoint.cpp
 *
 *  Created on: Jun 12, 2014
 *      Author: minh
 */

#include "checkpoint.h"

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
}

void Checkpoint::setFileName(string filename) {
	this->filename = filename;
}
void Checkpoint::load() {
	assert(filename != "");
    try {
        ifstream in;
        // set the failbit and badbit
        in.exceptions(ios::failbit | ios::badbit);
        in.open(filename.c_str());
        string line;
        getline(in, line);
        if (line != "Checkpoint file for IQ-TREE")
        	throw ("Invalid checkpoint file");
        // remove the failbit
        in.exceptions(ios::badbit);
        while (!in.eof()) {
        	getline(in, line);
        	size_t pos = line.find(" := ");
        	if (pos == string::npos)
        		throw "':=' is expected between key and value";
        	(*this)[line.substr(0, pos)] = line.substr(pos+3);
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

void Checkpoint::commit() {
	assert(filename != "");
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(filename.c_str());
        out << "Checkpoint file for IQ-TREE" << endl;
        for (iterator i = begin(); i != end(); i++)
        	out << i->first << " := " << i->second << endl;
        out.close();
    } catch (ios::failure &) {
        outError(ERR_WRITE_OUTPUT, filename.c_str());
    }
}

bool Checkpoint::containsKey(string key) {
	return (find(key) != end());
}

/**
 * series of get functions
 */

template<class T>
void Checkpoint::get(string key, T& value) {
	assert(containsKey(key));
	stringstream ss((*this)[key]);
	ss >> value;
}

bool Checkpoint::getBool(string key) {
	assert(containsKey(key));
	if ((*this)[key] == "1") return true;
	return false;
}

char Checkpoint::getChar(string key) {
	assert(containsKey(key));
	return (*this)[key][0];
}

double Checkpoint::getDouble(string key) {
	assert(containsKey(key));
	return convert_double((*this)[key].c_str());

}

int Checkpoint::getInt(string key) {
	assert(containsKey(key));
	return convert_int((*this)[key].c_str());

}

/**
 * series of put functions
 */

template<class T>
void Checkpoint::put(string key, T value) {
	stringstream ss;
	ss << value;
	(*this)[key] = ss.str();
}

template<class T>
void Checkpoint::putArray(string key, int num, T* value) {
	stringstream ss;
	for (int i = 0; i < num; i++) {
		if (i > 0) ss << ',';
		ss << value[i];
	}
	(*this)[key] = ss.str();
}


Checkpoint::~Checkpoint() {
}

