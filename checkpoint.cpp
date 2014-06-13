/*
 * checkpoint.cpp
 *
 *  Created on: Jun 12, 2014
 *      Author: minh
 */

#include "checkpoint.h"

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

