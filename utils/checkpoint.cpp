/*
 * checkpoint.cpp
 *
 *  Created on: Jun 12, 2014
 *      Author: minh
 */

#include "checkpoint.h"
#include "tools.h"
#include "timeutil.h"
#include "gzstream.h"
#include <cstdio>

const char* CKP_HEADER =     "--- # IQ-TREE Checkpoint ver >= 1.6";
const char* CKP_HEADER_OLD = "--- # IQ-TREE Checkpoint";

Checkpoint::Checkpoint() {
	filename = "";
    prev_dump_time = 0;
    dump_interval = 60; // dumping at most once per 60 seconds
    dump_count = 0;
    struct_name = "";
    compression = true;
    header = CKP_HEADER;
}


Checkpoint::~Checkpoint() {
}


void Checkpoint::setFileName(string filename) {
	this->filename = filename;
}


void Checkpoint::load(istream &in) {
    string line;
    string struct_name;
    size_t pos;
    int listid = 0;
    while (!in.eof()) {
        safeGetline(in, line);
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
            // mapping
            (*this)[struct_name + line.substr(0, pos)] = line.substr(pos+2);
        } else if (line[line.length()-1] == ':') {
            // start a new struct
            line.erase(line.length()-1);
            trimString(line);
            struct_name = line + CKP_SEP;
            listid = 0;
            continue;
        } else {
            // collection
            (*this)[struct_name + convertIntToString(listid)] = line;
            listid++;
        }
    }
}


bool Checkpoint::load() {
	ASSERT(filename != "");
    if (!fileExists(filename)) return false;
    try {
        igzstream in;
        // set the failbit and badbit
        in.exceptions(ios::failbit | ios::badbit);
        in.open(filename.c_str());
        // remove the failbit
        in.exceptions(ios::badbit);
        string line;
        if (!safeGetline(in, line)) {
            in.close();
            return false;
        }
        if (line == CKP_HEADER_OLD)
            throw "Incompatible checkpoint file from version 1.5.X or older.\nEither overwrite it with -redo option or run older version";
        if (line != header)
        	throw ("Invalid checkpoint file " + filename);
        // call load from the stream
        load(in);
        in.clear();
        // set the failbit again
        in.exceptions(ios::failbit | ios::badbit);
        in.close();
        return true;
    } catch (ios::failure &) {
        outError(ERR_READ_INPUT);
    } catch (const char *str) {
        outError(str);
    } catch (string &str) {
        outError(str);
    }
    return false;
}

void Checkpoint::setCompression(bool compression) {
    this->compression = compression;
}

/**
    set the header line to overwrite the default header
    @param header header line
*/
void Checkpoint::setHeader(string header) {
    this->header = "--- # " + header;
}

void Checkpoint::setDumpInterval(double interval) {
    dump_interval = interval;
}

void Checkpoint::dump(ostream &out) {
    string struct_name;
    size_t pos;
    int listid = 0;
    for (iterator i = begin(); i != end(); i++) {
        if ((pos = i->first.find(CKP_SEP)) != string::npos) {
            if (struct_name != i->first.substr(0, pos)) {
                struct_name = i->first.substr(0, pos);
                out << struct_name << ':' << endl;
                listid = 0;
            }
            // check if key is a collection
            out << ' ' << i->first.substr(pos+1) << ": " << i->second << endl;
        } else
            out << i->first << ": " << i->second << endl;
    }
}

void Checkpoint::dump(bool force) {
    if (filename == "")
        return;
        
    if (!force && getRealTime() < prev_dump_time + dump_interval) {
        return;
    }
    prev_dump_time = getRealTime();
    string filename_tmp = filename + ".tmp";
    if (fileExists(filename_tmp)) {
        outWarning("IQ-TREE was killed while writing temporary checkpoint file " + filename_tmp);
        outWarning("You should increase checkpoint interval from the default 60 seconds");
        outWarning("via -cptime option to avoid too frequent checkpoint for large datasets");
    }
    try {
        ostream *out;
        if (compression) 
            out = new ogzstream(filename_tmp.c_str());
        else
            out = new ofstream(filename_tmp.c_str());
        out->exceptions(ios::failbit | ios::badbit);
        *out << header << endl;
        // call dump stream
        dump(*out);
        if (compression)
            ((ogzstream*)out)->close();
        else
            ((ofstream*)out)->close();
        delete out;
//        cout << "Checkpoint dumped" << endl;
        if (fileExists(filename)) {
            if (std::remove(filename.c_str()) != 0)
                outError("Cannot remove file ", filename);
        }
        if (std::rename(filename_tmp.c_str(), filename.c_str()) != 0)
            outError("Cannot rename file ", filename_tmp);
    } catch (ios::failure &) {
        outError(ERR_WRITE_OUTPUT, filename.c_str());
    }
    if (Params::getInstance().print_all_checkpoints) {
        // Feature request by Nick Goldman
        dump_count++;
        filename_tmp = (string)Params::getInstance().out_prefix + "." + convertIntToString(dump_count) + ".ckp.gz";
        try {
            ostream *out;
            if (compression)
                out = new ogzstream(filename_tmp.c_str());
            else
                out = new ofstream(filename_tmp.c_str());
            out->exceptions(ios::failbit | ios::badbit);
            *out << header << endl;
            // call dump stream
            dump(*out);
            if (compression)
                ((ogzstream*)out)->close();
            else
                ((ofstream*)out)->close();
            delete out;
        } catch (ios::failure &) {
            outError(ERR_WRITE_OUTPUT, filename_tmp.c_str());
        }
    } else {
        // check that the dumping time is too long and increase dump_interval if necessary
        double dump_time = getRealTime() - prev_dump_time;
        if (dump_time*20 > dump_interval) {
            dump_interval = ceil(dump_time*20);
            cout << "NOTE: " << dump_time << " seconds to dump checkpoint file, increase to "
            << dump_interval << endl;
        }
    }
}

bool Checkpoint::hasKey(string key) {
	return (find(struct_name + key) != end());
}

bool Checkpoint::hasKeyPrefix(string key_prefix) {
    string prefix = key_prefix;
    if (!struct_name.empty()) {
        prefix = struct_name + key_prefix;
    }
	auto i = lower_bound(prefix);
    if (i != end()) {
        if (i->first.compare(0, prefix.size(), prefix) == 0)
            return true;
    }
    return false;
}

int Checkpoint::eraseKeyPrefix(string key_prefix) {
    int count = 0;
    iterator first_it = lower_bound(key_prefix);
    iterator i;
	for (i = first_it; i != end(); i++) {
        if (i->first.compare(0, key_prefix.size(), key_prefix) == 0)
            count++;
        else
            break;

    }
    if (count)
        erase(first_it, i);
    return count;
}

int Checkpoint::keepKeyPrefix(string key_prefix) {
    map<string,string> newckp;
    int count = 0;
    erase(begin(), lower_bound(key_prefix));
    
    for (iterator i = begin(); i != end(); i++) {
        if (i->first.compare(0, key_prefix.size(), key_prefix) == 0)
            count++;
        else {
            erase(i, end());
            break;
        }
        
    }
    return count;
}

/*-------------------------------------------------------------
 * series of get function to get value of a key
 *-------------------------------------------------------------*/

bool Checkpoint::getBool(string key, bool &ret) {
    string value;
    if (!get(key, value)) return false;
	if (value == "true") 
        ret = true;
    else if (value == "false") 
        ret = false;
    else
        outError("Invalid boolean value " + value + " for key " + key);
    return true;
}

bool Checkpoint::getBool(string key) {
    bool ret;
    if (!getBool(key, ret))
        return false;
    return ret;
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
    struct_name = struct_name + name + CKP_SEP;
}

/**
    end the current struct
*/
void Checkpoint::endStruct() {
    size_t pos = struct_name.find_last_of(CKP_SEP, struct_name.length()-2);
    if (pos == string::npos)
        struct_name = "";
    else
        struct_name.erase(pos+1);
}

void Checkpoint::startList(int nelem) {
    list_element.push_back(-1);
    if (nelem > 0)
        list_element_precision.push_back((int)ceil(log10(nelem)));
    else
        list_element_precision.push_back(0);
}

void Checkpoint::setListElement(int id) {
    list_element.back() = id;
    stringstream ss;
    ss << setw(list_element_precision.back()) << setfill('0') << list_element.back();
    struct_name += ss.str() + CKP_SEP;
}

void Checkpoint::addListElement() {
    list_element.back()++;
    if (list_element.back() > 0) {
        size_t pos = struct_name.find_last_of(CKP_SEP, struct_name.length()-2);
        ASSERT(pos != string::npos);
        struct_name.erase(pos+1);
    }
    stringstream ss;
    ss << setw(list_element_precision.back()) << setfill('0') << list_element.back();
//    ss << list_element.back();
    struct_name += ss.str() + CKP_SEP;
}

void Checkpoint::endList() {
    ASSERT(!list_element.empty());

    if (list_element.back() >= 0) {
        size_t pos = struct_name.find_last_of(CKP_SEP, struct_name.length()-2);
        ASSERT(pos != string::npos);
        struct_name.erase(pos+1);
    }

    list_element.pop_back();
    list_element_precision.pop_back();

}

void Checkpoint::getSubCheckpoint(Checkpoint *target, string partial_key) {
    int len = partial_key.length();
    for (auto it = lower_bound(partial_key); it != end() && it->first.substr(0, len) == partial_key; it++) {
        target->put(it->first.substr(len+1), it->second);
    }
}

void Checkpoint::putSubCheckpoint(Checkpoint *source, string partial_key) {
    if (!partial_key.empty())
        startStruct(partial_key);
    for (auto it = source->begin(); it != source->end(); it++) {
        put(it->first, it->second);
    }
    if (!partial_key.empty())
        endStruct();
}

void Checkpoint::transferSubCheckpoint(Checkpoint *target, string partial_key, bool overwrite) {
    int len = partial_key.length();
    for (auto it = lower_bound(partial_key); it != end() && it->first.substr(0, len) == partial_key; it++) {
        if (overwrite || !target->hasKey(it->first))
            target->put(it->first, it->second);
    }
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

void CheckpointFactory::startCheckpoint() {
    checkpoint->startStruct("CheckpointFactory");
}

void CheckpointFactory::saveCheckpoint() {
    // do nothing
}

void CheckpointFactory::restoreCheckpoint() {
    // do nothing
}

void CheckpointFactory::endCheckpoint() {
    checkpoint->endStruct();
}

