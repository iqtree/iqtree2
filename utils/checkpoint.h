/*
 * checkpoint.h
 *
 *  Created on: Jun 12, 2014
 *      Author: minh
 */

#ifndef CHECKPOINT_H_
#define CHECKPOINT_H_

#include <stdio.h>
#include <map>
#include <string>
#include <sstream>
#include <cassert>
#include <vector>
#include <typeinfo>
#include "tools.h"

using namespace std;

// several useful declaration to save to or restore from a checkpoint
#define CKP_SAVE(var) checkpoint->put(#var, var)
#define CKP_ARRAY_SAVE(num, arr) checkpoint->putArray(#arr, num, arr)
#define CKP_VECTOR_SAVE(arr) checkpoint->putVector(#arr, arr)

#define CKP_SAVE2(checkpoint, var) checkpoint->put(#var, var)
#define CKP_ARRAY_SAVE2(checkpoint, num, arr) checkpoint->putArray(#arr, num, arr)
#define CKP_VECTOR_SAVE2(checkpoint, arr) checkpoint->putVector(#arr, arr)

#define CKP_RESTORE(var) checkpoint->get(#var, var)
#define CKP_RESTORE_STRING(var) checkpoint->getString(#var, var)
#define CKP_ARRAY_RESTORE(num, arr) checkpoint->getArray(#arr, num, arr)
#define CKP_VECTOR_RESTORE(arr) checkpoint->getVector(#arr, arr)

#define CKP_RESTORE2(checkpoint, var) checkpoint->get(#var, var)
#define CKP_RESTORE_STRING2(checkpoint, var) checkpoint->getString(#var, var)
#define CKP_ARRAY_RESTORE2(checkpoint, num, arr) checkpoint->getArray(#arr, num, arr)
#define CKP_VECTOR_RESTORE2(checkpoint, arr) checkpoint->getVector(#arr, arr)

const char CKP_SEP = '!';

/** checkpoint stream */
class CkpStream : public stringstream {
public:
    explicit CkpStream (ios_base::openmode which = ios_base::in | ios_base::out) : stringstream(which) {}

    explicit CkpStream (const string& str, ios_base::openmode which = ios_base::in | ios_base::out) : 
        stringstream(str, which) {}

};

/* overload operators */
//ostream& operator<<(ostream& os, const T& obj) {
//        return os;
//}
//
//std::istream& operator>>(std::istream& is, T& obj) {
//    return is;
//}

/**
 * Checkpoint as map from key strings to value strings
 */
class Checkpoint : public map<string, string> {
public:

    /** constructor */
	Checkpoint();

    /** destructor */
	virtual ~Checkpoint();

	/**
	 * @param filename file name
	 */
	void setFileName(string filename);

    string &getFileName() { return filename; }

    /** 
        set compression for checkpoint file
        @param compression true to compress checkpoint file, or false: no compression 
    */
    void setCompression(bool compression);

    /**
        set the header line to overwrite the default header
        @param header header line
    */
    void setHeader(string header);

	/**
	 * load checkpoint information from an input stram
     * @param in input stream
	 */
	void load(istream &in);

	/**
	 * load checkpoint information from file
     * @return TRUE if loaded successfully, otherwise FALSE
	 */
	bool load();

	/**
	 * dump checkpoint information into an output stream
     * @param out output stream
	 */
	void dump(ostream &out);

	/**
	 * dump checkpoint information into file
	 * @param force TRUE to dump no matter if time interval exceeded or not
	 */
	void dump(bool force = false);

    /**
        set dumping interval in seconds
        @param interval dumping interval
    */
    void setDumpInterval(double interval);

	/**
	 * @return true if checkpoint contains the key
	 * @param key key to search for
	 */
	bool hasKey(string key);

	/**
	 * @return true if checkpoint contains the key prefix
	 * @param key_prefix key prefix to search for
	 */
	bool hasKeyPrefix(string key_prefix);

    /**
        erase all entries with a key prefix
        @param key_prefix key prefix
        @return number of entries removed
    */
    int eraseKeyPrefix(string key_prefix);

    /**
     erase all entries without a key prefix
     @param key_prefix key prefix
     @return number of entries kept
     */
    int keepKeyPrefix(string key_prefix);

    /*-------------------------------------------------------------
     * series of get function to get value of a key
     *-------------------------------------------------------------*/
     
	/**
        @param key key name
        @param[out] value value for key
        @return true if key exists, false otherwise
	 */
	template<class T>
    bool get(string key, T& value) {
        if (key.empty())
            key = struct_name.substr(0, struct_name.length()-1);
        else
            key = struct_name + key;
        iterator it = find(key);
        if (it == end())
            return false;
        CkpStream ss(it->second);
        ss >> value;
        return true;
    }

	/**
        @param key key name
        @param[out] value entire string
        @return true if key exists, false otherwise
	 */
    bool getString(string key, string &value) {
        if (key.empty())
            key = struct_name.substr(0, struct_name.length()-1);
        else
            key = struct_name + key;
        iterator it = find(key);
        if (it == end())
            return false;
        value = it->second;
        return true;
    }

    /**
        get an array from checkpoint
        @param key key name
        @param num number of elements
        @param[out] value value
    */
	template<class T>
    bool getArray(string key, int maxnum, T* value) {
        if (key.empty())
            key = struct_name.substr(0, struct_name.length()-1);
        else
            key = struct_name + key;
        iterator it = find(key);
        if (it == end())
            return false;
        size_t pos = 0, next_pos;
        for (int i = 0; i < maxnum; i++) {
        	next_pos = it->second.find(", ", pos);
            CkpStream ss(it->second.substr(pos, next_pos-pos));
        	if (!(ss >> value[i]))
                break;
            if (next_pos == string::npos) {
                ASSERT(i == maxnum-1);
                break;
            }
        	pos = next_pos+2;
        }
        ASSERT(next_pos == string::npos);
        return true;
    }

    /**
        get an array from checkpoint
        @param key key name
        @param num number of elements
        @param[out] value value
    */
	template<class T>
    bool getVector(string key, vector<T> &value) {
        if (key.empty())
            key = struct_name.substr(0, struct_name.length()-1);
        else
            key = struct_name + key;
        iterator it = find(key);
        if (it == end())
            return false;
        size_t pos = 0, next_pos;
        value.clear();
        for (int i = 0; ; i++) {
        	next_pos = it->second.find(", ", pos);
            CkpStream ss(it->second.substr(pos, next_pos-pos));
            T val;
            if (ss >> val) {
                value.push_back(val);
            }
            else {
                break;
            }
            if (next_pos == string::npos) {
                break;
            }
        	pos = next_pos+2;
        }
        return true;
    }

    /**
        get a vector in YAML syntax from checkpoint
        @param key key name
        @param num number of elements
        @param[out] value value
    */
    template<class T>
    bool ymlGetVector(string key, vector<T> &value) {
        if (key.empty())
            key = struct_name.substr(0, struct_name.length()-1);
        else
            key = struct_name + key;
        iterator it = find(key);
        if (it == end())
            return false;
        size_t pos = 0, next_pos;
        value.clear();
        if ((pos = it->second.find('[')) == string::npos)
            outError(key + " vector not starting with [");
        for (int i = 0; ; i++) {
            next_pos = it->second.find(", ", pos);
            CkpStream ss(it->second.substr(pos, next_pos-pos));
            T val;
            if (ss >> val)
                value.push_back(val);
            else
                break;
            if (next_pos == string::npos) break;
            pos = next_pos+2;
        }
        return true;
    }

    /** 
        @param key key name
        @return bool value for key
    */
	bool getBool(string key, bool &ret);
	bool getBool(string key);

//    /** 
//        @param key key name
//        @return double value for key
//    */
//	double getDouble(string key);
//
//    /** 
//        @param key key name
//        @return int value for key
//    */
//	int getInt(string key);


    /*-------------------------------------------------------------
     * series of put function to put pair of (key,value)
     *-------------------------------------------------------------*/

    /**
        put pair of (key,value) to checkpoint
        @param key key name
        @param value value
    */
	template<class T>
	void put(string key, T value) {
        if (key.empty())
            key = struct_name.substr(0, struct_name.length()-1);
        else
            key = struct_name + key;
        CkpStream ss;
        ss.precision(10);
        ss << value;
        (*this)[key] = ss.str();
    }
    
    /** 
        @param key key name
        @param value
    */
	void putBool(string key, bool value);

    /**
        put an array to checkpoint
        @param key key name
        @param num number of elements
        @param value value
    */
	template<class T>
	void putArray(string key, int num, T* value) {
        if (key.empty())
            key = struct_name.substr(0, struct_name.length()-1);
        else
            key = struct_name + key;
        CkpStream ss;
        ss.precision(10);
        for (int i = 0; i < num; i++) {
            if (i > 0) ss << ", ";
            ss << value[i];
        }
        (*this)[key] = ss.str();
    }

    /**
        put an STL vector to checkpoint
        @param key key name
        @param num number of elements
        @param value value
    */
	template<class T>
	void putVector(string key, vector<T> &value) {
        if (key.empty())
            key = struct_name.substr(0, struct_name.length()-1);
        else
            key = struct_name + key;
        CkpStream ss;
        ss.precision(10);
        for (int i = 0; i < value.size(); i++) {
            if (i > 0) ss << ", ";
            ss << value[i];
        }
        (*this)[key] = ss.str();
    }
    
    /*-------------------------------------------------------------
     * helper functions
     *-------------------------------------------------------------*/

    /**
        start a new struct
    */
    void startStruct(string name);

    /**
        end the current struct
    */
    void endStruct();

    /**
        @return struct_name
     */
    string getStructName() { return struct_name; }
    
    /**
        start a new list in the current scope
        @param nelem number of elements
    */
    void startList(int nelem);

    /**
        set the starting list element, should only be called right after startList
        @param id element ID
    */
    void setListElement(int id);

    /** 
        add an element to the current list
    */
    void addListElement();
    
    /**
        end the current list
    */
    void endList();

    /**
        get a subset of checkpoint where the key string contains a given substring
        @param[out] target checkpoint
        @param sub_key key substring to search for
    */
    void getSubCheckpoint(Checkpoint *target, string sub_key);

    /**
     put a checkpoint where the key string contains a given substring
     @param source checkpoint
     @param sub_key key substring to search for
     */
    void putSubCheckpoint(Checkpoint *source, string sub_key);

    /**
     transfer a checkpoint where the key string contains a given substring
     @param source checkpoint
     @param sub_key key substring to search for
     @param overwrite true to overwrite value even if key exists
     */
    void transferSubCheckpoint(Checkpoint *source, string sub_key, bool overwrite = false);

protected:

    /** filename to write checkpoint */
	string filename;
    
    /** previous dump time in seconds */
    double prev_dump_time;
    
    /** dumping time interval */
    double dump_interval;
    
    /** count number of dumping times */
    int dump_count;
    
    /** true (default) to compress checkpoint file, false: no compression */
    bool compression;
    
    /** header line of checkpoint file */
    string header;
    
private:

    /** name of the current nested key */
    string struct_name;

    /** current list element ID */
    vector<int> list_element;
    
    /** width to element ID for prefixing with '0' */
    vector<int> list_element_precision;
    
};



/**
    Root class handling all checkpoint facilities. Inherit this class
    if you want a class to be checkpointed.
*/
class CheckpointFactory {
public:

    /** constructor */
    CheckpointFactory();

    virtual ~CheckpointFactory() {}

    /**
        set checkpoint object
        @param checkpoint 
    */
    virtual void setCheckpoint(Checkpoint *checkpoint);

    /**
        get checkpoint object
        @return checkpoint 
    */
    Checkpoint *getCheckpoint();

    /**
        start structure for checkpointing
    */
    virtual void startCheckpoint();

    /**
        save object into the checkpoint
    */
    virtual void saveCheckpoint();

    /** 
        restore object from the checkpoint
    */
    virtual void restoreCheckpoint();

    /**
        end structure for checkpointing
    */
    virtual void endCheckpoint();

protected:

    Checkpoint *checkpoint;

};

#endif /* CHECKPOINT_H_ */
