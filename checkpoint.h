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

using namespace std;

// several useful declaration to save to or restore from a checkpoint
#define CKP_SAVE(var) checkpoint->put(#var, var)
#define CKP_ARRAY_SAVE(num, arr) checkpoint->putArray(#arr, num, arr)
#define CKP_VECTOR_SAVE(arr) checkpoint->putVector(#arr, arr)

#define CKP_RESTORE(var) checkpoint->get(#var, var)
#define CKP_ARRAY_RESTORE(num, arr) checkpoint->getArray(#arr, num, arr)
#define CKP_VECTOR_RESTORE(arr) checkpoint->getVector(#arr, arr)

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

	/**
	 * load checkpoint information from file
	 */
	void load();

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
	bool containsKey(string key);

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
        key = struct_name + key;
        iterator it = find(key);
        if (it == end())
            return false;
        CkpStream ss(it->second);
        ss >> value;
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
        key = struct_name + key;
        iterator it = find(key);
        if (it == end())
            return false;
        size_t pos = 0, next_pos;
        for (int i = 0; i < maxnum; i++) {
        	next_pos = it->second.find(", ", pos);
            CkpStream ss(it->second.substr(pos, next_pos-pos));
        	ss >> value[i];
        	if (next_pos == string::npos) break;
        	pos = next_pos+2;
        }
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
        	ss >> val;
            value.push_back(val);
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
        key = struct_name + key;
        CkpStream ss;
        if (typeid(T) == typeid(double))
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
        key = struct_name + key;
        CkpStream ss;
        if (typeid(T) == typeid(double))
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
        key = struct_name + key;
        CkpStream ss;
        if (typeid(T) == typeid(double))
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
        start a new list element in the current scope
    */
    void startListElement();
    
    /**
        end the current list element
    */
    void endListElement();

protected:

    /** filename to write checkpoint */
	string filename;
    
    /** previous dump time in seconds */
    double prev_dump_time;
    
    /** dumping time interval */
    double dump_interval;
    
private:

    /** name of the current nested key */
    string struct_name;

    /** current list element ID */
    int list_element;

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
        save object into the checkpoint
    */
    virtual void saveCheckpoint();

    /** 
        restore object from the checkpoint
    */
    virtual void restoreCheckpoint();

protected:

    Checkpoint *checkpoint;

};

#endif /* CHECKPOINT_H_ */
