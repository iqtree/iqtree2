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

using namespace std;

//#include "tools.h"

/**
 * Checkpoint as map from key strings to value strings
 */
class Checkpoint : public map<string, string> {
public:

    /** constructor */
	Checkpoint();

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
	 */
	void dump();

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
	 */
	template<class T>
    void get(string key, T& value) {
        assert(containsKey(key));
        stringstream ss((*this)[key]);
        ss >> value;
    }


    /** 
        @param key key name
        @return bool value for key
    */
	bool getBool(string key);

    /** 
        @param key key name
        @return char value for key
    */
	char getChar(string key);

    /** 
        @param key key name
        @return double value for key
    */
	double getDouble(string key);

    /** 
        @param key key name
        @return int value for key
    */
	int getInt(string key);


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
        stringstream ss;
        ss << value;
        (*this)[key] = ss.str();
    }
    

    /**
        put an array to checkpoint
        @param key key name
        @param num number of elements
        @param value value
    */
	template<class T>
	void putArray(string key, int num, T* value) {
        stringstream ss;
        for (int i = 0; i < num; i++) {
            if (i > 0) ss << ',';
            ss << value[i];
        }
        (*this)[key] = ss.str();
    }
    

    /** destructor */
	virtual ~Checkpoint();

protected:

    /** filename to write checkpoint */
	string filename;
    
    /** previous dump time in seconds */
    double prev_dump_time;
    
    /** dumping time interval */
    double dump_interval;
};


#endif /* CHECKPOINT_H_ */
