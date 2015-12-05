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
	 * commit checkpoint information into file
	 */
	void commit();

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
	void get(string key, T& value);

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
	void put(string key, T value);

    /**
        put an array to checkpoint
        @param key key name
        @param num number of elements
        @param value value
    */
	template<class T>
	void putArray(string key, int num, T* value);

    /** destructor */
	virtual ~Checkpoint();

    /** filename to write checkpoint */
	string filename;
};

#endif /* CHECKPOINT_H_ */
