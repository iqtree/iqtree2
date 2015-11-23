/*
 * checkpoint.h
 *
 *  Created on: Jun 12, 2014
 *      Author: minh
 */

#ifndef CHECKPOINT_H_
#define CHECKPOINT_H_

#include "tools.h"

/**
 * Checkpoint as map from key strings to value strings
 */
class Checkpoint : public map<string, string> {
public:
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

	/**
	 * series of get functions
	 */
	template<class T>
	void get(string key, T& value);

	bool getBool(string key);
	char getChar(string key);
	double getDouble(string key);
	int getInt(string key);


	/**
	 * series of put functions
	 */
	template<class T>
	void put(string key, T value);

	template<class T>
	void putArray(string key, int num, T* value);

	virtual ~Checkpoint();

	string filename;
};

#endif /* CHECKPOINT_H_ */
