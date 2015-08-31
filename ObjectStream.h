//
// Created by tung on 6/23/15.
//

#ifndef IQTREE_OBJECTSTREAM_H
#define IQTREE_OBJECTSTREAM_H
#include "TreeCollection.h"

/**
 *  This class is used to serialize object. It converts different object to byte stream
 *  and can also read in byte stream to reconstruct the object
 */
class ObjectStream {
public:

    /**
     * Constructor
     */
    ObjectStream(const char* data, size_t length);

    ObjectStream() {
        objectData = NULL;
    }

    virtual ~ObjectStream() {
        delete [] objectData;
    }

    /**
     *  Convert a tree collection into the internal byte stream
     *  @param[IN] trees
     */
    void writeObject(TreeCollection& trees);

    /**
     *  Reconstruct TreeCollection from a byte stream
     *  @param[OUT] trees
     */
    void readObject(TreeCollection& trees);


public:
    size_t getDataLength() const {
        return objectDataSize;
    }

public:
    char *getObjectData() const {
        return objectData;
    }

private:
    /**
     *  Byte stream representing the object
     */
    char* objectData;

    size_t objectDataSize;


    /**
     *  Convert vector of strings to array of chars
     *  @param [IN] strings the vector strings
     *  @param [OUT] the char array
     *  @return size of the char array
     */
    size_t serializeStringVector(vector<string> &strings, char*& data);

    /**
     *  Convert array of chars to vector of strings
     *  @param [IN] data byte stream representing vector<string>
     *  @param [IN] length size of data
     *  @param [OUT] strings the reconstructed vector<string>
     */
    void deserializeStringVector(char* data, size_t length, vector<string>& strings);

};
#endif // IQTREE_OBJECTSTREAM_H


