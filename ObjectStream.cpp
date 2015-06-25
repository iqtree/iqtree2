//
// Created by tung on 6/23/15.
//

#include "ObjectStream.h"

ObjectStream::ObjectStream(const char *data, size_t length) {
    objectData = new char[length];
    memcpy(objectData, data, length);
    objectDataSize = length;
}

void ObjectStream::writeObject(TreeCollection &trees) {
    vector<string> treeStrings = trees.getTreeStrings();
    vector<double> scores = trees.getScores();
    size_t numTrees = trees.getNumTrees();
    if (objectData != NULL)
        delete [] objectData;
    char* stringData;
    size_t stringDataSize = serializeStringVector(treeStrings, stringData);
    size_t doubleDataSize = scores.size() * sizeof(double);

    objectDataSize = sizeof(size_t) * 2 + stringDataSize + doubleDataSize;
    char* pos = objectData;

    // Copy the size of the string block and double block into the beginning of objectData
    memcpy(pos, &stringDataSize, sizeof(size_t));
    pos = pos + sizeof(size_t);
    memcpy(pos, &doubleDataSize, sizeof(size_t));
    pos = pos + sizeof(size_t);

    // Add string block and double block afterwards
    memcpy(pos, stringData, stringDataSize);
    pos = pos + stringDataSize;
    memcpy(pos, scores.data(), doubleDataSize);

    delete [] stringData;
}

void ObjectStream::readObject(TreeCollection &trees) {
    trees.clear();
    size_t metaInfo[2];
    memcpy(metaInfo, objectData, sizeof(size_t) * 2);
    size_t stringDataSize = metaInfo[0];
    size_t doubleDataSize = metaInfo[1];
    size_t numTrees = doubleDataSize / sizeof(double);
    vector<string> treeStrings;
    vector<double> scores;
    deserializeStringVector(objectData + sizeof(size_t) * 2, stringDataSize, treeStrings);
    assert(treeStrings.size() == numTrees);
    double scoreArr[numTrees];
    memcpy(scoreArr, objectData + sizeof(size_t) * 2 + stringDataSize, doubleDataSize);
    for (int i = 0; i < numTrees; ++i)
        scores.push_back(scoreArr[i]);
    trees.setTreeStrings(treeStrings);
    trees.setScores(scores);
}


size_t ObjectStream::serializeStringVector(vector<string> &strings, char *&data) {
    size_t numStrings = strings.size();
    size_t totalSize = 0;
    // Determine the total bytes required
    for (int i = 0; i < numStrings; i++) {
        totalSize += strings[i].length() + 1;
    }
    data = new char[totalSize];
    char* pos = data;
    for (int i = 0; i < numStrings; i++) {
        size_t length = strings[i].length();
        const char* cString = strings[i].c_str();
        strncpy(pos, cString, length + 1);
        pos = pos + length + 1;
    }
    return totalSize;
}

void ObjectStream::deserializeStringVector(char *data, size_t length, vector<string> &strings) {
    strings.clear();
    stringstream ss;
    ss.str("");
    for (int i = 0; i < length; i++) {
        if (data[i] == '\0') {
            strings.push_back(ss.str());
            ss.str("");
        } else {
            ss << data[i];
        }
    }
}
