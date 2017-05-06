//
// Created by tung on 6/23/15.
//

#include "ObjectStream.h"

ObjectStream::ObjectStream(const char *data, size_t length) {
    objectData = new char[length];
    memcpy(objectData, data, length);
    objectDataSize = length;
}

ObjectStream::ObjectStream(TreeCollection &trees) {
    objectData = NULL;
    objectDataSize = 0;
    initFromTreeCollection(trees);
}

void ObjectStream::initFromTreeCollection(TreeCollection &trees) {
    vector<string> treeStrings = trees.getTreeStrings();
    vector<double> scores = trees.getScores();
    vector<int> sourceProcID = trees.getSourceProcID();

    char* stringData;
    size_t stringDataSize = serializeStrings(treeStrings, stringData);
    size_t doubleDataSize = scores.size() * sizeof(double);
    size_t intDataSize = sourceProcID.size() * sizeof(int);

    objectDataSize = sizeof(size_t) * 3 + stringDataSize + doubleDataSize + intDataSize;

    if (objectData != NULL) {
        delete[] objectData;
    }
    objectData = new char[objectDataSize];

    char* pos = objectData;
    // Copy the size of the string block and double block into the beginning of objectData
    memcpy(pos, &stringDataSize, sizeof(size_t));
    pos = pos + sizeof(size_t);
    memcpy(pos, &doubleDataSize, sizeof(size_t));
    pos = pos + sizeof(size_t);
    memcpy(pos, &intDataSize, sizeof(size_t));
    pos = pos + sizeof(size_t);

    // Add string block and double block afterwards
    memcpy(pos, stringData, stringDataSize);
    pos = pos + stringDataSize;
    
    memcpy(pos, scores.data(), doubleDataSize);
    pos = pos + doubleDataSize;
    
    memcpy(pos, sourceProcID.data(), intDataSize);

    delete [] stringData;
}

TreeCollection ObjectStream::getTreeCollection() {
    size_t metaInfo[3];
    memcpy(metaInfo, objectData, sizeof(size_t) * 3);
    size_t stringDataSize = metaInfo[0];
    size_t doubleDataSize = metaInfo[1];
    size_t intDataSize = metaInfo[2];
    size_t numTrees = doubleDataSize / sizeof(double);
    vector<string> treeStrings;
    deserializeStrings(objectData + sizeof(size_t) * 3, stringDataSize, treeStrings);
    ASSERT(treeStrings.size() == numTrees);

    double scoreArr[numTrees];
    memcpy(scoreArr, objectData + sizeof(size_t) * 3 + stringDataSize, doubleDataSize);
    vector<double> scores(scoreArr, scoreArr + sizeof(scoreArr) / sizeof(scoreArr[0]));

    int sourceProcIDArr[numTrees];
    memcpy(sourceProcIDArr, objectData + sizeof(size_t) * 3 + stringDataSize + doubleDataSize, intDataSize);
    vector<int> sourceProcID(sourceProcIDArr, sourceProcIDArr + sizeof(sourceProcIDArr) / sizeof(sourceProcIDArr[0]));

    TreeCollection decodedTrees(treeStrings, scores, sourceProcID);
    return decodedTrees;
}


size_t ObjectStream::serializeStrings(vector<string> &strings, char *&data) {
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

void ObjectStream::deserializeStrings(char *data, size_t length, vector<string> &strings) {
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


