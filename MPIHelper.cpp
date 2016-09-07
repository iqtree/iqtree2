//
// Created by tung on 6/18/15.
//

#include "MPIHelper.h"
#include "timeutil.h"

/**
 *  Initialize the single getInstance of MPIHelper
 */

MPIHelper& MPIHelper::getInstance() {
    static MPIHelper instance;
#ifndef _IQTREE_MPI
    instance.setProcessID(0);
    instance.setNumProcesses(1);
#endif
    return instance;
}

void MPIHelper::distributeTrees(vector<string> treeStrings, vector<double> scores, int tag) {
    if (getNumProcesses() == 1)
        return;
#ifdef _IQTREE_MPI
    vector<int> sourceProcID;
    sourceProcID.insert(sourceProcID.end(), scores.size(), getProcessID());
    TreeCollection outTrees(treeStrings, scores, sourceProcID);
    cleanUpMessages();
    for (int i = 0; i < getNumProcesses(); i++) {
        if (i != getProcessID()) {
            MPI_Request *request = new MPI_Request;
            ObjectStream *os = new ObjectStream(outTrees);
            MPI_Isend(os->getObjectData(), os->getDataLength(), MPI_CHAR, i, tag, MPI_COMM_WORLD, request);
            sentMessages.push_back(make_pair(request, os));
        }
    }
    //numTreeSent += treeStrings.size();
#endif
}

void MPIHelper::distributeTree(string treeString, double score, int tag) {
    if (getNumProcesses() == 1)
        return;
#ifdef _IQTREE_MPI
    double start = getRealTime();
    vector<string> trees;
    vector<double> scores;
    trees.push_back(treeString);
    scores.push_back(score);
    distributeTrees(trees, scores, tag);
    if (verbose_mode >= VB_MED)
        cout << "Sent tree to other processes in " << getRealTime() - start << " seconds" << endl;
    numTreeSent++;
#endif
}

void MPIHelper::sendTrees(int dest, vector<string> treeStrings, vector<double> scores, int tag) {
    if (getNumProcesses() == 1 || dest == getProcessID())
        return;
#ifdef _IQTREE_MPI
    vector<int> sourceProcID;
    sourceProcID.insert(sourceProcID.end(), scores.size(), getProcessID());
    TreeCollection outTrees(treeStrings, scores, sourceProcID);
    cleanUpMessages();
    MPI_Request *request = new MPI_Request;
    ObjectStream *os = new ObjectStream(outTrees);
    MPI_Isend(os->getObjectData(), os->getDataLength(), MPI_CHAR, dest, tag, MPI_COMM_WORLD, request);
    sentMessages.push_back(make_pair(request, os));
    numTreeSent += treeStrings.size();
#endif
}

void MPIHelper::sendTree(int dest, string treeString, double score, int tag) {
    if (getNumProcesses() == 1 || dest == getProcessID())
        return;
#ifdef _IQTREE_MPI
    StrVector treeStrings;
    treeStrings.push_back(treeString);
    DoubleVector scores;
    scores.push_back(score);
    sendTrees(dest, treeStrings, scores, tag);
#endif
}

void MPIHelper::sendMsg(int tag, string msg) {
    if (getNumProcesses() == 1)
        return;
#ifdef _IQTREE_MPI
    if (tag == STOP_TAG)
        cleanUpMessages();
    for (int i = 0; i < getNumProcesses(); i++) {
        if (i != getProcessID()) {
            MPI_Request *request = new MPI_Request;
            ObjectStream *os = new ObjectStream(msg.c_str(), msg.size()+1);
            MPI_Isend(os->getObjectData(), os->getDataLength(), MPI_CHAR, i, tag, MPI_COMM_WORLD, request);
            sentMessages.push_back(make_pair(request, os));
        }
    }
#endif
}

bool MPIHelper::checkMsg(int tag, string &msg) {
    if (getNumProcesses() == 1)
        return true;
#ifdef _IQTREE_MPI
    int flag=0;
    MPI_Status status;
    char *recvBuffer;
    int numBytes;
    // Check for incoming messages
    MPI_Iprobe(PROC_MASTER, tag, MPI_COMM_WORLD, &flag, &status);
    // flag == true if there is a message
    if (flag) {
        MPI_Get_count(&status, MPI_CHAR, &numBytes);
        recvBuffer = new char[numBytes];
        MPI_Recv(recvBuffer, numBytes, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
        msg = recvBuffer;
        delete[] recvBuffer;
        return true;
    }
#endif
    return false;
}

bool MPIHelper::checkMsg(int tag) {
    if (getNumProcesses() == 1) {
        return false;
    }
#ifdef _IQTREE_MPI
    string msg;
    if (checkMsg(tag, msg)) {
        cout << "Worker " << getProcessID() << " gets message " << msg << endl;
        return true;
    }
#endif
    return false;
}


void MPIHelper::receiveTrees(bool fromAll, int maxNumTrees, TreeCollection &trees, int tag) {
    if (getNumProcesses() == 1) {
        return;
    }
#ifdef _IQTREE_MPI
    int flag = 0;
    int minNumTrees = 0;
    bool nodes[getNumProcesses()];
    if (fromAll)
        minNumTrees = getNumProcesses() - 1;
    for (int i = 0; i < getNumProcesses(); i++)
        nodes[i] = false;
    nodes[getProcessID()] = true;
    // Process all pending messages
    MPI_Status status;
    do {
        char* recvBuffer;
        int numBytes;
        // Check for incoming messages
        MPI_Iprobe(MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &flag, &status);
        // flag == true if there is a message
        if (flag) {
            //cout << "Getting messages from node " << status.MPI_SOURCE << endl;
            MPI_Get_count(&status, MPI_CHAR, &numBytes);
            recvBuffer = new char[numBytes];
            MPI_Recv(recvBuffer, numBytes, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
            ObjectStream os(recvBuffer, numBytes);
            if (status.MPI_TAG == STOP_TAG) {
                stringstream stopMsg;
                stopMsg << os.getObjectData();
                cout << stopMsg << endl;
                MPI_Finalize();
                exit(0);
            }
            TreeCollection curTrees = os.getTreeCollection();
            trees.addTrees(curTrees);
            if (trees.getNumTrees() >= maxNumTrees) {
                break;
            }
            if (fromAll && !nodes[status.MPI_SOURCE]) {
                nodes[status.MPI_SOURCE] = true;
                minNumTrees--;
            }
            delete [] recvBuffer;
        }
    } while (minNumTrees > 0 || flag);
    numTreeReceived += trees.getNumTrees();
#endif
}

int MPIHelper::receiveTrees(TreeCollection &trees, int tag) {
    if (getNumProcesses() == 1) {
        return -1;
    }
#ifdef _IQTREE_MPI
    int flag = 0;
    // Process all pending messages
    MPI_Status status;
    char* recvBuffer;
    int numBytes;
    // Check for incoming messages
    MPI_Iprobe(MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &flag, &status);
    // flag == true if there is a message
    if (flag) {
        //cout << "Getting messages from node " << status.MPI_SOURCE << endl;
        MPI_Get_count(&status, MPI_CHAR, &numBytes);
        recvBuffer = new char[numBytes];
        MPI_Recv(recvBuffer, numBytes, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
        ObjectStream os(recvBuffer, numBytes);
        TreeCollection curTrees = os.getTreeCollection();
        trees.addTrees(curTrees);
        delete [] recvBuffer;
        return status.MPI_SOURCE;
    }
#endif
    return -1;
}

int MPIHelper::cleanUpMessages() {
#ifdef _IQTREE_MPI
    int numMsgCleaned = 0;
    int flag = 0;
    MPI_Status status;
    vector< pair<MPI_Request*, ObjectStream*> >::iterator it;

    for (it = sentMessages.begin(); it != sentMessages.end(); ) {
        MPI_Test(it->first, &flag, &status);
        if (flag) {
            delete it->first;
            delete it->second;
            numMsgCleaned++;
            it = sentMessages.erase(it);
        } else {
            ++it;
        }
    }
    return numMsgCleaned;
#else
    return 0;
#endif
}


