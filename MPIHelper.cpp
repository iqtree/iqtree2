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

#ifdef _IQTREE_MPI
void MPIHelper::distributeTrees(vector<string> treeStrings, vector<double> scores, int tag) {
    TreeCollection outTrees(treeStrings, scores);
    if (getNumProcesses() == 1)
        return;
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
}

void MPIHelper::distributeTree(string treeString, double score, int tag) {
    double start = getRealTime();
    vector<string> trees;
    vector<double> scores;
    trees.push_back(treeString);
    scores.push_back(score);
//    cout << "Sent tree to other processes in ";
    MPIHelper::getInstance().sendTreesToOthers(trees, scores, tag);
//    cout << getRealTime() - st    art << " seconds" << endl;
    numTreeSent++;
}

void MPIHelper::sendStopMsg(string msg) {
    if (getNumProcesses() == 1)
        return;
    cleanUpMessages();
    for (int i = 0; i < getNumProcesses(); i++) {
        if (i != getProcessID()) {
            MPI_Request *request = new MPI_Request;
            ObjectStream *os = new ObjectStream(msg.c_str(), msg.size()+1);
            MPI_Isend(os->getObjectData(), os->getDataLength(), MPI_CHAR, i, STOP_TAG, MPI_COMM_WORLD, request);
            sentMessages.push_back(make_pair(request, os));
        }
    }
}

bool MPIHelper::checkStopMsg(string &msg) {
    int flag;
    MPI_Status status;
    char *recvBuffer;
    int numBytes;
    // Check for incoming messages
    MPI_Iprobe(MPI_ANY_SOURCE, STOP_TAG, MPI_COMM_WORLD, &flag, &status);
    // flag == true if there is a message
    if (flag) {
        //cout << "Getting messages from node " << status.MPI_SOURCE << endl;
        MPI_Get_count(&status, MPI_CHAR, &numBytes);
        recvBuffer = new char[numBytes];
        MPI_Recv(recvBuffer, numBytes, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, NULL);
        ObjectStream os(recvBuffer, numBytes);
        string strMsg(os.getObjectData());
        cout << strMsg << endl;
        msg = strMsg;
        delete[] recvBuffer;
        return true;
    }
}

void MPIHelper::receiveTrees(bool fromAll, int maxNumTrees, TreeCollection &trees) {
    if (getNumProcesses() == 1) {
        return;
    }
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
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
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
}

int MPIHelper::cleanUpMessages() {
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
}
#endif


