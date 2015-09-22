//
// Created by tung on 6/18/15.
//

#include "MPIHelper.h"

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
void MPIHelper::sendTreesToOthers(TreeCollection trees, int tag) {
    if (getNumProcesses() == 1)
        return;
    cleanUpMessages();
    for (int i = 0; i < getNumProcesses(); i++) {
        if (i != getProcessID()) {
            MPI_Request *request = new MPI_Request;
            ObjectStream *os = new ObjectStream(trees);
            MPI_Isend(os->getObjectData(), os->getDataLength(), MPI_CHAR, i, tag, MPI_COMM_WORLD, request);
            sentMessages.push_back(make_pair(request, os));
        }
    }
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

bool MPIHelper::receiveTrees(bool fromAll, TreeCollection &trees) {
    if (getNumProcesses() == 1) {
        return false;
    }
    bool stop = false;
    int flag, totalMsg;
    bool nodes[getNumProcesses()];
    if (fromAll)
        totalMsg = getNumProcesses() - 1;
    else
        totalMsg = 0;
    for (int i = 0; i < getNumProcesses(); i++)
        nodes[i] = false;
    nodes[getProcessID()] = true;
    // Process all pending messages
    do {
        MPI_Status status;
        char* recvBuffer;
        int numBytes;
        // Check for incoming messages
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
        // flag == true if there is a message
        if (flag) {
            //cout << "Getting messages from node " << status.MPI_SOURCE << endl;
            MPI_Get_count(&status, MPI_CHAR, &numBytes);
            if (status.MPI_SOURCE == MASTER && status.MPI_TAG == STOP_TAG) {
                stop = true;
                break;
            }
            recvBuffer = new char[numBytes];
            MPI_Recv(recvBuffer, numBytes, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, NULL);
            ObjectStream os(recvBuffer, numBytes);
            TreeCollection curTrees = os.getTreeCollection();
            trees.addTrees(curTrees);
            if (fromAll && !nodes[status.MPI_SOURCE]) {
                nodes[status.MPI_SOURCE] = true;
                totalMsg--;
            }
            delete [] recvBuffer;
        }
    } while (totalMsg != 0 || flag);
    return stop;
}

int MPIHelper::cleanUpMessages() {
    int numMsgCleaned = 0;
    int flag = 0;
    MPI_Status status;
    list< pair<MPI_Request*, ObjectStream*> >::iterator it;
    it = sentMessages.begin();
    while ( it != sentMessages.end() ) {
        MPI_Test(it->first, &flag, &status);
        if (flag) {
            delete it->first;
            delete it->second;
            numMsgCleaned++;
            sentMessages.erase(it++);
        } else {
            ++it;
        }
    }
    return numMsgCleaned;
}
#endif


