//
// Created by tung on 6/18/15.
//

#include "MPIHelper.h"

/**
 *  Initialize the single getInstance of MPIHelper
 */

MPIHelper& MPIHelper::getInstance() {
    static MPIHelper instance;
    return instance;
}

void MPIHelper::sendTreesToOthers(TreeCollection &trees, int tag) {
    cleanUpMessages();
    for (int i = 0; i < getNumProcesses(); i++) {
        if (i != getProcessID()) {
            MPI_Request *request = new MPI_Request;
            ObjectStream *os = new ObjectStream(trees);
            MPI_Isend(os->getObjectData(), os->getDataLength(), MPI_CHAR, i, tag, MPI_COMM_WORLD, request);
            messages.push_back(make_pair(request, os));
        }
    }
}

TreeCollection MPIHelper::getTreesForMe(bool fromAll) {
    TreeCollection allTrees;
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
                cout << "Stop signal received. Process is now being ended." << endl;
                MPI_Finalize();
                exit(0);
            }
            recvBuffer = new char[numBytes];
            MPI_Recv(recvBuffer, numBytes, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, NULL);
            ObjectStream os(recvBuffer, numBytes);
            TreeCollection curTrees = os.getTreeCollection();
            allTrees.addTrees(curTrees);
            if (fromAll && !nodes[status.MPI_SOURCE]) {
                nodes[status.MPI_SOURCE] = true;
                totalMsg--;
            }
            delete [] recvBuffer;
        }
    } while (totalMsg != 0 || flag);
    return allTrees;
}

int MPIHelper::cleanUpMessages() {
    int numMsgCleaned = 0;
    int flag = 0;
    MPI_Status status;
    list< pair<MPI_Request*, ObjectStream*> >::iterator it;
    it = messages.begin();
    while ( it != messages.end() ) {
        MPI_Test(it->first, &flag, &status);
        if (flag) {
            delete it->first;
            delete it->second;
            numMsgCleaned++;
            messages.erase(it++);
        } else {
            ++it;
        }
    }
    return numMsgCleaned;
}


