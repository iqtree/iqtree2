//
// Created by tung on 6/18/15.
//

#include "MPIHelper.h"
#include "ObjectStream.h"

/**
 *  Initialize the single getInstance of MPIHelper
 */

MPIHelper& MPIHelper::getInstance() {
    static MPIHelper instance;
    return instance;
}

void MPIHelper::sendTreesToAll(CandidateSet& candidateTrees, bool blocking) {

}

void MPIHelper::sendTreesToNode(int nodeID, CandidateSet& candidateTrees, bool blocking) {

}

CandidateSet MPIHelper::receiveTreesFromAll() {
    MPI_Status status;
    char* recvBuffer;
    int MPISource, flag;
    int numBytes;
    TreeCollection allTrees;
    // Process all pending messages
    do {
        // Check for incoming messages
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
        // flag == true if there is a message
        if (flag) {
            MPISource = status.MPI_SOURCE;
            MPI_Get_count(&status, MPI_CHAR, &numBytes);
            recvBuffer = new char[numBytes];
            MPI_Recv(recvBuffer, numBytes, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, NULL);
            ObjectStream os(recvBuffer, numBytes);
            TreeCollection curTrees;
            os.readObject(curTrees);
            allTrees.addTrees(curTrees);
            delete [] recvBuffer;
        }
    } while (flag);
    return allTrees;
}


void MPIHelper::sendTreesToNode(int nodeID, TreeCollection &trees) {

}
