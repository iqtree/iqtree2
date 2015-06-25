//
// Created by tung on 6/18/15.
//

#include "mpiHelper.h"
#include "ObjectStream.h"

/**
 *  Initialize the single instance of MPIHelper
 */
MPIHelper* MPIHelper::MPIHelperInstance = NULL;

MPIHelper* MPIHelper::instance() {
    if (!MPIHelperInstance) {
        MPIHelperInstance = new MPIHelper();
    }
    return MPIHelperInstance;
}

void MPIHelper::sendTreesToAll(TreeCollection& trees) {

}

TreeCollection MPIHelper::receiveTreesFromAnySource() {
    MPI_Status status;
    char* recvBuffer;
    int MPISource, flag;
    int numBytes;
    TreeCollection allTrees;
    // Process all pending messages
    do {
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
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






