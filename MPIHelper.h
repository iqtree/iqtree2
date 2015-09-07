//
// Created by tung on 6/18/15.
//

#ifndef IQTREE_MPIHELPER_H
#define IQTREE_MPIHELPER_H
#include <string>
#include <vector>
#include <mpi.h>
#include "tools.h"
#include "TreeCollection.h"
#include "ObjectStream.h"

#define MASTER 0

using namespace std;

class MPIHelper {

    /**
     *  Singleton method: get one and only one getInstance of the class
     */
public:
    static MPIHelper& getInstance();

    /**
     *  Get all the trees that have been sent to the current node
     *  @param fromAll wait for messages from all other processes
     */
    TreeCollection getTreesForMe(bool fromAll);


    /**
     *  Send a set of candidate trees to all remaining nodes
     */
    void sendTreesToOthers(TreeCollection &trees);

    /**
     *  Remove the buffers for finished messages
     */
    int cleanUpMessages();

    int getNumProcesses() const {
        return numProcesses;
    }

    void setNumProcesses(int numProcesses) {
        MPIHelper::numProcesses = numProcesses;
    }

    int getProcessID() const {
        return processID;
    }

    void setProcessID(int processID) {
        MPIHelper::processID = processID;
    }

private:
    MPIHelper() {}; // Disable constructor
    MPIHelper (MPIHelper const&) {}; // Disable copy constructor
    void operator=(MPIHelper const&) {}; // Disable assignment

    int processID;

    int numProcesses;

    list< pair<MPI_Request*, ObjectStream*> > messages;
};

#endif //IQTREE_MPIHELPER_H
