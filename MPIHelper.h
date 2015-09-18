//
// Created by tung on 6/18/15.
//

#ifndef MPIHELPER_H
#define MPIHELPER_H
#include <string>
#include <vector>
#include "tools.h"
#include "TreeCollection.h"
#include "ObjectStream.h"
#ifdef _IQTREE_MPI
#include <mpi.h>
#endif

#define MASTER 0
#define CNT_TAG 0
#define STOP_TAG 1

using namespace std;

class MPIHelper {

    /**
     *  Singleton method: get one and only one getInstance of the class
     */
public:

    static MPIHelper& getInstance();

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

#ifdef _IQTREE_MPI

    /**
     *  Get all the trees that have been sent to the current node
     *  @param fromAll wait for messages from all other processes
     */
    TreeCollection getTreesForMe(bool fromAll);


    /**
     *  Send a set of trees with scores to other nodes
     *  @param trees a TreeCollection object
     *  @param tag tag used with the message
     */
    void sendTreesToOthers(TreeCollection &trees, int tag);


    /**
     *  Remove the buffers for finished messages
     */
    int cleanUpMessages();
#endif


private:
    MPIHelper() {}; // Disable constructor
    MPIHelper (MPIHelper const&) {}; // Disable copy constructor
    void operator=(MPIHelper const&) {}; // Disable assignment

    int processID;

    int numProcesses;

#ifdef _IQTREE_MPI
    list< pair<MPI_Request*, ObjectStream*> > messages;
#endif

};

#endif
