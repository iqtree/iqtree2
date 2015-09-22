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
#define TREE_TAG 1 // Message contain trees
#define STOP_TAG 2 // Stop message

using namespace std;

class MPIHelper {
public:
    /**
    *  Singleton method: get one and only one getInstance of the class
    */
    static MPIHelper &getInstance();

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
     *  Receive trees that sent to the current process
     *
     *  @param fromAll
     *      wait for messages from all other processes
     *  @param trees[OUT]
     *      Trees received from other processes
     *  @return
     *      Whether or not the search is finished (only relevant for non-master processes)
     */
    bool receiveTrees(bool fromAll, TreeCollection &trees);


    /**
     *  Send a set of trees with scores to other nodes
     *
     *  @param trees a TreeCollection object
     *  @param tag tag used with the message
     */
    void sendTreesToOthers(TreeCollection trees, int tag);

    /**
     *  Send a stop message to other process
     *
     *  @param msg
     *      Some useful message about the reason for stopping
     */
    void sendStopMsg(string msg);

    /**
     *  Remove the buffers for finished messages
     */
    int cleanUpMessages();

#endif


        private:
        MPIHelper()
        { }; // Disable constructor
        MPIHelper(MPIHelper const&) { }; // Disable copy constructor
        void operator=(MPIHelper const &) { }; // Disable assignment

        int processID;

        int numProcesses;

#ifdef _IQTREE_MPI
        list<pair<MPI_Request *, ObjectStream *> > sentMessages;
#endif

    };

#endif
