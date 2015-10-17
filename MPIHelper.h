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
#define NOTREE_TAG 3 // Message does not contain a tree but a score only

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
     *      wait until at least one tree from each remaining process has been received
     *  @param maxNumTrees
     *      Only received up to maxNumTrees to prevent the function to block because it can constantly receive
     *      new trees
     *  @param trees[OUT]
     *      Trees received from other processes
     */
    void receiveTrees(bool fromAll, int maxNumTrees, TreeCollection &trees);

    /**
     *   @param treeStrings
     *   @param scores
     *   @param tag message tag
     */
    void sendTreesToOthers(vector<string> treeStrings, vector<double> scores, int tag = TREE_TAG);

    /**
    *   @param treeStrings
    *   @param scores
    *   @param tag message tag
    */
    void sendTreeToOthers(string treeString, double score, int tag);

    /**
     *  Send a stop message to other process
     *
     *  @param msg
     *      Some useful message about the reason for stopping
     */
    void sendStopMsg(string msg);

    bool checkStopMsg(string &msg);

private:
    /**
    *  Remove the buffers for finished messages
    */
    int cleanUpMessages();

#endif

private:
    MPIHelper() { }; // Disable constructor
    MPIHelper(MPIHelper const &) { }; // Disable copy constructor
    void operator=(MPIHelper const &) { }; // Disable assignment

    int processID;

    int numProcesses;

public:
    int getNumTreeReceived() const {
        return numTreeReceived;
    }

    void setNumTreeReceived(int numTreeReceived) {
        MPIHelper::numTreeReceived = numTreeReceived;
    }

    int getNumTreeSent() const {
        return numTreeSent;
    }

    void setNumTreeSent(int numTreeSent) {
        MPIHelper::numTreeSent = numTreeSent;
    }

private:
    int numTreeSent;

    int numTreeReceived;

public:
    int getNumNNISearch() const {
        return numNNISearch;
    }

    void setNumNNISearch(int numNNISearch) {
        MPIHelper::numNNISearch = numNNISearch;
    }

private:
    int numNNISearch;

#ifdef _IQTREE_MPI
    list<pair<MPI_Request *, ObjectStream *> > sentMessages;
#endif


};

#endif
