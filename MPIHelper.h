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

#define MASTER 0

using namespace std;

class MPIHelper {

    /**
     *  Singleton method: get one and only one getInstance of the class
     */
public:
    static MPIHelper& getInstance();

    CandidateSet receiveTrees(int &MPISource);

    /**
     *  Get all the trees that have been sent to the current node
     */
    TreeCollection getTreesFromOtherNodes(int numNodes = 1);

    /**
     *  Send a set of candidate trees to all remaining nodes
     */
    void sendTreesToAllNodes(TreeCollection &trees, bool blocking);

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
};

#endif //IQTREE_MPIHELPER_H
