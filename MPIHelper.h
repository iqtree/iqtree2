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

    TreeCollection receiveTrees(int &MPISource);

    TreeCollection receiveTreesFromAnySource();

    void sendTreesToAll(TreeCollection& trees);

    void sendTreesToNode(int nodeID, TreeCollection& trees);

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
