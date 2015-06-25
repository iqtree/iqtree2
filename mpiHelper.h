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
     *  Singleton method: get one and only one instance of the class
     */
public:
    static MPIHelper* instance();

    TreeCollection receiveTrees(int &MPISource);

    TreeCollection receiveTreesFromAnySource();

    void sendTreesToAll(TreeCollection& trees);

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
    int processID;

    int numProcesses;

    /**
     *  Private constructor
     */
    MPIHelper() {};

    /**
     *  Disable copy constructor from outside
     */
    MPIHelper(MPIHelper const&) {};

    /**
     *  Disable assignment constructor
     */
    MPIHelper& operator=(MPIHelper const&) {};

    static MPIHelper* MPIHelperInstance;

};

#endif //IQTREE_MPIHELPER_H
