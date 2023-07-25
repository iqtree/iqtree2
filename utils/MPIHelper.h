/***************************************************************************
 *   Copyright (C) 2015 by                                                 *
 *   Lam-Tung Nguyen <nltung@gmail.com>                                    *
 *                                                                         *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef MPIHELPER_H
#define MPIHELPER_H

#include <string>
#include <vector>
#include "utils/tools.h"
#include "utils/checkpoint.h"

#ifdef _IQTREE_MPI
#include <mpi.h>
#endif

#define PROC_MASTER 0
#define TREE_TAG 1 // Message contain trees
#define STOP_TAG 2 // Stop message
#define BOOT_TAG 3 // Message to please send bootstrap trees
#define BOOT_TREE_TAG 4 // bootstrap tree tag
#define LOGL_CUTOFF_TAG 5 // send logl_cutoff for ultrafast bootstrap

using namespace std;

class MPIHelper {
public:
    /**
    *  Singleton method: get one and only one getInstance of the class
    */
    static MPIHelper &getInstance();

    /** initialize MPI */
    void init(int argc, char *argv[]);
    
    /** finalize MPI */
    void finalize();
    
    /**
        destructor
    */
    ~MPIHelper();

    int getNumProcesses() const {
        return numProcesses;
    }

    void setNumProcesses(int numProcesses) {
        MPIHelper::numProcesses = numProcesses;
    }

    int getProcessID() const {
        return processID;
    }

    bool isMaster() const {
        return processID == PROC_MASTER;
    }

    bool isWorker() const {
        return processID != PROC_MASTER;
    }

    void setProcessID(int processID) {
        MPIHelper::processID = processID;
    }

    void barrier() {
        #ifdef _IQTREE_MPI
            MPI_Barrier(MPI_COMM_WORLD);
        #endif
    }

    /** synchronize random seed from master to all workers */
    void syncRandomSeed();
    
    /** count the number of host with the same name as the current host */
    int countSameHost();

    /** @return true if got any message from another process */
    bool gotMessage();


    /** wrapper for MPI_Send a string
        @param str string to send
        @param dest destination process
        @param tag message tag
    */

#ifdef _IQTREE_MPI
    void sendString(string &str, int dest, int tag);

    /** wrapper for MPI_Recv a string
        @param[out] str string received
        @param src source process
        @param tag message tag
        @return the source process that sent the message
    */
    int recvString(string &str, int src = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG);

    /** wrapper for MPI_Send an entire Checkpoint object
        @param ckp Checkpoint object to send
        @param dest destination process
    */
    void sendCheckpoint(Checkpoint *ckp, int dest);

    /** wrapper for MPI_Recv an entire Checkpoint object
        @param[out] ckp Checkpoint object received
        @param src source process
        @param tag message tag
        @return the source process that sent the message
    */
    int recvCheckpoint(Checkpoint *ckp, int src = MPI_ANY_SOURCE);

    /**
        wrapper for MPI_Bcast to broadcast checkpoint from Master to all Workers
        @param ckp Checkpoint object
    */
    void broadcastCheckpoint(Checkpoint *ckp);

    /**
        wrapper for MPI_Gather to gather all checkpoints into Master
        @param ckp Checkpoint object
    */
    void gatherCheckpoint(Checkpoint *ckp);
#endif

    void increaseTreeSent(int inc = 1) {
        numTreeSent += inc;
    }

    void increaseTreeReceived(int inc = 1) {
        numTreeReceived += inc;
    }

private:
    /**
    *  Remove the buffers for finished messages
    */
    int cleanUpMessages();

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
    
    void resetNumbers() {
//        numTreeSent = 0;
//        numTreeReceived = 0;
//        numNNISearch = 0;
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


};

#endif
