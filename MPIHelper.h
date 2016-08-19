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
#define BOOT_TREE_TAG 4 // bootstrap tree tag

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

    bool isMaster() const {
        return processID == MASTER;
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
     *  @param tag MPI tag
     */
    void receiveTrees(bool fromAll, int maxNumTrees, TreeCollection &trees, int tag);


    /**
     *  Receive trees that sent to the current process
     *
     *  @param trees[OUT]
     *      Trees received from other processes
     *  @param tag MPI tag
     *  @return source process ID
     */
    int receiveTrees(TreeCollection &trees, int tag);

    /**
     *   Send trees to all other processes
     *   @param treeStrings vector of trees
     *   @param scores vector containing scores of the trees with same order as in treeStrings
     *   @param tag used to classified the message
     */
    void distributeTrees(vector<string> treeStrings, vector<double> scores, int tag = TREE_TAG);

    /**
    *   Similar to distributeTrees but only 1 tree is sent
    *   @param treeString
    *   @param score
    *   @param tag
    */
    void distributeTree(string treeString, double score, int tag);

    /**
     *   Send trees to a dest process
     *   @param dest MPI rank of destination process
     *   @param treeStrings vector of trees
     *   @param scores vector containing scores of the trees with same order as in treeStrings
     *   @param tag used to classified the message
     */
    void sendTrees(int dest, vector<string> treeStrings, vector<double> scores, int tag);

    /**
     *  Send a stop message to other process
     *
     *  @param msg
     *      Some useful message about the reason for stopping
     */
    void sendStopMsg();

    bool checkStopMsg();

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
    
    void resetNumbers() {
        numTreeSent = 0;
        numTreeReceived = 0;
        numNNISearch = 0;
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
    // A list storing messages and the corresponding requests that have been sent from the current process.
    // When a message has been successfully received, it will be deleted from the list
    vector< pair<MPI_Request *, ObjectStream *> > sentMessages;
#endif


};

#endif
