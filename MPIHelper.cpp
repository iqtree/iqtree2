//
// Created by tung on 6/18/15.
//

#include "MPIHelper.h"
#include "timeutil.h"

/**
 *  Initialize the single getInstance of MPIHelper
 */

MPIHelper& MPIHelper::getInstance() {
    static MPIHelper instance;
#ifndef _IQTREE_MPI
    instance.setProcessID(0);
    instance.setNumProcesses(1);
#endif
    return instance;
}

void MPIHelper::distributeTrees(vector<string> &treeStrings, vector<double> &scores, int tag) {
    if (getNumProcesses() == 1)
        return;
#ifdef _IQTREE_MPI
    vector<int> sourceProcID;
    sourceProcID.insert(sourceProcID.end(), scores.size(), getProcessID());
    TreeCollection outTrees(treeStrings, scores, sourceProcID);
    cleanUpMessages();
    for (int i = 0; i < getNumProcesses(); i++) {
        if (i != getProcessID()) {
            MPI_Request *request = new MPI_Request;
            ObjectStream *os = new ObjectStream(outTrees);
            MPI_Isend(os->getObjectData(), os->getDataLength(), MPI_CHAR, i, tag, MPI_COMM_WORLD, request);
            sentMessages.push_back(make_pair(request, os));
            int flag = 0;
            MPI_Status status;
            MPI_Test(request, &flag, &status);
        }
    }
    //numTreeSent += treeStrings.size();
#endif
}

void MPIHelper::distributeTree(string treeString, double score, int tag) {
    if (getNumProcesses() == 1)
        return;
#ifdef _IQTREE_MPI
    double start = getRealTime();
    vector<string> trees;
    vector<double> scores;
    trees.push_back(treeString);
    scores.push_back(score);
    distributeTrees(trees, scores, tag);
    if (verbose_mode >= VB_MED)
        cout << "Sent tree to other processes in " << getRealTime() - start << " seconds" << endl;
    numTreeSent++;
#endif
}

void MPIHelper::sendTrees(int dest, vector<string> &treeStrings, vector<double> &scores, int tag) {
    if (getNumProcesses() == 1 || dest == getProcessID())
        return;
#ifdef _IQTREE_MPI
    vector<int> sourceProcID;
    sourceProcID.insert(sourceProcID.end(), scores.size(), getProcessID());
    TreeCollection outTrees(treeStrings, scores, sourceProcID);
    cleanUpMessages();
    MPI_Request *request = new MPI_Request;
    ObjectStream *os = new ObjectStream(outTrees);
    MPI_Isend(os->getObjectData(), os->getDataLength(), MPI_CHAR, dest, tag, MPI_COMM_WORLD, request);
    sentMessages.push_back(make_pair(request, os));
    numTreeSent += treeStrings.size();

    int flag = 0;
    MPI_Status status;
    MPI_Test(request, &flag, &status);
#endif
}

void MPIHelper::sendTree(int dest, string treeString, double score, int tag) {
    if (getNumProcesses() == 1 || dest == getProcessID())
        return;
#ifdef _IQTREE_MPI
    StrVector treeStrings;
    treeStrings.push_back(treeString);
    DoubleVector scores;
    scores.push_back(score);
    sendTrees(dest, treeStrings, scores, tag);
#endif
}

int MPIHelper::sendRecvTrees(int dest, vector<string> &treeStrings, vector<double> &scores, int tag) {
    if (getNumProcesses() == 1 || dest == getProcessID())
        return tag;
#ifdef _IQTREE_MPI
    double beginTime = getRealTime();
    // prepare message
    vector<int> sourceProcID;
    sourceProcID.insert(sourceProcID.end(), scores.size(), getProcessID());
    TreeCollection outTrees(treeStrings, scores, sourceProcID);
    ObjectStream *os = new ObjectStream(outTrees);

    // blocking send
    MPI_Send(os->getObjectData(), os->getDataLength(), MPI_CHAR, dest, tag, MPI_COMM_WORLD);
    numTreeSent += treeStrings.size();
    delete os;

    // blocking probe
    MPI_Status status;
    MPI_Probe(dest, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    int msgCount;
    MPI_Get_count(&status, MPI_CHAR, &msgCount);

    // receive the message
    char *recvBuffer = new char[msgCount];
    MPI_Recv(recvBuffer, msgCount, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
    treeStrings.clear();
    scores.clear();

    if (status.MPI_TAG != STOP_TAG) {
        os = new ObjectStream(recvBuffer, msgCount);
        TreeCollection curTrees = os->getTreeCollection();
        treeStrings = curTrees.getTreeStrings();
        scores = curTrees.getScores();
        numTreeReceived += treeStrings.size();
    }
    delete [] recvBuffer;

    double endTime = getRealTime();
    cout << "INFO: " << endTime - beginTime << " seconds for " << __func__ << endl;

    return status.MPI_TAG;
#endif
}

int MPIHelper::recvSendTrees(vector<string> &treeStrings, vector<double> &scores, vector<bool> &should_send, int tag) {
    if (getNumProcesses() == 1)
        return 0;
#ifdef _IQTREE_MPI
    double beginTime = getRealTime();
    // blocking probe
    MPI_Status status;
    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    int msgCount;
    MPI_Get_count(&status, MPI_CHAR, &msgCount);
    int dest = status.MPI_SOURCE;

    // receive the message
    char *recvBuffer = new char[msgCount];
    MPI_Recv(recvBuffer, msgCount, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);

    // now send message
    if (!should_send[dest]) {
        treeStrings.resize(1, "notree");
        scores.resize(1, -DBL_MAX);
    }
    IntVector sourceProcID;
    sourceProcID.insert(sourceProcID.end(), scores.size(), getProcessID());
    TreeCollection outTrees(treeStrings, scores, sourceProcID);
    ObjectStream *os = new ObjectStream(outTrees);

    // blocking send
    MPI_Send(os->getObjectData(), os->getDataLength(), MPI_CHAR, dest, tag, MPI_COMM_WORLD);
    numTreeSent += treeStrings.size();
    delete os;

    // now extract trees from received buffer
    treeStrings.clear();
    scores.clear();
    os = new ObjectStream(recvBuffer, msgCount);
    TreeCollection curTrees = os->getTreeCollection();
    treeStrings = curTrees.getTreeStrings();
    scores = curTrees.getScores();
    delete [] recvBuffer;
    numTreeReceived += treeStrings.size();
    
    should_send[dest] = false;

    double endTime = getRealTime();
    if (endTime - beginTime > 1)
        cout << "WARNING: " << endTime - beginTime << " seconds for " << __func__ << endl;

    return dest;
#endif
}

void MPIHelper::gatherTrees(TreeCollection &trees) {
    if (getNumProcesses() == 1)
        return;
#ifdef _IQTREE_MPI
    double beginTime = getRealTime();

    if (isMaster()) {
        trees.clear();
        // Master: receive from all Workers
        for (int w = 1; w < getNumProcesses(); w++) {
            // blocking probe
            MPI_Status status;
            MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            int msgCount;
            MPI_Get_count(&status, MPI_CHAR, &msgCount);
            // receive the message
            char *recvBuffer = new char[msgCount];
            MPI_Recv(recvBuffer, msgCount, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
            ObjectStream *os = new ObjectStream(recvBuffer, msgCount);
            TreeCollection curTrees = os->getTreeCollection();
            trees.addTrees(curTrees);
            numTreeReceived += curTrees.getNumTrees();
            delete [] recvBuffer;
        }
        cout << trees.getNumTrees() << " trees gathered from workers in ";
    } else {
        // Worker: send trees to Master
        ObjectStream *os = new ObjectStream(trees);
        // blocking send
        MPI_Send(os->getObjectData(), os->getDataLength(), MPI_CHAR, PROC_MASTER, TREE_TAG, MPI_COMM_WORLD);
        numTreeSent += trees.getNumTrees();
        delete os;
        cout << trees.getNumTrees() << " trees sent to master in ";
    }

    double endTime = getRealTime();
    cout << endTime - beginTime << " seconds" << endl;
#endif
}

void MPIHelper::broadcastTrees(TreeCollection &trees) {
    if (getNumProcesses() == 1)
        return;
#ifdef _IQTREE_MPI
    double beginTime = getRealTime();

    // prepare data from Master
    ObjectStream *os;
    int msgCount = 0;
    if (isMaster()) {
        os = new ObjectStream(trees);
        msgCount = os->getDataLength();
    }

    // broadcast the count for workers
    MPI_Bcast(&msgCount, 1, MPI_INT, PROC_MASTER, MPI_COMM_WORLD);

    char *recvBuffer = new char[msgCount];
    if (isMaster())
        memcpy(recvBuffer, os->getObjectData(), msgCount);

    // broadcast trees to workers
    MPI_Bcast(recvBuffer, msgCount, MPI_CHAR, PROC_MASTER, MPI_COMM_WORLD);

    if (isWorker()) {
        os = new ObjectStream(recvBuffer, msgCount);
        trees = os->getTreeCollection();
    }
    delete os;
    delete [] recvBuffer;

    double endTime = getRealTime();
    cout << trees.getNumTrees() << " trees broadcasted to workers in " << endTime - beginTime << " seconds" << endl;

#endif
}


bool MPIHelper::gotMessage() {
    // Check for incoming messages
    if (getNumProcesses() == 1)
        return false;
#ifdef _IQTREE_MPI
    int flag = 0;
    MPI_Status status;
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
    if (flag)
        return true;
    else
        return false;
#endif
}

void MPIHelper::sendMsg(int tag, string msg) {
    if (getNumProcesses() == 1)
        return;
#ifdef _IQTREE_MPI
    if (tag == STOP_TAG)
        cleanUpMessages();
    for (int i = 0; i < getNumProcesses(); i++) {
        if (i != getProcessID()) {
            MPI_Request *request = new MPI_Request;
            ObjectStream *os = new ObjectStream(msg.c_str(), msg.size()+1);
            MPI_Isend(os->getObjectData(), os->getDataLength(), MPI_CHAR, i, tag, MPI_COMM_WORLD, request);
            sentMessages.push_back(make_pair(request, os));
            int flag = 0;
            MPI_Status status;
            MPI_Test(request, &flag, &status);
        }
    }
#endif
}

bool MPIHelper::checkMsg(int tag, string &msg) {
    if (getNumProcesses() == 1)
        return true;
#ifdef _IQTREE_MPI
    int flag=0;
    MPI_Status status;
    char *recvBuffer;
    int numBytes;
    // Check for incoming messages
    MPI_Iprobe(PROC_MASTER, tag, MPI_COMM_WORLD, &flag, &status);
    // flag == true if there is a message
    if (flag) {
        MPI_Get_count(&status, MPI_CHAR, &numBytes);
        recvBuffer = new char[numBytes];
        MPI_Recv(recvBuffer, numBytes, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
        msg = recvBuffer;
        delete[] recvBuffer;
        return true;
    }
#endif
    return false;
}

bool MPIHelper::checkMsg(int tag) {
    if (getNumProcesses() == 1) {
        return false;
    }
#ifdef _IQTREE_MPI
    string msg;
    if (checkMsg(tag, msg)) {
        cout << "Worker " << getProcessID() << " gets message " << msg << endl;
        return true;
    }
#endif
    return false;
}


void MPIHelper::receiveTrees(bool fromAll, int maxNumTrees, TreeCollection &trees, int tag) {
    if (getNumProcesses() == 1) {
        return;
    }
#ifdef _IQTREE_MPI
    int flag = 0;
    int minNumTrees = 0;
    bool nodes[getNumProcesses()];
    if (fromAll)
        minNumTrees = getNumProcesses() - 1;
    for (int i = 0; i < getNumProcesses(); i++)
        nodes[i] = false;
    nodes[getProcessID()] = true;
    // Process all pending messages
    MPI_Status status;
    size_t totalMsgSize = 0;
    do {
        char* recvBuffer;
        int numBytes;
        flag = 0;
        // Check for incoming messages
        MPI_Iprobe(MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &flag, &status);
        // flag == true if there is a message
        if (flag) {
            //cout << "Getting messages from node " << status.MPI_SOURCE << endl;
            MPI_Get_count(&status, MPI_CHAR, &numBytes);
            totalMsgSize += numBytes;
            recvBuffer = new char[numBytes];
            MPI_Recv(recvBuffer, numBytes, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
            ObjectStream os(recvBuffer, numBytes);
            if (status.MPI_TAG == STOP_TAG) {
                stringstream stopMsg;
                stopMsg << os.getObjectData();
                cout << stopMsg << endl;
                MPI_Finalize();
                exit(0);
            }
            TreeCollection curTrees = os.getTreeCollection();
            trees.addTrees(curTrees);
            if (trees.getNumTrees() >= maxNumTrees) {
                break;
            }
            if (fromAll && !nodes[status.MPI_SOURCE]) {
                nodes[status.MPI_SOURCE] = true;
                minNumTrees--;
            }
            delete [] recvBuffer;
        }
    } while (minNumTrees > 0 || flag);
    numTreeReceived += trees.getNumTrees();
    if (trees.getNumTrees() > 0) {
        cout << "Proc " << getProcessID() << ": " << trees.getNumTrees() << " trees received from other processes (" << totalMsgSize << " bytes)" << endl;
    }
#endif
}

int MPIHelper::receiveTrees(TreeCollection &trees, int tag) {
    if (getNumProcesses() == 1) {
        return -1;
    }
#ifdef _IQTREE_MPI
    int flag = 0;
    // Process all pending messages
    MPI_Status status;
    char* recvBuffer;
    int numBytes;
    // Check for incoming messages
    MPI_Iprobe(MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &flag, &status);
    // flag == true if there is a message
    if (flag) {
        //cout << "Getting messages from node " << status.MPI_SOURCE << endl;
        MPI_Get_count(&status, MPI_CHAR, &numBytes);
        recvBuffer = new char[numBytes];
        MPI_Recv(recvBuffer, numBytes, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
        ObjectStream os(recvBuffer, numBytes);
        TreeCollection curTrees = os.getTreeCollection();
        trees.addTrees(curTrees);
        delete [] recvBuffer;
        return status.MPI_SOURCE;
    }
#endif
    return -1;
}

int MPIHelper::cleanUpMessages() {
#ifdef _IQTREE_MPI
    int numMsgCleaned = 0;
    // change iterator to index to avoid iterator being invalidated after erase()
    for (int i = 0; i < sentMessages.size(); ) {
        int flag = 0;
        MPI_Status status;
        MPI_Test(sentMessages[i].first, &flag, &status);
        if (flag) {
            delete sentMessages[i].first;
            delete sentMessages[i].second;
            numMsgCleaned++;
            sentMessages.erase(sentMessages.begin()+i);
        } else {
            i++;
        }
    }
    if (verbose_mode >= VB_MED && numMsgCleaned)
        cout << numMsgCleaned << " messages sent and cleaned up" << endl;
    return numMsgCleaned;
#else
    return 0;
#endif
}


MPIHelper::~MPIHelper() {
//    cleanUpMessages();
}

