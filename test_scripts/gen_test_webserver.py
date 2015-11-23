#!/usr/bin/env python
'''
Created on Feb. 01, 2015
This script collects all commands submited by users of the IQ-TREE web service 
that were crashed
Default location: /project/web-iqtree/user-data 

@author: Tung Nguyen
@email: nltung@gmai.com
'''
import sys, os, time, multiprocessing, optparse, fnmatch 
import subprocess, logging, datetime
import cmd
import shutil
from operator import itemgetter
from cmd import Cmd

def collect_logs(log_dir):
    bugLogs = []
    numLog = 0
    for root, dirnames, filenames in os.walk(log_dir):
        for filename in fnmatch.filter(filenames, '*.log'):
            numLog = numLog + 1
            logFile = os.path.join(root,filename)
            if 'CRASH' in open(logFile).read():
                bugLogs.append(logFile)
    return (bugLogs, numLog)

def collect_cmds(logFiles, out_dir):
    if not os.path.exists(out_dir):
        #shutil.rmtree(options.out_dir)
        os.makedirs(out_dir)
    runs = []
    for log in logFiles:
        id = log.split('/')[-2]
        email = log.split('/')[-3]
        with open(log) as f:            
            for line in f:
                if line.startswith('Command:'):
                    cmd = " ".join(line.split()[2:])
                    aln = line.split()[3]
                    run_dir = os.path.abspath(os.path.join(log, os.pardir))
                    shutil.copy2(os.path.join(run_dir, aln), out_dir)
                    if line.find('-sp') != -1:
                        partitionFile = line.split()[5]
                        shutil.copy2(os.path.join(run_dir, partitionFile), out_dir)
                if line.startswith('Seed:'):
                    seed = line.split()[1]
                    runs.append((int(id), email,seed, aln, cmd))
                    #print run_dir
                    break
    return runs

def create_test_cmds(runs, iqtree_binary, filename):
    outfile = open(filename, "wb")
    for run in runs:
        run_id = run[1] + "_" + str((run[0]))
        seed = run[2]
        args = run[4]
        if not " -b " in args:
            cmd = run_id + ' ' + iqtree_binary + ' ' + args + ' -seed ' + seed + ' -pre ' + run_id
            print >> outfile, cmd
    outfile.close()
                                        
if __name__ == '__main__':
    usage = "USAGE: %prog [options]"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-d', '--indir', dest="in_dir", 
                      help='Path to user directory of the IQ-TREE web server [default: %default]', default="/project/web-iqtree/user-data/")
    parser.add_option('-o', '--outdir', dest="out_dir", 
                      help='Directory containing alignments [default: %default]', default="webserver_alignments")
    parser.add_option('-b', '--iqtree_bin', dest='iqtree_bin', help='Path to IQ-Tree binary')
    (options, args) = parser.parse_args()
    if not options.iqtree_bin:
        print "Please specify the path to your IQ-TREE binary"
        parser.print_help()
        exit(0)
    print "Collecting buggy runs from " + options.in_dir
    (bugLogs, numLog) = collect_logs(options.in_dir)
    print ("Found %d job submissions, %d of them caused bugs" % (numLog, len(bugLogs)))
    runs = collect_cmds(bugLogs, options.out_dir)
    runs.sort(key=lambda tup: tup[0])
    create_test_cmds(runs, options.iqtree_bin, os.path.basename(options.iqtree_bin) + '_test_webserver_cmds.txt')

    
