#!/usr/bin/env python
'''
Created on Jan. 26, 2015

@author: Tung Nguyen <nltung@gmail.com>
'''
import sys, os, time, multiprocessing, optparse
import subprocess, logging, datetime

def parse_config(config_file):
  singleAln, partitionAln, partOpts, genericOpts = [], [], [], []
  with open(config_file) as f:
    #lines = f.readlines()
    lines = [line.strip() for line in f if line.strip()]
  readSingleAln = False
  readPartAln = False
  partOpt = False
  genericOpt = False
  for line in lines:
    #print line
    if line == 'START_SINGLE_ALN':
      readSingleAln = True
      continue
    if line == 'END_SINGLE_ALN':
      readSingleAln = False
      continue
    if readSingleAln:
      singleAln.append(line)
    if line == 'START_PARTITION_ALN':
      readPartAln = True
      continue
    if line == 'END_PARTITION_ALN':
      readPartAln = False
      continue
    if readPartAln:
      partitionAln.append(line.split())
    if line == 'START_PARTITION_OPTIONS':
      partOpt = True
      continue
    if line == 'END_PARTITION_OPTIONS':
      partOpt = False
      continue
    if line == 'START_GENERIC_OPTIONS':
      genericOpt = True
      continue
    if line == 'END_GENERIC_OPTIONS':
      genericOpt = False
      continue
    if partOpt:
      partOpts.append(line)
    if genericOpt:
      genericOpts.append(line)
  return (singleAln, partitionAln, genericOpts, partOpts)


if __name__ == '__main__':
  usage = "USAGE: %prog [options]"
  parser = optparse.OptionParser(usage=usage)
  parser.add_option('-b','--binary', dest="iqtree_bin", help='Path to your IQ-TREE binary')
  parser.add_option('-c','--config', dest="config_file", help='Path to test configuration file')
  parser.add_option('-o', '--output', dest="outFile", help='Output file for test cases')
  parser.add_option('-f', '--flags', dest="flags", help='Additional flags for IQ-TREE')
  (options, args) = parser.parse_args()
  if not options.iqtree_bin or not options.config_file:
    parser.print_help()
    exit(0)
  (singleAln, partitionAln, genericOpts, partOpts) = parse_config(options.config_file)
  testCmds = []
  # Generate test commands for single model
  for aln in singleAln:
    for opt in genericOpts:
      cmd = '-s ' + aln + ' -redo ' + opt
      if options.flags:
        cmd = cmd + ' ' + options.flags
        testCmds.append(cmd)
  # Generate test commands for partition model
  for aln in partitionAln:
    for opt in genericOpts:
      for partOpt in partOpts:
        cmd = '-s ' + aln[0] + ' -redo ' + opt + ' ' + partOpt + ' ' + aln[1]
        if options.flags:
            cmd = cmd + ' ' + options.flags
        testCmds.append(cmd)

  testNr = 1
  jobs = []
  for cmd in testCmds:
    testIDRel = os.path.basename(options.iqtree_bin) + "_TEST_" + str(testNr)
    testCMD = testIDRel + " " + os.path.abspath(options.iqtree_bin) + " -pre " + testIDRel + " " + cmd
    testNr = testNr + 1
    jobs.append(testCMD)
#  print "\n".join(jobs)
  outfile = open(options.outFile, "wb")
  for job in jobs:
    print >> outfile, job
  outfile.close()
