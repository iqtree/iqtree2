#!/usr/bin/env python
'''
Created on Jan. 26, 2015

@author: Tung Nguyen <nltung@gmail.com>
'''
import sys, os, time, multiprocessing, optparse
import subprocess, logging, datetime

def parse_config(config_file):
  singleAln = []
  partitionAln = []
  partOpts = []
  genericOpts = []
  preliminaryTests = []
  universalOptions = []
  with open(config_file) as f:
    #lines = f.readlines()
    lines = [line.strip() for line in f if line.strip()]
  readSingleAln = False
  readPartAln = False
  partOpt = False
  genericOpt = False
  preliminaryTest = False
  universalOpt = False
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
    if partOpt:
        partOpts.append(line)

    if line == 'START_GENERIC_OPTIONS':
      genericOpt = True
      continue
    if line == 'END_GENERIC_OPTIONS':
      genericOpt = False
      continue
    if genericOpt:
      genericOpts.append(line)

    if line == 'START_PRELIMINARY_TESTS':
      preliminaryTest = True
      continue
    if line == 'END_PRELIMINARY_TESTS':
      preliminaryTest = False
      continue
    if preliminaryTest:
      preliminaryTests.append(line)
      
    if line == 'START_UNIVERSAL_OPTIONS':
      universalOpt = True
      continue
    if line == 'END_UNIVERSAL_OPTIONS':
      universalOpt = False
      continue
    if universalOpt:
      universalOptions.append(line)

  return (singleAln, partitionAln, genericOpts, partOpts, preliminaryTests, universalOptions)


if __name__ == '__main__':
  usage = "USAGE: %prog [options]"
  parser = optparse.OptionParser(usage=usage)
  parser.add_option('-b', '--binary', dest="iqtree_bin", help='Path to your IQ-TREE binary')
  parser.add_option('-c', '--config', dest="config_file", help='Path to test configuration file')
  parser.add_option('-o', '--output', dest="outFile", help='Output file for test cases')
  parser.add_option('-f', '--flags',  dest="flags", help='Additional flags for IQ-TREE')
  (options, args) = parser.parse_args()
  if not options.iqtree_bin or not options.config_file:
    parser.print_help()
    exit(0)
  (singleAln, partitionAln, genericOpts, partOpts, preliminaryTests, universalOpts) = parse_config(options.config_file)
  testCmds = []
  
  #Generate suffix parameters (used on every command execution)
  suffix = ''
  for opt in universalOpts:
    suffix += ' ' + opt
  
  # Generate preliminary tests
  for prelim in preliminaryTests:
    cmd = prelim
    testCmds.append(cmd)
  
  # Generate test commands for single model
  for aln in singleAln:
    for opt in genericOpts:
      cmd = '-s ' + aln + ' ' + opt
      if options.flags:
        cmd = cmd + ' ' + options.flags
        testCmds.append(cmd)
        
  # Generate test commands for partition model
  for aln in partitionAln:
    for opt in genericOpts:
      for partOpt in partOpts:
        cmd = '-s ' + aln[0] + ' ' + opt + ' ' + partOpt + ' ' + aln[1]
        if options.flags:
            cmd = cmd + ' ' + options.flags
        testCmds.append(cmd)

  testNr = 1
  jobs = []
  for cmd in testCmds:
    testIDRel = os.path.basename(options.iqtree_bin) + "_TEST_" + str(testNr)
    testCMD = testIDRel + " " + os.path.abspath(options.iqtree_bin) + " -pre " + testIDRel + " " + cmd + suffix
    testNr = testNr + 1
    jobs.append(testCMD)
#  print "\n".join(jobs)
  if options.outFile == "STDOUT":
    outfile = sys.stdout
  else:
    outfile = open(options.outFile, "w")
  for job in jobs:
    outfile.write(job + "\n")
    # Or, Python 2 only: print >> outfile, job
    # Or, Python 3 only: print(job, file=outfile)
  outfile.close()
