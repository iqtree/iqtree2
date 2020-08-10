#!/usr/bin/env python
'''
Created on 10-Aug-2020, by James Barbetti
'''

import subprocess
import re
from datetime import datetime

def grep_output(stringPat, filePat):
    command = 'grep ' + stringPat + ' test_output/' + filePat
    bytes   = subprocess.run(['bash', '-c', command], stdout=subprocess.PIPE).stdout
    return bytes.decode('UTF8')
    
def get_prefix(line):
    if len(line)==0:
        return ''
    return re.search('/(.*?)\\.', line).group(1)

class Run(object):
    pass

if __name__ == '__main__':
    runs = {}
    param_lines = grep_output('^Command', '*_*.log')
    for line in param_lines.split('\n'):
        prefix = ''
        words = line.split(' ')
        run = Run()
        for i, word in enumerate(words):
            if i+1<len(words):
                nextWord = words[i+1]
                if word=='-pre':
                    run.prefix = nextWord
                if word=='-s':
                    run.alignment = nextWord
                if word=='-sp':
                    run.partition = nextWord
                if word=='-m':
                    run.model = nextWord
                if word=='-nt':
                    run.numberOfThreads = int(nextWord)
                if word=='-t' or word=='-starttree':
                    run.startingTree = nextWord
        if hasattr(run,'prefix'):
            runs[run.prefix] = run
            run.errorCode = ''

    align_lines = grep_output('"Alignment has"', '*_*.log')
    for line in align_lines.split('\n'):
        prefix = get_prefix(line)
        if 0<len(line) and prefix in runs:
            words    = line.split(' ')
            lastWord = ''
            for i, word in enumerate(words):
                if 0<i:
                    if word=='sequences':
                        runs[prefix].sequences = int(lastWord)
                    if word=='columns,':
                        runs[prefix].sites = int(lastWord)
                    if word=='distinct':
                        runs[prefix].patterns = int(lastWord)
                lastWord = word

    parsimony_lines = grep_output('"constant sites"', '*_*.log')
    for line in parsimony_lines.split('\n'):
        prefix = get_prefix(line)
        if 0<len(line) and prefix in runs:
            words    = line.split(' ')
            lastWord = ''
            for i, word in enumerate(words):
                if 0<i:
                    if word=='constant':
                        runs[prefix].constant_sites = int(lastWord)
                lastWord = word
                
    likelihood_lines = grep_output('"Log-likelihood"', '*_*.iqtree')
    for line in likelihood_lines.split('\n'):
        prefix = get_prefix(line)
        if 0<len(line) and prefix in runs:
            words = line.split(' ')
            for i, word in enumerate(words):
                if i+1<len(words):
                    nextWord = words[i+1]
                    if word=='tree:':
                        runs[prefix].tree_likelihood = float(nextWord)

    # Problem: we want only the most recent test_commands log file.  But
    # I haven't figured out a good way to get that yet.
    log_lines = grep_output("INFO", "test_commands.*.log")
    for line in log_lines.split('\n'):
        #Times in the format, 2020-08-07 11:50:50,520
        #There doesn't seem to be a "milliseconds" format specifier,
        #so I fake microseconds (format specifier %f) by tacking on three zeroes (!).
        if " - " in line:
            time = datetime.strptime( line.split(' - ')[0] + '000', '%Y-%m-%d %H:%M:%S,%f')
            if "job " in line or "Job " in line:
                prefix = re.search('[Jj]ob (.*?)[: ]', line).group(1)
                if prefix in runs:
                    if "Executing" in line:
                        runs[prefix].start_time = time
                    if "finished" in line:
                        runs[prefix].finish_time  = time
                        runs[prefix].elapsed_time = (time - runs[prefix].start_time).total_seconds()
                        if "ERROR CODE" in line:
                            runs[prefix].errorCode = line.split("ERROR CODE")[1]
                        
    for run in runs.values():
        if not hasattr(run,'tree_likelihood'):
            run.tree_likelihood = ''
        if hasattr(run,'elapsed_time'):
            print(run.prefix + ',\t' + str(run.elapsed_time) + ',\t' + run.errorCode + ',\t' + str(run.tree_likelihood))
        #print(run.prefix + ' ' + run.errorCode)
        
