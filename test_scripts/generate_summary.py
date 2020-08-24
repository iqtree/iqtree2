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
    
def last_file_matching(folder, file_pattern):
    command = 'ls -t1 ' + folder + '/' + file_pattern + ' | head -1'
    bytes   = subprocess.run(['bash', '-c', command], stdout=subprocess.PIPE).stdout
    return bytes.decode('UTF8').split('\n')[0]
    
def all_files_matching(folder, file_pattern):
    command = 'ls -1 ' + folder + '/' + file_pattern
    bytes   = subprocess.run(['bash', '-c', command], stdout=subprocess.PIPE).stdout
    list    = bytes.decode('UTF8').split('\n')
    return list[0:len(list)-2]
    
def content_of_file(filePath):
    bytes   = subprocess.run(['cat', filePath], stdout=subprocess.PIPE).stdout
    return bytes.decode('UTF8')

class Run(object):
    pass

if __name__ == '__main__':
    runs = {}
    param_lines = grep_output('^Command', '*_*.log')
    for line in param_lines.split('\n'):
        prefix = ''
        words = line.split(' ')
        run = Run()
        run.startingTree = ''
        run.numberOfThreads = 0
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
            # The following will be updated if they're
            # found in the output of the test
            run.errorCode = '0'
            run.sequences = 0
            run.sites = 0
            run.patterns = 0
            run.numberOfThreads = 1

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

    command_lines = grep_output('^Command\:', '*_*.log')
    for line in command_lines.split('\n'):
        prefix = get_prefix(line)
        if 0<len(line) and prefix in runs:
            words    = line.split(' ')
            #... e.g. Command: [path-to-iqtree2] -pre iqtree2-mpi_TEST_1 -s example.phy
            #         -redo -m TEST -sp example.nex
            #... We want the fifth and subsequent words
            params = ''
            for i, word in enumerate(words):
                if (4<=i):
                    if (5<=i):
                        params += ' '
                    params += word
            runs[prefix].params = params

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
        prefix   = get_prefix(line)
        if 0<len(line) and prefix in runs:
            prevWord = ''
            words = line.split(' ')
            read_tree_likelihood = 0
            for i, word in enumerate(words):
                if 0<i:
                    prevWord = words[i-1]
                if i+1<len(words):
                    nextWord = words[i+1]
                    if word=='tree:':
                        if prevWord=='consensus':
                            runs[prefix].consensus_tree_likelihood = float(nextWord)
                        else:
                            runs[prefix].tree_likelihood = float(nextWord)
                            read_tree_likelihood = 1
                    if word=='(s.e.' and 0<read_tree_likelihood:
                        if nextWord[len(nextWord)-1] == ')':
                            nextWord = nextWord[0: (len(nextWord)-1)]
                        runs[prefix].tree_standard_error = float(nextWord)
                        

    latest_test_logfile = last_file_matching('test_output', 'test_commands.*.log')
    log_lines = content_of_file(latest_test_logfile)
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
                            
    scheme_file_paths = all_files_matching('test_output', '*best_scheme.nex')
    for scheme_file_path in scheme_file_paths:
        reading_scheme = 0
        prefix = re.search('\/(.*?)\.', scheme_file_path).group(1)
        if prefix in runs:
            runs[prefix].schemes = []
            for line_no, line in enumerate(content_of_file(scheme_file_path).split('\n')):
                if "charpartition" in line:
                    reading_scheme = 1
                    continue
                if "end" in line:
                    reading_schem = 0
                elif (0<len(line) and 0<reading_scheme):
                    runs[prefix].schemes.append ( line.strip() )

    model_file_paths = all_files_matching('test_output', '*best_model.nex')
    for model_file_path in model_file_paths:
        reading_model = 0
        prefix = re.search('\/(.*?)\.', model_file_path).group(1)
        if prefix in runs:
            runs[prefix].models = []
            for line_no, line in enumerate(content_of_file(model_file_path).split('\n')):
                if "charpartition" in line:
                    reading_model = 1
                    continue
                if "end" in line:
                    reading_model = 0
                elif (0<len(line) and 0<reading_model):
                    runs[prefix].models.append ( line.strip() )
            
    summary_path = latest_test_logfile.replace('test_commands', 'test_summary')
    summary_path = summary_path.replace('.log', '.csv')
    summary_file = open(summary_path, 'w')
    print(summary_path)

    count_failed    = 0
    count_succeeded = 0
    line_number     = 0
    for run in runs.values():
        line_number += 1
        if not hasattr(run,'params'):
            run.params = ''
        if not hasattr(run,'tree_likelihood'):
            run.tree_likelihood = ''
        if not hasattr(run,'consensus_tree_likelihood'):
            run.consensus_tree_likelihood = ''
        if not hasattr(run, 'elapsed_time'):
            run.elapsed_time = ''
        if not hasattr(run, 'models'):
            run.models = []
        if not hasattr(run, 'schemes'):
            run.schemes = []
        if not hasattr(run, 'partition'):
            run.partition = ''
        if run.errorCode != '0':
            count_failed += 1
        else:
            count_succeeded += 1
        if line_number == 1:
            summary_file.write(
                'Params,Prafix,Time,Error,Alignment,Partition' +
                ',Sequences,Sites,Patterns,Tree,Threads' +
                ',Likelihood,Consensus_Likelihood' +
                ',Schemes,Models\n' )
        summary_file.write( ''
            + '"'    + run.params + '"'
            + ',\t'  + run.prefix
            + ',\t'  + str(run.elapsed_time)
            + ',\t'  + run.errorCode
            + ',\t'  + run.alignment
            + ',\t'  + run.partition
            + ',\t'  + str(run.sequences)
            + ',\t'  + str(run.sites)
            + ',\t'  + str(run.patterns)
            + ',\t'  + run.startingTree
            + ',\t'  + str(run.numberOfThreads)
            + ',\t'  + str(run.tree_likelihood)
            + ',\t'  + str(run.consensus_tree_likelihood)
            + ',\t"' + str(run.schemes) + '"'
            + ',\t"' + str(run.models) + '"\n')

    print("Wrote details of " + str(count_failed) + " failed, and "
        + str(count_succeeded) + " successful tests to " + summary_path)
