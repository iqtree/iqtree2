#!/usr/bin/env python
'''
Created on 22-Feb-2021, by James Barbetti
'''

import subprocess

def grep_output(stringPat, filePat):
    command = 'grep ' + stringPat + ' test_output/' + filePat
    bytes   = subprocess.run(['bash', '-c', command], stdout=subprocess.PIPE).stdout
    return bytes.decode('UTF8')
    
class Run(object):
    pass

#' -parsimony-spr 10 -lazy-spr -spr-radius' #(no good, large inputs)
#' -parsimony-tbr 10 -lazy-tbr -tbr-radius' #(ditto)
#' -parsimony-nni'                          #(not really commensurate)
for radius in range(1,51):
  suffixes = [ ' -parsimony-spr 10 -spr-radius',
               ' -parsimony-tbr 10 -tbr-radius' ]
  parsimony = ""
  time = ""
  for run in range(2):       
    command = '../build/iqtree2-mpi -s ../build/old/clean_mar_aligned.fasta.gz '
    command += ' -pre crud/last_spr_and_tbr_run -no-ml-dist -nt 12'
    command += ' -fast -redo -m JC -t NJ-R -experimental -v -n 0 -fixbr -seed 1'
    command +=  suffixes[run] + ' ' + str(radius)
    bytes  = subprocess.run(['bash', '-c', command], stdout=subprocess.PIPE).stdout
    data = bytes.decode('UTF8')
    for line in data.split('\n'):
      if "parsimony now" in line:
        if "last iteration" in line:
          parsimony += '\t' + line.split(')')[2].split(' ')[7]
      if "Looking for parsimony" in line:
         time += '\t' + line.split(' ')[6]
  print (str(radius) + parsimony + time)  

