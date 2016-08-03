#!/bin/bash - 
#===============================================================================
#
#          FILE: submit_jobs.sh
# 
#         USAGE: ./submit_jobs.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Tung Nguyen 
#  ORGANIZATION: 
#       CREATED: 02/06/2015 02:21:40 PM CET
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error

if [ $# -ne 4 ]
then
  echo "USAGE: $0 <number_of_threads> <cmd_file> <aln_dir> <out_dir>" 
  exit 1
fi
numThreads=$1
cmd_file=$2
aln_dir=$3
out_dir=$4

if [ -d $out_dir ]
then
  rm -rf $out_dir
fi
mkdir $out_dir
cp ${aln_dir}/* $out_dir
cp $cmd_file $out_dir
cd $out_dir
submitCMD="submit2sge -N iqtree_system_test -s $numThreads \"../jobmanager.py -f $cmd_file -c $numThreads\""
#echo "../jobmanager.py -f $cmd_file -c $numThreads" | qsub -V -S /bin/bash -cwd -j y -r y -N iqtree_system_test -l zuseX -l cluster -pe threads 16 -q q.norm@zuse02  
$submitCMD
cd ..
