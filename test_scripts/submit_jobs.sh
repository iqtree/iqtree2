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

if [ $# -ne 5 ]
then
  echo "USAGE: $0 <number_of_threads> <cmd_file> <aln_dir> <out_dir> <binary_dir>"
  exit 1
fi
numThreads=$1
cmd_file=$2
aln_dir=$3
out_dir=$4
binary_dir=$5

if [ -d $out_dir ]
then
  rm -rf $out_dir
else
  mkdir $out_dir
fi
cp ${aln_dir}/* $out_dir
cp $cmd_file $out_dir
cp -r ${binary_dir} ${out_dir}/
cd $out_dir
submitCMD="submit2sge -N iqtree_system_test -q cluster -r zuseX -s $numThreads \"../jobmanager.py -f $cmd_file -c $numThreads\""
#echo "../jobmanager.py -f $cmd_file -c $numThreads" | qsub -V -S /bin/bash -cwd -j y -r y -N iqtree_system_test -l zuseX -l cluster -pe threads 16 -q q.norm@zuse02  
$submitCMD
cd ..
