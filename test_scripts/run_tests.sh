#!/bin/bash -
#===============================================================================
#
#          FILE: run_tests.sh
#
#         USAGE: ./run_tests.sh
#
#   DESCRIPTION:
#
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Tung Nguyen (nltung@gmail.com),
#  ORGANIZATION:
#       CREATED: 2016-08-12 16:43:54 CEST
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error

if [ "$#" -lt 1 ]
then
    echo "Please enter the name of the local branch you want to compile"
    echo "USAGE: $0 <branch_name> [<iqtree_flags_in_quotes>]" >&2
    exit 1
fi

branchName=$1
flags=$2

#Compile the specified branch
source compile.sh ${branchName} "$flags"

#Generate test cases
echo "Generate test cases for the sequential version"
./gen_test_standard.py -b ${buildDir}/iqtree -c test_configs.txt
