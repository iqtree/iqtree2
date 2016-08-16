#!/bin/bash -
#===============================================================================
#
#          FILE: compile.sh
#
#         USAGE: ./compile.sh
#
#   DESCRIPTION: This script checkouts and compile the specified branch of IQ-TREE
#
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Tung Nguyen (nltung@gmail.com)
#  ORGANIZATION:
#       CREATED: 2015-01-26 13:02:57 CET
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error

require_clean_work_tree () {
    # Update the index
    git update-index -q --ignore-submodules --refresh
    err=0

    # Disallow unstaged changes in the working tree
    if ! git diff-files --quiet --ignore-submodules --
    then
        echo >&2 "cannot $0: you have unstaged changes."
        git diff-files --name-status -r --ignore-submodules -- >&2
        err=1
    fi

    # Disallow uncommitted changes in the index
    if ! git diff-index --cached --quiet HEAD --ignore-submodules --
    then
        echo >&2 "cannot $0: your index contains uncommitted changes."
        git diff-index --cached --name-status -r --ignore-submodules HEAD -- >&2
        err=1
    fi

    if [ $err = 1 ]
    then
        echo >&2 "Please commit or stash them."
        exit 1
    fi
}

#Check whether the git work tree is clean
#require_clean_work_tree

if [ "$#" -lt 1 ]
then
    echo "Please enter the name of the local branch you want to compile"
    echo "USAGE: $0 <branch_name> [<iqtree_flags>]" >&2
    exit 1
fi


#Determine hash code of current branch
#branch=`git status | grep "On branch" | awk '{print $3}'`
branch=$1
flags=$2
flagOMP="${flags} omp" # flags used to compile OpenMP version of IQ-TREE
echo "COMPILING BRANCH ${branch} USING FLAGS ${flags}"
#Take the first 6 characters of the current head commit
commit_cur=`git log | head -n1 | awk '{print $2}' | cut -c 1-6`

#Assign names to build and binary directories
flagSuffix=`echo ${flags} | sed 's/ /-/g'`
buildDir="build-${branch}-${flagSuffix}"
buildDirOMP="build-${branch}-${flagSuffix}-omp"
binaryName="iqtree-${branch}"
binaryNameOMP="${binaryName}-omp"
binDir="iqtree-${branch}-bin"

#Create the build directory
if [[ ! -e $buildDir ]]
then
  mkdir $buildDir
fi
if [[ ! -e $buildDirOMP ]]
then
  mkdir $buildDirOMP
fi

#Create binary directory
if [[ ! -e $binDir ]]
then
    mkdir $binDir
fi

#Fetch changes from server
git fetch
curBranch=`git status | grep 'On branch' | awk '{print $3}'`
if [[ ${curBranch} != ${branch} ]]
then
    echo "Switch to branch ${branch} and pull code from the server ... "
    git stash
    echo "Current changes stashed."
    git checkout $branch
    git pull
    #git submodule update
fi

#Build the selected

echo -e "\nGENERATING MAKEFILE FOR SEQUENTIAL VERSION OF IQ-TREE FOR BRANCH ${branch}\n"
cmake -B${buildDir} -H.. -DIQTREE_FLAGS="${flags}"
echo -e "\nBUILDING SEQUENTIAL VERSION OF IQ-TREE FOR BRANCH ${branch}\n"
make -C ${buildDir} -j4

echo -e "\nGENERATING MAKEFILE FOR OPENMP VERSION OF IQ-TREE FOR BRANCH ${branch}\n"
echo ${flagOMP}
echo ${buildDirOMP}
cmake -B${buildDirOMP} -H.. -DIQTREE_FLAGS="${flagOMP}"
echo -e "\nBUILDING OPENMP VERSION OF IQ-TREE FOR BRANCH ${branch}\n"
make -C ${buildDirOMP} -j4

#cp ${buildDir}/iqtree- ${binDir}/${binaryName}

#Clean up
#rm -rf $buildDir

#echo -e "Binaries of IQ-TREE for branch ${branch} are stored in $binDir"
#rm -rf $release_build
