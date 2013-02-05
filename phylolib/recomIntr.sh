#!/bin/sh

rm *TEST*
rm  *.o
rm raxmlLight*
make -f Makefile.SSE3.gcc
make -f Makefile.SSE3.PTHREADS.gcc

NUM_T=2
DIR=testdata
DATASET=${DIR}/20.dna.singlegene.binary
TREE=${DIR}/20.tree

for MODEL in PSR GAMMA
do
  for ALGORI in fd i5 fb 
  do
    PARAMS="-${ALGORI} -s ${DATASET} -t ${TREE} -m ${MODEL}"
    NAME=TEST_${MODEL}_${ALGORI}
    echo ${PARAMS}
    ./raxmlLight ${PARAMS} -n ${NAME} > /dev/null
    for REC in 0.4 0.9 #0.5 0.7 
    do
      echo $REC
      ./raxmlLight-PTHREADS -T ${NUM_T} ${PARAMS} -n ${NAME}_${REC}_T${NUM_T} -L ${REC} > /dev/null
      ./raxmlLight ${PARAMS} -n ${NAME}_${REC} -L ${REC} > /dev/null
    done
    if [ $ALGORI = "fb" ]; then
      grep "tree" RAxML_info.${NAME}* | grep -v testdata
      grep "lh:" RAxML_info.${NAME}*
    else
      grep "Likelihood" RAxML_info.${NAME}*
    fi
  done
done
