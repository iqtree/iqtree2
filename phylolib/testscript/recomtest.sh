#! /bin/sh
# This is the working directory you want to run the test from
WORKINGDIR=`pwd`
DATADIR=`pwd`/../testdata
# Set the number of threads you want to use
NUMPTHREADS=4
# NO NEED TO EDIT BEYOND THIS POINT
SSE3_GCC="testscript_recom.SSE3.gcc"
SSE3_PTHREADS_GCC="testscript.SSE3.PTHREADS.gcc"

# KEEP A LOG OF STDOUT OF RAxML
LOGFILE=RAxML_log.txt
ERRLOGFILE=RAxML_ERR_log.txt

run_SSE3_GCC()
{
  sh ${WORKINGDIR}/${SSE3_GCC} ${OPTIONS}
  RETVAL=$?
  if [ $RETVAL -ne 0 ] ; then
    echo " Error: ${SSE3_GCC} failed. Check the messages above ...";
    echo "\n"
    echo " Supertestscript did not finish successfully."
    echo "\n"
    exit 1
  fi
  echo "\n"
}

run_SSE3_PTHREADS_GCC()
{
  sh ${WORKINGDIR}/${SSE3_PTHREADS_GCC} ${OPTIONS} -T ${NUMPTHREADS}
  RETVAL=$?
  if [ $RETVAL -ne 0 ] ; then
    echo " Error: ${SSE3_PTHREADS_GCC} failed. Check the messages above ...";
    echo "\n"
    echo " Supertestscript did not finish successfully."
    echo "\n"
    exit 1
  fi
  echo "\n"
}


if [ ! -f ${WORKINGDIR}/${SSE3_GCC} ] ; then
  echo " Error: file ${SSE3_GCC} is missing from ${WORKINGDIR}, exiting ...";
  exit 1
fi
if [ ! -f ${WORKINGDIR}/${SSE3_PTHREADS_GCC} ] ; then
  echo " Error: file ${SSE3_PTHREADS_GCC} is missing from ${WORKINGDIR}, exiting ...";
  exit 1
fi

if [ $# -eq 0 ] ; then
  echo "\n"
  echo " usage: sh recomtest.sh [ [ 1-3 ] | [raxml options] ]"
  echo " options: "
  echo "        [1]  small size test"
  echo "        [3]  large size test"
  echo "        [raxml options]  specific test"
  echo "\n"
  exit 0
fi

if [ $# -eq 1 ] ; then
  if [ $1 -ne 0 ] && [ $1 -ne 1 ] && [ $1 -ne 2 ] && [ $1 -ne 3 ] ; then
    echo "\n"
    echo " usage: sh supertestscript.sh [ [ 0-3 ] | [raxml options] ]"
    echo " options: "
    echo "        [1]  small size test"
    echo "        [3]  large size test"
    echo "        [raxml options]  specific test"
    echo "\n"
    exit 1
  fi
  if [ $1 -eq 1 ] ; then
    DESC="small"
    GAPPY_DATASET=50
  fi
  if [ $1 -eq 3 ] ; then
    DESC="large"
    GAPPY_DATASET=1288
  fi

  #SIMPLE=""                # this will be ignored in the loop
  #BL_PARTITION="-M"
  #SIMPLE_GAPPY="-S"
  #BL_PARTITION_GAPPY="-M -S" # !! this will be taken as 2
  #SIMPLE_RF_CONV="-D"

  echo "Starting recom custom test `date`, Errors" > $ERRLOGFILE 
  echo "Starting recom custom test `date`, Log" > $LOGFILE 

  for VERSION in SSE3_GCC #SSE3_PTHREADS_GCC 
  do 
    echo $VERSION
    for MODEL in PSR #GAMMA 
    do
      echo $MODEL
      for ALPHABET in dna #aa 
      do
        TREE="${DATADIR}/${DESC}.startingTree.${ALPHABET}.tree"
        for PARTITION_LABEL in singlegene.binary #binary 
        do
          echo "   "
          ALIGNMENT="${DATADIR}/${DESC}.${ALPHABET}.${PARTITION_LABEL}"
          BASE_OPTIONS="-m ${MODEL} -s ${ALIGNMENT} -t ${TREE}"
          echo ${BASE_OPTIONS} 

          #Normal
          OPTIONS="${BASE_OPTIONS}"
          echo "$VERSION with $OPTIONS" | tee -a $ERRLOGFILE $LOGFILE
          (run_${VERSION} 2>> $ERRLOGFILE) >> $LOGFILE

          # RF convergence
          OPTIONS="${BASE_OPTIONS} -D"
          echo "$VERSION with $OPTIONS" | tee -a $ERRLOGFILE $LOGFILE
          (run_${VERSION} 2>> $ERRLOGFILE) >> $LOGFILE

        done
        # Now partitioned with branch lengths
        ALIGNMENT="${DATADIR}/${DESC}.${ALPHABET}.binary"
        OPTIONS="-m ${MODEL} -s ${ALIGNMENT} -t ${TREE} -M"
        echo "$VERSION with $OPTIONS" | tee -a $ERRLOGFILE $LOGFILE
        (run_${VERSION} 2>> $ERRLOGFILE) >> $LOGFILE

      done
      # Now check the special datasets (gappy stuff)
      ALIGNMENT=${DATADIR}/gappy/${GAPPY_DATASET}_gappy.binary # this is a single gene
      TREE=${DATADIR}/gappy/RAxML_parsimonyTree.${GAPPY_DATASET}
      OPTIONS="-m ${MODEL} -s ${ALIGNMENT} -t ${TREE} -S"
      echo "$VERSION with $OPTIONS" | tee -a $ERRLOGFILE $LOGFILE
      (run_${VERSION} 2>> $ERRLOGFILE) >> $LOGFILE

    done
  done
else
  OPTIONS="$@"
  for VERSION in SSE3_GCC SSE3_PTHREADS_GCC 
  do 
    run_${VERSION}
  done
fi

echo "\n"
echo " recom test script finished successfully. Check file ${LOGFILE} to ensure that everything is ok!"
echo "\n"
