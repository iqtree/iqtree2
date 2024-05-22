#!/bin/bash
#
# Creates one static library from several.
#
# Usage: copy all your libs to a single directory and call this script.
#
if [[ $# -ne 2 ]]; then
  echo "Usage: unify-static-libs.sh output-name path-to-libs"
  exit 2
fi
# Inputs
LIBNAME=$1
LIBSDIR=$2
# Tmp dir
OBJDIR=/tmp/unify-static-libs
# Command to evaluate
cmd="ar -crs $LIBNAME"
mkdir -p ${OBJDIR}
# Extract .o
echo "Extracting objects to ${OBJDIR}..."
for i in ${LIBSDIR}/*.a
do
    echo $i
    mkdir -p ${OBJDIR}/$i
    ar --output ${OBJDIR}/$i -x $i
    cmd="${cmd} ${OBJDIR}/$i/*.o"
done
# Link objects into a single lib
# echo "Creating $LIBNAME from objects..."
eval "$cmd"
# Clean up
rm -rf ${OBJDIR}
echo "Done."
