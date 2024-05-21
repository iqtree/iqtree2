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
mkdir -p ${OBJDIR}
# Extract .o
echo "Extracting objects to ${OBJDIR}..."
for i in ${LIBSDIR}/*.a
do
    echo $i
    ar --output $OBJDIR -x $i
done
# Link objects into a single lib
echo "Creating $LIBNAME from objects..."
ar -crs $LIBNAME $OBJDIR/*.o
# Clean up
rm -rf ${OBJDIR}
echo "Done."
