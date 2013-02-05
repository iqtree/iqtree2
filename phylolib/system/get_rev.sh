#!/bin/bash
which hg &> /dev/null
if [ $? -eq 0 ] ; then 
	hg heads | awk '{if(NR==1) print $2;}'
else
	echo "no_rev"
fi
