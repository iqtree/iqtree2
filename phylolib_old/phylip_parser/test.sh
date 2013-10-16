#!/bin/sh

for f in testdata/il/*
do
  ./parser $f PHYLIP_INTERLEAVED | grep -i correct
done

for f in testdata/seq/*
do
  ./parser $f PHYLIP_SEQUENTIAL | grep -i correct
done
