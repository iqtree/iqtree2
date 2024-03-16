# Legacy tests to be ported to new testing infrastructure

## Comparison tests of log likelihood between iqtree1 and iqtree2 

### Single alignments

-  -s data/example.phy -redo -m TEST -cmin 2 -pre output
-  -s data/example.phy -redo -m TEST -mtree -pre output
-  -s data/example.phy -redo -m TEST -ft -pre output
-  -s data/example.phy -redo -m TEST -fs -pre output
-  -s data/example.phy -redo -m TEST -fmax -pre output
-  -s data/example.phy -redo -m TEST -allnni -pre output
-  -s data/example.phy -redo -m TEST -djc -pre output
-  -s data/example.phy -redo -m TEST -ntop 10 -pre output
-  -s data/example.phy -redo -m TEST -nbest 5 -pre output
-  -s data/example.phy -redo -m TEST -nstop 100 -pre output
-  -s data/example.phy -redo -m TEST -pers 0.5 -pre output
-  -s data/example.phy -redo -m TEST -sprrad 6 -pre output
-  -s data/example.phy -redo -m TEST -bb 1000 -pre output
-  -s data/example.phy -redo -m TEST -bcor 0.5 -pre output
-  -s data/example.phy -redo -m TEST -beps 0.5 -pre output
-  -s data/example.phy -redo -m TEST -bnni -pre output

-  -s data/d59_8.phy -redo -m TEST -cmin 2 -pre output
-  -s data/d59_8.phy -redo -m TEST -mtree -pre output
-  -s data/d59_8.phy -redo -m TEST -ft -pre output
-  -s data/d59_8.phy -redo -m TEST -fs -pre output
-  -s data/d59_8.phy -redo -m TEST -fmax -pre output
-  -s data/d59_8.phy -redo -m TEST -allnni -pre output
-  -s data/d59_8.phy -redo -m TEST -djc -pre output
-  -s data/d59_8.phy -redo -m TEST -ntop 10 -pre output
-  -s data/d59_8.phy -redo -m TEST -nbest 5 -pre output
-  -s data/d59_8.phy -redo -m TEST -nstop 100 -pre output
-  -s data/d59_8.phy -redo -m TEST -pers 0.5 -pre output
-  -s data/d59_8.phy -redo -m TEST -sprrad 6 -pre output
-  -s data/d59_8.phy -redo -m TEST -bb 1000 -pre output
-  -s data/d59_8.phy -redo -m TEST -bcor 0.5 -pre output
-  -s data/d59_8.phy -redo -m TEST -beps 0.5 -pre output
-  -s data/d59_8.phy -redo -m TEST -bnni -pre output

### Partition alignments

#### Edge-unlinked partition model
-  -s data/example.phy -redo -sp data/example.nex -m TEST -cmin 2 -pre output
-  -s data/example.phy -redo -sp data/example.nex -m TEST -mtree -pre output
-  -s data/example.phy -redo -sp data/example.nex -m TEST -ft -pre output
-  -s data/example.phy -redo -sp data/example.nex -m TEST -fs -pre output
-  -s data/example.phy -redo -sp data/example.nex -m TEST -fmax -pre output
-  -s data/example.phy -redo -sp data/example.nex -m TEST -allnni -pre output
-  -s data/example.phy -redo -sp data/example.nex -m TEST -djc -pre output
-  -s data/example.phy -redo -sp data/example.nex -m TEST -ntop 10 -pre output
-  -s data/example.phy -redo -sp data/example.nex -m TEST -nbest 5 -pre output
-  -s data/example.phy -redo -sp data/example.nex -m TEST -nstop 100 -pre output
-  -s data/example.phy -redo -sp data/example.nex -m TEST -pers 0.5 -pre output
-  -s data/example.phy -redo -sp data/example.nex -m TEST -sprrad 6 -pre output
-  -s data/example.phy -redo -sp data/example.nex -m TEST -bb 1000 -pre output
-  -s data/example.phy -redo -sp data/example.nex -m TEST -bcor 0.5 -pre output
-  -s data/example.phy -redo -sp data/example.nex -m TEST -beps 0.5 -pre output
-  -s data/example.phy -redo -sp data/example.nex -m TEST -bnni -pre output

-  -s data/d59_8.phy -redo -sp data/d58_9.nex -m TEST -cmin 2 -pre output
-  -s data/d59_8.phy -redo -sp data/d58_9.nex -m TEST -mtree -pre output
-  -s data/d59_8.phy -redo -sp data/d58_9.nex -m TEST -ft -pre output
-  -s data/d59_8.phy -redo -sp data/d58_9.nex -m TEST -fs -pre output
-  -s data/d59_8.phy -redo -sp data/d58_9.nex -m TEST -fmax -pre output
-  -s data/d59_8.phy -redo -sp data/d58_9.nex -m TEST -allnni -pre output
-  -s data/d59_8.phy -redo -sp data/d58_9.nex -m TEST -djc -pre output
-  -s data/d59_8.phy -redo -sp data/d58_9.nex -m TEST -ntop 10 -pre output
-  -s data/d59_8.phy -redo -sp data/d58_9.nex -m TEST -nbest 5 -pre output
-  -s data/d59_8.phy -redo -sp data/d58_9.nex -m TEST -nstop 100 -pre output
-  -s data/d59_8.phy -redo -sp data/d58_9.nex -m TEST -pers 0.5 -pre output
-  -s data/d59_8.phy -redo -sp data/d58_9.nex -m TEST -sprrad 6 -pre output
-  -s data/d59_8.phy -redo -sp data/d58_9.nex -m TEST -bb 1000 -pre output
-  -s data/d59_8.phy -redo -sp data/d58_9.nex -m TEST -bcor 0.5 -pre output
-  -s data/d59_8.phy -redo -sp data/d58_9.nex -m TEST -beps 0.5 -pre output
-  -s data/d59_8.phy -redo -sp data/d58_9.nex -m TEST -bnni -pre output

#### partition-specific rates

-  -s data/example.phy -redo -spp data/example.nex -m TEST -cmin 2 -pre output
-  -s data/example.phy -redo -spp data/example.nex -m TEST -mtree -pre output
-  -s data/example.phy -redo -spp data/example.nex -m TEST -ft -pre output
-  -s data/example.phy -redo -spp data/example.nex -m TEST -fs -pre output
-  -s data/example.phy -redo -spp data/example.nex -m TEST -fmax -pre output
-  -s data/example.phy -redo -spp data/example.nex -m TEST -allnni -pre output
-  -s data/example.phy -redo -spp data/example.nex -m TEST -djc -pre output
-  -s data/example.phy -redo -spp data/example.nex -m TEST -ntop 10 -pre output
-  -s data/example.phy -redo -spp data/example.nex -m TEST -nbest 5 -pre output
-  -s data/example.phy -redo -spp data/example.nex -m TEST -nstop 100 -pre output
-  -s data/example.phy -redo -spp data/example.nex -m TEST -pers 0.5 -pre output
-  -s data/example.phy -redo -spp data/example.nex -m TEST -sprrad 6 -pre output
-  -s data/example.phy -redo -spp data/example.nex -m TEST -bb 1000 -pre output
-  -s data/example.phy -redo -spp data/example.nex -m TEST -bcor 0.5 -pre output
-  -s data/example.phy -redo -spp data/example.nex -m TEST -beps 0.5 -pre output
-  -s data/example.phy -redo -spp data/example.nex -m TEST -bnni -pre output

-  -s data/d59_8.phy -redo -spp data/d58_9.nex -m TEST -cmin 2 -pre output
-  -s data/d59_8.phy -redo -spp data/d58_9.nex -m TEST -mtree -pre output
-  -s data/d59_8.phy -redo -spp data/d58_9.nex -m TEST -ft -pre output
-  -s data/d59_8.phy -redo -spp data/d58_9.nex -m TEST -fs -pre output
-  -s data/d59_8.phy -redo -spp data/d58_9.nex -m TEST -fmax -pre output
-  -s data/d59_8.phy -redo -spp data/d58_9.nex -m TEST -allnni -pre output
-  -s data/d59_8.phy -redo -spp data/d58_9.nex -m TEST -djc -pre output
-  -s data/d59_8.phy -redo -spp data/d58_9.nex -m TEST -ntop 10 -pre output
-  -s data/d59_8.phy -redo -spp data/d58_9.nex -m TEST -nbest 5 -pre output
-  -s data/d59_8.phy -redo -spp data/d58_9.nex -m TEST -nstop 100 -pre output
-  -s data/d59_8.phy -redo -spp data/d58_9.nex -m TEST -pers 0.5 -pre output
-  -s data/d59_8.phy -redo -spp data/d58_9.nex -m TEST -sprrad 6 -pre output
-  -s data/d59_8.phy -redo -spp data/d58_9.nex -m TEST -bb 1000 -pre output
-  -s data/d59_8.phy -redo -spp data/d58_9.nex -m TEST -bcor 0.5 -pre output
-  -s data/d59_8.phy -redo -spp data/d58_9.nex -m TEST -beps 0.5 -pre output
-  -s data/d59_8.phy -redo -spp data/d58_9.nex -m TEST -bnni -pre output

#### Edge-linked partition model

-  -s data/example.phy -redo -q data/example.nex -m TEST -cmin 2 -pre output
-  -s data/example.phy -redo -q data/example.nex -m TEST -mtree -pre output
-  -s data/example.phy -redo -q data/example.nex -m TEST -ft -pre output
-  -s data/example.phy -redo -q data/example.nex -m TEST -fs -pre output
-  -s data/example.phy -redo -q data/example.nex -m TEST -fmax -pre output
-  -s data/example.phy -redo -q data/example.nex -m TEST -allnni -pre output
-  -s data/example.phy -redo -q data/example.nex -m TEST -djc -pre output
-  -s data/example.phy -redo -q data/example.nex -m TEST -ntop 10 -pre output
-  -s data/example.phy -redo -q data/example.nex -m TEST -nbest 5 -pre output
-  -s data/example.phy -redo -q data/example.nex -m TEST -nstop 100 -pre output
-  -s data/example.phy -redo -q data/example.nex -m TEST -pers 0.5 -pre output
-  -s data/example.phy -redo -q data/example.nex -m TEST -sprrad 6 -pre output
-  -s data/example.phy -redo -q data/example.nex -m TEST -bb 1000 -pre output
-  -s data/example.phy -redo -q data/example.nex -m TEST -bcor 0.5 -pre output
-  -s data/example.phy -redo -q data/example.nex -m TEST -beps 0.5 -pre output
-  -s data/example.phy -redo -q data/example.nex -m TEST -bnni -pre output

-  -s data/d59_8.phy -redo -q data/d58_9.nex -m TEST -cmin 2 -pre output
-  -s data/d59_8.phy -redo -q data/d58_9.nex -m TEST -mtree -pre output
-  -s data/d59_8.phy -redo -q data/d58_9.nex -m TEST -ft -pre output
-  -s data/d59_8.phy -redo -q data/d58_9.nex -m TEST -fs -pre output
-  -s data/d59_8.phy -redo -q data/d58_9.nex -m TEST -fmax -pre output
-  -s data/d59_8.phy -redo -q data/d58_9.nex -m TEST -allnni -pre output
-  -s data/d59_8.phy -redo -q data/d58_9.nex -m TEST -djc -pre output
-  -s data/d59_8.phy -redo -q data/d58_9.nex -m TEST -ntop 10 -pre output
-  -s data/d59_8.phy -redo -q data/d58_9.nex -m TEST -nbest 5 -pre output
-  -s data/d59_8.phy -redo -q data/d58_9.nex -m TEST -nstop 100 -pre output
-  -s data/d59_8.phy -redo -q data/d58_9.nex -m TEST -pers 0.5 -pre output
-  -s data/d59_8.phy -redo -q data/d58_9.nex -m TEST -sprrad 6 -pre output
-  -s data/d59_8.phy -redo -q data/d58_9.nex -m TEST -bb 1000 -pre output
-  -s data/d59_8.phy -redo -q data/d58_9.nex -m TEST -bcor 0.5 -pre output
-  -s data/d59_8.phy -redo -q data/d58_9.nex -m TEST -beps 0.5 -pre output
-  -s data/d59_8.phy -redo -q data/d58_9.nex -m TEST -bnni -pre output
