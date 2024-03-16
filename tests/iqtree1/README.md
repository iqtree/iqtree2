# IQ-TREE1 

This directory contains the legacy IQ-TREE version 1 binary. This binary is used for regression testing and is not included in the main IQ-TREE2 repository.

To build the IQ-TREE1 binary, clone the IQ-TREE1 repository as follows

```bash
git clone git@github.com:iqtree/iqtree1.git
```

Then, build the binary using the following commands

```bash
cd iqtree1
mkdir build
cd build
cmake ..
make
```

Then copy the binary to the `iqtree2/tests/iqtree` directory

```bash
cp /iqtree1/iqtree /iqtree2/tests/iqtree1/iqtree1
```

# Running the tests

To run the tests, follow the instructions in the [Running the Tests](../docs/running_tests.md) document.

# Caching the IQ-TREE1 results

IQ-Tree1 regression tests compare the results from IQ-Tree2 with the results from IQ-Tree1 using the same input parameters.  As IQ-Tree1 is no longer under active development, the results from IQ-Tree1 are not expected to change, so the results from IQ-Tree1 are cached in the directory `tests/iqtree/cached`.

