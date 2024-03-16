# IQ-TREE Testing Directory

This directory contains tests for the IQ-TREE project. The tests are written in Python using the `pytest` framework and are located in the `/tests` path. Each test script contains one or more test functions, and are named `test_<feature being tested>.py`. The `pytest` framework is used to execute these test functions and report the results.

## Directory Structure

- `tests/`: Contains the Python test scripts.
- `tests/data`: Contains sample data files used for testing.
  - `tests/data/iqtree?.cache.json`: are cached results stored in this directory when tests are run. 
- `tests/iqtree1`: An IQ-TREE 1.x binary should be built and copied into this directory for regression testing.  This binary will not be stored in the repository as it should be specific to the OS platform on which the tests will be run.
- `tests/test_utils` : Contains utility functions used by the test scripts to run iqtree binaries, and cache their results.
- `tests/conftest.py`: Contains fixtures for the pytest framework, that include the creation of temporary directories for each test, and the copying of test specific `tests/data` data files to a `data` subdirectory under the tests temporary directory .
- `requirements.txt`: Contains the Python dependencies required for testing.

## Getting Started

- [How to write tests](docs/writing_tests.md)
- [How to run tests](docs/running_tests.md)
- [Regression tests to be ported](docs/legacy_tests.md)

## Contributing

Contributions to the test suite are welcome. If you are implementing new features or fixing bugs in the IQ-TREE code, please consider adding tests that verify the correctness of your changes.