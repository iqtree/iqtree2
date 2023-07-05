# IQ-TREE Testing Directory

This directory contains tests for the IQ-TREE project. The tests are written in Python using the `pytest` framework and are located in the `/tests` path. Each test script contains one or more test functions, and are named `test_<feature being tested>.py`. The `pytest` framework is used to execute these test functions and report the results.

## Directory Structure

- `tests/`: Contains the Python test scripts.
- `tests/data`: Contains sample data files used for testing.
- `tests/bin`: Contains binaries including the iqtree version 1 binary.
- `tests/conftest.py`: Contains fixtures for the pytest framework. Fixtures are reusable pieces of code that can be used to set up resources (like temporary directories), as well as pass data or configurations to tests.
- `requirements.txt`: Contains the Python dependencies required for testing.

## Getting Started

- [How to write tests](docs/writing_tests.md)
- [How to run tests](docs/running_tests.md)
- [Regression tests to be ported](docs/legacy_tests.md)

## Contributing

Contributions to the test suite are welcome. If you are implementing new features or fixing bugs in the IQ-TREE code, please consider adding tests that verify the correctness of your changes.