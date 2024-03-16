This document provides instructions on how to set up the environment and run the tests.

## Prerequisites

- Python 3.11 or higher
- C++ Compiler (e.g., GCC, Clang)
- CMake

### Building the iqtree2 Binary

Before running the tests, you need to build the `iqtree2` binary. Navigate to the root of the IQ-TREE repository and run the following commands:

```sh
rm -rf build
mkdir -p build
cd build
cmake ..
make
```

This will create the `iqtree2` binary in the `build` directory.

## Setting Up the Environment

### Installing Python 3.11

1. Go to the official Python downloads page: [Python Downloads](https://www.python.org/downloads/).

2. Download Python 3.11 for your operating system (Windows, macOS, or Linux).

3. Run the installer and follow the instructions to install Python 3.11 on your system.

4. Verify the installation by opening a terminal or command prompt and running:

    ```sh
    python --version
    ```

    This should display the version of Python that you installed.

### Setting Up a Virtual Environment

1. Navigate to the `tests` directory.
    
        ```sh
        cd tests
        ```
2. It is recommended to create a virtual environment to isolate the dependencies. Navigate to the project directory and run:

    ```sh
    python -m venv venv
    ```

3. Activate the virtual environment.

    - On Windows:

        ```sh
        .\venv\Scripts\activate
        ```

    - On macOS and Linux:

        ```sh
        source venv/bin/activate
        ```

4. Install the Python dependencies from the `requirements.txt` file.

    ```sh
    pip install -r requirements.txt
    ```

## Running the Tests

1. Navigate to the `tests` directory.
    
        ```sh
        cd tests
        ```
2. Activate the virtual environment.

    - On Windows:

        ```sh
        .\venv\Scripts\activate
        ```

    - On macOS and Linux:

        ```sh
        source venv/bin/activate
        ```
3. Run the tests using `pytest`.

    ```sh
    pytest
    ```

This document provides instructions on how to set up the environment and run the tests.

## Prerequisites

- Python 3.11 or higher
- C++ Compiler (e.g., GCC, Clang)
- CMake

### Building the iqtree2 Binary

Before running the tests, you need to build the `iqtree2` binary. Navigate to the root of the IQ-TREE repository and run the following commands:

```sh
rm -rf build
mkdir -p build
cd build
cmake ..
make
```

This will create the `iqtree2` binary in the `build` directory.

## Setting Up the Environment

### Installing Python 3.11

1. Go to the official Python downloads page: [Python Downloads](https://www.python.org/downloads/).

2. Download Python 3.11 for your operating system (Windows, macOS, or Linux).

3. Run the installer and follow the instructions to install Python 3.11 on your system.

4. Verify the installation by opening a terminal or command prompt and running:

    ```sh
    python --version
    ```

    This should display the version of Python that you installed.

### Setting Up a Virtual Environment

1. Navigate to the `tests` directory.
    
        ```sh
        cd tests
        ```
2. It is recommended to create a virtual environment to isolate the dependencies. Navigate to the project directory and run:

    ```sh
    python -m venv venv
    ```

3. Activate the virtual environment.

    - On Windows:

        ```sh
        .\venv\Scripts\activate
        ```

    - On macOS and Linux:

        ```sh
        source venv/bin/activate
        ```

4. Install the Python dependencies from the `requirements.txt` file.

    ```sh
    pip install -r requirements.txt
    ```

## Running the Tests

1. Navigate to the `/tests` directory.
    
2. Activate the virtual environment (see above).

3. Run all the tests by executing `pytest`.