import os
import shutil
import tempfile
from pathlib import Path

import pytest

from .utils import iqtree1_dir, iqtree2_dir, repo_root, tests_data, tests_root


@pytest.fixture(scope="function")
def temp_dir(request, data_files):
    """
    Create a temporary directory with a data subdirectory, change to it for running tests, and clean up.

    This fixture creates a temporary directory and a data subdirectory within it,
    and changes the current working directory to this temporary directory before running tests.
    After the tests are complete, it changes back to the original directory and removes the
    temporary directory.

    Yields
    ------
    dict
        The updated repo_paths dictionary with the paths to the temporary directory and its data subdirectory.
    """

    # Save the current working directory
    original_dir = os.getcwd()

    # Create a temporary directory and change to it
    temp_path = Path(tempfile.mkdtemp(dir=tests_root))
    os.chdir(temp_path)

    # Create a data subdirectory
    data_subdir = temp_path / "data"
    data_subdir.mkdir()

    # Copy the files specified in data_files to the data subdirectory
    for data_file in data_files:
        shutil.copy(tests_data / data_file, data_subdir)

    # Yield control to the test function
    yield temp_path

    # Change back to the original directory and remove the temporary directory
    os.chdir(original_dir)
    shutil.rmtree(temp_path)


@pytest.fixture(scope="function")
def data_files(request):
    """
    Fixture to pass data files as parameters to the test function.

    This fixture is used in conjunction with pytest.mark.parametrize to pass
    data files as parameters to the test function. The data files should be
    specified as a list of filenames.  Setting indirect=True in the arguments
    of the pytest.mark.parametrize decorator indicates that the data_files will
    be passed through to the temp_dir fixture for copying to the temporary data
    directory.


    Parameters
    ----------
    request : _pytest.fixtures.SubRequest
        The pytest request object.

    Returns
    -------
    list
        The list of data files passed as parameters, or an empty list if no data files are specified.
    """
    return request.param if hasattr(request, "param") else []
