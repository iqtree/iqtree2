# remove from git repo
import os
import pathlib
from pathlib import Path

import pytest

from .test_utils import (iqtree1_dir, iqtree2_dir, repo_root, tests_data,
                         tests_root)


@pytest.mark.parametrize(
    "data_files", [(["sample_file1.txt", "sample_file2.txt"])], indirect=True
)
def test_fixtures(temp_dir, data_files):
    """
    Test the functionality of the custom fixtures.

    This test function asserts that:
    - The current directory is a temporary directory.
    - The specified data_files are in the <temp>/data subdirectory.
    - The repo_paths build_dir, tests_dir, and data_dir are all ultimately in the path below the repo_root.

    The test function is parameterized with a list of data files using the pytest.mark.parametrize
    decorator. The 'indirect=True' argument indicates that the 'data_files' argument should be passed
    through a fixture for processing.

    Parameters
    ----------
    temp_dir : Path
        Path object representing the temporary directory created by the temp_dir fixture.
    data_files : list
        List of data files processed by the data_files fixture.
    repo_paths : dict
        Dictionary containing paths to the repository root, build directory, tests directory,
        and data directory. It is created by the repo_paths fixture.
    """

    # Assert that the current directory is a subdirectory of the /tests directory
    # and starts with the prefix 'tmp'
    current_dir = Path(os.getcwd())
    assert current_dir.parent == tests_root
    assert current_dir.name.startswith("tmp")

    # Assert that the data_files specified are in the <temp>/data subdirectory
    for data_file in data_files:
        data_file_path = temp_dir / "data" / data_file
        assert data_file_path.is_file()

    # Assert that repo_paths iqtree2_dir, tests_dir, and data_dir are all ultimately in the path below the repo_root
    assert iqtree2_dir.parts[: len(repo_root.parts)] == repo_root.parts
    assert tests_root.parts[: len(repo_root.parts)] == repo_root.parts
    assert tests_data.parts[: len(repo_root.parts)] == repo_root.parts
    assert iqtree1_dir.parts[: len(repo_root.parts)] == repo_root.parts
