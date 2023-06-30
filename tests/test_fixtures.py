import pytest
import os
from pathlib import Path
from utils import run_test_using

@pytest.mark.parametrize(
    "data_files", [(["sample_file1.txt", "sample_file2.txt"])], indirect=True
)
def test_fixtures(temp_dir, data_files, repo_paths):
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

    Examples
    --------
    To run this test, simply execute pytest from the repository root or the tests directory:

        pytest test_fixtures.py
    """

    # Assert that the current directory is a subdirectory of the /tests directory
    # and starts with the prefix 'tmp'
    current_dir = Path(os.getcwd())
    tests_root = repo_paths["tests_root"]
    assert current_dir.parent == tests_root
    assert current_dir.name.startswith("tmp")

    # Assert that the data_files specified are in the <temp>/data subdirectory
    tests_data = repo_paths["tests_data"]
    for data_file in data_files:
        assert (tests_data / data_file).is_file()

    # Assert that repo_paths build_dir, tests_dir, and data_dir are all ultimately in the path below the repo_root
    repo_root = repo_paths["repo_root"]
    assert repo_paths["build_dir"].parts[: len(repo_root.parts)] == repo_root.parts
    assert repo_paths["tests_root"].parts[: len(repo_root.parts)] == repo_root.parts
    assert repo_paths["tests_data"].parts[: len(repo_root.parts)] == repo_root.parts


if __name__ == "__main__":
    import subprocess
    import sys

    # Get the path to the current file
    current_path = Path(__file__).resolve()

    # Find the root directory by removing everything after the first 'tests'
    parts = current_path.parts
    test_index = parts.index("tests")
    root_dir = Path(*parts[: test_index + 1])

    # Run pytest for the current file
    subprocess.check_call(
        [sys.executable, "-m", "pytest", "--rootdir", str(root_dir), str(current_path)]
    )

if __name__ == "__main__":
    run_test_using(__file__)
