import os
from pathlib import Path
import pytest
from utils import run_test_using 

@pytest.mark.parametrize(
    "data_files", [(["sample_file1.txt", "sample_file2.txt"])], indirect=True
)
def test_sample_feature(temp_dir, data_files, repo_paths):
    """
    Sample test function to demonstrate the usage of fixtures.

    This test function can be used as a template for writing new tests.
    It demonstrates how to use the temp_dir, data_files, and repo_paths fixtures.

    Parameters
    ----------
    temp_dir : Path
        Path object representing the temporary directory created by the temp_dir fixture.
    data_files : list
        List of data files processed by the data_files fixture.
    repo_paths : dict
        Dictionary containing paths to the repo_root, build_dir,tests_dir, tests_data.
        It is created by the repo_paths fixture.
    """

    data_dir = Path(temp_dir) / "data"

    # Iterate through each file in data_files and check if it exists in data_dir
    for file_name in data_files:
        file_path = data_dir / file_name
        assert file_path.is_file(), f"{file_name} does not exist in {data_dir}"
    
    # Add your test logic here
    assert True

if __name__ == "__main__":
    run_test_using(__file__)
