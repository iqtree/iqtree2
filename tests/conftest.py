import pytest
import tempfile
import shutil
import os
from pathlib import Path

@pytest.fixture(scope="session")
def repo_paths():
    """
    Determine the repository root and other important paths.
    
    This fixture determines the repository root by searching for the .git directory
    and defines other important paths relative to the repository root.
    
    Returns
    -------
    dict
        A dictionary containing paths to the repository root, build directory,
        tests directory, and data directory. The keys are 'repo_root', 'build_dir',
        'tests_dir', and 'data_dir'.
    """

   
    # Find the repository root by searching for the .git directory
    current_dir = Path(__file__).parent
    repo_root = current_dir
    while repo_root != Path('/') and not (repo_root / '.git').is_dir():
        repo_root = repo_root.parent
    
    # Define other important paths
    build_dir = repo_root / 'build'
    tests_dir = repo_root / 'tests'
    data_dir = tests_dir / 'data'
    
    # Return the paths in a dictionary
    return {
        'repo_root': repo_root,
        'build_dir': build_dir,
        'tests_dir': tests_dir,
        'data_dir': data_dir
    }

@pytest.fixture(scope="function")
def temp_dir(request, repo_paths):
    """
    Create a temporary directory with a data subdirectory for running tests.
    
    This fixture creates a temporary directory and a data subdirectory within it
    inside the /tests directory.
    
    If the test has declared specific data files, they are copied to the data subdirectory.

    Yields
    ------
    pathlib.Path
        The file path of the temporary directory.
    """
    
    # Save the current working directory
    original_dir = os.getcwd()
    
    # Create a temporary directory inside the /tests directory
    tests_dir = repo_paths['tests_dir']
    temp_path = Path(tempfile.mkdtemp(dir=tests_dir))

    # Change the current working directory to the temporary directory
    os.chdir(temp_path)

    # Create a data subdirectory
    data_subdir = temp_path / "data"
    data_subdir.mkdir()

    # If the test has declared specific data files, copy them to the data subdirectory
    if hasattr(request, "param"):
        data_dir = repo_paths['data_dir']

        for data_file in request.param:
            shutil.copy(data_dir / data_file, data_subdir)

    try:
        # Yield control to the test function
        yield temp_path
    finally:
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
        The list of data files passed as parameters.
    """
    
    return request.param
