import subprocess
import sys
from pathlib import Path

def run_test_using(filename):
    """
    Utility function to run pytest for the given test file.

    This function is intended to be used in the __main__ block of a test file.
    It runs pytest for the given test file.

    Parameters
    ----------
    filename : str
        The __file__ attribute of the test script.
    """
    
    # Get the path to the current file
    current_path = Path(filename).resolve()

    # Find the root directory by removing everything after the first 'tests'
    parts = current_path.parts
    test_index = parts.index("tests")
    root_dir = Path(*parts[: test_index + 1])

    # Run pytest for the current file
    subprocess.check_call(
        [sys.executable, "-m", "pytest", "--rootdir", str(root_dir), str(current_path)]
    )
