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


def iqtree2_log_liklihood(from_yml_file):
    from yaml import load as load_yml, Loader
    from cogent3 import open_

    data = load_yml(open_(from_yml_file), Loader=Loader)
    # getting the key for the best tree
    best = min(data["CandidateSet"])
    # extract log-likelihood
    lnL = float(data["CandidateSet"][best].split()[0])
    return lnL

def iqtree1_log_liklihood(iqtree1_dir,iqtree1_params, iqtree1_output):
    from yaml import load as load_yml, Loader
    from cogent3 import open_
    
    iqtree1_binary = iqtree1_dir / "iqtree1"
    _ = exec_command(str(iqtree1_binary)+" "+iqtree1_params)
    
    data = load_yml(open_(iqtree1_output), Loader=Loader)
    # getting the key for the best tree
    best = min(data["CandidateSet"])
    # extract log-likelihood
    lnL = float(data["CandidateSet"][best].split()[0])
    return lnL


def exec_command(cmnd, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
    """executes shell command and returns stdout if completes exit code 0

    Parameters
    ----------

    cmnd : str
      shell command to be executed
    stdout, stderr : streams
      Default value (PIPE) intercepts process output, setting to None
      blocks this."""
    proc = subprocess.Popen(cmnd, shell=True, stdout=stdout, stderr=stderr)
    out, err = proc.communicate()
    if proc.returncode != 0:
        msg = err
        raise RuntimeError(f"FAILED: {cmnd}\n{msg}")
    return out.decode("utf8") if out is not None else None
