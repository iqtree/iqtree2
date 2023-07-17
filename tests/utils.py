import json
import os
import pathlib
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

class Iqtree:
    def __init__(self, iqtree_binary):
        self.iqtree_binary = iqtree_binary
        self.checkpointYAML = None
        self.invariantparams = " -redo "
        
    def exec(self, alignment_file, parameters):
        _ = exec_command(str(self.iqtree_binary)+" -s "+str(alignment_file)+" "+parameters+ " " + self.invariantparams)
        from yaml import load as load_yml, Loader
        from cogent3 import open_
        self.checkpointYAML = load_yml(open_(str(alignment_file)+".ckp.gz"), Loader=Loader)
        return self
        
    def get_log_likelihood(self):
        assert self.checkpointYAML is not None, "No checkpoint file found"
        # getting the key for the best tree
        best = min(int(key) for key in self.checkpointYAML["CandidateSet"].keys())
        # get section of checkpoint file with best tree
        try:
            section = self.checkpointYAML["CandidateSet"][best]
        except KeyError:
            try:
                section = self.checkpointYAML["CandidateSet"][str(best)]
            except KeyError:
                assert False, "No best tree found"
        
        # extract log-likelihood
        lnL = float(section.split()[0])
        return lnL

class Iqtree2(Iqtree):
    def __init__(self):
        super().__init__(iqtree2_dir / "iqtree2")
        self.invariantparams += "  "

class Iqtree1(Iqtree):
    def __init__(self):
        super().__init__(iqtree1_dir / "iqtree")
        self.invariantparams += "  "
        self.cache_file = iqtree1_dir / "cache.json"  # The path to the cache file
        # Load the cache from the file if it exists, otherwise start with an empty cache
        if os.path.exists(self.cache_file):
            with open(self.cache_file, 'r') as f:
                self.cache = json.load(f)
        else:
            self.cache = {}
            
    def exec(self, alignment_file, parameters):
        # Use the base name of the alignment file as part of the key
        alignment_file_name = pathlib.Path(alignment_file).name
        key = f" -s {alignment_file_name} {parameters}"

        # Check if the results are in the cache
        if key in self.cache:
            # Load the results from the cache
            self.checkpointYAML = self.cache[key]
        else:
            # Execute the command and load the checkpoint file
            super().exec(alignment_file, parameters)
            # Save the results to the cache
            self.cache[key] = self.checkpointYAML
            # Save the cache to the file
            with open(self.cache_file, 'w') as f:
                json.dump(self.cache, f)

        return self        
    
def iqtree2_log_liklihood(from_yml_file):
    from yaml import load as load_yml, Loader
    from cogent3 import open_

    data = load_yml(open_(from_yml_file), Loader=Loader)
    # getting the key for the best tree
    best = min(int(key) for key in data["CandidateSet"].keys())
    # extract log-likelihood
    lnL = float(data["CandidateSet"][best].split()[0])
    return lnL

def iqtree1_log_liklihood(iqtree1_dir,iqtree1_params, iqtree1_output):
    from yaml import load as load_yml, Loader
    from cogent3 import open_
    
    iqtree1_binary = iqtree1_dir / "iqtree"
    _ = exec_command(str(iqtree1_binary)+" "+iqtree1_params)
    
    data = load_yml(open_(iqtree1_output), Loader=Loader)
    # getting the key for the best tree
    
    best = min(int(key) for key in data["CandidateSet"].keys())
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

# set global variables to point to important paths 
current_dir = Path(__file__).parent
repo_root = current_dir
while repo_root != Path("/") and not (repo_root / ".git").is_dir():
    repo_root = repo_root.parent

# Define other important paths
iqtree2_dir = repo_root / "build"
tests_root = repo_root / "tests"
tests_data = tests_root / "data"
iqtree1_dir = tests_root / "iqtree1"