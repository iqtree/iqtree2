from json import dump, load
from pathlib import Path
from re import sub
from subprocess import PIPE, Popen, check_call
from sys import executable
from typing import Any, Dict

from cogent3 import open_
from filelock import FileLock
from yaml import Loader
from yaml import load as load_yml

# Set global variables to point to important paths
current_dir = Path(__file__).parent
repo_root = current_dir
while repo_root != Path("/") and not (repo_root / ".git").is_dir():
    repo_root = repo_root.parent

# Define other important paths
iqtree2_dir = repo_root / "build"
tests_root = repo_root / "tests"
tests_data = tests_root / "data"
iqtree1_dir = tests_root / "iqtree1"


def exec_command(cmnd: str, stdout=PIPE, stderr=PIPE) -> str:
    """
    Executes a command and returns stdout if completes exit code 0.

    Parameters
    ----------
    cmnd : str
        The command to be executed.
    stdout, stderr : streams
        Default value (PIPE) intercepts process output, setting to None
        blocks this.

    Returns
    -------
    str
        The output of the command execution.
    """
    proc = Popen(cmnd, shell=True, stdout=stdout, stderr=stderr)
    out, err = proc.communicate()
    if proc.returncode != 0:
        msg = err
        raise RuntimeError(f"FAILED: {cmnd}\n{msg}")
    return out.decode("utf8") if out is not None else None


class Iqtree:
    """
    Base class for executing IQ-TREE (1 or 2).

    Parameters
    ----------
    iqtree_binary : pathlib.Path
        Path to the IQ-TREE binary.

    Attributes
    ----------
    checkpointYAML : Dict[str, Any]
        The checkpoint file loaded as a dictionary.
    """

    def __init__(self, iqtree_binary: Path):
        self.iqtree_binary = iqtree_binary
        self.checkpointYAML: Dict[str, Any] = None

    def exec(self, alignment_file: str, parameters: str) -> "Iqtree":
        """
        Executes IQ-TREE and loads the checkpoint file.

        Parameters
        ----------
        alignment_file : str
            The path to the alignment file.
        parameters : str
            The parameters for the IQ-TREE command.

        Returns
        -------
        Iqtree
            The current instance.
        """
        _ = exec_command(
            str(self.iqtree_binary)
            + " -s "
            + str(alignment_file)
            + " "
            + parameters
        )
        self.checkpointYAML = load_yml(
            open_(str(alignment_file) + ".ckp.gz"), Loader=Loader
        )
        return self

    def get_log_likelihood(self) -> float:
        """
        Gets the log likelihood from the checkpoint file.

        Returns
        -------
        float
            The log likelihood.
        """
        assert self.checkpointYAML is not None, "No checkpoint file found"
        # getting the key for the best tree
        best = min(int(key) for key in self.checkpointYAML["CandidateSet"].keys())
        # extract log-likelihood
        lnL = float(self.checkpointYAML["CandidateSet"][str(best)].split()[0])
        return lnL


class Iqtree1(Iqtree):
    """
    Class for running IQ-TREE 1 with caching and file locking.

    Attributes
    ----------
    cache_file : str
        The path to the cache file.
    cache : Dict[str, Any]
        The cache loaded from the cache file.
    """

    def __init__(self):
        super().__init__(iqtree1_dir / "iqtree")
        self.cache_file = tests_data / "cache1.json"
        self.cache: Dict[str, Any] = {}

    def exec(self, alignment_file: str, parameters: str) -> "Iqtree1":
        """
        If this parameter set is in the cache, loads the previously generated
        checkpoint file. If this parameter set has not been cached executes
        IQ-TREE 1, loads the resulting checkpoint file, and adds the parameters
        and the result to the cache.

        Parameters
        ----------
        alignment_file : str
            The path to the alignment file.
        parameters : str
            The parameters for the IQ-TREE 1.

        Returns
        -------
        Iqtree1
            The current instance.
        """
        alignment_file_name = Path(alignment_file).name
        parameters = sub(r"-s\s+\S+", f"-s {alignment_file_name}", parameters)
        cache_key = f"{parameters}"

        lock = FileLock(str(self.cache_file) + ".lock")

        with lock:
            if Path(self.cache_file).exists():
                with open(self.cache_file, "r") as f:
                    self.cache = load(f)
            else:
                self.cache = {}

            if cache_key in self.cache:
                self.checkpointYAML = self.cache[cache_key]
            else:
                super().exec(alignment_file, parameters)
                self.cache[cache_key] = self.checkpointYAML
                with open(self.cache_file, "w") as f:
                    dump(self.cache, f)
                # reread the checkpoint from the cache file as dump converts integer keys to strings
                with open(self.cache_file, "r") as f:
                    self.cache = load(f)
                self.checkpointYAML = self.cache[cache_key]
                    

        return self


class Iqtree2(Iqtree):
    """
    Class for running IQ-TREE 2 with caching, cache invalidation when the
    binary changes, and file locking.

    Attributes
    ----------
    cache_file : str
        The path to the cache file.
    cache : Dict[str, Any]
        The cache loaded from the cache file.
    """

    def __init__(self):
        super().__init__(iqtree2_dir / "iqtree2")
        self.cache_file = tests_data / "cache2.json"

    def exec(self, alignment_file: str, parameters: str) -> "Iqtree2":
        """
        If this parameter set is in the cache, and the binary has not been
        updated since that result was stored, loads the previously generated
        checkpoint file. If this parameter set has not been cached, or the
        cache has been invalidated, executes IQ-TREE 2, loads the resulting
        checkpoint file, and adds the parameters and the result to the cache.

        Parameters
        ----------
        alignment_file : str
            The path to the alignment file.
        parameters : str
            The parameters for the IQ-TREE 2.

        Returns
        -------
        Iqtree2
            The current instance.
        """
        alignment_file_name = Path(alignment_file).name
        parameters = sub(r"-s\s+\S+", f"-s {alignment_file_name}", parameters)
        key = f"{parameters}"
        iqtree2_mod_time = Path(self.iqtree_binary).stat().st_mtime

        lock = FileLock(str(self.cache_file) + ".lock")
        cache : Dict[str, Any] = {}


        with lock:
            if Path(self.cache_file).exists():
                with open(self.cache_file, "r") as f:
                    cache = load(f)
            else:
                cache = {"mod_time": iqtree2_mod_time}

            if cache.get("mod_time") != iqtree2_mod_time:
                cache = {"mod_time": iqtree2_mod_time}

            if key in cache:
                self.checkpointYAML = cache[key]
            else:
                super().exec(alignment_file, parameters)
                cache[key] = self.checkpointYAML
                with open(self.cache_file, "w") as f:
                    dump(cache, f)
                # reread the checkpoint from the cache file as dump converts integer keys to strings
                with open(self.cache_file, "r") as f:
                    cache = load(f)
                self.checkpointYAML = cache[key]

        return self
