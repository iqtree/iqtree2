from json import dump, load
from pathlib import Path
from re import sub
from typing import Any, Dict
from filelock import FileLock

from .log_file import Log
from .checkpoint_file import Checkpoint
from .iqtree import Iqtree
from .paths import iqtree1_dir, tests_data


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
        self.cache_file = tests_data / "iqtree1.cache.json"
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
                    cache = load(f)
            else:
                cache = {}

            if cache_key in cache:
                self.results["checkpoint"] = Checkpoint.from_cache(cache[cache_key])
                self.results["log"] = Log.from_cache(cache[cache_key])
            else:
                super().exec(alignment_file, parameters)
                cache[cache_key] = self.results
                with open(self.cache_file, "w") as f:
                    dump(cache, f)
                # reread the checkpoint from the cache file as dump converts integer keys to strings
                with open(self.cache_file, "r") as f:
                    cache = load(f)
                self.results["checkpoint"] = Checkpoint.from_cache(cache[cache_key])
                self.results["log"] = Log.from_cache(cache[cache_key])
        return self
