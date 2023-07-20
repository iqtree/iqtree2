from json import dump, load
from pathlib import Path
from re import sub
from typing import Any, Dict

from filelock import FileLock

from .checkpoint_file import Checkpoint
from .iqtree import Iqtree
from .paths import iqtree1_dir, tests_data, iqtree2_dir


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
        self.cache_file = tests_data / "iqtree2.cache"

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
        cache_key = f"{parameters}"
        iqtree2_mod_time = Path(self.iqtree_binary).stat().st_mtime

        lock = FileLock(str(self.cache_file) + ".lock")
        cache: Dict[str, Any] = {}

        with lock:
            if Path(self.cache_file).exists():
                with open(self.cache_file, "r") as f:
                    cache = load(f)
            else:
                cache = {"mod_time": iqtree2_mod_time}

            if cache.get("mod_time") != iqtree2_mod_time:
                cache = {"mod_time": iqtree2_mod_time}

            if cache_key in cache:
                self.checkpoint = Checkpoint.from_cache(cache[cache_key])
            else:
                super().exec(alignment_file, parameters)
                cache[cache_key] = self.checkpoint.data
                with open(self.cache_file, "w") as f:
                    dump(cache, f)
                # reread the checkpoint from the cache file as dump converts integer keys to strings
                with open(self.cache_file, "r") as f:
                    cache = load(f)
                self.checkpoint = Checkpoint.from_cache(cache[cache_key])

        return self
