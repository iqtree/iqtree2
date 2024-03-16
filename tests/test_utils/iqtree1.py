from pathlib import Path
from typing import Any, Dict

from .cache import Cache
from .iqtree import Iqtree
from .paths import iqtree1_dir, tests_data


class Iqtree1(Iqtree):
    """
    Class for running IQ-TREE 1 with caching of previously calculated results

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

    def process(self, alignment_file: str, parameters: str) -> Dict[str, Any]:
        """
        Run IQ-TREE 1 with caching

        Parameters
        ----------
        alignment_file : str
            The path to the alignment file.
        parameters : str
            The parameters for the IQ-TREE run.

        Returns
        -------
        Iqtree1
            The instance so that the process method can be chained.
        """
        cache = Cache(self.cache_file)
        self.results = cache.get(
            alignment_file,
            parameters,
            lambda: self.invoke_iqtree(alignment_file, parameters),
        )
        return self
