from pathlib import Path
from typing import Any, Dict

from .cache import Cache
from .iqtree import Iqtree
from .paths import iqtree2_dir, tests_data


class Iqtree2(Iqtree):
    """
    Class for running IQ-TREE 2 with caching of previously calculated results
    The cache is invalidated when the binary is changed

    Attributes
    ----------
    cache_file : str
        The path to the cache file.
    cache : Dict[str, Any]
        The cache loaded from the cache file.
    """

    def __init__(self):
        super().__init__(iqtree2_dir / "iqtree2")
        self.cache_file = tests_data / "iqtree2.cache.json"

    def process(self, alignment_file: str, parameters: str) -> Dict[str, Any]:
        """
        Run IQ-TREE 2 with caching, cache invalidation when the binary changes.

        Parameters
        ----------
        alignment_file : str
            The path to the alignment file.
        parameters : str
            The parameters for the IQ-TREE run.

        Returns
        -------
        Iqtree2
            The instance so that the process method can be chained.
        """
        cache = Cache(self.cache_file, Path(self.iqtree_binary).stat().st_mtime)
        self.results = cache.get(
            alignment_file,
            parameters,
            lambda: self.invoke_iqtree(alignment_file, parameters),
        )
        return self
