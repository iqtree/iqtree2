from pathlib import Path
from typing import Any, Dict
from .cache import Cache
from .iqtree import Iqtree
from .paths import tests_data, iqtree2_dir

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
        self.cache_file = tests_data / "iqtree2.cache.json"

    def get(self, alignment_file: str, parameters: str) -> Dict[str, Any]:
        cache = Cache(self.cache_file)
        cache.get(alignment_file, parameters, self.exec_iqtree_binary, Path(self.iqtree_binary).stat().st_mtime)
