from json import dump, load
from pathlib import Path

from filelock import FileLock


class Cache:
    """
    Class for managing the cache file.

    Parameters
    ----------
    cache_file_name : str
        The name of the cache file.
    timestamp : float, optional
        The timestamp of the binary used to create results that have been cached

    Attributes
    ----------
    file_name : str
        The name of the cache file.
    lock : FileLock
        The file lock for the cache file.
    """

    def __init__(self, cache_file_name, timestamp=None):
        self.file_name = cache_file_name
        self.lock = FileLock(str(cache_file_name) + ".lock")
        self.timestamp = timestamp

    def get(
        self,
        alignment_file: str,
        parameters: str,
        generate_new_result,
    ) -> dict[str, any]:
        alignment_file_name = Path(alignment_file).name
        cache_key = f"{alignment_file_name} {parameters}"
        cache = {"mod_time": self.timestamp}
        with self.lock:
            try:
                with open(self.file_name, "r") as f:
                    cache = load(f)
            except FileNotFoundError:
                cache = {"mod_time": self.timestamp}
            # invalidate the cache if the timestamp has changed
            if cache.get("mod_time") != self.timestamp:
                cache = {"mod_time": self.timestamp}
            result = cache.get(cache_key)
            if result is None:
                result = generate_new_result()
                # Convert the result into a dictionary
                result = {
                    key: value.to_dict() if hasattr(value, "to_dict") else value
                    for key, value in result.items()
                }
                cache[cache_key] = result
                with open(self.file_name, "w") as f:
                    dump(cache, f)
            return result
