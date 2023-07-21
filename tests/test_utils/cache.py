from json import dump, load
from pathlib import Path
from filelock import FileLock

class Cache:
    def __init__(self, cache_file_name, time_stamp=None):
        self.file_name = cache_file_name
        self.lock = FileLock(str(cache_file_name) + ".lock")
    
    def get(self, alignment_file: str, parameters: str, generate_new_result, timestamp : float = None)->dict[str,any]:
        alignment_file_name = Path(alignment_file).name
        cache_key = f"{alignment_file_name} {parameters}"
        cache = {"mod_time" : timestamp}
        with self.lock:
            try:
                with open(self.file_name, "r") as f:
                    cache = load(f)
            except FileNotFoundError:
                cache = {"mod_time" : timestamp}
            # invalidate the cache if the timestamp has changed
            if cache.get("mod_time") != timestamp:
                cache = {"mod_time" : timestamp}
            result = cache.get(cache_key)
            if result is None:
                result = generate_new_result()
                cache[cache_key] = result    
                with open(self.file_name, "w") as f:
                    dump(cache, f)
            return result    