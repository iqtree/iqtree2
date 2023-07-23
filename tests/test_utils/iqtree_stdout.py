from typing import Any, Dict, List


class Iqtree_stdout:
    def __init__(self, iqtree_output: List[str] = None):
        # If iqtree_output is not None, load the lines into data
        # Otherwise, set data to None
        if iqtree_output is not None:
            self.data = {}
            self.data["lines"] = iqtree_output
        else:
            self.data = None

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "Iqtree_stdout":
        instance = cls.__new__(cls)
        instance.data = data
        return instance

    def to_dict(self) -> Dict[str, Any]:
        return self.data
