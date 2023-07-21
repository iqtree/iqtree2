from pathlib import Path
from typing import Any, Dict
from .checkpoint_file import Checkpoint
from .execute import exec_command
from .log_file import Log

class Iqtree:
    """
    Base class for running generic IQ-TREE (1 or 2) binaries.

    Parameters
    ----------
    iqtree_binary : pathlib.Path
        Path to the IQ-TREE binary.

    Attributes
    ----------
    results : Dict[str, Union[Checkpoint, Logfile]]
        A dictionary to store the results of the query, 
        including the checkpoint and log file.
    """

    def __init__(self, iqtree_binary: Path):
        self.iqtree_binary = iqtree_binary
        self.results = {}

    def calculate(self, alignment_file: str, parameters: str) -> "Iqtree":
        self.results = self.get_result(alignment_file, parameters)
        return self

    def get_result(self, alignment_file: str, parameters: str) -> Dict[str, Any]:
        return self.exec_iqtree_binary(alignment_file, parameters)
    
    def exec_iqtree_binary(self, alignment_file: str, parameters: str)-> Dict[str, Any]:
        _ = exec_command(
            str(self.iqtree_binary) + " -s " + str(alignment_file) + " " + parameters + " -redo -nt AUTO"
        )
        checkpoint = Checkpoint(str(alignment_file) + ".ckp.gz")
        log = Log(str(alignment_file) + ".log")
        return {"checkpoint": checkpoint, "log": log}

    @property
    def checkpoint(self) -> Checkpoint:
        return self.results.get("checkpoint") 
        
    @property
    def log(self) -> Log:
        return self.results.get("log")