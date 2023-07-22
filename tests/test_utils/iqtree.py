from pathlib import Path
from typing import Any, Dict
from .checkpoint_file import Checkpoint
from .execute import exec_command
from .log_file import Log

class Iqtree:
    def __init__(self, iqtree_binary: Path):
        self.iqtree_binary = iqtree_binary
        self.results = {}

    def process(self, alignment_file: str, parameters: str) -> Dict[str, Any]:
        self.results = self.exec_iqtree_binary(alignment_file, parameters)
        return self
    
    def exec_iqtree_binary(self, alignment_file: str, parameters: str)-> Dict[str, Any]:
        stdout = exec_command(
            str(self.iqtree_binary) + " -s " + str(alignment_file) + " " + parameters + " -redo -nt AUTO"
        )
        checkpoint = Checkpoint(str(alignment_file) + ".ckp.gz")
        log = Log(str(alignment_file) + ".log")
        return {"checkpoint": checkpoint, "log": log}

    @property
    def checkpoint(self) -> Checkpoint:
        return Checkpoint.from_dict(self.results.get("checkpoint")) 
        
    @property
    def log(self) -> Log:
        return Log.from_dict(self.results.get("log"))