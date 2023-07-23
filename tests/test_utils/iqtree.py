from pathlib import Path
from typing import Any, Dict

from .iqtree_stdout import Iqtree_stdout
from .checkpoint_file import Checkpoint
from .execute import exec_command
from .log_file import Log

class Iqtree:
    """
    Class for running generic IQ-TREE and processing the results.

    Parameters
    ----------
    iqtree_binary : Path
        The path to the IQ-TREE binary.

    Attributes
    ----------
    results : Dict[str, Any]
        The results of the IQ-TREE run.
    """

    def __init__(self, iqtree_binary: Path):
        self.iqtree_binary = iqtree_binary
        self.results = {}

    def process(self, alignment_file: str, parameters: str) -> Dict[str, Any]:
        """
        Run IQ-TREE and process the results.

        Parameters
        ----------
        alignment_file : str
            The path to the alignment file.
        parameters : str
            The parameters for the IQ-TREE run.

        Returns
        -------
        Dict[str, Any]
            The results of the IQ-TREE run.
        """
        self.results = self.exec_iqtree_binary(alignment_file, parameters)
        return self
    
    def exec_iqtree_binary(self, alignment_file: str, parameters: str)-> Dict[str, Any]:
        """
        Execute the IQ-TREE binary and return the results.

        Parameters
        ----------
        alignment_file : str
            The path to the alignment file.
        parameters : str
            The parameters for the IQ-TREE run.

        Returns
        -------
        Dict[str, Any]
            The results of the IQ-TREE run.
        """
        stdout = exec_command(
            str(self.iqtree_binary) + " -s " + str(alignment_file) + " " + parameters + " -redo -nt AUTO"
        )
        output = Iqtree_stdout(stdout.split("\n"))
        checkpoint = Checkpoint(str(alignment_file) + ".ckp.gz")
        log = Log(str(alignment_file) + ".log")
        return {"checkpoint": checkpoint, "log": log, "output": output}

    @property
    def checkpoint(self) -> Checkpoint:
        """
        Get the checkpoint from the results.

        Returns
        -------
        Checkpoint
            The checkpoint from the results.
        """
        return Checkpoint.from_dict(self.results.get("checkpoint")) 
        
    @property
    def log(self) -> Log:
        """
        Get the log from the results.

        Returns
        -------
        Log
            The log from the results.
        """
        return Log.from_dict(self.results.get("log"))

    @property
    def stdout(self) -> Iqtree_stdout:
        """
        Get the console output of the execution.

        Returns
        -------
        Iqtree_stdout
            The console output of the execution.
        """
        return Iqtree_stdout.from_dict(self.results.get("output"))
