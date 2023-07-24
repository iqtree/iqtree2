from pathlib import Path
from typing import Any, Dict

from .checkpoint_file import Checkpoint
from .execute import exec_command
from .iqtree_stdout import Iqtree_stdout
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

        Stub method that should be overridden by subclasses.

        Returns
        -------
        Iqtree
            The instance so that the process method can be chained.
        """
        self.results = self.invoke_iqtree(alignment_file, parameters)
        return self

    def invoke_iqtree(self, alignment_file: str, parameters: str) -> Dict[str, Any]:
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
            The results of the IQ-TREE invokation.
        """
        stdout = exec_command(
            str(self.iqtree_binary)
            + " -s "
            + str(alignment_file)
            + " "
            + parameters
            + " -redo -nt AUTO"  # invariant parameters for all tests: overwrite existing results, use all available threads
        )
        output = Iqtree_stdout(stdout.split("\n"))
        checkpoint = Checkpoint(str(alignment_file) + ".ckp.gz")
        log = Log(str(alignment_file) + ".log")
        return {"checkpoint": checkpoint, "log": log, "output": output}

    @property
    def checkpoint(self) -> Checkpoint:
        """
        Get the checkpoint from the results.
        """
        return Checkpoint.from_dict(self.results.get("checkpoint"))

    @property
    def log(self) -> Log:
        """
        Get the log from the results.
        """
        return Log.from_dict(self.results.get("log"))

    @property
    def stdout(self) -> Iqtree_stdout:
        """
        Get the console output of the execution.
        """
        return Iqtree_stdout.from_dict(self.results.get("output"))
