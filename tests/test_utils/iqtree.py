from pathlib import Path
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
        A dictionary to store the results of the execution, including the checkpoint and log file.
    """

    def __init__(self, iqtree_binary: Path):
        self.iqtree_binary = iqtree_binary
        self.results = {
            "checkpoint": Checkpoint(),
            "log": Log()
        }
            
    @property
    def checkpoint(self) -> Checkpoint:
        return self.results["checkpoint"]
    
    @property
    def log(self) -> Log:
        return self.results["log"]

    def exec(self, alignment_file: str, parameters: str) -> "Iqtree":
        """
        Executes IQ-TREE and loads the checkpoint file.

        Parameters
        ----------
        alignment_file : str
            The path to the alignment file.
        parameters : str
            The parameters for the IQ-TREE command.

        Returns
        -------
        Iqtree
            The current instance.
        """
        _ = exec_command(
            str(self.iqtree_binary) + " -s " + str(alignment_file) + " " + parameters + " -redo -nt AUTO"
        )
        self.results["checkpoint"] = Checkpoint(str(alignment_file) + ".ckp.gz")
        self.results["log"] = Log(str(alignment_file) + ".log")
        return self
