from typing import Any, Dict
from cogent3 import open_
from yaml import Loader
from yaml import load as load_yml


class Checkpoint:
    """
    Class for loading and accessing data from an IQ-TREE checkpoint file.

    Parameters
    ----------
    file_path : str
        The path to the checkpoint file.

    Attributes
    ----------
    data : Dict[str, Any]
        The checkpoint file loaded as a dictionary.
    """

    def __init__(self, file_path: str = None):
        if file_path is not None:
            # load the checkpoint file using cogent3's open_ function to handle gzipped files
            self.data = load_yml(open_(file_path), Loader=Loader)
        else:
            self.data = None
            
    @classmethod
    def from_cache(cls, data: Dict[str, Any]) -> "Checkpoint":
        """
        Creates a Checkpoint instance from a cached instance.

        Parameters
        ----------
        data : Dict[str, Any]
            The checkpoint data as a dictionary.

        Returns
        -------
        Checkpoint
            The new Checkpoint instance.
        """
        instance = cls.__new__(cls)
        instance.data = data
        return instance

    @property
    def log_likelihood(self) -> float:
        """
        Gets the log likelihood from the checkpoint file.

        Returns
        -------
        float
            The log likelihood.
        """
        assert self.data is not None, "No checkpoint data found"
        # getting the key for the best tree - this has to be cast to int as the keys are strings and sometimes include a leading zero
        best = min(int(key) for key in self.data["CandidateSet"].keys())
        # extract log-likelihood
        lnL = float(self.data["CandidateSet"][str(best)].split()[0])
        return lnL
