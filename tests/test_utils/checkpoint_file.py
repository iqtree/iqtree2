from typing import Any, Dict

from cogent3 import open_
from yaml import Loader
from yaml import load as load_yml


class Checkpoint:
    """
    Class for loading and accessing data from an IQ-TREE checkpoint file.
    This class provides methods to load the checkpoint file and extract the log likelihood.

    Parameters
    ----------
    file_path : str
        The path to the checkpoint file. If None, the data attribute is set to None.

    Attributes
    ----------
    data : Dict[str, Any]
        The checkpoint file loaded as a dictionary.
    """

    def __init__(self, file_path: str = None):
        # If file_path is not None, load the checkpoint file using cogent3's open_ function to handle gzipped files
        # Otherwise, set data to None
        if file_path is not None:
            self.data = load_yml(open_(file_path), Loader=Loader)
        else:
            self.data = None

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "Checkpoint":
        """
        Create a Checkpoint instance from a dictionary.
        """
        instance = cls.__new__(cls)
        instance.data = data
        return instance

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert the Checkpoint instance to a dictionary.
        """
        return self.data

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
        # Getting the key for the best tree. The key is cast to int as the keys are strings and sometimes include a leading zero.
        # The best tree is the one with the smallest key.
        best = min(int(key) for key in self.data["CandidateSet"].keys())
        # Note: the CandidateSet key in the YAML returned by IQTREE is an integer,
        # but the key once it's been put into the cache as JSON is a string so we have to check both
        best_CandidateSet = self.data["CandidateSet"].get(
            str(best), self.data["CandidateSet"].get(best)
        )
        # Extract log-likelihood from the best CandidateSet
        lnL = float(best_CandidateSet.split()[0])
        return lnL
