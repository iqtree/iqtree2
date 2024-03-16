from typing import Any, Dict


class Log:
    """
    Class for loading and accessing data from a log file.

    Parameters
    ----------
    file_path : str, optional
        The path to the log file. If None, the data attribute is set to None.

    Attributes
    ----------
    data : Dict[str, Any]
        The log file loaded as a dictionary containing a list of lines.
    """

    def __init__(self, file_path: str = None):
        # If file_path is not None, load the log file
        # Otherwise, set data to None
        if file_path is not None:
            self.data = self._load_log_file(file_path)
        else:
            self.data = None

    @staticmethod
    def _load_log_file(file_path: str) -> Dict[str, Any]:
        """
        Load a log file into a Dict containing a list of lines.
        """
        data = {}
        with open(file_path, "r") as file:
            data["lines"] = file.readlines()
        return data

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "Log":
        """
        Create a Log instance from a dictionary.
        """
        instance = cls.__new__(cls)
        instance.data = data
        return instance

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert the Log instance to a dictionary.
        """
        return self.data

    @property
    def best_score(self) -> float:
        """
        Get the best score from the log file.

        Returns
        -------
        float
            The best score.
        """
        for line in self.data.get("lines", []):
            if "BEST SCORE FOUND" in line:
                _, score = line.split(":")
                return float(score.strip())
        # Raise an error if no best score is found
        raise ValueError("No best score found in .log file")
