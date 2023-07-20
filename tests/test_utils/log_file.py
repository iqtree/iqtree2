from typing import Any, Dict

class Log:
    """
    Class for loading and accessing data from an IQ-TREE log file.

    Parameters
    ----------
    file_path : str
        The path to the log file.

    Attributes
    ----------
    data : Dict[str, Any]
        The log file loaded as a dictionary.
    """

    def __init__(self, file_path: str = None):
        if file_path is not None:
            self.data = self._load_log_file(file_path)
        else:
            self.data = None

    @staticmethod
    def _load_log_file(file_path: str) -> Dict[str, Any]:
        data = {}
        with open(file_path, 'r') as file:
            data["lines"] = file.readlines()
        return data

    @classmethod
    def from_cache(cls, data: Dict[str, Any]) -> "Log":
        """
        Creates a LogFile instance from a cached instance.

        Parameters
        ----------
        data : Dict[str, Any]
            The log file data as a dictionary.

        Returns
        -------
        LogFile
            The new LogFile instance.
        """
        instance = cls.__new__(cls)
        instance.data = data
        return instance

    @property
    def best_score(self) -> float:
        """
        Gets the best score from the log file.

        Returns
        -------
        float
            The best score.
        """
        for line in self.data.get("lines", []):
            if "BEST SCORE FOUND" in line:
                _, score = line.split(":")
                return float(score.strip())
        raise ValueError("No best score found in .log file")
