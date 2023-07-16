import os
from pathlib import Path
import pytest
from .utils import run_test_using


@pytest.mark.parametrize("options", ["-cmin 2", "-mtree"])
@pytest.mark.parametrize("data_files", [(["example.phy"])], indirect=True)
def test_single_alignment(temp_dir, data_files, repo_paths, options):
    """
    single alignment comparison between iqtree1 and iqtree2

    Parameters
    ----------
    options : list
        Sample list of options from legacy testing framework's config.yaml
    """

    # Add your test logic here
    assert True


if __name__ == "__main__":
    run_test_using(__file__)
