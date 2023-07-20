import os
from pathlib import Path

import pytest
from numpy.testing import assert_allclose

from .utils import (
    Iqtree1,
    Iqtree2,
    iqtree1_dir,
    iqtree2_dir,
    repo_root,
    tests_data,
    tests_root,
)


# @pytest.mark.parametrize("options", ["-cmin 2", "-nbest 5"])
@pytest.mark.xfail
@pytest.mark.parametrize("options", ["-cmin 2", "-nbest 5"])
@pytest.mark.parametrize("data_files", [(["example.phy"])], indirect=True)
def test_single_alignment(temp_dir, data_files, options):
    """
    single alignment comparison between iqtree1 and iqtree2

    Parameters
    ----------
    options : list
        Sample list of options from legacy testing framework's config.yaml
    """
    iqtree2_binary = iqtree2_dir / "iqtree2"
    assert iqtree2_binary.is_file(), "IQ-Tree2 binary not found"
    assert options is not None, "No options specified"
    assert len(data_files) == 1, "Only one alignment file should be specified per test"
    alignment_file = temp_dir / "data" / data_files[0]
    iqtree_params = " " + options + " -m TEST"

    lnL1 = Iqtree1().exec(alignment_file, iqtree_params).checkpoint.log_likelihood
    lnL2 = Iqtree2().exec(alignment_file, iqtree_params).checkpoint.log_likelihood
    # hope they're the same!
    assert_allclose(lnL1, lnL2)
