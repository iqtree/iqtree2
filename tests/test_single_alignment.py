import os
from pathlib import Path

import pytest
from numpy import allclose
from numpy.testing import assert_allclose

from .test_utils import (Iqtree1, Iqtree2, iqtree1_dir, iqtree2_dir, repo_root,
                         tests_data, tests_root)


@pytest.mark.parametrize("options", ["-cmin 2", "-nbest 5"])
@pytest.mark.parametrize(
    "data_files", [(["example.phy"]), (["d59_8.phy"])], indirect=True
)
def test_single_alignment_via_checkpoint(temp_dir, data_files, options):
    """
    single alignment comparison between iqtree1 and iqtree2, comparing log likelihoods from checkpoint files

    Parameters
    ----------
    options : list
        Sample list of options from legacy testing framework's config.yaml
    data_files : list[tuple[list[str]]]
        Sample data files from legacy testing framework
    """
    iqtree2_binary = iqtree2_dir / "iqtree2"
    assert iqtree2_binary.is_file(), "IQ-Tree2 binary not found"
    assert options is not None, "No options specified"
    assert len(data_files) == 1, "Only one alignment file should be specified per test"
    alignment_file = temp_dir / "data" / data_files[0]
    iqtree_params = " " + options + " -m TEST"

    lnL1 = Iqtree1().process(alignment_file, iqtree_params).checkpoint.log_likelihood
    lnL2 = Iqtree2().process(alignment_file, iqtree_params).checkpoint.log_likelihood

    # hope they're the same!
    assert_allclose(
        lnL1,
        lnL2,
        rtol=1e-6,
        err_msg=f"Log-likelihood results from IQ-Tree1 ({lnL1}) and IQ-Tree2 ({lnL2}) differ with {options}",
    )


@pytest.mark.parametrize("options", ["-cmin 2", "-nbest 5"])
@pytest.mark.parametrize(
    "data_files", [(["example.phy"]), (["d59_8.phy"])], indirect=True
)
def test_single_alignment_via_log(temp_dir, data_files, options):
    """
    single alignment comparison between iqtree1 and iqtree2, comparing BEST SCORE FOUND from log files

    Parameters
    ----------
    options : list
        Sample list of options from legacy testing framework's config.yaml
    data_files : list[tuple[list[str]]]
        Sample data files from legacy testing framework
    """
    iqtree2_binary = iqtree2_dir / "iqtree2"
    assert iqtree2_binary.is_file(), "IQ-Tree2 binary not found"
    assert options is not None, "No options specified"
    assert len(data_files) == 1, "Only one alignment file should be specified per test"
    alignment_file = temp_dir / "data" / data_files[0]
    iqtree_params = " " + options + " -m TEST"

    lnL1 = Iqtree1().process(alignment_file, iqtree_params).log.best_score
    lnL2 = Iqtree2().process(alignment_file, iqtree_params).log.best_score

    # hope they're the same!
    assert_allclose(
        lnL1,
        lnL2,
        rtol=1e-6,
        err_msg=f"Log-likelihood results from IQ-Tree1 ({lnL1}) and IQ-Tree2 ({lnL2}) differ",
    )
