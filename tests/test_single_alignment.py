import os
from pathlib import Path
import pytest
from .utils import exec_command, iqtree1_log_liklihood, iqtree2_log_liklihood, run_test_using
from numpy.testing import assert_allclose


# @pytest.mark.parametrize("options", ["-cmin 2", "-mtree"])
@pytest.mark.parametrize("options", ["-cmin 2"])
@pytest.mark.parametrize("data_files", [(["example.phy"])], indirect=True)
def test_single_alignment(temp_dir, data_files, repo_paths, options):
    """
    single alignment comparison between iqtree1 and iqtree2

    Parameters
    ----------
    options : list
        Sample list of options from legacy testing framework's config.yaml
    """
    iqtree2_binary = repo_paths["build_dir"] / "iqtree2"
    assert iqtree2_binary.is_file(), "IQ-Tree2 binary not found"
    assert options is not None, "No options specified"
    assert len(data_files) == 1, "Only one alignment file should be specified per test"
    alignment_file = temp_dir / "data" / data_files[0]
    iqtree_params = " -s "+ str(alignment_file) + " " + options + " -m TEST"
    
    lnL1 = iqtree1_log_liklihood(repo_paths["iqtree1_dir"],iqtree_params, str(alignment_file)+".ckp.gz")
    
    _ = exec_command(str(iqtree2_binary)+" "+iqtree_params)
    lnL2 = iqtree2_log_liklihood(str(alignment_file)+".ckp.gz")    
    # hope they're the same!
    assert_allclose(lnL1, lnL2)

if __name__ == "__main__":
    run_test_using(__file__)
