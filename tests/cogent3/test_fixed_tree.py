import subprocess
import pathlib
from numpy.testing import assert_allclose
import pytest
from .cogent3_utils import get_cogent3_result
from tests.utils import iqtree2_log_liklihood, exec_command, run_test_using


@pytest.mark.parametrize("data_files", [(["three-ungapped.fa"])], indirect=True)
def test_iqtree_simple(data_files, temp_dir, repo_paths):
    outpath = pathlib.Path("three-ungapped.fa.ckp.gz")
    # need better determination of path
    build_path = repo_paths["build_dir"]
    cmnd = build_path / "iqtree2 -s data/three-ungapped.fa -m HKY -redo"
    _ = exec_command(cmnd)
    lnL = iqtree2_log_liklihood(outpath)
    # cogent3 value
    c3_lf = get_cogent3_result()
    # hope they're the same!
    assert_allclose(lnL, c3_lf.lnL)


if __name__ == "__main__":
    run_test_using(__file__)
