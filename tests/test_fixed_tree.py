import os
import subprocess
import pathlib
from numpy.testing import assert_allclose
import pytest
from .utils import iqtree2_log_liklihood, exec_command, run_test_using
from cogent3 import load_aligned_seqs, get_model, make_tree, open_

def get_cogent3_result(alignment_file):
    aln = load_aligned_seqs(alignment_file, moltype="dna")
    tree = make_tree("(Human,Rhesus,Mouse)")
    # this is using average nuc freqs, which means it will match -iqtree -m HKY
    sm = get_model("HKY85")
    lf = sm.make_likelihood_function(tree)
    lf.set_alignment(aln)
    lf.optimise(show_progress=False)
    return lf

@pytest.mark.parametrize(
    "data_files", [(["three-ungapped.fa"])], indirect=True
)
def test_iqtree_simple(data_files, temp_dir, repo_paths):
    # need better determination of path
    build_path = repo_paths["build_dir"]
    cmnd = build_path / "iqtree2 -s data/three-ungapped.fa -m HKY -redo"
    current = os.getcwd()
    _ = exec_command(str(cmnd))
    lnL = iqtree2_log_liklihood("data/three-ungapped.fa.ckp.gz")
    # cogent3 value
    c3_lf = get_cogent3_result("data/three-ungapped.fa")
    # hope they're the same!
    assert_allclose(lnL, c3_lf.lnL)


if __name__ == "__main__":
    run_test_using(__file__)
