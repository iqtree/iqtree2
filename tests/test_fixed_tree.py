import os
import subprocess
import pathlib
from numpy.testing import assert_allclose
import pytest
from .utils import Iqtree2, iqtree2_log_liklihood, exec_command, run_test_using, repo_root, tests_root, tests_data, iqtree2_dir, iqtree1_dir
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
def test_iqtree_simple(data_files, temp_dir):
    assert len(data_files) == 1, "Only one alignment file should be specified per test"
    iqtree2_binary = iqtree2_dir / "iqtree2"
    assert iqtree2_binary.is_file(), "IQ-Tree2 binary not found"
    alignment_file = temp_dir / "data" / data_files[0]
    iqtree2_params = " -s "+ str(alignment_file) +" -m HKY -redo"
    lnL = Iqtree2().exec(alignment_file, iqtree2_params).get_log_likelihood()

    # cogent3 log_likelihood from the same alignment
    c3_lf = get_cogent3_result(alignment_file)
    # hope they're the same!
    assert_allclose(lnL, c3_lf.lnL)

if __name__ == "__main__":
    run_test_using(__file__)
