import os
import pathlib
import subprocess

import pytest
import utils

from cogent3 import get_model, load_aligned_seqs, make_tree, open_
from numpy.testing import assert_allclose

from .test_utils import Iqtree2, iqtree2_dir, repo_root, tests_data, tests_root


def get_cogent3_result(alignment_file):
    aln = load_aligned_seqs(alignment_file, moltype="dna")
    tree = make_tree("(Human,Rhesus,Mouse)")
    # this is using average nuc freqs, which means it will match -iqtree -m HKY
    sm = get_model("HKY85")
    lf = sm.make_likelihood_function(tree)
    lf.set_alignment(aln)
    lf.optimise(show_progress=False)
    return lf


@pytest.mark.parametrize("data_files", [(["three-ungapped.fa"])], indirect=True)
def test_against_cogent3(data_files, temp_dir):
    assert len(data_files) == 1, "Only one alignment file should be specified per test"
    iqtree2_binary = iqtree2_dir / "iqtree2"
    assert iqtree2_binary.is_file(), "IQ-Tree2 binary not found"
    alignment_file = temp_dir / "data" / data_files[0]
    iqtree2_params = " -m HKY -redo"
    lnL = Iqtree2().process(alignment_file, iqtree2_params).checkpoint.log_likelihood

    # cogent3 log_likelihood from the same alignment
    c3_lf = get_cogent3_result(alignment_file)
    # hope they're the same!
    assert_allclose(lnL, c3_lf.lnL)
