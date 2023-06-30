from cogent3 import load_aligned_seqs, get_model, make_tree, open_


def get_cogent3_result():
    aln = load_aligned_seqs("three-ungapped.fa", moltype="dna")
    tree = make_tree("(Human,Rhesus,Mouse)")
    # this is using average nuc freqs, which means it will match -iqtree -m HKY
    sm = get_model("HKY85")
    lf = sm.make_likelihood_function(tree)
    lf.set_alignment(aln)
    lf.optimise(show_progress=False)
    return lf
