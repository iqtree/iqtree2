IQ-TREE
=======

Efficient and versatile phylogenomic software by maximum likelihood <http://www.iqtree.org>

Introduction
------------

The IQ-TREE software was created as the successor of IQPNNI and [TREE-PUZZLE](http://www.tree-puzzle.de) (thus the name IQ-TREE). IQ-TREE was motivated by the rapid accumulation of phylogenomic data, leading to a need for efficient phylogenomic software that can handle a large amount of data and provide more complex models of sequence evolution. To this end, IQ-TREE can utilize multicore computers and distributed parallel computing to speed up the analysis. IQ-TREE automatically performs checkpointing to resume an interrupted analysis.

As input IQ-TREE accepts all common sequence alignment formats including PHYLIP, FASTA, Nexus, Clustal and MSF. As output IQ-TREE will write a self-readable report file (name suffix `.iqtree`), a NEWICK tree file (`.treefile`)  which can be visualized by tree viewer programs such as [FigTree](http://tree.bio.ed.ac.uk/software/figtree/), [Dendroscope](http://dendroscope.org) or [iTOL](http://itol.embl.de).


Key features of IQ-TREE
-----------------------

* __Efficient search algorithm__: Fast and effective stochastic algorithm to reconstruct phylogenetic trees by maximum likelihood. IQ-TREE compares favorably to RAxML and PhyML in terms of likelihood while requiring similar amount of computing time ([Nguyen et al., 2015]).
* __Ultrafast bootstrap__: An ultrafast bootstrap approximation (UFBoot) to assess branch supports. UFBoot is 10 to 40 times faster than RAxML rapid bootstrap and obtains less biased support values ([Minh et al., 2013]).
* __Ultrafast model selection__: An ultrafast and automatic model selection (ModelFinder) which is 10 to 100 times faster than jModelTest and ProtTest. ModelFinder also finds best-fit partitioning scheme like PartitionFinder.
* __Phylogenetic testing__: Several fast branch tests like SH-aLRT and aBayes test ([Anisimova et al., 2011]) and tree topology tests like the approximately unbiased (AU) test ([Shimodaira, 2002]).


The strength of IQ-TREE is the availability of a wide variety of phylogenetic models:

* __Common models__: All [common substitution models](http://www.iqtree.org/doc/Substitution-Models) for DNA, protein, codon, binary and morphological data with [rate heterogeneity among sites](http://www.iqtree.org/doc/Substitution-Models/#rate-heterogeneity-across-sites) and [ascertainment bias correction](http://www.iqtree.org/doc/Substitution-Models/#ascertainment-bias-correction) for e.g. SNP data.
* __[Partition models](http://www.iqtree.org/doc/Complex-Models/#partition-models)__: Allowing individual models for different genomic loci (e.g. genes or codon positions), mixed data types, mixed rate heterogeneity types, linked or unlinked branch lengths between partitions.
* __Mixture Models__: [fully customizable mixture models](http://www.iqtree.org/doc/Complex-Models/#mixture-models) and [empirical protein mixture models](http://www.iqtree.org/doc/Substitution-Models/#protein-models) and.

IQ-TREE web service
-------------------

For a quick start you can also try the IQ-TREE web server, which performs online computation using a dedicated computing cluster. It is very easy to use with as few as just 3 clicks! Try it out at

<http://iqtree.cibiv.univie.ac.at>


User support
------------

Please refer to the [user documentation](http://www.iqtree.org/doc/) and [frequently asked questions](http://www.iqtree.org/doc/Frequently-Asked-Questions). If you have further questions, feedback, feature requests, and bug reports, please sign up the following Google group (if not done yet) and post a topic to the 

<https://groups.google.com/d/forum/iqtree>

_The average response time is one working day._

Citations
---------

To cite IQ-TREE please use:

* L.-T. Nguyen, H.A. Schmidt, A. von Haeseler, and B.Q. Minh (2015) IQ-TREE: A fast and effective stochastic algorithm for estimating maximum likelihood phylogenies. *Mol. Biol. Evol.*, 32, 268-274. [DOI: 10.1093/molbev/msu300](http://dx.doi.org/10.1093/molbev/msu300)

For the ultrafast bootstrap (UFBoot) please cite:

* B.Q. Minh, M.A.T. Nguyen, and A. von Haeseler (2013) Ultrafast approximation for phylogenetic bootstrap. *Mol. Biol. Evol.*, 30:1188-1195. [DOI: 10.1093/molbev/mst024](http://dx.doi.org/10.1093/molbev/mst024)

#### Credits and Acknowledgements

Some parts of the code were taken from the following packages/libraries: [Phylogenetic likelihood library](http://www.libpll.org), [TREE-PUZZLE](http://www.tree-puzzle.de), 
[BIONJ](http://dx.doi.org/10.1093/oxfordjournals.molbev.a025808), [Nexus Class Libary](http://dx.doi.org/10.1093/bioinformatics/btg319), [Eigen library](http://eigen.tuxfamily.org/),
[SPRNG library](http://www.sprng.org), [Zlib library](http://www.zlib.net), gzstream library, [vectorclass library](http://www.agner.org/optimize/), [GNU scientific library](https://www.gnu.org/software/gsl/).


IQ-TREE was partially funded by the [Austrian Science Fund - FWF](http://www.fwf.ac.at/) (grant no. I760-B17 from 2012-2015 and and I 2508-B29 from 2016-2019) and the [University of Vienna](https://www.univie.ac.at/) (Initiativkolleg I059-N).


[Anisimova et al., 2011]: http://dx.doi.org/10.1093/sysbio/syr041
[Guindon et al., 2010]: http://dx.doi.org/10.1093/sysbio/syq010
[Minh et al., 2013]: http://dx.doi.org/10.1093/molbev/mst024
[Nguyen et al., 2015]: http://dx.doi.org/10.1093/molbev/msu300
[Shimodaira, 2002]: http://dx.doi.org/10.1080/10635150290069913
