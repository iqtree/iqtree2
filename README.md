IQ-TREE
=======

[![Github IQ-TREE 1 Releases](https://img.shields.io/github/downloads/Cibiv/IQ-TREE/total.svg?style=social&logo=github&label=iqtree1%20download)](https://github.com/Cibiv/IQ-TREE/releases)
[![Github IQ-TREE 2 Releases](https://img.shields.io/github/downloads/iqtree/iqtree2/total.svg?style=social&logo=github&label=iqtree2%20download)](https://github.com/iqtree/iqtree2/releases)
[![BioConda downloads](https://img.shields.io/conda/dn/bioconda/iqtree.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/iqtree)
[![BioConda platforms](https://img.shields.io/conda/pn/bioconda/iqtree?style=flag)](https://github.com/bioconda/bioconda-recipes/tree/master/recipes/iqtree)
[![Build Status](https://github.com/iqtree/iqtree2/workflows/Build/badge.svg)](https://github.com/iqtree/iqtree2/actions)
[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

Efficient and versatile phylogenomic software by maximum likelihood <http://www.iqtree.org>

Introduction
------------

The IQ-TREE software was created as the successor of IQPNNI and [TREE-PUZZLE](http://www.tree-puzzle.de) (thus the name IQ-TREE). IQ-TREE was motivated by the rapid accumulation of phylogenomic data, leading to a need for efficient phylogenomic software that can handle a large amount of data and provide more complex models of sequence evolution. To this end, IQ-TREE can utilize multicore computers and distributed parallel computing to speed up the analysis. IQ-TREE automatically performs checkpointing to resume an interrupted analysis.

As input IQ-TREE accepts all common sequence alignment formats including PHYLIP, FASTA, Nexus, Clustal and MSF. As output IQ-TREE will write a self-readable report file (name suffix `.iqtree`), a NEWICK tree file (`.treefile`)  which can be visualized by tree viewer programs such as [FigTree](http://tree.bio.ed.ac.uk/software/figtree/), [Dendroscope](http://dendroscope.org) or [iTOL](http://itol.embl.de).


Key features of IQ-TREE
-----------------------

* __Efficient search algorithm__: Fast and effective stochastic algorithm to reconstruct phylogenetic trees by maximum likelihood. IQ-TREE compares favorably to RAxML and PhyML in terms of likelihood while requiring similar amount of computing time ([Nguyen et al., 2015]).
* __Ultrafast bootstrap__: An ultrafast bootstrap approximation (UFBoot) to assess branch supports. UFBoot is 10 to 40 times faster than RAxML rapid bootstrap and obtains less biased support values ([Minh et al., 2013]).
* __Accurate model selection__: An ultrafast and automatic model selection (ModelFinder) which is 10 to 100 times faster than jModelTest and ProtTest. ModelFinder also finds best-fit partitioning scheme like PartitionFinder ([Kalyaanamoorthy et al., 2017]).
* __Alignment simulation__: A flexible simulator (AliSim) which
    allows to simulate sequence alignments under more realistic models than
    Seq-Gen and INDELible ([Ly-Trong et al., 2023]).
* __Phylogenetic testing__: Several fast branch tests like SH-aLRT and aBayes test ([Anisimova et al., 2011]) and tree topology tests like the approximately unbiased (AU) test ([Shimodaira, 2002]).


The strength of IQ-TREE is the availability of a wide variety of phylogenetic models:

* __Common models__: All [common substitution models](http://www.iqtree.org/doc/Substitution-Models) for DNA, protein, codon, binary and morphological data with [rate heterogeneity among sites](http://www.iqtree.org/doc/Substitution-Models/#rate-heterogeneity-across-sites) and [ascertainment bias correction](http://www.iqtree.org/doc/Substitution-Models/#ascertainment-bias-correction) for e.g. SNP data.
* __[Partition models](http://www.iqtree.org/doc/Complex-Models/#partition-models)__: Allowing individual models for different genomic loci (e.g. genes or codon positions), mixed data types, mixed rate heterogeneity types, linked or unlinked branch lengths between partitions.
* __Mixture Models__: [fully customizable mixture models](http://www.iqtree.org/doc/Complex-Models/#mixture-models) and [empirical protein mixture models](http://www.iqtree.org/doc/Substitution-Models/#protein-models) and.
* __Polymorphism-aware models (PoMo)__: <http://www.iqtree.org/doc/Polymorphism-Aware-Models>


IQ-TREE web service
-------------------

For a quick start you can also try the IQ-TREE web server, which performs 
online computation using a dedicated computing cluster. It is very easy to use 
with as few as just 3 clicks! Try it out at:

* Vienna Bioinformatics Cluster: <http://iqtree.cibiv.univie.ac.at>
* CIPRES Gateway: <https://www.phylo.org>
* Los Alamos Laboratories: <https://www.hiv.lanl.gov/content/sequence/IQTREE/iqtree.html>

User support
------------

Please refer to the [user documentation](http://www.iqtree.org/doc/) and 
[frequently asked questions](http://www.iqtree.org/doc/Frequently-Asked-Questions). 
If you have further questions and feedback, please create a topic at 
[Github discussions](https://github.com/iqtree/iqtree2/discussions).
For feature requests bug reports please post a topic at
[Github issues](https://github.com/iqtree/iqtree2/issues).

Citations
---------

General citation for IQ-TREE 2:

* B.Q. Minh, H.A. Schmidt, O. Chernomor, D. Schrempf, M.D. Woodhams, A. von Haeseler, R. Lanfear (2020) 
  IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era.
  *Mol. Biol. Evol.*, 37:1530-1534. <https://doi.org/10.1093/molbev/msaa015>

Moreover, there are other papers associated with notable features in IQ-TREE,
which are normally mentioned in the corresponding documentation. We ask that you also
cite these papers, which are important for us to obtain fundings to 
continuously maintain the code of IQ-TREE. These papers are also listed below.

When using tree mixture models (MAST) please cite:

* T.K.F. Wong, C. Cherryh, A.G. Rodrigo, M.W. Hahn, B.Q. Minh, R. Lanfear (2024)
  MAST: Phylogenetic Inference with Mixtures Across Sites and Trees.
  _Syst. Biol._, in press. <https://doi.org/10.1093/sysbio/syae008>

When computing concordance factors please cite:

* Y.K. Mo, R. Lanfear, M.W. Hahn, B.Q. Minh (2023)
  Updated site concordance factors minimize effects of homoplasy and taxon sampling.
  _Bioinformatics_, 39:btac741. <https://doi.org/10.1093/bioinformatics/btac741>

When using AliSim to simulate alignments please cite:

* N. Ly-Trong, G.M.J. Barca, B.Q. Minh (2023)
  AliSim-HPC: parallel sequence simulator for phylogenetics.
  *Bioinformatics*, 39:btad540. <https://doi.org/10.1093/bioinformatics/btad540>

When estimating amino-acid Q matrix please cite:

* B.Q. Minh, C. Cao Dang, L.S. Vinh, R. Lanfear (2021)
  QMaker: Fast and accurate method to estimate empirical models of protein evolution.
  _Syst. Biol._, 70:1046–1060. <https://doi.org/10.1093/sysbio/syab010>

When using the heterotachy GHOST model "+H" please cite:

* S.M. Crotty, B.Q. Minh, N.G. Bean, B.R. Holland, J. Tuke, L.S. Jermiin, A. von Haeseler (2020)
  GHOST: Recovering Historical Signal from Heterotachously Evolved Sequence Alignments.
  _Syst. Biol._, 69:249-264. <https://doi.org/10.1093/sysbio/syz051>

When using the tests of symmetry please cite:

* S. Naser-Khdour, B.Q. Minh, W. Zhang, E.A. Stone, R. Lanfear (2019) 
  The Prevalence and Impact of Model Violations in Phylogenetic Analysis. 
  *Genome Biol. Evol.*, 11:3341-3352. <https://doi.org/10.1093/gbe/evz193>

When using polymorphism-aware models please cite:

* D. Schrempf, B.Q. Minh, A. von Haeseler, C. Kosiol (2019) 
  Polymorphism-aware species trees with advanced mutation models, bootstrap, and rate heterogeneity. 
  *Mol. Biol. Evol.*, 36:1294–1301. <https://doi.org/10.1093/molbev/msz043>

For the ultrafast bootstrap (UFBoot) please cite:

* D.T. Hoang, O. Chernomor, A. von Haeseler, B.Q. Minh, and L.S. Vinh (2018) 
  UFBoot2: Improving the ultrafast bootstrap approximation. 
  *Mol. Biol. Evol.*, 35:518–522. <https://doi.org/10.1093/molbev/msx281>

When using posterior mean site frequency model (PMSF) please cite:

* H.C. Wang, B.Q. Minh, S. Susko, A.J. Roger (2018) 
  Modeling site heterogeneity with posterior mean site frequency profiles 
  accelerates accurate phylogenomic estimation. 
  *Syst. Biol.*, 67:216–235. <https://doi.org/10.1093/sysbio/syx068>

When using ModelFinder please cite:

* S. Kalyaanamoorthy, B.Q. Minh, T.K.F. Wong, A. von Haeseler, L.S. Jermiin (2017) 
  ModelFinder: Fast model selection for accurate phylogenetic estimates. 
  *Nat. Methods*, 14:587-589. <https://doi.org/10.1038/nmeth.4285>

When using partition models please cite:

* O. Chernomor, A. von Haeseler, B.Q. Minh (2016) 
  Terrace aware data structure for phylogenomic inference from supermatrices. 
  *Syst. Biol.*, 65:997-1008. <https://doi.org/10.1093/sysbio/syw037>

When using IQ-TREE web server please cite:

* J. Trifinopoulos, L.-T. Nguyen, A. von Haeseler, B.Q. Minh (2016) 
  W-IQ-TREE: a fast online phylogenetic tool for maximum likelihood analysis.
  *Nucleic Acids Res.*, 44:W232-W235. <https://doi.org/10.1093/nar/gkw256>

When using IQ-TREE version 1 please cite:

* L. Nguyen, H.A. Schmidt, A. von Haeseler, B.Q. Minh (2015)
  IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating Maximum-Likelihood Phylogenies.
  _Mol. Biol. and Evol._, 32:268-274. <https://doi.org/10.1093/molbev/msu300>

#### Credits and Acknowledgements

Some parts of the code were taken from the following packages/libraries: [Phylogenetic likelihood library](http://www.libpll.org), [TREE-PUZZLE](http://www.tree-puzzle.de), 
[BIONJ](http://dx.doi.org/10.1093/oxfordjournals.molbev.a025808), [Nexus Class Libary](http://dx.doi.org/10.1093/bioinformatics/btg319), [Eigen library](http://eigen.tuxfamily.org/),
[SPRNG library](http://www.sprng.org), [Zlib library](http://www.zlib.net), [gzstream library](http://www.cs.unc.edu/Research/compgeom/gzstream/), [vectorclass library](http://www.agner.org/optimize/), [GNU scientific library](https://www.gnu.org/software/gsl/).


IQ-TREE was funded by the [Austrian Science Fund - FWF](http://www.fwf.ac.at/) 
(grant no. I 760-B17 from 2012-2015 and and I 2508-B29 from 2016-2019),
the [University of Vienna](https://www.univie.ac.at/) (Initiativkolleg I059-N),
the [Australian National University](https://www.anu.edu.au),
[Chan-Zuckerberg Initiative](https://chanzuckerberg.com) (open source software for science grants),
[Simons Foundation](https://www.simonsfoundation.org), [Moore Foundation](https://www.moore.org),
and [Australian Research Council](https://www.arc.gov.au).


[Anisimova et al., 2011]: http://dx.doi.org/10.1093/sysbio/syr041
[Guindon et al., 2010]: http://dx.doi.org/10.1093/sysbio/syq010
[Kalyaanamoorthy et al., 2017]: https://doi.org/10.1038/nmeth.4285
[Ly-Trong et al., 2023]: https://doi.org/10.1093/bioinformatics/btad540
[Minh et al., 2013]: http://dx.doi.org/10.1093/molbev/mst024
[Nguyen et al., 2015]: http://dx.doi.org/10.1093/molbev/msu300
[Shimodaira, 2002]: http://dx.doi.org/10.1080/10635150290069913

