IQ-TREE
-------

Efficient phylogenetic software by maximum likelihood

Please see our github wiki for more information: <https://github.com/Cibiv/IQ-TREE/wiki>

IQ-TREE PoMo
------------

IQ-TREE+PoMo is still under development.  Please check out

    iqtree --help

Especially, the section titled `POLYMORPHISM AWARE MODELS (PoMo)`.

```
POLYMORPHISM AWARE MODELS (PoMo):
PoMo is run when
- a Counts File is used as input file, and/or when
- it is specified in the model string (see below).
  -st C[FR] or C[FR]ps Counts File (automatically detected).
                       Useful to customize the virtual population size `ps`
                       3 <= ps <= 19; ps has to be an odd number, 2 or 10.
                       F: Sum over partial likelihoods at the tip of the tree (weighted).
                       R: Random binomial sampling of PoMo states from data (sampled).
                       Default is `CF9`.
  -m <sm>+<pm>+<ft>    Default: `HKY+rP+FO`.
                 <sm>: Substitution model.
                  DNA: HKY (default), JC, F81, K2P, K3P, K81uf, TN/TrN, TNef,
                       TIM, TIMef, TVM, TVMef, SYM, GTR, or a 6-digit model
                       specification (e.g., 010010 = HKY).
                 <pm>: PoMo model.
                       - rP (default; reversible PoMo with tree inference).
                 <ft>: Frequency type (optional; default: +F, counted).
                       F or +FO or +FU or +FQ.
                       Counted, optimized, user-defined, equal state frequency.
                       This overwrites the specifications of the DNA model.
  The default model string is: -m HKY+rP+F.
  Until now, only DNA models work with PoMo.
  Model testing and rate heterogeneity do not work with PoMo yet.
```
