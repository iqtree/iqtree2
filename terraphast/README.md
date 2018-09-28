Documentation Home                         {#mainpage}
==================

C++ tool to check for and enumerate terraces in phylogenetic tree space.

-----

**Usage**: `terraces/build/release/app <nwk file> <gene/site file>`

Terraphast takes a .nkw file in Newick format and a genes/sites file, which denotes whether (1) or not (0) gene i is present in species j.

Program output states some imput data properties, the species whose leaf edge is used as a new tree root, and the resulting supertree in compressed newick format.



**Compressed Newick Format**: The resulting supertree representation cann be plain Newick, but can also contain the following two notation enhancements:
- `{a,b,c}` represents any conceivable binary subtree comprising the taxa a, b, and c.
- `(A|B,C|D)` represents any conceivable binary subtree comprising either subtrees A or B on the left, and either subtrees C or D on the right branch.

Both enhancements were chosen such that the result is standard newick format if there's only one possible supertree.



## The Terrace Phenomenon and Problem
In recent years, it has become common practice to infer phylogenies on so-called multi-gene datasets. Concatenated multi-gene datasets usually exhibit holes, that is, sequence data for some species might not be available for some genes Gi in our concatenated dataset. This can be due to a plethora of reasons, for instance, a specific species might simply not have a specific gene G i or the specific gene has simply not been sequenced for some of the species. After concatenating genes (partitions) we therefore end up with an alignment that contains patches of missing data:

```
index       0123

Species 1   AC--
Species 2   AG--
Species 3   ACTT
Species 4   --AG
Species 5   --GG
```

Under the likelihood model conditions that generate terraces, the log likelihood LnL(T) of a tree T can be computed as follows: LnL(T) = LnL(T|G1) + LnL(T|G2) where T|Gi denotes the tree topology induced by T for the species/sequences in partition i for which we have sequence data. In our example, the trees induced by G1 and G2 contain only three taxa. We know that there's only one tree topology with three taxa. On the other hand, there are 15 possible topologys for 5-taxa trees. So all 15 possible 5-taxon trees for our example dataset will induce the same per-gene/partition trees and therefore span a terrace of size 15.
This example dataset is bad: It does not contain any signal for disentangling
the phylogenetic history of these 5 species, since they are only connected via species 3.

**Terraces**: two distinct comprehensive (containing all n species) trees are on a terrace if all induced per-partition subtrees of the two trees are identical. This phenomenon was named and described in [SMS11].

Knowing about the phenomenon of terraces, researchers might want to know
(i) if a given tree is on a terrace, 
(ii) how many trees there are on that terrace, and 
(iii) how the trees on that terrace look like.

## The Basic Approach

TO PUT IN HERE:
-  Take NWK and DATA file and re-root the tree so that the comprehensive taxon is a leaf under the root, and the rest of the tree is a subtree under the root.
-  Extract constraints according to Constantinescu's algorithm (Only reference and very short outline)
-  Using the constraints, generate trees according to C's algorithm
-  No guarantees for completeness (simply unknown)

## A Short Guide to the Code

This can be found [here](documentation/walkthrough.md).

## Improvements and Optimizations to the basic approach

### Implemented:

- We introduced an optional, **compressed tree output format**. This format makes printing to terminal faster, since not all possible trees are listed in full detail. See section *Enhanced Newick Format* above. [https://git.scc.kit.edu/bioinfo2017/terraces/issues/8]
- **Memory allocation in large blocks**, and managing them with free lists. [https://git.scc.kit.edu/bioinfo2017/terraces/issues/37, https://git.scc.kit.edu/bioinfo2017/terraces/issues/13]
- **Deletion of unnecessary constraints**. [https://git.scc.kit.edu/bioinfo2017/terraces/issues/23]
- **Improved data structures**: We replaced index vectors representing the current leaves by bitvectors with rank support, thus improving space requirements and the efficiency of constraint filtering. The union-find data structure could be improved by storing the set ranks in out-of-bounds indices, thus halving the storage. [https://git.scc.kit.edu/bioinfo2017/terraces/issues/29]
- **Remap constraints**: By removing inner nodes from the constraint numbering, we were able to halve the space requirements of most of our data structures. [https://git.scc.kit.edu/bioinfo2017/terraces/issues/21]
- **Use specialized bit manipulation instructions**: Bipartition iteration and bitvector operations (bit iteration and rank computation) were improved significantly by using specialized CPU instructions supported by Compiler intrinsics.
- **Provide a fast terrace check**: We discovered that checking for the existence of multiple trees on a terrace can be done without explicitly building any of these trees, thus decreasing the runtime even further.
- **Implemented validation methods**: For checking whether the trees generated by our algorithm are indeed distinct and equivalent to the input tree with respect to the missing data, we implemented a fast isomorphy check operating directly on our data structures.
- Since we wanted to implement different versions of the algorithm, we used a **generic enumerator** which relies on callback methods to implement the concrete version of the algorithm (terrace checking, tree counting or multitree construction). This method also allows us to attach **logging and status update decorators** to the algorithm to check the internal computations or monitor the progress of the algorithm.
- Short of support for arbitrary-precision math, our implementation is **fully compatible with Visual C++** in addition to the normal gcc/clang support.

### Planned:

- Enumerate subtrees in **parallel**. One challenge would be separation of the workload so that multiple threads have "enough to do". Another challenge would be the merging of the individual threads' results. [https://git.scc.kit.edu/bioinfo2017/terraces/issues/6]
- Finding **good heuristics** for choosing a subtree into which we want to descend first. Ideally, we'd have a nice heuristic that tells us which subtrees and associated constraints probably give us several options to construct a supertree, in which case we can safely answer "yes" to the question "are we on a terrace?". **Ideas**: "Smallest subtree first", "Least constraints first". [https://git.scc.kit.edu/bioinfo2017/terraces/issues/3]
- ...

## References
[SMS11] Michael J Sanderson, Michelle M McMahon, and Mike Steel. Terraces in phylogenetic
tree space. Science, 333(6041):448–450, 2011.