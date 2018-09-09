## A Short Guide to the Codebase

This is a short and rather shallow guide to the codebase. Starting from main.cpp, you'll see the crucial steps and data structures involved in solving the terrace detection and enumeration problem.

### 0) Interpreting the CLI Arguments
If two file paths are provided, it is assumed that the first path points to the tree file (.nwk) and the second path points to the data file (.data).

If one file path Y/X is provided, it is assumed that the files X.nwk and X.data are present in the directory Y.

### 1) Parsing the Input Files
This is done by calling terraces::parse_nwk and terraces::parse_bitmatrix. Both throw exceptions, if the input files are not in the right format. If the .data file does not contain a species that possesses all gene sites, this is denoted by a terraces::none value in the std::pair returned by terraces::parse_bitmatrix. In this case, the current course of action is to exit with error code 1.

### 2) Re-Rooting the Input Tree
A call to terraces::reroot_at_taxon_inplace re-roots the given tree at the given species. This is done in-place by traversing the input tree from the given species to the original root, and adjusting all edges so that the parent reference of every node points to the node that was traversed before. This procedure has O(tree-height) time complexity.

### 3) Extracting Subtrees
Using an occurence bitmatrix (is gene i present in species j?) and the re-rooted tree, all subtrees of this tree are extracted.

### 4) Computing Constraints
By collecting the left- and rightmost nodes x and y for every inner node i, we establish lowest-common-ancestor relationships of the form lca(x, y) = i. With this information, we then obtain our constraints of the form lca(x, y) < lca(x, z) which mean that the lowest common ancestor node od x and y is a descendant of the lca of x and z.

### 5) Deduplicating Constaints
Constraint calculation can result in constraints being found multiple times. terraces::deduplicate_constraints ensures that the given constraints vector contains only unique elements by removing all duplicates. 

### 6) Supertree-Assembly and Output
Now that we have extracted the constraints from the re-rooted version of our tree, we assemble our supertrees according to Constantinescu's Algorithm [CS95] which merges sets of nodes so that for each constraint lca(i, j) < lca(k, l), the sets of nodes i and j are merged. This is done for every constraint. In the end, if at least two node sets remain, there's at least one supertree. This gives us the leaves under the two subtrees under our new root's non-comprehensive-leaf child. In these two subtrees, we then use the constraints to further narrow down the tree's overall structure. Doing so, we can end up with subtress with more than two subtrees and without constraints. In this case, we enumerate all possible subtree topologies.


### References
[CS95] Mariana Constantinescu and David Sankoff. An efficient algorithm for supertrees. Journal of Classification, 12(1):101–112, 1995.