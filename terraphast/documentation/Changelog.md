Since 1st meeting
=================

Optimizations
-------------

* implement bitvector with rank support
* use bitvector instead of index vectors for subsets
* inline lots of often-used methods
* deduplicate constraints
* remap constraints (removing inner nodes from the numbering)
* halved union-find storage by out-of-bounds parent trick
* avoid allocations by reusing old storage
* some micro-optimizations

Fixes and Features
------------------

* fast terrace check without traversing the tree to the end
* fixed subtree computation (problematic input data caused an assertion to fail)
* fix counting bug (communication problem)
* don't count incorrectly rooted trees
* extracted all computations to callback methods
* implemented logging + stack state decorators
* implemented isomorphy check
* Visual C++ compatibility checked with appveyor
* extracted intrinsics to stay portable
