add_library(tree
constrainttree.cpp
constrainttree.h
candidateset.cpp candidateset.h
iqtree.cpp
iqtree.h
matree.cpp
matree.h
memslot.cpp memslot.h
mexttree.cpp
mexttree.h
mtree.cpp
mtree.h
mtreeset.cpp
mtreeset.h
ncbitree.cpp
ncbitree.h
node.cpp
node.h
phylokernel.h
phylokernelnew.h
phylokernelnonrev.h
phylonode.cpp
phylonode.h
phylonodemixlen.cpp
phylonodemixlen.h
phylosupertree.cpp
phylosupertree.h
phylosupertreeplen.cpp
phylosupertreeplen.h
phylosupertreeunlinked.cpp
phylosupertreeunlinked.h
phylotree.cpp phylotreethreadingcontext.cpp
phylotree.h   phylotreethreadingcontext.h
phylotreemixlen.cpp
phylotreemixlen.h
phylotreepars.cpp parsimonyjoining.cpp
parsimonysearchparameters.h parsimonysearch.h
parsimonymove.cpp   parsimonymove.h
parsimonynni.cpp    parsimonynni.h
parsimonyspr.cpp    parsimonyspr.h
parsimonytbr.cpp    parsimonytbr.h
parsimonyhybrid.cpp parsimonyhybrid.h
likelihoodspr.cpp
phylotreesse.cpp
placement.cpp
quartet.cpp
supernode.cpp
supernode.h
tinatree.cpp
tinatree.h
parstree.cpp
parstree.h
discordance.cpp
likelihoodbufferset.cpp
likelihoodbufferset.h
alignedalloc.h
nnimove.cpp
nnimove.h
)

target_link_libraries(tree pll model alignment)
