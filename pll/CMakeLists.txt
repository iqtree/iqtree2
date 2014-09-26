if (IQTREE_FLAGS MATCHES "avx" AND IQTREE_FLAGS MATCHES "omp")
  message("PLL Kernel    : AVX/Pthreads")
#  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D__AVX -mavx")
  add_library(pll 
  avxLikelihood.c
  alignment.c
  bipartitionList.c
  evaluateGenericSpecial.c
  evaluatePartialGenericSpecial.c
  fastDNAparsimony.c
  hardware.c
  hash.c
  lexer.c
  makenewzGenericSpecial.c
  models.c
  newick.c
  newviewGenericSpecial.c
  genericParallelization.c
  optimizeModel.c
  parsePartition.c
  queue.c
  randomTree.c
  recom.c
  restartHashTable.c
  searchAlgo.c
  ssort.c
  stack.c
  topologies.c
  trash.c
  treeIO.c
  utils.c)
elseif (IQTREE_FLAGS MATCHES "omp")
  message("PLL Kernel    : SSE3/Pthreads")
#  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse3 -D__SSE3")
  add_library(pll 
  alignment.c
  bipartitionList.c
  evaluateGenericSpecial.c
  evaluatePartialGenericSpecial.c
  fastDNAparsimony.c
  hardware.c
  hash.c
  lexer.c
  makenewzGenericSpecial.c
  models.c
  newick.c
  genericParallelization.c
  newviewGenericSpecial.c
  optimizeModel.c
  parsePartition.c
  queue.c
  randomTree.c
  recom.c
  restartHashTable.c
  searchAlgo.c
  ssort.c
  stack.c
  topologies.c
  trash.c
  treeIO.c
  utils.c)
elseif (IQTREE_FLAGS MATCHES "avx")
  message("PLL Kernel    : AVX")  
#  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D__AVX -mavx")
  add_library(pll 
    avxLikelihood.c
    alignment.c
    bipartitionList.c
    evaluateGenericSpecial.c
    evaluatePartialGenericSpecial.c
    fastDNAparsimony.c
    hardware.c
    hash.c
    lexer.c
    makenewzGenericSpecial.c
    models.c
    newick.c
    newviewGenericSpecial.c
    optimizeModel.c
    parsePartition.c
    queue.c
    randomTree.c
    recom.c
    restartHashTable.c
    searchAlgo.c
    ssort.c
    stack.c
    topologies.c
    trash.c
    treeIO.c
    utils.c)
else ()
  message("PLL Kernel    : SSE3")
#  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse3 -D__SSE3")
  add_library(pll 
    alignment.c
    bipartitionList.c
    evaluateGenericSpecial.c
    evaluatePartialGenericSpecial.c
    fastDNAparsimony.c
    hardware.c
    hash.c
    lexer.c
    makenewzGenericSpecial.c
    models.c
    newick.c
    newviewGenericSpecial.c
    optimizeModel.c
    parsePartition.c
    queue.c
    randomTree.c
    recom.c
    restartHashTable.c
    searchAlgo.c
    ssort.c
    stack.c
    topologies.c
    trash.c
    treeIO.c
    utils.c)       
endif()
