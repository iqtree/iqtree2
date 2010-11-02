bin_PROGRAMS = pda 

# set the include path found by configure
INCLUDES = $(all_includes)

# the library search path.
pda_LDFLAGS = $(all_libraries) -ldl
pda_SOURCES = alignment.cpp alignmentpairwise.cpp circularnetwork.cpp \
	eigendecomposition.cpp greedy.cpp gtrmodel.cpp gurobiwrapper.cpp hashsplitset.cpp iqptree.cpp \
	lpwrapper.c matree.cpp mexttree.cpp modeldna.cpp modelfactory.cpp modelprotein.cpp \
	modeltest_wrapper.c mpdablock.cpp msetsblock.cpp msplitsblock.cpp mtree.cpp mtreeset.cpp \
	node.cpp optimization.cpp pattern.cpp pda.cpp pdnetwork.cpp pdtree.cpp \
	pdtreeset.cpp phyloanalysis.cpp phylonode.cpp phylotree.cpp pruning.cpp rategamma.cpp \
	rategammainvar.cpp rateheterogeneity.cpp rateinvar.cpp ratemeyerhaeseler.cpp split.cpp \
	splitgraph.cpp splitset.cpp stoprule.cpp substmodel.cpp tools.cpp
noinst_HEADERS = alignment.h matree.h alignmentpairwise.h bionj.h circularnetwork.h \
	eigendecomposition.h greedy.h gtrmodel.h gurobiwrapper.h hashsplitset.h iqptree.h \
	lpwrapper.h mexttree.h modeldna.h modelfactory.h modelprotein.h modeltest_wrapper.h \
	mpdablock.h msetsblock.h msplitsblock.h mtree.h mtreeset.h myreader.h node.h \
	optimization.h pattern.h pdnetwork.h pdtree.h pdtreeset.h phyloanalysis.h phylonode.h \
	phylotree.h pruning.h rategamma.h rategammainvar.h rateheterogeneity.h rateinvar.h \
	ratemeyerhaeseler.h split.h splitgraph.h splitset.h stoprule.h substmodel.h tools.h
SUBDIRS = lp modeltest ncl
pda_LDADD = $(top_builddir)/src/lp/liblpsolve.la \
	$(top_builddir)/src/modeltest/libmodeltest.la $(top_builddir)/src/ncl/libNCL.la $(top_srcdir)/lib/libacml_mv.a
AM_CXXFLAGS = -Wno-deprecated
