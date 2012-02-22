bin_PROGRAMS = iqtree

# set the include path found by configure
INCLUDES = $(all_includes)

# the library search path.
iqtree_LDFLAGS = $(all_libraries) -ldl -static
iqtree_SOURCES = alignment.cpp alignmentpairwise.cpp circularnetwork.cpp \
	eigendecomposition.cpp greedy.cpp gss.cpp gtrmodel.cpp guidedbootstrap.cpp gurobiwrapper.cpp \
	hashsplitset.cpp iqptree.cpp lpwrapper.c maalignment.cpp matree.cpp mexttree.cpp \
	modelbin.cpp modeldna.cpp modelfactory.cpp modelnonrev.cpp modelprotein.cpp \
	modeltest_wrapper.c mpdablock.cpp msetsblock.cpp msplitsblock.cpp mtree.cpp mtreeset.cpp \
	ngs.cpp node.cpp optimization.cpp parsmultistate.cpp partitionmodel.cpp \
	pattern.cpp pda.cpp pdnetwork.cpp pdtree.cpp pdtreeset.cpp phyloanalysis.cpp \
	phylonode.cpp phylosupertree.cpp phylotree.cpp phylotreesse.cpp pruning.cpp \
	rategamma.cpp rategammainvar.cpp rateheterogeneity.cpp rateinvar.cpp \
	ratemeyerdiscrete.cpp ratemeyerhaeseler.cpp split.cpp splitgraph.cpp splitset.cpp \
	stoprule.cpp substmodel.cpp superalignment.cpp superalignmentpairwise.cpp \
	supernode.cpp tinatree.cpp tools.cpp whtest_wrapper.cpp
noinst_HEADERS = alignment.h alignmentpairwise.h bionj.h circularnetwork.h \
	eigendecomposition.h greedy.h gss.h gtrmodel.h gurobiwrapper.h hashsplitset.h iqptree.h \
	lpwrapper.h maalignment.h matree.h mexttree.h modelbin.h modeldna.h modelfactory.h modelnonrev.h \
	modelprotein.h modeltest_wrapper.h mpdablock.h msetsblock.h msplitsblock.h mtree.h \
	mtreeset.h myreader.h ngs.h node.h optimization.h parsmultistate.h \
	partitionmodel.h pattern.h pdnetwork.h pdtree.h pdtreeset.h phyloanalysis.h phylonode.h \
	phylosupertree.h phylotree.h pruning.h rategamma.h rategammainvar.h rateheterogeneity.h \
	rateinvar.h ratemeyerdiscrete.h ratemeyerhaeseler.h split.h splitgraph.h splitset.h \
	stoprule.h substmodel.h superalignment.h superalignmentpairwise.h supernode.h \
	tinatree.h tools.h whtest_wrapper.h
SUBDIRS = kmeans lp modeltest ncl whtest
iqtree_LDADD = $(top_builddir)/src/kmeans/libkmeans.a \
	$(top_builddir)/src/lp/liblpsolve.a $(top_builddir)/src/modeltest/libmodeltest.a \
	$(top_builddir)/src/ncl/libNCL.a $(top_builddir)/src/whtest/libwhtest.a
AM_CXXFLAGS = -Wno-deprecated -DIQ_TREE -O3
AM_CFLAGS = -DIQ_TREE -O3

EXTRA_DIST = Eigen/Core Eigen/Dense Eigen/Eigen \
Eigen/src/Core/ArrayBase.h \
Eigen/src/Core/Array.h \
Eigen/src/Core/ArrayWrapper.h \
Eigen/src/Core/Assign.h \
Eigen/src/Core/BandMatrix.h \
Eigen/src/Core/Block.h \
Eigen/src/Core/BooleanRedux.h \
Eigen/src/Core/CommaInitializer.h \
Eigen/src/Core/CwiseBinaryOp.h \
Eigen/src/Core/CwiseNullaryOp.h \
Eigen/src/Core/CwiseUnaryOp.h \
Eigen/src/Core/CwiseUnaryView.h \
Eigen/src/Core/DenseBase.h \
Eigen/src/Core/DenseCoeffsBase.h \
Eigen/src/Core/DenseStorageBase.h \
Eigen/src/Core/Diagonal.h \
Eigen/src/Core/DiagonalMatrix.h \
Eigen/src/Core/DiagonalProduct.h \
Eigen/src/Core/Dot.h \
Eigen/src/Core/EigenBase.h \
Eigen/src/Core/Flagged.h \
Eigen/src/Core/ForceAlignedAccess.h \
Eigen/src/Core/Functors.h \
Eigen/src/Core/Fuzzy.h \
Eigen/src/Core/GenericPacketMath.h \
Eigen/src/Core/GlobalFunctions.h \
Eigen/src/Core/IO.h \
Eigen/src/Core/MapBase.h \
Eigen/src/Core/Map.h \
Eigen/src/Core/MathFunctions.h \
Eigen/src/Core/MatrixBase.h \
Eigen/src/Core/Matrix.h \
Eigen/src/Core/MatrixStorage.h \
Eigen/src/Core/NestByValue.h \
Eigen/src/Core/NoAlias.h \
Eigen/src/Core/NumTraits.h \
Eigen/src/Core/PermutationMatrix.h \
Eigen/src/Core/ProductBase.h \
Eigen/src/Core/Product.h \
Eigen/src/Core/Random.h \
Eigen/src/Core/Redux.h \
Eigen/src/Core/Replicate.h \
Eigen/src/Core/ReturnByValue.h \
Eigen/src/Core/Reverse.h \
Eigen/src/Core/Select.h \
Eigen/src/Core/SelfAdjointView.h \
Eigen/src/Core/SelfCwiseBinaryOp.h \
Eigen/src/Core/SolveTriangular.h \
Eigen/src/Core/StableNorm.h \
Eigen/src/Core/Stride.h \
Eigen/src/Core/Swap.h \
Eigen/src/Core/Transpose.h \
Eigen/src/Core/Transpositions.h \
Eigen/src/Core/TriangularMatrix.h \
Eigen/src/Core/VectorBlock.h \
Eigen/src/Core/VectorwiseOp.h \
Eigen/src/Core/Visitor.h \
Eigen/src/Core/arch/AltiVec/Complex.h \
Eigen/src/Core/arch/AltiVec/PacketMath.h \
Eigen/src/Core/arch/Default/Settings.h \
Eigen/src/Core/arch/NEON/Complex.h \
Eigen/src/Core/arch/NEON/PacketMath.h \
Eigen/src/Core/arch/SSE/Complex.h \
Eigen/src/Core/arch/SSE/MathFunctions.h \
Eigen/src/Core/arch/SSE/PacketMath.h \
Eigen/src/Core/products/CoeffBasedProduct.h \
Eigen/src/Core/products/GeneralBlockPanelKernel.h \
Eigen/src/Core/products/GeneralMatrixMatrix.h \
Eigen/src/Core/products/GeneralMatrixVector.h \
Eigen/src/Core/products/Parallelizer.h \
Eigen/src/Core/products/SelfadjointMatrixMatrix.h \
Eigen/src/Core/products/SelfadjointMatrixVector.h \
Eigen/src/Core/products/SelfadjointProduct.h \
Eigen/src/Core/products/SelfadjointRank2Update.h \
Eigen/src/Core/products/TriangularMatrixMatrix.h \
Eigen/src/Core/products/TriangularMatrixVector.h \
Eigen/src/Core/products/TriangularSolverMatrix.h \
Eigen/src/Core/util/BlasUtil.h \
Eigen/src/Core/util/Constants.h \
Eigen/src/Core/util/DisableMSVCWarnings.h \
Eigen/src/Core/util/EnableMSVCWarnings.h \
Eigen/src/Core/util/ForwardDeclarations.h \
Eigen/src/Core/util/Macros.h \
Eigen/src/Core/util/Memory.h \
Eigen/src/Core/util/Meta.h \
Eigen/src/Core/util/StaticAssert.h \
Eigen/src/Core/util/XprHelper.h \
Eigen/src/plugins/ArrayCwiseBinaryOps.h \
Eigen/src/plugins/ArrayCwiseUnaryOps.h \
Eigen/src/plugins/BlockMethods.h \
Eigen/src/plugins/CommonCwiseBinaryOps.h \
Eigen/src/plugins/CommonCwiseUnaryOps.h \
Eigen/src/plugins/MatrixCwiseBinaryOps.h \
Eigen/src/plugins/MatrixCwiseUnaryOps.h 
