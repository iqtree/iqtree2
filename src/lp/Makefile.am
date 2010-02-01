INCLUDES = -I$(top_srcdir)/src/lp/shared -I$(top_srcdir)/src/lp/bfp \
	-I$(top_srcdir)/src/lp/bfp/bfp_LUSOL -I$(top_srcdir)/src/lp/bfp/bfp_LUSOL/LUSOL \
	-I$(top_srcdir)/src/lp/colamd $(all_includes)
noinst_LTLIBRARIES = liblpsolve.la
liblpsolve_la_SOURCES = colamd.c commonlib.c ini.c lp_Hash.c lp_LUSOL.c \
	lp_MDO.c lp_MPS.c lp_SOS.c lp_crash.c lp_lib.c lp_matrix.c lp_mipbb.c \
	lp_params.c lp_presolve.c lp_price.c lp_pricePSE.c lp_report.c lp_rlp.c lp_scale.c \
	lp_simplex.c lp_utils.c lp_wlp.c lusol.c mmio.c myblas.c yacc_read.c
AM_CFLAGS = -DINVERSE_ACTIVE=INVERSE_LUSOL -DPARSER_LP \
	-DRoleIsExternalInvEngine -DYY_NEVER_INTERACTIVE
liblpsolve_la_LDFLAGS = -ldl
SUBDIRS = shared bfp colamd
noinst_HEADERS = declare.h fortify.h ini.h lp_Hash.h lp_MDO.h lp_MPS.h lp_SOS.h \
	lp_crash.h lp_explicit.h lp_fortify.h lp_lib.h lp_matrix.h lp_mipbb.h \
	lp_presolve.h lp_price.h lp_pricePSE.h lp_report.h lp_rlp.h lp_scale.h lp_simplex.h \
	lp_solveDLL.h lp_types.h lp_utils.h lp_wlp.h lpkit.h lpsolve.h ufortify.h yacc_read.h
