add_library(lp colamd/colamd.c shared/commonlib.c ini.c lp_Hash.c bfp/bfp_LUSOL/lp_LUSOL.c
        lp_MDO.c lp_MPS.c lp_SOS.c lp_crash.c lp_lib.c lp_matrix.c lp_mipbb.c 
        lp_params.c lp_presolve.c lp_price.c lp_pricePSE.c lp_report.c lp_rlp.c lp_scale.c 
        lp_simplex.c lp_utils.c lp_wlp.c bfp/bfp_LUSOL/LUSOL/lusol.c shared/mmio.c shared/myblas.c yacc_read.c)

#add_subdirectory(shared)
#add_subdirectory(bfp)
#add_subdirectory(colamd)

include_directories("${PROJECT_SOURCE_DIR}/lp/shared"
"${PROJECT_SOURCE_DIR}/lp/bfp"
"${PROJECT_SOURCE_DIR}/lp/bfp/bfp_LUSOL"
"${PROJECT_SOURCE_DIR}/lp/bfp/bfp_LUSOL/LUSOL"
"${PROJECT_SOURCE_DIR}/lp/colamd"
"${PROJECT_SOURCE_DIR}/lp"
)

add_definitions(-DINVERSE_ACTIVE=INVERSE_LUSOL -DPARSER_LP 
        -DRoleIsExternalInvEngine -DYY_NEVER_INTERACTIVE)
