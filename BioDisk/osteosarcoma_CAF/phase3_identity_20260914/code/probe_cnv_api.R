# SYNTAX_CHECK: Rscript --vanilla -e 'invisible(parse(file=commandArgs(TRUE)[1]))' probe_cnv_api.R
# RUN_COMMAND: Rscript --vanilla probe_cnv_api.R
# EXPECTED_SUCCESS_MARKER: CNV_API_CHECKED
root <- '/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase3_identity_20260914'
.libPaths(c(file.path(root,'Rlib'),.libPaths(),'/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/lib/R/library'))
print(formals(infercnv::run))
if(requireNamespace('copykat',quietly=TRUE)) {
 print(formals(copykat::copykat))
 writeLines(deparse(copykat::copykat),file.path(root,'resources','copykat_installed_function.R'))
}
cat('CNV_API_CHECKED\n')
