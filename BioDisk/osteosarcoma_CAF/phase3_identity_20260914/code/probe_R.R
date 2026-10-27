# SYNTAX_CHECK: Rscript --vanilla -e 'parse(file="probe_R.R")'
# RUN_COMMAND: Rscript --vanilla probe_R.R
# EXPECTED_SUCCESS_MARKER: R_PROBE_COMPLETE
cat(R.version.string,"\n")
print(.libPaths())
pk <- c('infercnv','copykat','Matrix','data.table','jsonlite','dlm','parallelDist','mixtools','MCMCpack','KernSmooth','RColorBrewer','gplots','foreach','doParallel')
for(p in pk) {
  ok <- requireNamespace(p,quietly=TRUE)
  cat(p,ok,if(ok) as.character(packageVersion(p)) else 'ABSENT','\n')
}
cat('R_PROBE_COMPLETE\n')
