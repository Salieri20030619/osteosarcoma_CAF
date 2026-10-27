args_all<-commandArgs(FALSE)
f<-Sys.getenv('CAF_FROZEN_R_SCRIPT');stopifnot(nzchar(f));parse(file=f);source(f,echo=FALSE)
