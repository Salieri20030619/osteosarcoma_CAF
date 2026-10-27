P <- '/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF'
Q <- file.path(P,'phase48_clinical_programs_delivery_20261001')
O <- file.path(P,'phase46_spatial_functional_20261001')
readtsv <- function(p) read.delim(p,check.names=FALSE,stringsAsFactors=FALSE)
checkbh <- function(s) {
 q<-rep(NA_real_,nrow(s))
 for(f in unique(s$family)) {ix<-which(s$family==f);v<-s$P[ix];ok<-is.finite(v);v[!ok]<-1;qq<-p.adjust(v,'BH');qq[!ok]<-NA_real_;q[ix]<-qq}
 stopifnot(identical(is.na(q),is.na(s$FDR)));stopifnot(max(c(0,abs(q-s$FDR)),na.rm=TRUE)<1e-12)
}
s<-readtsv(file.path(Q,'tables/clinical_program_statistics.tsv'));r<-readtsv(file.path(Q,'tables/independent_patient_roster.tsv'));p<-readtsv(file.path(Q,'tables/paired_patient_program_deltas.tsv'));models<-readtsv(file.path(Q,'tables/HC3_model_designs.tsv'));keys<-c('layer','dataset','class_id','contrast','endpoint');errors<-list();cache<-list();set.seed(20261001)
matchrow <- function(t,row) {ok<-rep(TRUE,nrow(t));for(k in keys)ok<-ok&t[[k]]==row[[k]];t[ok,,drop=FALSE]}
for(i in seq_len(nrow(s))) {
 row<-s[i,];expected<-NA_real_;tol<-1e-10
 if(row$design=='paired') {
  t<-matchrow(p,row);d<-t$effect;n<-length(d);stopifnot(n==row$N)
  if(n>=2&&any(d!=0)) {key<-as.character(n);if(is.null(cache[[key]]))cache[[key]]<-as.matrix(expand.grid(rep(list(c(-1,1)),n)));expected<-mean(abs(cache[[key]]%*%d/n)>=abs(mean(d))-1e-12)}
 } else if(row$design=='independent') {
  t<-matchrow(r,row);a<-t$value[t$group0_or_1==0];b<-t$value[t$group0_or_1==1];stopifnot(length(a)==row$N0,length(b)==row$N1)
  if(length(a)>=2&&length(b)>=2&&length(unique(c(a,b)))>1) {
   z<-c(a,b);n<-length(a);obs<-abs(mean(b)-mean(a));total<-choose(length(z),n)
   if(total<=200000) {ii<-combn(seq_along(z),n);av<-colMeans(matrix(z[ii],nrow=n));expected<-mean(abs((sum(z)-n*av)/length(b)-av)>=obs-1e-12)} else {
    null<-replicate(19999,{av<-mean(sample(z,n));abs((sum(z)-n*av)/length(b)-av)});expected<-(1+sum(null>=obs-1e-12))/20000;tol<-max(.002,5*sqrt(max(expected*(1-expected),row$P*(1-row$P))/20000*2))
   }
  }
 } else if(is.finite(row$P)) {
  id<-paste(unlist(row[keys]),collapse='|');t<-models[models$model_id==id,];cols<-names(t)[grepl('^X[0-9]+$',names(t))];x<-as.matrix(t[,cols,drop=FALSE]);x<-x[,colSums(is.finite(x))==nrow(x),drop=FALSE];y<-t$value;fit<-lm.fit(x,y);v<-solve(crossprod(x));hat<-rowSums((x%*%v)*x);u<-fit$residuals/(1-hat);cov<-v%*%crossprod(x*u)%*%v;se<-sqrt(cov[2,2]);df<-nrow(x)-ncol(x);expected<-2*pt(-abs(fit$coefficients[2]/se),df=df);stopifnot(abs(fit$coefficients[2]-row$mean_effect)<1e-8,abs(se-row$HC3_SE)<1e-8,df==row$residual_df)
 }
 stopifnot(is.na(expected)==is.na(row$P));err<-if(is.finite(expected))abs(expected-row$P) else 0;stopifnot(err<tol)
 errors[[i]]<-data.frame(row=i,design=row$design,P_error=err,tolerance=tol,passed=TRUE)
}
checkbh(s)
ss<-readtsv(file.path(O,'tables/spatial_model_sensitivity_statistics.tsv'));sp<-readtsv(file.path(O,'tables/spatial_sensitivity_patient_effects.tsv'));stopifnot(all(ss$N==19));mat<-as.matrix(expand.grid(rep(list(c(-1,1)),19)))
for(i in seq_len(nrow(ss))) {d<-sp$effect[sp$program==ss$program[i]&sp$model==ss$model[i]];d<-d[is.finite(d)];expected<-mean(abs(mat%*%d/19)>=abs(mean(d))-1e-12);stopifnot(abs(expected-ss$P[i])<1e-12)}
checkbh(ss);write.table(do.call(rbind,errors),file.path(Q,'tables/independent_R_clinical_program_QA.tsv'),sep='\t',quote=FALSE,row.names=FALSE)

fs<-readtsv(file.path(Q,'tables/D_fraction_HC3_current_statistics.tsv'));fm<-readtsv(file.path(Q,'tables/D_fraction_HC3_valid_model_designs.tsv'))
for(i in which(is.finite(fs$P))) {
 row<-fs[i,];t<-fm[fm$model_id==as.character(i-1),];cols<-names(t)[grepl('^X[0-9]+$',names(t))];x<-as.matrix(t[,cols,drop=FALSE]);x<-x[,colSums(is.finite(x))==nrow(x),drop=FALSE];fit<-lm.fit(x,t$value);v<-solve(crossprod(x));hat<-rowSums((x%*%v)*x);stopifnot(max(hat)<1-1e-12);u<-fit$residuals/(1-hat);cov<-v%*%crossprod(x*u)%*%v;se<-sqrt(cov[2,2]);ts<-fit$coefficients[2]/se;df<-nrow(x)-ncol(x);pp<-2*pnorm(-abs(ts));ptv<-2*pt(-abs(ts),df=df);stopifnot(abs(se-row$SE_HC3)<1e-10,abs(pp-row$P)<1e-10,abs(ptv-row$P_t)<1e-10)
}
checkbh(fs);ft<-fs;ft$P<-ft$P_t;ft$FDR<-ft$FDR_t;checkbh(ft)
writeLines('SUCCESS',file.path(Q,'INDEPENDENT_R_QA_COMPLETE'));writeLines('SUCCESS',file.path(O,'MODEL_SENSITIVITY_REPAIR_QA_COMPLETE'));cat('CLINICAL_SPATIAL_AND_FRACTION_HC3_QA_SUCCESS\n')
