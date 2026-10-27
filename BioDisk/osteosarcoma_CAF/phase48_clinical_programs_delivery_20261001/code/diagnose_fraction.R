Q <- '/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase48_clinical_programs_delivery_20261001'
s<-read.delim(file.path(Q,'tables/D_fraction_HC3_current_statistics.tsv'));m<-read.delim(file.path(Q,'tables/D_fraction_HC3_valid_model_designs.tsv'));out<-list()
for(i in which(is.finite(s$P))) {
 r<-s[i,];t<-m[m$model_id==as.character(i-1),];cols<-names(t)[grepl('^X[0-9]+$',names(t))];x<-as.matrix(t[,cols,drop=FALSE]);x<-x[,colSums(is.finite(x))==nrow(x),drop=FALSE];fit<-lm.fit(x,t$value);v<-solve(crossprod(x));hat<-rowSums((x%*%v)*x);u<-fit$residuals/(1-hat);cov<-v%*%crossprod(x*u)%*%v;se<-sqrt(cov[2,2]);ts<-fit$coefficients[2]/se;pp<-2*pnorm(-abs(ts));pv<-2*pt(-abs(ts),nrow(x)-ncol(x));out[[length(out)+1]]<-data.frame(row=i,feature=r$feature,den=r$denominator,PythonCoef=r$target_effect,Rcoef=fit$coefficients[2],PythonSE=r$SE_HC3,Rse=se,PythonP=r$P,RP=pp,PythonPt=r$P_t,Rpt=pv,P_error=abs(pp-r$P),maxhat=max(hat))
}
write.table(do.call(rbind,out),file.path(Q,'tables/fraction_HC3_R_numeric_diagnostic.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
