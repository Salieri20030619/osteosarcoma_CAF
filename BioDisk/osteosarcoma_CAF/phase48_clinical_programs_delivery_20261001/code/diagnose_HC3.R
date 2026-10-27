Q <- '/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase48_clinical_programs_delivery_20261001'
s<-read.delim(file.path(Q,'tables/clinical_program_statistics.tsv'));m<-read.delim(file.path(Q,'tables/HC3_model_designs.tsv'));keys<-c('layer','dataset','class_id','contrast','endpoint');out<-list()
for(i in which(s$design=='adjusted_HC3_t'&is.finite(s$P))) {
 row<-s[i,];id<-paste(unlist(row[keys]),collapse='|');t<-m[m$model_id==id,];cols<-names(t)[grepl('^X[0-9]+$',names(t))];x<-as.matrix(t[,cols,drop=FALSE]);x<-x[,colSums(is.finite(x))==nrow(x),drop=FALSE];fit<-lm.fit(x,t$value);v<-solve(crossprod(x));hat<-rowSums((x%*%v)*x);u<-fit$residuals/(1-hat);cov<-v%*%crossprod(x*u)%*%v;se<-sqrt(cov[2,2]);out[[length(out)+1]]<-data.frame(row=i,model_id=id,max_hat=max(hat),SE_Python=row$HC3_SE,SE_R=se,SE_error=abs(se-row$HC3_SE))
}
write.table(do.call(rbind,out),file.path(Q,'tables/HC3_numerical_diagnostic.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
