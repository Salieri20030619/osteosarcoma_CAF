p <- "/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase98_reference_thresholds_20261020"
read <- function(n) read.delim(file.path(p,"tables",n),check.names=FALSE,stringsAsFactors=FALSE,fileEncoding="UTF-8")
tum <- read("original_tumor_patient_expression.tsv")
bone <- read("bone_reference_with_coverage_sensitivity.tsv")
lung <- read("lung_reference_pseudobulk.tsv")
res <- read("normal_background_all_results.tsv")
checks <- list()
for(i in seq_len(nrow(res))) {
 r <- res[i,]
 if(is.na(r$P)) next
 a <- subset(tum,source==r$source & gene==r$gene & site=="Lung")$expression
 b <- subset(tum,source==r$source & gene==r$gene & site=="Bone")$expression
 c <- subset(lung,gene==r$gene & specification==r$lung_reference)$expression
 d <- subset(bone,gene==r$gene & specification==r$bone_reference)$expression
 vals <- list(TL=a,TB=b,NL=c[is.finite(c)],NB=d[is.finite(d)])
 dat <- data.frame(y=unlist(vals,use.names=FALSE),g=factor(rep(names(vals),lengths(vals)),levels=names(vals)))
 fit <- lm(y~0+g,data=dat)
 X <- model.matrix(fit); A <- solve(crossprod(X))
 meat <- crossprod(X*(residuals(fit)/sqrt(1-hatvalues(fit))))
 V <- A %*% meat %*% A
 C <- c(1,-1,-1,1)
 effect <- sum(C*coef(fit)); se <- sqrt(drop(t(C)%*%V%*%C))
 v <- sapply(vals,function(z)var(z)/length(z))
 df <- sum(v)^2/sum(v^2/(lengths(vals)-1))
 P <- 2*pt(-abs(effect/se),df)
 checks[[length(checks)+1]] <- data.frame(source=r$source,gene=r$gene,specification=r$specification,effect_error=abs(effect-r$effect),SE_error=abs(se-r$SE),P_error=abs(P-r$P))
}
stopifnot(length(checks)>0,length(checks)==sum(!is.na(res$P)))
z <- do.call(rbind,checks)
write.table(z,file.path(p,"tables/independent_R_HC2_checks.tsv"),sep="\t",fileEncoding="UTF-8",row.names=FALSE,quote=FALSE)
stopifnot(max(z$effect_error)<1e-10,max(z$SE_error)<1e-10,max(z$P_error)<1e-10)
print(as.data.frame(sapply(z[,c("effect_error","SE_error","P_error")],max)))
cat("INDEPENDENT_R_HC2_PASS",nrow(z),"\n")
