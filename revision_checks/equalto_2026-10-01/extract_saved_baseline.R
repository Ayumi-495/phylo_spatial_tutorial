suppressPackageStartupMessages(library(glmmTMB))
f<-readRDS('revision_checks/visualisation_interval_outputs/moura_glmmTMB_bm.rds')
v<-c('beta:(Intercept)'=fixef(f)$cond[[1]],iid_variance=sigma(f)^2,'variance:study.id'=VarCorr(f)$cond$study.id[1,1],'variance:species.id'=VarCorr(f)$cond$species.id[1,1],'variance:g.1'=VarCorr(f)$cond[['g.1']][1,1],logLik=as.numeric(logLik(f)),AIC=AIC(f))
write.csv(data.frame(parameter=names(v),value=unname(v),extraction_version=as.character(packageVersion('glmmTMB'))),'revision_checks/equalto_2026-10-01/saved_moura_baseline.csv',row.names=FALSE)
