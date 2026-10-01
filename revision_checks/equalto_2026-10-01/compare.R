.libPaths(c('/private/tmp/ayumi-equalto-library',.libPaths()))
suppressPackageStartupMessages(library(glmmTMB))
out <- 'revision_checks/equalto_2026-10-01'
n <- read.csv(file.path(out,'new_results.csv'));stopifnot(all(n$version=='1.1.15.2'),all(n$pdHess),all(n$convergence==0L));o<-read.csv(file.path(out,'old_results.csv'))
saved <- read.csv(file.path(out,'saved_moura_baseline.csv'))
v <- setNames(saved$value,saved$parameter)
b <- rbind(data.frame(model=o$model,parameter=o$parameter,baseline=o$value,provenance='preserved tutorial rerun with glmmTMB 1.1.15',tolerance=1e-5),data.frame(model='moura_ma',parameter=names(v),baseline=unname(v),provenance='pre-existing saved glmmTMB fit',tolerance=1e-5))
# Moura MR has no saved glmmTMB baseline. Compare to its displayed rounded values,
# retaining different precision for variance obtained by squaring printed SD.
v<-c('beta:(Intercept)'=.35618152,'beta:temporally.pooledyes'=.03953359,iid_variance=.12036710^2,'variance:study.id'=.13930061^2,'variance:species.id.phy'=.23227256^2,'variance:g.1'=.05198478)
b<-rbind(b,data.frame(model='moura_mr',parameter=names(v),baseline=unname(v),provenance='preserved tutorial printed estimates (rounded)',tolerance=5e-7))
c<-merge(b,n[c('model','parameter','value')]);c$delta<-c$value-c$baseline;c$within_tolerance<-abs(c$delta)<=c$tolerance
write.csv(c,file.path(out,'comparison.csv'),row.names=FALSE)
stopifnot(all(c$within_tolerance),length(unique(c$model))==5)
cat('FIVE_MODELS_BASELINE_COMPARISON_PASSED\n')
print(aggregate(abs(delta)~model,c,max))
