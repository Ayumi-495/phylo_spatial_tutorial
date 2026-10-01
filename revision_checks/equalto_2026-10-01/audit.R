#!/usr/bin/env Rscript
# Execute the current and preserved tutorial blocks in separate R sessions.
args <- commandArgs(TRUE)
mode <- args[[1]]
if (mode == 'new') .libPaths(c('/private/tmp/ayumi-equalto-library', .libPaths()))
suppressPackageStartupMessages({library(glmmTMB); library(metadat); library(metafor); library(ape); library(gtools); library(dplyr); library(sf); library(geosphere); library(here)})
if(mode == 'new') stopifnot(packageVersion('glmmTMB') == '1.1.15.2')
out <- 'revision_checks/equalto_2026-10-01'
if(length(args)>1 && args[[2]]=='verify') {
 r <- read.csv(file.path(out,paste0(mode,'_results.csv')))
 stopifnot(length(unique(r$model))==5L, all(r$version=='1.1.15.2'), all(r$pdHess), all(r$convergence==0L))
 for(n in unique(r$model)) {f<-readRDS(file.path(out,paste0(mode,'_',n,'.rds')));stopifnot(isTRUE(f$sdr$pdHess), f$fit$convergence==0L, is.finite(as.numeric(logLik(f))))}
 writeLines(character(),file.path(out,'new_errors.txt'))
 message('FIVE_CURRENT_MODELS_RERUN_PASSED');quit(status=0)
}
src <- readLines(if(mode == 'new') 'tutorial_v2.qmd' else file.path(out,'tutorial_before.qmd'))
starts <- grep('^```\\{r',src); ends <- vapply(starts,function(i) i+which(src[(i+1):length(src)]=='```')[[1]],integer(1))
blocks <- mapply(function(a,b) paste(src[(a+1):(b-1)],collapse='\n'),starts,ends)
run_block <- function(code,env) {
 for(expr in parse(text=code)) {
  # Run preparation and model fit; avoid printing summaries/confint and plotting.
  txt <- paste(deparse(expr),collapse=' ')
  if(grepl('^(summary|head|confint|class|sigma|VarCorr)\\(',txt)) next
  eval(expr, env)
 }
}
rows <- list(); errors <- list()
record <- function(fit,name,warnings) {
 saveRDS(fit,file.path(out,paste0(mode,'_',name,'.rds')))
 beta <- fixef(fit)$cond; vc <- VarCorr(fit)$cond
 values <- c(setNames(as.numeric(beta),paste0('beta:',names(beta))), iid_variance=sigma(fit)^2,
             setNames(vapply(vc,function(x) as.numeric(x[1,1]),numeric(1)),ifelse(names(vc) %in% c('g','const'),'known_sampling_variance_first',paste0('variance:',names(vc)))),
             logLik=as.numeric(logLik(fit)), AIC=AIC(fit))
 if(name=='spain') {theta<-fit$fit$par[names(fit$fit$par)=='theta']; stopifnot(length(theta)==2); values<-c(values,rho_km=unname(exp(theta[2])))}
 rows[[name]] <<- data.frame(model=name,parameter=names(values),value=unname(values),version=as.character(packageVersion('glmmTMB')),pdHess=isTRUE(fit$sdr$pdHess),convergence=fit$fit$convergence,warnings=paste(unique(warnings),collapse=' | '))
 write.csv(do.call(rbind,rows),file.path(out,paste0(mode,'_results.csv')),row.names=FALSE)
}
fit_case <- function(name,env,codes,fitname) {
 warnings <- character(); message('START ',mode,' ',name)
 tryCatch({withCallingHandlers({for(code in codes)run_block(code,env)},warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}); record(get(fitname,env),name,warnings);message('DONE ',mode,' ',name)},error=function(e){errors[[name]]<<-conditionMessage(e);message('ERROR ',name,': ',conditionMessage(e))})
}
e<-new.env(); e$dat_moura2021<-metafor::escalc(measure='ZCOR',ri=ri,ni=ni,data=metadat::dat.moura2021$dat)
e$dat_moura2021$effect.size.id<-factor(seq_len(nrow(e$dat_moura2021)));e$dat_moura2021$species.id.phy<-e$dat_moura2021$species.id
e$A<-ape::vcv(ape::compute.brlen(metadat::dat.moura2021$tree),corr=TRUE)
e$dat_moura2021$species.id.phy<-factor(e$dat_moura2021$species.id.phy,levels=rownames(e$A))
ma <- blocks[grepl('phylo_eg1_tmb <- glmmTMB',blocks,fixed=TRUE)];mr<-blocks[grepl('phylo_eg1.1_tmb <- glmmTMB',blocks,fixed=TRUE)]
prep <- blocks[grepl('VCV <- diag',blocks,fixed=TRUE)&!grepl('glmmTMB(',blocks,fixed=TRUE)]
fit_case('moura_ma',e,c('dat_moura2021$g <- factor("all")',prep,ma),'phylo_eg1_tmb')
fit_case('moura_mr',e,mr,'phylo_eg1.1_tmb')
e<-new.env();e$dat_lim2014<-metafor::escalc(measure='ZCOR',ri=ri,ni=ni,data=metadat::dat.lim2014$o_o_unadj)
e$dat_lim2014$id<-seq_len(nrow(e$dat_lim2014));e$dat_lim2014$phy<-e$dat_lim2014$species
e$A<-ape::vcv(ape::compute.brlen(metadat::dat.lim2014$o_o_unadj_tree),corr=TRUE)
fit_case('lim_ma',e,blocks[grepl('phylo_eg2_tmb_ma <- glmmTMB',blocks,fixed=TRUE)],'phylo_eg2_tmb_ma')
fit_case('lim_mr',e,blocks[grepl('phylo_eg2_tmb_mr <- glmmTMB',blocks,fixed=TRUE)],'phylo_eg2_tmb_mr')
e<-new.env()
fit_case('spain',e,c(blocks[grepl('dat_spain <- read.csv',blocks,fixed=TRUE)],blocks[grepl('spain_glmmTMB <- glmmTMB::glmmTMB',blocks,fixed=TRUE)]),'spain_glmmTMB')
writeLines(c(paste('Mode',mode),capture.output(sessionInfo())),file.path(out,paste0(mode,'_sessionInfo.txt')))
writeLines(as.character(unlist(lapply(names(errors),function(n)paste(n,errors[[n]],sep=': ')))),file.path(out,paste0(mode,'_errors.txt')))
if(mode=='new') {stopifnot(length(rows)==5L,length(errors)==0L);message('FIVE_CURRENT_MODELS_RERUN_PASSED')}
