# Create compact output records from completed Totoro fits; never refits models.
root <- normalizePath(".")
in_dir <- file.path(root, "revision_checks", "reviewer18_influential_effects_outputs")
out_dir <- file.path(root, "revision_checks", "cleaned_gaussian_kernel_outputs")
p <- readRDS(file.path(in_dir, "published_cleaned_prepared.rds"))
d <- p$dat
D <- p$distance_km
stopifnot(nrow(d) == 2355L, nlevels(d$study_id) == 390L, nlevels(d$site_id) == 380L,
          identical(rownames(D), levels(d$site_id)), identical(colnames(D), levels(d$site_id)))
unstr <- readRDS(file.path(in_dir, "unstructured_only.rds"))
stopifnot(isTRUE(all.equal(as.numeric(unstr$yi), as.numeric(d$d_Hedges), tolerance = 0)),
          isTRUE(all.equal(as.numeric(unstr$vi), as.numeric(d$var_Hedges), tolerance = 0)), ncol(unstr$X) == 1L)
write.csv(data.frame(n_effects=2355L,n_studies=390L,n_sites=380L,n_const_levels=1L,
  distance_units="WGS84 ellipsoidal geodesic kilometres",distance_rows_match_site_levels=TRUE,
  distance_cols_match_site_levels=TRUE,same_yi_V_X_as_saved_cleaned_SPEXP=TRUE),
  file.path(out_dir,"cleaned_gaussian_input_validation.csv"),row.names=FALSE)
row_for <- function(x,label,model,start=c(NA,NA,NA,NA),elapsed=NA_real_) data.frame(
  model=model,label=label,kernel="SPGAU",n_effects=2355L,n_studies=390L,n_sites=380L,distance_units="km",
  correlation="exp(-d^2/rho^2)",start_effect_variance=start[1],start_study_variance=start[2],
  start_spatial_variance=start[3],start_rho_km=start[4],fixed_spatial_variance=NA_real_,fixed_rho_km=NA_real_,
  pooled_mean=as.numeric(x$b[1]),ci_lb=as.numeric(x$ci.lb[1]),ci_ub=as.numeric(x$ci.ub[1]),
  iid_effect_variance=as.numeric(x$sigma2[1]),study_variance=if(length(x$sigma2)>1)as.numeric(x$sigma2[2]) else NA_real_,
  spatial_variance=as.numeric(x$tau2[1]),rho_km=as.numeric(x$rho[1]),logLik_REML=as.numeric(x$fit.stats["ll","REML"]),
  AIC_REML=as.numeric(x$fit.stats["AIC","REML"]),convergence_status="completed_no_explicit_optimizer_status",
  elapsed_seconds=elapsed,warnings="",stringsAsFactors=FALSE)
# The sequential Totoro completion timestamps give these elapsed times to the
# nearest few seconds; the original R process ended before it could compile them.
sp <- row_for(readRDS(file.path(out_dir,"spatial_only_spgau.rds")),"primary","spatial_only",elapsed=108)
write.csv(sp,file.path(out_dir,"spatial_only_spgau_result.csv"),row.names=FALSE)
starts <- list(short_312_km=c(.751,1.140,.080,312),intermediate_800_km=c(.751,1.160,.060,800),long_3000_km=c(.751,1.210,.040,3000))
elapsed <- c(short_312_km=1209,intermediate_800_km=987,long_3000_km=328)
multi <- do.call(rbind,lapply(names(starts),function(label) row_for(readRDS(file.path(out_dir,paste0("combined_spgau_",label,".rds"))),label,"combined",starts[[label]],elapsed[[label]])))
multi$delta_logLik_from_best <- max(multi$logLik_REML)-multi$logLik_REML
multi$delta_AIC_from_best <- multi$AIC_REML-min(multi$AIC_REML)
write.csv(multi,file.path(out_dir,"combined_spgau_targeted_multistart.csv"),row.names=FALSE)
best <- which.max(multi$logLik_REML)[1]
best_label <- multi$label[best]
writeLines(best_label,file.path(out_dir,"combined_spgau_best_label.txt"))
tau0 <- data.frame(model="combined_tau2_zero",label="tau2_fixed_zero",kernel="SPGAU",n_effects=2355L,n_studies=390L,n_sites=380L,distance_units="km",correlation="exp(-d^2/rho^2)",start_effect_variance=NA_real_,start_study_variance=NA_real_,start_spatial_variance=NA_real_,start_rho_km=NA_real_,fixed_spatial_variance=0,fixed_rho_km=NA_real_,pooled_mean=as.numeric(unstr$b[1]),ci_lb=as.numeric(unstr$ci.lb[1]),ci_ub=as.numeric(unstr$ci.ub[1]),iid_effect_variance=as.numeric(unstr$sigma2[1]),study_variance=as.numeric(unstr$sigma2[2]),spatial_variance=0,rho_km=NA_real_,logLik_REML=as.numeric(unstr$fit.stats["ll","REML"]),AIC_REML=as.numeric(unstr$fit.stats["AIC","REML"]),convergence_status="reused_exact_unstructured_restriction",elapsed_seconds=0,warnings="No refit: tau2=0 makes rho unidentified and yields the saved cleaned unstructured model.")
write.csv(tau0,file.path(out_dir,"combined_spgau_tau2_zero_result.csv"),row.names=FALSE)
write.csv(data.frame(best_combined_label=best_label,combined_logLik_spread=max(multi$logLik_REML)-min(multi$logLik_REML),best_vs_tau2_zero_logLik_loss=multi$logLik_REML[best]-tau0$logLik_REML,tau_zero_vs_saved_unstructured_logLik_difference=0),file.path(out_dir,"cleaned_gaussian_identification_summary.csv"),row.names=FALSE)
v <- .110937235977794
i2 <- rbind(data.frame(model="spatial_only",solution="primary",component=c("iid_effect","spatial","total"),variance=c(sp$iid_effect_variance,sp$spatial_variance,sp$iid_effect_variance+sp$spatial_variance)),data.frame(model="combined",solution=best_label,component=c("iid_effect","study","spatial","total"),variance=c(multi$iid_effect_variance[best],multi$study_variance[best],multi$spatial_variance[best],sum(multi[best,c("iid_effect_variance","study_variance","spatial_variance")]))) )
i2$v_tilde <- v
i2$I2_percent <- NA_real_
for (idx in split(seq_len(nrow(i2)),interaction(i2$model,i2$solution,drop=TRUE))) { denom <- i2$variance[idx[i2$component[idx]=="total"]]+v; i2$I2_percent[idx] <- 100*i2$variance[idx]/denom }
write.csv(i2,file.path(out_dir,"cleaned_gaussian_generalized_i2.csv"),row.names=FALSE)
cat("CLEANED_GAUSSIAN_FINALIZED\n")
