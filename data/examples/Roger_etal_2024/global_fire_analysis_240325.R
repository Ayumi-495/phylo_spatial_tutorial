rm(list=ls(all=TRUE))

library(metafor)

reg <- read.csv("global_fire_DB_240325.csv")


# Preliminary analyses - outliers ----

# main
reg.ab <- reg[reg$response=="abundance",]
reg.di <- reg[reg$response=="diversity",]
reg.fit <- reg[reg$response=="fitness",]

myrma.ab <- rma.mv(d_Hedges, var_Hedges,
                   random = list(~1|study_num/ES_num),
                   method="ML",
                   data = reg.ab)
myrma.di <- rma.mv(d_Hedges, var_Hedges,
                   random = list(~1|study_num/ES_num),
                   method="ML",
                   data = reg.di)
myrma.fi <- rma.mv(d_Hedges, var_Hedges,
                   random = list(~1|study_num/ES_num),
                   method="ML",
                   data = reg.fi)

# Cooks distance on complete models
# ! Takes a few hours to compute
# cooksd.ab <- cooks.distance.rma.mv(myrma.ab,
#                                    progbar = TRUE, parallel = "multicore",
#                                    ncpus = 4, reestimate = FALSE)
# plot(cooksd.ab)
# reg.ab[cooksd.ab>0.015,] # Ngugi_2022-2, Gagnon_2015-2, Moris_2017-1

reg.ab$d_Hedges[reg.ab$study_id=="Ngugi_2022"&reg.ab$ES_num==2] <- NA
reg.ab$d_Hedges[reg.ab$study_id=="Gagnon_2015"&reg.ab$ES_num==2] <- NA
reg.ab$d_Hedges[reg.ab$study_id=="Moris_2017"&reg.ab$ES_num==1] <- NA

# cooksd.di <- cooks.distance.rma.mv(myrma.di,
#                                    progbar = TRUE, parallel = "multicore",
#                                    ncpus = 4, reestimate = FALSE)
# plot(cooksd.di)
# reg.di[cooksd.di>0.015,] #Schwilk_1997-1, Silveira_2016-4

reg.di$d_Hedges[reg.di$study_id=="Schwilk_1997"&reg.di$ES_num==1] <- NA
reg.di$d_Hedges[reg.di$study_id=="Silveira_2016"&reg.di$ES_num==4] <- NA

# cooksd.fi <- cooks.distance.rma.mv(myrma.fi,
#                                    progbar = TRUE, parallel = "multicore",
#                                    ncpus = 4, reestimate = FALSE)
# plot(cooksd.fi)
# reg.fi[cooksd.fi>0.06,] #Launonen_1999-1, Ansley_2015-1

reg.fit$d_Hedges[reg.fit$study_id=="Launonen_1999"&reg.fit$ES_num==1] <- NA
reg.fit$d_Hedges[reg.fit$study_id=="Ansley_2015"&reg.fit$ES_num==1] <- NA


# Datasets ----

# Main

reg.ab <- reg.ab[!is.na(reg.ab$d_Hedges),]
reg.di <- reg.di[!is.na(reg.di$d_Hedges),]
reg.fit <- reg.fit[!is.na(reg.fit$d_Hedges),]


# Subsetted (for testing moderator categories against other moderator categories)

# abundance

reg.ab.fi1 <- reg.ab[!is.na(reg.ab$fire_type)&reg.ab$fire_type=="prescribed",]
reg.ab.fi2 <- reg.ab[!is.na(reg.ab$fire_type)&reg.ab$fire_type=="wildfire",]

reg.ab.tt1 <- reg.ab[reg.ab$FR_comp=="frequency",]
reg.ab.tt2 <- reg.ab[reg.ab$FR_comp=="severity",]

reg.ab.ft1 <- reg.ab[!is.na(reg.ab$TSF_factor)&reg.ab$TSF_factor=="short",]
reg.ab.ft2 <- reg.ab[!is.na(reg.ab$TSF_factor)&reg.ab$TSF_factor=="long",]

reg.ab.rh1 <- reg.ab[!is.na(reg.ab$FR_hist)&reg.ab$FR_hist=="no-fire",]
reg.ab.rh2 <- reg.ab[!is.na(reg.ab$FR_hist)&reg.ab$FR_hist=="surface",]
reg.ab.rh3 <- reg.ab[!is.na(reg.ab$FR_hist)&reg.ab$FR_hist=="crown",]

reg.ab.ha1 <- reg.ab[!is.na(reg.ab$habitat)&reg.ab$habitat=="forest-broadleaf",]
reg.ab.ha2 <- reg.ab[!is.na(reg.ab$habitat)&reg.ab$habitat=="forest-conifer",]
reg.ab.ha3 <- reg.ab[!is.na(reg.ab$habitat)&reg.ab$habitat=="forest-mixed",]
reg.ab.ha4 <- reg.ab[!is.na(reg.ab$habitat)&reg.ab$habitat=="grassland",]
reg.ab.ha5 <- reg.ab[!is.na(reg.ab$habitat)&reg.ab$habitat=="shrubland",]
reg.ab.ha7 <- reg.ab[!is.na(reg.ab$habitat)&reg.ab$habitat=="woodland",]

reg.ab.pfg1 <- reg.ab[!is.na(reg.ab$PLF)&reg.ab$PLF=="bryophyte",]
reg.ab.pfg2 <- reg.ab[!is.na(reg.ab$PLF)&reg.ab$PLF=="herb",]
reg.ab.pfg3 <- reg.ab[!is.na(reg.ab$PLF)&reg.ab$PLF=="woody",]

reg.ab.cl1 <- reg.ab[!is.na(reg.ab$climate)&reg.ab$climate=="arid",]
reg.ab.cl2 <- reg.ab[!is.na(reg.ab$climate)&reg.ab$climate=="cold",]
reg.ab.cl3 <- reg.ab[!is.na(reg.ab$climate)&reg.ab$climate=="temperate-dry",]
reg.ab.cl4 <- reg.ab[!is.na(reg.ab$climate)&reg.ab$climate=="temperate-nodry",]
reg.ab.cl5 <- reg.ab[!is.na(reg.ab$climate)&reg.ab$climate=="tropical",]

# Diversity

reg.di.fi1 <- reg.di[!is.na(reg.di$fire_type)&reg.di$fire_type=="prescribed",]
reg.di.fi2 <- reg.di[!is.na(reg.di$fire_type)&reg.di$fire_type=="wildfire",]

reg.di.tt1 <- reg.di[reg.di$FR_comp=="frequency",]
reg.di.tt2 <- reg.di[reg.di$FR_comp=="severity",]

reg.di.rh1 <- reg.di[!is.na(reg.di$FR_hist)&reg.di$FR_hist=="no-fire",]
reg.di.rh2 <- reg.di[!is.na(reg.di$FR_hist)&reg.di$FR_hist=="surface",]
reg.di.rh3 <- reg.di[!is.na(reg.di$FR_hist)&reg.di$FR_hist=="crown",]

reg.di.ha1 <- reg.di[!is.na(reg.di$habitat)&reg.di$habitat=="forest-broadleaf",]
reg.di.ha2 <- reg.di[!is.na(reg.di$habitat)&reg.di$habitat=="forest-conifer",]
reg.di.ha3 <- reg.di[!is.na(reg.di$habitat)&reg.di$habitat=="forest-mixed",]
reg.di.ha4 <- reg.di[!is.na(reg.di$habitat)&reg.di$habitat=="grassland",]
reg.di.ha5 <- reg.di[!is.na(reg.di$habitat)&reg.di$habitat=="shrubland",]
reg.di.ha7 <- reg.di[!is.na(reg.di$habitat)&reg.di$habitat=="woodland",]

reg.di.pfg1 <- reg.di[!is.na(reg.di$PLF)&reg.di$PLF=="bryophyte",]
reg.di.pfg2 <- reg.di[!is.na(reg.di$PLF)&reg.di$PLF=="herb",]
reg.di.pfg3 <- reg.di[!is.na(reg.di$PLF)&reg.di$PLF=="woody",]

reg.di.cl1 <- reg.di[!is.na(reg.di$climate)&reg.di$climate=="arid",]
reg.di.cl2 <- reg.di[!is.na(reg.di$climate)&reg.di$climate=="cold",]
reg.di.cl3 <- reg.di[!is.na(reg.di$climate)&reg.di$climate=="temperate-dry",]
reg.di.cl4 <- reg.di[!is.na(reg.di$climate)&reg.di$climate=="temperate-nodry",]
reg.di.cl5 <- reg.di[!is.na(reg.di$climate)&reg.di$climate=="tropical",]

# fitness

reg.fit.fi1 <- reg.fit[!is.na(reg.fit$fire_type)&reg.fit$fire_type=="prescribed",]
reg.fit.fi2 <- reg.fit[!is.na(reg.fit$fire_type)&reg.fit$fire_type=="wildfire",]

reg.fit.tt1 <- reg.fit[!is.na(reg.fit$FR_comp)&reg.fit$FR_comp=="frequency",]
reg.fit.tt2 <- reg.fit[!is.na(reg.fit$FR_comp)&reg.fit$FR_comp=="severity",]

reg.fit.ft1 <- reg.fit[!is.na(reg.fit$TSF_factor)&reg.fit$TSF_factor=="short",]
reg.fit.ft2 <- reg.fit[!is.na(reg.fit$TSF_factor)&reg.fit$TSF_factor=="long",]

reg.fit.rh1 <- reg.fit[!is.na(reg.fit$FR_hist)&reg.fit$FR_hist=="no-fire",]
reg.fit.rh2 <- reg.fit[!is.na(reg.fit$FR_hist)&reg.fit$FR_hist=="surface",]
reg.fit.rh3 <- reg.fit[!is.na(reg.fit$FR_hist)&reg.fit$FR_hist=="crown",]

reg.fit.ha1 <- reg.fit[!is.na(reg.fit$habitat)&reg.fit$habitat=="forest-broadleaf",]
reg.fit.ha2 <- reg.fit[!is.na(reg.fit$habitat)&reg.fit$habitat=="forest-conifer",]
reg.fit.ha3 <- reg.fit[!is.na(reg.fit$habitat)&reg.fit$habitat=="forest-mixed",]
reg.fit.ha4 <- reg.fit[!is.na(reg.fit$habitat)&reg.fit$habitat=="grassland",]
reg.fit.ha5 <- reg.fit[!is.na(reg.fit$habitat)&reg.fit$habitat=="shrubland",]
reg.fit.ha7 <- reg.fit[!is.na(reg.fit$habitat)&reg.fit$habitat=="woodland",]

reg.fit.pfg1 <- reg.fit[!is.na(reg.fit$PLF)&reg.fit$PLF=="bryophyte",]
reg.fit.pfg2 <- reg.fit[!is.na(reg.fit$PLF)&reg.fit$PLF=="herb",]
reg.fit.pfg3 <- reg.fit[!is.na(reg.fit$PLF)&reg.fit$PLF=="woody",]

reg.fit.cl1 <- reg.fit[!is.na(reg.fit$climate)&reg.fit$climate=="arid",]
reg.fit.cl2 <- reg.fit[!is.na(reg.fit$climate)&reg.fit$climate=="cold",]
reg.fit.cl3 <- reg.fit[!is.na(reg.fit$climate)&reg.fit$climate=="temperate-dry",]
reg.fit.cl4 <- reg.fit[!is.na(reg.fit$climate)&reg.fit$climate=="temperate-nodry",]
reg.fit.cl5 <- reg.fit[!is.na(reg.fit$climate)&reg.fit$climate=="tropical",]


# RMA.MV analyses main ----

## Abundance ----

myrma.ab.fr <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                      random = list(~1|study_num/ES_num),
                      data = reg.ab)

myrma.ab.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                      random = list(~1|study_num/ES_num),
                      data = reg.ab)

myrma.ab.ti <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                      random = list(~1|study_num/ES_num),
                      data = reg.ab)

myrma.ab.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                      random = list(~1|study_num/ES_num),
                      data = reg.ab)

myrma.ab.PLF <- rma.mv(d_Hedges ~ PLF - 1, var_Hedges,
                        random = list(~1|study_num/ES_num),
                        data = reg.ab)

myrma.ab.ha <- rma.mv(d_Hedges ~ habitat - 1, var_Hedges,
                      random = list(~1|study_num/ES_num),
                      data = reg.ab)

myrma.ab.cl <- rma.mv(d_Hedges ~ climate - 1, var_Hedges,
                       random = list(~1|study_num/ES_num),
                       data = reg.ab)

myrma.ab <- rma.mv(d_Hedges, var_Hedges, 
                   random = list(~1|study_num/ES_num),
                   data = reg.ab)

## Diversity ----

myrma.di.fr <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                      random = list(~1|study_num/ES_num),
                      data = reg.di)

myrma.di.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                      random = list(~1|study_num/ES_num),
                      data = reg.di)

myrma.di.ti <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                      random = list(~1|study_num/ES_num),
                      data = reg.di)

myrma.di.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                      random = list(~1|study_num/ES_num),
                      data = reg.di)

myrma.di.PLF <- rma.mv(d_Hedges ~ PLF - 1, var_Hedges,
                       random = list(~1|study_num/ES_num),
                       data = reg.di)

myrma.di.ha <- rma.mv(d_Hedges ~ habitat - 1, var_Hedges,
                      random = list(~1|study_num/ES_num),
                      data = reg.di)

myrma.di.cl <- rma.mv(d_Hedges ~ climate - 1, var_Hedges,
                      random = list(~1|study_num/ES_num),
                      data = reg.di)

myrma.di <- rma.mv(d_Hedges, var_Hedges, 
                   random = list(~1|study_num/ES_num),
                   data = reg.di)


## Fitness ----

myrma.fit.fr <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                      random = list(~1|study_num/ES_num),
                      data = reg.fit)

myrma.fit.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                      random = list(~1|study_num/ES_num),
                      data = reg.fit)

myrma.fit.ti <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                      random = list(~1|study_num/ES_num),
                      data = reg.fit)

myrma.fit.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                      random = list(~1|study_num/ES_num),
                      data = reg.fit)

myrma.fit.PLF <- rma.mv(d_Hedges ~ PLF - 1, var_Hedges,
                       random = list(~1|study_num/ES_num),
                       data = reg.fit)

myrma.fit.ha <- rma.mv(d_Hedges ~ habitat - 1, var_Hedges,
                      random = list(~1|study_num/ES_num),
                      data = reg.fit)

myrma.fit.cl <- rma.mv(d_Hedges ~ climate - 1, var_Hedges,
                      random = list(~1|study_num/ES_num),
                      data = reg.fit)

myrma.fit <- rma.mv(d_Hedges, var_Hedges, 
                   random = list(~1|study_num/ES_num),
                   data = reg.fit)



# RMA.MV by subgroups ----

## Fire regime driver ----

# Abundance

myrma.ab.fi1.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.fi1)
myrma.ab.fi2.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.fi2)

myrma.ab.rh1.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.ab.rh1)
myrma.ab.rh2.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.ab.rh2)
myrma.ab.rh3.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.ab.rh3)

myrma.ab.pfg1.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.ab.pfg1)
myrma.ab.pfg2.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.ab.pfg2)
myrma.ab.pfg3.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.ab.pfg3)

myrma.ab.ha1.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.ha1)
myrma.ab.ha2.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.ha2)
myrma.ab.ha3.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.ha3)
myrma.ab.ha4.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.ha4)
myrma.ab.ha5.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.ha5)
myrma.ab.ha7.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.ha7)

myrma.ab.cl1.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.cl1)
myrma.ab.cl2.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.cl2)
myrma.ab.cl3.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.cl3)
myrma.ab.cl4.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.cl4)
myrma.ab.cl5.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.cl5)

# Diversity

myrma.di.fi1.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.fi1)
myrma.di.fi2.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.fi2)

myrma.di.rh1.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.di.rh1)
myrma.di.rh2.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.di.rh2)
myrma.di.rh3.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.di.rh3)

myrma.di.pfg1.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.di.pfg1)
myrma.di.pfg2.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.di.pfg2)
myrma.di.pfg3.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.di.pfg3)

myrma.di.ha1.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.ha1)
myrma.di.ha2.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.ha2)
myrma.di.ha3.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.ha3)
myrma.di.ha4.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.ha4)
myrma.di.ha5.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.ha5)
myrma.di.ha7.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.ha7)

myrma.di.cl1.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.cl1)
myrma.di.cl2.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.cl2)
myrma.di.cl3.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.cl3)
myrma.di.cl4.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.cl4)
myrma.di.cl5.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.cl5)

# Fitness

myrma.fit.fi1.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.fi1)
myrma.fit.fi2.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.fi2)

myrma.fit.rh1.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                            random = list(~1|study_num/ES_num),
                            data = reg.fit.rh1)
myrma.fit.rh2.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                            random = list(~1|study_num/ES_num),
                            data = reg.fit.rh2)
myrma.fit.rh3.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                            random = list(~1|study_num/ES_num),
                            data = reg.fit.rh3)

myrma.fit.pfg1.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                            random = list(~1|study_num/ES_num),
                            data = reg.fit.pfg1)
myrma.fit.pfg2.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                            random = list(~1|study_num/ES_num),
                            data = reg.fit.pfg2)
myrma.fit.pfg3.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                            random = list(~1|study_num/ES_num),
                            data = reg.fit.pfg3)

myrma.fit.ha1.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.ha1)
myrma.fit.ha2.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.ha2)
myrma.fit.ha3.tt <- rma.mv(d_Hedges, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.ha3[!is.na(reg.fit.ha3$FR_comp),]) #all severity
myrma.fit.ha4.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.ha4)
myrma.fit.ha5.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.ha5)
myrma.fit.ha7.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.ha7)

myrma.fit.cl1.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.cl1)
myrma.fit.cl2.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.cl2)
myrma.fit.cl3.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.cl3)
myrma.fit.cl4.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.cl4)
myrma.fit.cl5.tt <- rma.mv(d_Hedges ~ FR_comp - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.cl5)


## Fire type ----

# Abundance

myrma.ab.rh1.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.ab.rh1)
myrma.ab.rh2.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.ab.rh2)
myrma.ab.rh3.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.ab.rh3)

myrma.ab.pfg1.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.ab.pfg1)
myrma.ab.pfg2.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.ab.pfg2)
myrma.ab.pfg3.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.ab.pfg3)

myrma.ab.ha1.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.ha1)
myrma.ab.ha2.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.ha2)
myrma.ab.ha3.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.ha3)
myrma.ab.ha4.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.ha4)
myrma.ab.ha5.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.ha5)
myrma.ab.ha7.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.ha7)

myrma.ab.cl1.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.cl1)
myrma.ab.cl2.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.cl2)
myrma.ab.cl3.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.cl3)
myrma.ab.cl4.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.cl4)
myrma.ab.cl5.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.cl5)

# Diversity

myrma.di.rh1.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.di.rh1)
myrma.di.rh2.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.di.rh2)
myrma.di.rh3.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.di.rh3)

myrma.di.pfg1.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.di.pfg1)
myrma.di.pfg2.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.di.pfg2)
myrma.di.pfg3.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.di.pfg3)

myrma.di.ha1.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.ha1)
myrma.di.ha2.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.ha2)
myrma.di.ha3.fi <- rma.mv(d_Hedges, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.ha3[!is.na(reg.di.ha3$fire_type),]) #all wildfire
myrma.di.ha4.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.ha4)
myrma.di.ha5.fi <- rma.mv(d_Hedges, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.ha5[!is.na(reg.di.ha5$fire_type),]) #all wildfire
myrma.di.ha7.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.ha7)

myrma.di.cl1.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.cl1)
myrma.di.cl2.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.cl2)
myrma.di.cl3.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.cl3)
myrma.di.cl4.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.cl4)
myrma.di.cl5.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.cl5)

# Fitness

myrma.fit.rh1.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                            random = list(~1|study_num/ES_num),
                            data = reg.fit.rh1)
myrma.fit.rh2.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                            random = list(~1|study_num/ES_num),
                            data = reg.fit.rh2)
myrma.fit.rh3.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                            random = list(~1|study_num/ES_num),
                            data = reg.fit.rh3)

myrma.fit.pfg1.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                            random = list(~1|study_num/ES_num),
                            data = reg.fit.pfg1)
myrma.fit.pfg2.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                            random = list(~1|study_num/ES_num),
                            data = reg.fit.pfg2)
myrma.fit.pfg3.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                            random = list(~1|study_num/ES_num),
                            data = reg.fit.pfg3)

myrma.fit.ha1.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.ha1)
myrma.fit.ha2.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.ha2)
myrma.fit.ha3.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.ha3)
myrma.fit.ha4.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.ha4)
myrma.fit.ha5.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.ha5)
myrma.fit.ha7.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.ha7)

myrma.fit.cl1.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.cl1)
myrma.fit.cl2.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.cl2)
myrma.fit.cl3.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.cl3)
myrma.fit.cl4.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.cl4)
myrma.fit.cl5.fi <- rma.mv(d_Hedges ~ fire_type - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.cl5)


## Time since fire ----

# Abundance

myrma.ab.fi1.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.fi1)
myrma.ab.fi2.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.fi2)

myrma.ab.rh1.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.ab.rh1)
myrma.ab.rh2.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.ab.rh2)
myrma.ab.rh3.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.ab.rh3)

myrma.ab.pfg1.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.ab.pfg1)
myrma.ab.pfg2.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.ab.pfg2)
myrma.ab.pfg3.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.ab.pfg3)

myrma.ab.ha1.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.ha1)
myrma.ab.ha2.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.ha2)
myrma.ab.ha3.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.ha3)
myrma.ab.ha4.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.ha4)
myrma.ab.ha5.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.ha5)
myrma.ab.ha7.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.ha7)

myrma.ab.cl1.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.cl1)
myrma.ab.cl2.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.cl2)
myrma.ab.cl3.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.cl3)
myrma.ab.cl4.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.cl4)
myrma.ab.cl5.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.cl5)

# Diversity

myrma.di.fi1.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.fi1)
myrma.di.fi2.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.fi2)

myrma.di.rh1.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.di.rh1)
myrma.di.rh2.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.di.rh2)
myrma.di.rh3.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.di.rh3)

myrma.di.pfg1.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.di.pfg1)
myrma.di.pfg2.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.di.pfg2)
myrma.di.pfg3.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.di.pfg3)

myrma.di.ha1.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.ha1)
myrma.di.ha2.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.ha2)
myrma.di.ha3.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.ha3)
myrma.di.ha4.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.ha4)
myrma.di.ha5.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.ha5)
myrma.di.ha7.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.ha7)

myrma.di.cl1.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.cl1)
myrma.di.cl2.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.cl2)
myrma.di.cl3.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.cl3)
myrma.di.cl4.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.cl4)
myrma.di.cl5.ft <- rma.mv(d_Hedges, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.cl5[!is.na(reg.di.cl5$TSF_factor),])

# Fitness

myrma.fit.fi1.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.fi1)
myrma.fit.fi2.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.fi2)

myrma.fit.rh1.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                            random = list(~1|study_num/ES_num),
                            data = reg.fit.rh1)
myrma.fit.rh2.ft <- rma.mv(d_Hedges, var_Hedges,
                            random = list(~1|study_num/ES_num),
                            data = reg.fit.rh2)
myrma.fit.rh3.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                            random = list(~1|study_num/ES_num),
                            data = reg.fit.rh3)

myrma.fit.pfg1.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                            random = list(~1|study_num/ES_num),
                            data = reg.fit.pfg1)
myrma.fit.pfg2.ft <- rma.mv(d_Hedges, var_Hedges,
                            random = list(~1|study_num/ES_num),
                            data = reg.fit.pfg2)
myrma.fit.pfg3.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                            random = list(~1|study_num/ES_num),
                            data = reg.fit.pfg3)

myrma.fit.ha1.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.ha1)
myrma.fit.ha2.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.ha2)
myrma.fit.ha3.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.ha3)
myrma.fit.ha4.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.ha4)
myrma.fit.ha5.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.ha5)
myrma.fit.ha7.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.ha7)

myrma.fit.cl1.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.cl1)
myrma.fit.cl2.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.cl2)
myrma.fit.cl3.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.cl3)
myrma.fit.cl4.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.cl4)
myrma.fit.cl5.ft <- rma.mv(d_Hedges ~ TSF_factor - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit.cl5)


## Historical fire regime ----

# Abundance

myrma.ab.pfg1.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.ab.pfg1)
myrma.ab.pfg2.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.ab.pfg2)
myrma.ab.pfg3.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.ab.pfg3)

myrma.ab.ha1.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.ha1)
myrma.ab.ha2.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.ha2)
myrma.ab.ha3.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.ha3)
myrma.ab.ha4.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.ha4)
myrma.ab.ha5.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.ha5)
myrma.ab.ha7.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.ha7)

myrma.ab.cl1.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.cl1)
myrma.ab.cl2.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.cl2)
myrma.ab.cl3.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.cl3)
myrma.ab.cl4.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.cl4)
myrma.ab.cl5.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.cl5)

# Diversity

myrma.di.pfg1.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.di.pfg1)
myrma.di.pfg2.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.di.pfg2)
myrma.di.pfg3.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.di.pfg3)

myrma.di.ha1.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.ha1)
myrma.di.ha2.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.ha2)
myrma.di.ha3.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.ha3)
myrma.di.ha4.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.ha4) 
myrma.di.ha5.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.ha5) 
myrma.di.ha7.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.ha7)

myrma.di.cl1.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.cl1)
myrma.di.cl2.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.cl2)
myrma.di.cl3.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.cl3)
myrma.di.cl4.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.cl4)
myrma.di.cl5.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.cl5)

# Fitness

myrma.fit.pfg1.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                            random = list(~1|study_num/ES_num),
                            data = reg.fit3.pfg1)
myrma.fit.pfg2.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                            random = list(~1|study_num/ES_num),
                            data = reg.fit3.pfg2)
myrma.fit.pfg3.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                            random = list(~1|study_num/ES_num),
                            data = reg.fit3.pfg3)

myrma.fit.ha1.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit3.ha1)
myrma.fit.ha2.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit3.ha2)
myrma.fit.ha3.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit3.ha3)
myrma.fit.ha4.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit3.ha4)
myrma.fit.ha5.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit3.ha5)
myrma.fit.ha7.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit3.ha7)

myrma.fit.cl1.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit3.cl1)
myrma.fit.cl2.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit3.cl2)
myrma.fit.cl3.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit3.cl3)
myrma.fit.cl4.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit3.cl4)
myrma.fit.cl5.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fit3.cl5)


# RMA.MV surface-to-crown fires ----

# datasets

# arid, tropical, wetland, have too few observations

reg.crs <- subset(reg, surface_crown=="surface-crown")
reg.crs <- reg.crs[!reg.crs$PLF%in%c("bryophyte")&
                     !reg.crs$habitat%in%c("shrubland","grassland","wetland","woodland")&
                     !reg.crs$climate%in%c("arid","tropical")&
                     !reg.crs$FR_hist%in%c("no-fire"),]
reg.ab.crs <- reg.crs[reg.crs$response=="abundance",]
reg.di.crs <- reg.crs[reg.crs$response=="diversity",]
reg.fi.crs <- reg.crs[reg.crs$response=="fitness",]

reg.ab.crs <- droplevels(reg.ab.crs)
reg.di.crs <- droplevels(reg.di.crs)
reg.fi.crs <- droplevels(reg.fi.crs)


# abundance

myrma.ab.crs.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.ab.crs)

myrma.ab.crs.pfg <- rma.mv(d_Hedges ~ PLF - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.ab.crs)

myrma.ab.crs.ha <- rma.mv(d_Hedges ~ habitat - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.crs)

myrma.ab.crs.cl <- rma.mv(d_Hedges ~ climate - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.ab.crs)

myrma.ab.crs <- rma.mv(d_Hedges, var_Hedges,
                       random = list(~1|study_num/ES_num),
                       data = reg.ab.crs)

# diversity

myrma.di.crs.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.di.crs)
myrma.di.crs.pfg <- rma.mv(d_Hedges ~ PLF - 1, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.di.crs)
myrma.di.crs.ha <- rma.mv(d_Hedges ~ habitat - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.crs)
myrma.di.crs.cl <- rma.mv(d_Hedges ~ climate - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.di.crs)
myrma.di.crs <- rma.mv(d_Hedges, var_Hedges,
                       random = list(~1|study_num/ES_num),
                       data = reg.di.crs)

# fitness

myrma.fi.crs.rh <- rma.mv(d_Hedges ~ FR_hist - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.fi.crs)
myrma.fi.crs.pfg <- rma.mv(d_Hedges, var_Hedges,
                           random = list(~1|study_num/ES_num),
                           data = reg.fi.crs)
myrma.fi.crs.ha <- rma.mv(d_Hedges ~ habitat - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.fi.crs)
myrma.fi.crs.cl <- rma.mv(d_Hedges ~ climate - 1, var_Hedges,
                          random = list(~1|study_num/ES_num),
                          data = reg.fi.crs)
myrma.fi.crs <- rma.mv(d_Hedges, var_Hedges,
                       random = list(~1|study_num/ES_num),
                       data = reg.fi.crs)

# Publication bias ----

bias.ab <- lm(residuals(myrma.ab) ~ reg.ab$eff_n)
summary(bias.ab)

bias.di <- lm(residuals(myrma.di) ~ reg.di$eff_n)
summary(bias.di)

bias.fi <- lm(residuals(myrma.fit) ~ reg.fit$eff_n)
summary(bias.fi)