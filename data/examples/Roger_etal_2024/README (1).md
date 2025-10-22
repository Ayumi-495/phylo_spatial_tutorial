**Description of the data and the file structure**

The R code for carrying out the meta-analysis is provided in file "global_fire_analysis_240325.R".

The database is provided in the file "global_fire_DB_240325.csv". The columns in the file are described below:

study_id: study identification, composed of the name of the first author and the year of publication.

study_num: study number.

ES_num: number of effect size within study.

country: country where the study was carried out.

latitude/longitude: coordinates indicating where the study was carried out.

FR_hist: type of historical fire regime (surface fire, crown fire, or non fire prone).

FR_comp: fire regime component (frequency or severity).

fire_type: type of fire (wildfire or prescribed fire).

habitat: type of habitat (broadleaf forest, conifer forest, mixed forest, grassland, shrubland, or woodland).

response: broad response variable measured (abundance, diversity, or fitness).

response_detailed: more detailed response variable measured (e.g., size or biomass, density or frequency, Shannon or Simpson indices, gamma diversity, survival).

TSF_months: time since fire in months.

TSF_factor: factor variable of time since fire (short, 24 months or less, or long, more than 24 months).

PLF: plant life form studied (woody plant, herb, or bryophyte).

climate: climate prevalent in the study region (arid, cold, temperate with a dry season, temperate without a dry season, or tropical).

Xtreat1_type: type of crossed factor number one.

Xtreat1_fact: category for crossed factor number one.

Xtreat2_type: type of crossed factor number two.

Xtreat2_fact: category for crossed factor number two.

surface-crown: type of change in fire behaviour with intensification of fire regime (surface to crown, surface to surface, or crown to crown).

eff_n: effective sample size.

d_Hedges: Hedges’ d.

var_Hedges: variance of Hedges’ d.

imputed: whether the effect size was imputed or not.

source: source of the data used to calculate the effect size (e.g., Table1, Figure2, main text).

'NA' indicate missing values, e.g., because they were not reported in the source study, or because the variable does not apply to that study.
