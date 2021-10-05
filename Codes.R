###IVW, MR-Egger, weighted median methods
library(MRInstruments)
library(TwoSampleMR)

exposure_dat<-clump_data(exposure_dat, clump_r2 = 0.001,clump_kb=5000)#Clumping the exposure files

outcome_dat <-extract_outcome_data(snps=exposure_dat$SNP, outcomes="ieu-a-966", proxies = FALSE,
 maf_threshold = 0.01,) #extracting outcome data from the MRbase
 
mydata <- harmonise_data(
       exposure_dat=exposure_dat,
       outcome_dat=outcome_dat,
       action= 2
       )

res <- mr(mydata, method_list = c("mr_wald_ratio","mr_egger_regression","mr_weighted_median","mr_ivw","mr_ivw_fe"))
res
generate_odds_ratios(res)

pleio <- mr_pleiotropy_test(mydata)
pleio

het <- mr_heterogeneity(mydata)
het

single <- mr_leaveoneout(mydata)
mr_leaveoneout_plot(single)

mr_scatter_plot(res,mydata)

res_single <- mr_singlesnp(mydata)
mr_forest_plot(res_single)

mr_funnel_plot(res_single)

###MR-PRESSO methods
library(MRPRESSO)
mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", 
SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = mydata,
 NbDistribution = 2000,  SignifThreshold = 0.05)












