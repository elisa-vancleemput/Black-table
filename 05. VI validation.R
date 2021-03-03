#######################################################
###                 05. VI validation              ###
#######################################################

### The objective of this script is to relate vegetaion index (VI) 'optical traits' estimated from hyperspectral data with lab data
#   Different relations might exist between estimated and real trait values (see e.g. Liu et al. 2017)
#   - linear
#   - exponential
#   - logarithmic

# The obtained optical traits are in this script compared to the lab measurements of the same observations.
# In other words: the validation/test set = lab measurements.

# The input for this validation procedure are:
#   - Paired lab measurements: created in '01b. Paired Field-Table-Lab data'
#   - Paired spectral measurements: created in '01b. Paired Field-Table-Lab data'
#       --> the same VI's as calculated in '04. Interspecific VI comparisons' are quantified here

# Script by Elisa Van Cleemput, 2018-19
#--------------------------------------------------------------------------------

# Clean workspace

rm(list=ls())

#--------------------------------------------------------------------------------

# Calculate the same VI's as used for the main analyses (script '04. Optical trait estimation VI')
VI_calc_lab <- function (wd1, jumpcorr) {
  setwd(wd1)
  datasets<-readRDS(paste("Spectra_paired_",".rds",sep=jumpcorr))
  
  indices<-c("(R925-R710)/(R925+R710)","TCARI2/OSAVI2","MCARI2/OSAVI2","MTCI",
             "REP_LE","PRI","(R860-R720)/(R860+R720)",
             "R539/R490","R807/R490",
             "(R860-R1240)/(R860+R1240)","(R858-R2130)/(R858+R2130)","PRI_norm",
             "(R1662-R1732)/(R1662+R1732)","(R1550-R1750)/(R1550+R1750)","(R1649-R1722)/(R1649+R1722)")
  
  # All individual measurements
  FG_indices <- vegindex(datasets$library_FG_paired,indices)
  TG_raw_indices <- vegindex(datasets$library_TG_raw_paired,indices)
  TG_unmixed_indices <- vegindex(datasets$library_TG_unmixed_paired,indices)
  # Mean spectrum per species for each Site_Date
  FG_indices.mean <- vegindex(datasets$library_FG.mean_paired,indices)
  TG_raw_indices.mean <- vegindex(datasets$library_TG_raw.mean_paired,indices)
  TG_unmixed_indices.mean <- vegindex(datasets$library_TG_unmixed.mean_paired,indices)
  
  colnames(FG_indices) <- c("ND_Chl","TC/OS","MC/OS","MTCI",
                            "REP","PRI","ND_N_P",
                            "SR_Car1","SR_Car2",
                            "NDWI","NDWI_SWIR","PRI_norm",
                            "ND_LMA1","ND_LMA2","NDMI")
  colnames(TG_raw_indices) <-colnames(FG_indices) 
  colnames(TG_unmixed_indices) <-colnames(FG_indices)
  colnames(FG_indices.mean) <-colnames(FG_indices)
  colnames(TG_raw_indices.mean) <- colnames(FG_indices)
  colnames(TG_unmixed_indices.mean) <- colnames(FG_indices)
  
  VI_output <- list ("VI_FG" = FG_indices,
                     "VI_TG_raw" = TG_raw_indices,
                     "VI_TG_unmixed" = TG_unmixed_indices,
                     "VI_FG.mean" = FG_indices.mean,
                     "VI_TG_raw.mean" = TG_raw_indices.mean,
                     "VI_TG_unmixed.mean" = TG_unmixed_indices.mean)
  return(VI_output)
}

#--------------------------------------------------------------------------------

# How well do the predicted results mirror the real lab measurements?
### -------- Calculate R² and RMSE of many VI's 
FitFunction <- function(PFT,VI_values,trait){
  # PFT <- data$PFT
  
  fit.stat<-matrix(names(VI_values)[1:ncol(VI_values)],ncol=1)
  fit.stat<-cbind(fit.stat,array(rep(0),dim=c(length(names(VI_values)[1:ncol(VI_values)]),6)))
  colnames(fit.stat)<-c("VI","R2_lin","nRMSE_lin","R2_exp","nRMSE_exp","R2_log","nRMSE_log")
  
  for (i in 1 : ncol(VI_values)){
    fit.lin<-lm(PFT[,trait]~VI_values[,i]) # y = a+b*x
    fit.exp<-lm(log(PFT[,trait])~VI_values[,i]) # back transform to y = a*exp(b*x)  
    # https://stat.ethz.ch/pipermail/r-help/2010-January/224131.html
    # https://www.theanalysisfactor.com/r-tutorial-5/
    fit.log<-tryCatch(lm(PFT[,trait]~log(VI_values[,i])) , warning = function(w) NA) # equivalent to y = a+b*log(x) 
    
    # --------- R²
    fit.stat[i,"R2_lin"]<-summary(fit.lin)$r.squared # idem to: cor(VI_values[,i],PFT[,c("Cab")])^2
    fit.stat[i,"R2_exp"]<-summary(fit.exp)$r.squared
    if (is.na(fit.log[1])){ fit.stat[i,"R2_log"]<-NA
    } else { fit.stat[i,"R2_log"]<-summary(fit.log)$r.squared }
    
    # --------- RMSE: Different calculation methods
    fit.stat[i,"nRMSE_lin"]<-sqrt(c(crossprod(fit.lin$residuals))/length(fit.lin$residuals))/mean(PFT[,trait])
    # library(qpcR) 
    # RMSE(fit.lin)
    # scores<-predict(lm(PFT[,trait]~VI_values[,i]))
    # sqrt(mean((PFT[,trait]-scores)^2))
    scores_exp<-exp(predict(lm(log(PFT[,trait])~VI_values[,i])))
    fit.stat[i,"nRMSE_exp"]<-sqrt(mean((PFT[,trait]-scores_exp)^2))/mean(PFT[,trait])
    if (is.na(fit.log[1])){  fit.stat[i,"nRMSE_log"]<-NA
    } else {  fit.stat[i,"nRMSE_log"]<-sqrt(c(crossprod(fit.log$residuals))/length(fit.log$residuals))/mean(PFT[,trait])}
    # scores_log<-predict(lmPFT[,trait]~log(VI_values[,i])))
    # sqrt(mean((PFT[,trait]-scores_log)^2))
  }
  
  return(fit.stat)
}
### -------- Only retain VI rendering max R² and/or min RMSE
BestFunction <- function(trait,fit.input,best.output){
  best.output[trait,"VI_R2"]<-fit.input[which(fit.input == max(as.numeric(fit.input[,c("R2_lin","R2_exp","R2_log")]),na.rm=T), arr.ind = T)[1,1],"VI"]
  best.output[trait,"R2"]<-round(max(as.numeric(fit.input[,c("R2_lin","R2_exp","R2_log")]),na.rm=T),digits=2)
  best.output[trait,"relation_R2"]<-colnames(fit.input)[which(fit.input == max(as.numeric(fit.input[,c("R2_lin","R2_exp","R2_log")]),na.rm=T), arr.ind = T)[1,2]]
  best.output[trait,"VI_nRMSE"]<-fit.input[which(fit.input == min(as.numeric(fit.input[,c("nRMSE_lin","nRMSE_exp","nRMSE_log")]),na.rm=T), arr.ind = T)[1,1],"VI"]
  best.output[trait,"nRMSE"]<-round(min(as.numeric(fit.input[,c("nRMSE_lin","nRMSE_exp","nRMSE_log")]),na.rm=T),digits=2)
  best.output[trait,"relation_nRMSE"]<-colnames(fit.input)[which(fit.input == min(as.numeric(fit.input[,c("nRMSE_lin","nRMSE_exp","nRMSE_log")]),na.rm=T), arr.ind = T)[1,2]]
  
  return(best.output)
}
BestFunction_perRel <- function(trait,fit.input,best.output){
  best.output[trait,"VI_R2_lin"]<-fit.input[which(fit.input == max(as.numeric(fit.input[,c("R2_lin")]),na.rm=T), arr.ind = T)[1,1],"VI"]
  best.output[trait,"R2_lin"]<-round(max(as.numeric(fit.input[,c("R2_lin")]),na.rm=T),digits=2)
  best.output[trait,"VI_nRMSE_lin"]<-fit.input[which(fit.input == min(as.numeric(fit.input[,c("nRMSE_lin")]),na.rm=T), arr.ind = T)[1,1],"VI"]
  best.output[trait,"nRMSE_lin"]<-round(min(as.numeric(fit.input[,c("nRMSE_lin")]),na.rm=T),digits=2)
  
  best.output[trait,"VI_R2_exp"]<-fit.input[which(fit.input == max(as.numeric(fit.input[,c("R2_exp")]),na.rm=T), arr.ind = T)[1,1],"VI"]
  best.output[trait,"R2_exp"]<-round(max(as.numeric(fit.input[,c("R2_exp")]),na.rm=T),digits=2)
  best.output[trait,"VI_nRMSE_exp"]<-fit.input[which(fit.input == min(as.numeric(fit.input[,c("nRMSE_exp")]),na.rm=T), arr.ind = T)[1,1],"VI"]
  best.output[trait,"nRMSE_exp"]<-round(min(as.numeric(fit.input[,c("nRMSE_exp")]),na.rm=T),digits=2)
  
  best.output[trait,"VI_R2_log"]<-fit.input[which(fit.input == max(as.numeric(fit.input[,c("R2_log")]),na.rm=T), arr.ind = T)[1,1],"VI"]
  best.output[trait,"R2_log"]<-round(max(as.numeric(fit.input[,c("R2_log")]),na.rm=T),digits=2)
  best.output[trait,"VI_nRMSE_log"]<-fit.input[which(fit.input == min(as.numeric(fit.input[,c("nRMSE_log")]),na.rm=T), arr.ind = T)[1,1],"VI"]
  best.output[trait,"nRMSE_log"]<-round(min(as.numeric(fit.input[,c("nRMSE_log")]),na.rm=T),digits=2)
  
  return(best.output)
}
BestFunction_perVI <- function(trait,fit.input,indices){
  best.output<-array(rep(0),dim=c(length(indices),4))
  rownames(best.output)<-indices
  colnames(best.output)<-c("R2","relation_R2","nRMSE","relation_nRMSE")
  
  for (i in 1:nrow(best.output)){
    # best.output[trait,"VI_R2"]<-fit.input[which(fit.input == max(as.numeric(fit.input[,c("R2_lin","R2_exp","R2_log")]),na.rm=T), arr.ind = T)[1,1],"VI"]
    best.output[i,"R2"]<-round(max(as.numeric(fit.input[i,c("R2_lin","R2_exp","R2_log")]),na.rm=T),digits=2)
    best.output[i,"relation_R2"]<-colnames(fit.input)[which(fit.input == max(as.numeric(fit.input[i,c("R2_lin","R2_exp","R2_log")]),na.rm=T), arr.ind = T)[1,2]]
    # best.output[trait,"VI_nRMSE"]<-fit.input[which(fit.input == min(as.numeric(fit.input[,c("nRMSE_lin","nRMSE_exp","nRMSE_log")]),na.rm=T), arr.ind = T)[1,1],"VI"]
    best.output[i,"nRMSE"]<-round(min(as.numeric(fit.input[i,c("nRMSE_lin","nRMSE_exp","nRMSE_log")]),na.rm=T),digits=2)
    best.output[i,"relation_nRMSE"]<-colnames(fit.input)[which(fit.input == min(as.numeric(fit.input[i,c("nRMSE_lin","nRMSE_exp","nRMSE_log")]),na.rm=T), arr.ind = T)[1,2]]
  }
  
  return(best.output)
}

### -------- Apply FitFunction and BestFunction to a list of traits
Best_VI_per_Trait <- function(data, VI_values,VI_values_string, trait_list) {
  fit.stat.best<-array(rep(0),dim=c(length(trait_list),6))
  rownames(fit.stat.best)<-trait_list
  colnames(fit.stat.best)<-c("VI_R2","R2","relation_R2","VI_nRMSE","nRMSE","relation_nRMSE")
  fit.stat.best.predefinedVI<-fit.stat.best
  for (i in 1:length(trait_list)) {
    fit.stat <- FitFunction(data, VI_values, trait_list[i])
    fit.stat.best<-BestFunction(trait_list[i],fit.stat,fit.stat.best)
    # fit.stat.best.predefinedVI<-BestFunction(trait_list[i],fit.stat[1:15,],fit.stat.best.predefinedVI)
  }
  
  fit.stat.best.perRel<-array(rep(0),dim=c(length(trait_list),12))
  rownames(fit.stat.best.perRel)<-trait_list
  colnames(fit.stat.best.perRel)<-c("VI_R2_lin","R2_lin","VI_nRMSE_lin","nRMSE_lin","VI_R2_exp","R2_exp","VI_nRMSE_exp","nRMSE_exp","VI_R2_log","R2_log","VI_nRMSE_log","nRMSE_log")
  fit.stat.best.perRel.predefinedVI <- fit.stat.best.perRel
  for (i in 1:length(trait_list)) {
    fit.stat <- FitFunction(data, VI_values, trait_list[i])
    fit.stat.best.perRel<-BestFunction_perRel(trait_list[i],fit.stat,fit.stat.best.perRel)
    # fit.stat.best.perRel.predefinedVI<-BestFunction_perRel(trait_list[i],fit.stat[1:15,],fit.stat.best.perRel.predefinedVI)
  }
  
  output <- list("fit.stat.best" = fit.stat.best,
                 # "fit.stat.best.predefinedVI" = fit.stat.best.predefinedVI,
                 "fit.stat.best.perRel" =  fit.stat.best.perRel)
                 # "fit.stat.best.perRel.predefinedVI" =  fit.stat.best.perRel.predefinedVI)
  return(output)
  
  
}

VI_validation <- function(PFT_traits, VI, datasetname){
  indices_chl <- c("ND_Chl","TC/OS","MC/OS","MTCI","PRI")
  indices_carot <- c("SR_Car1","SR_Car2")
  indices_LNC <- c("MC/OS","MTCI","REP","PRI","ND_N_P")
  indices_LPC <- c("PRI", "ND_N_P")
  indices_Water <- c("NDWI", "NDWI_SWIR", "PRI_norm")
  indices_LMA <- c("ND_LMA1", "ND_LMA2", "NDMI")
  
  fit_Chl_area_Nocorr<-FitFunction(PFT_traits,VI[[datasetname]][,indices_chl],"Chltotal_area")
  fit_Chl_mass_Nocorr<-FitFunction(PFT_traits,VI[[datasetname]][,indices_chl],"Chltotal_mass")
  fit_Car_area_Nocorr<-FitFunction(PFT_traits,VI[[datasetname]][,indices_carot],"Carot_area")
  fit_Car_mass_Nocorr<-FitFunction(PFT_traits,VI[[datasetname]][,indices_carot],"Carot_mass")
  fit_LNC_area_Nocorr<-FitFunction(PFT_traits,VI[[datasetname]][,indices_LNC],"LNC_area")
  fit_LNC_mass_Nocorr<-FitFunction(PFT_traits,VI[[datasetname]][,indices_LNC],"LNC_mass")
  fit_LPC_area_Nocorr<-FitFunction(PFT_traits,VI[[datasetname]][,indices_LPC],"LPC_area")
  fit_LPC_mass_Nocorr<-FitFunction(PFT_traits,VI[[datasetname]][,indices_LPC],"LPC_mass")
  fit_Water_area_Nocorr<-FitFunction(PFT_traits,VI[[datasetname]][,indices_Water],"Water_area")
  fit_Water_mass_Nocorr<-FitFunction(PFT_traits,VI[[datasetname]][,indices_Water],"Water_mass")
  fit_SLA_Nocorr<-FitFunction(PFT_traits,VI[[datasetname]][,indices_LMA],"SLA")
  fit_LMA_Nocorr<-FitFunction(PFT_traits,VI[[datasetname]][,indices_LMA],"LMA")
  
  chl_area <- BestFunction_perVI("Chltotal_area", fit_Chl_area_Nocorr, indices_chl)
  chl_mass <- BestFunction_perVI("Chltotal_mass", fit_Chl_mass_Nocorr, indices_chl)
  car_area <- BestFunction_perVI("Carot_area", fit_Car_area_Nocorr, indices_carot)
  car_mass <- BestFunction_perVI("Carot_mass", fit_Car_mass_Nocorr, indices_carot)
  LNC_area <- BestFunction_perVI("LNC_area", fit_LNC_area_Nocorr, indices_LNC)
  LNC_mass <- BestFunction_perVI("LNC_mass", fit_LNC_mass_Nocorr, indices_LNC)
  LPC_area <- BestFunction_perVI("LPC_area", fit_LPC_area_Nocorr, indices_LPC)
  LPC_mass <- BestFunction_perVI("LPC_mass", fit_LPC_mass_Nocorr, indices_LPC)
  Water_area <- BestFunction_perVI("Water_area", fit_Water_area_Nocorr, indices_Water)
  Water_mass <- BestFunction_perVI("Water_mass", fit_Water_mass_Nocorr, indices_Water)
  SLA <- BestFunction_perVI("SLA", fit_SLA_Nocorr, indices_LMA)
  LMA <- BestFunction_perVI("LMA", fit_LMA_Nocorr, indices_LMA)
  
  output <- list("chl_area" = chl_area,
                 "chl_mass" = chl_mass,
                 "car_area" = car_area,
                 "car_mass" = car_mass,
                 "LNC_area" = LNC_area,
                 "LNC_mass" = LNC_mass,
                 "LPC_area" = LPC_area,
                 "LPC_mass" = LPC_mass,
                 "Water_area" = Water_area,
                 "Water_mass" = Water_mass,
                 "SLA" = SLA,
                 "LMA"= LMA)
}

#--------------------------------------------------------------------------------
# Load trait data 

setwd("C:/Users/u0091812/Box Sync/03. Black table procedure/Data/Paired Field-Table-Lab list")
PFT<-read.csv("PFT_paired.csv",sep=';',header=T)

PFT <- PFT[order(PFT$ID),] # order the observations according to ID
target <- which(names(PFT) == 'ID')[1]
PFT_meta<-PFT[,1:target]
PFT_traits<-PFT[,-c(1:target)]

# Traits can be estimated on area or mass basis: 
# Chtotal, Chla, Chlb: mg/g
# Chltotal_area, Chla_area, Chlb_area: µg/cm²
PFT_traits$Chltotal_mass <- PFT_traits$Chltotal
PFT_traits$Carot_mass <- PFT_traits$Carot
PFT_traits$LNC_mass <- 10 * PFT_traits$LNC # % dry mass to mg/g
PFT_traits$LNC_area <- PFT_traits$LNC_mass / PFT_traits$SLA * 10^-3 # mg/mm²
PFT_traits$LPC_mass <- 10 * PFT_traits$P # % dry mass to mg/g
PFT_traits$LPC_area <- PFT_traits$LPC_mass / PFT_traits$SLA * 10^-3 # mg/mm²
PFT_traits$LMA <- 1/PFT_traits$SLA
PFT_traits$Water_area<-PFT$EWT # cm = cm³/cm² = 10²/(SLA*LDMC)
PFT_traits$Water_mass<-(1/(PFT$LDMC*10^-3))-1 #  g/g


#--------------------------------------------------------------------------------
#  APPLY FUNCTIONS
wd1 = "C:/Users/u0091812/Box Sync/03. Black table procedure/Data/Paired Field-Table-Lab list"
VI_Nocorr_lab <- VI_calc_lab(wd1,"Nocorr")


trait_list <- c("Chltotal_area","Chltotal_mass","Carot_area","Carot_mass",
                "LNC_area","LNC_mass", "LPC_area","LPC_mass",
                "Water_area","Water_mass","SLA","LMA","LDMC")
VI_perfomance_FG_Nocorr <- Best_VI_per_Trait(PFT_traits, VI_Nocorr_lab$VI_FG.mean, "test", trait_list)
VI_perfomance_TG_raw_Nocorr <- Best_VI_per_Trait(PFT_traits, VI_Nocorr_lab$VI_TG_raw.mean, "test", trait_list)
VI_perfomance_TG_unmixed_Nocorr <- Best_VI_per_Trait(PFT_traits, VI_Nocorr_lab$VI_TG_unmixed.mean, "test", trait_list)

VI_perfomance_FG_Nocorr$fit.stat.best
VI_perfomance_TG_raw_Nocorr$fit.stat.best
VI_perfomance_TG_unmixed_Nocorr$fit.stat.best

# The results we are really interested in:
VI_FG.mean_validation <- VI_validation(PFT_traits,VI_Nocorr_lab, datasetname="VI_FG.mean")
VI_TG_unmixed.mean_validation <- VI_validation(PFT_traits,VI_Nocorr_lab,"VI_TG_unmixed.mean")
