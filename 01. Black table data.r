#######################################################
###              01. Black table data              ###
#######################################################

### The objective of all 'black table scripts' is to evaluate the performance of the black table
# In this specific script we prepare the data needed in following scripts for
#   - plotting spectral signatures (02. Spectral signature plots)
#   - calculating (dis)similarity indices within and between species (03. Inter- and intraspecific similarities)

### The input of these steps are the average spectra of the bouquets/patches of either
# - the raw spectra (created in "04. Average bouquets") --> average field patches and raw table measurements
# - the unmixed spectra (created in "05. Signal unmixing") --> average unmixed table measurements

# ! These signatures still need following preprocessing: smoothing, removal of atmospheric noise windows
# asd and sv are smoothed with a 51 and 101 window respectively

# The script contains 4 parts:
#   1) Create spectral library of both raw and unmixed datasets
#   2) Subset the libraries according to field/table and flower/green --> needed for visualisation in script "02. Spectral signature plots"
#   3) Create average libraries of field and table spectra --> needed for (dis)similarity index calculations in script "03. Inter- and intraspecific similarities"
#   4) Create 'green' libraries of field and table spectra --> needed for optical traits and ANOVA in "04. Optical trait comparisons"

# Script by Elisa Van Cleemput, 2018
#--------------------------------------------------------------------------------

# Clean workspace

rm(list=ls())

#--------------------------------------------------------------------------------

# Load libraries
library(hsdar)
library(colorRamps)

#--------------------------------------------------------------------------------
##################################################################################
# Create function that subsets data accordign to the needs for different analyses
Black_data <- function (wd1,wd2,wd3,wd4,jumpcorr,fracperc,remove) {
  
  # Raw spectra
  setwd(wd1)
  # dataset_All<-read.csv("CompleteDataset_All_Average.csv",sep=";",header=TRUE)
  dataset_Veg<-read.csv(paste("CompleteDataset_Veg_Average_",".csv",sep=jumpcorr),sep=";",header=TRUE)
  dataset_Veg<-dataset_Veg[!dataset_Veg$Code %in% c("20158015_Papen0028.sv","20158015_Papen0051.sv","20160817_Btop00014.asd","20180628_L00018.asd") ,] #We remove the four measurements which were also removed from the unmixed dataset due to unreal values
  # dataset_Black<-read.csv(paste("CompleteDataset_Black_Average_",".csv",sep=jumpcorr),sep=";",header=TRUE)
  # dataset_White<-read.csv(paste("CompleteDataset_White_Average_",".csv",sep=jumpcorr),sep=";",header=TRUE)
  levels(dataset_Veg$Site)[levels(dataset_Veg$Site) == "BattNT" | levels(dataset_Veg$Site) == "BattNM"] <- "BattN"
  
  # Signal unmixed spectra
  setwd(wd2)
  dataset_Veg_unmixed<-read.csv(paste(paste("CompleteDataset_Veg_Average_Unmixed_",fracperc,sep=jumpcorr),".csv",sep=""),sep=";",header=TRUE)
  levels(dataset_Veg_unmixed$Site)[levels(dataset_Veg_unmixed$Site) == "BattNT" | levels(dataset_Veg_unmixed$Site) == "BattNM"] <- "BattN"

  # Remove CIRSP measurements because several Cirsium species were measured during field work
  if (remove == "1"){
    dataset_Veg<-dataset_Veg[!dataset_Veg$Species == "CIRSP",] 
    dataset_Veg_unmixed<-dataset_Veg_unmixed[!dataset_Veg_unmixed$Species == "CIRSP",] 
  }
  
#--------------------------------------------------------------------------------

#####   1) Create spectral library of both raw and unmixed dataset

# Raw data (field and table spectra)
target_raw <- which(names(dataset_Veg) == 'nm350')[1]
spectra_raw<-dataset_Veg[,target_raw:ncol(dataset_Veg)]
meta_raw<-dataset_Veg[,2:target_raw-1]
rownames(spectra_raw)=dataset_Veg[,1]
rownames(meta_raw)=dataset_Veg[,1]
# Due to the correction made above we have to update some other columns
meta_raw$Site_Date<-as.factor(paste(meta_raw$Site, meta_raw$Date))
meta_raw$Species_Nr_Site_Date<-as.factor(paste(meta_raw$Species_Nr, meta_raw$Site_Date))
meta_raw$Species_Cond_Site_Date<-as.factor(paste(meta_raw$Species, meta_raw$Cond , meta_raw$Site_Date))
# Extra metacolumns
meta_raw$Species_Site_Date<-as.factor(paste(meta_raw$Species_Nr, meta_raw$Site_Date))
meta_raw$Sp_SiteDate_Phase<-as.factor(paste(meta_raw$Species,meta_raw$Site_Date,meta_raw$DvpPhase)) 
meta_raw$Sp_Phase<-as.factor(paste(meta_raw$Species,meta_raw$DvpPhase))
meta_raw$ID<-paste(meta_raw$Species, meta_raw$Site_Date)


# unmixed data (only table spectra)
target_unmixed <- which(names(dataset_Veg_unmixed) == 'nm350')[1]
spectra_unmixed<-dataset_Veg_unmixed[,target_unmixed:ncol(dataset_Veg_unmixed)]
meta_unmixed<-dataset_Veg_unmixed[,2:target_unmixed-1]
rownames(spectra_unmixed)=dataset_Veg_unmixed[,1]
rownames(meta_unmixed)=dataset_Veg_unmixed[,1]
# Due to the correction made above we have to update some other columns
meta_unmixed$Site_Date<-as.factor(paste(meta_unmixed$Site, meta_unmixed$Date))
meta_unmixed$Species_Nr_Site_Date<-as.factor(paste(meta_unmixed$Species_Nr, meta_unmixed$Site_Date))
meta_unmixed$Species_Cond_Site_Date<-as.factor(paste(meta_unmixed$Species, meta_unmixed$Cond , meta_unmixed$Site_Date))
# Extra metacolumns
meta_unmixed$Species_Site_Date<-as.factor(paste(meta_unmixed$Species_Nr, meta_unmixed$Site_Date))
meta_unmixed$Sp_SiteDate_Phase<-as.factor(paste(meta_unmixed$Species,meta_unmixed$Site_Date,meta_unmixed$DvpPhase)) 
meta_unmixed$Sp_Phase<-as.factor(paste(meta_unmixed$Species,meta_unmixed$DvpPhase))
meta_unmixed$ID<-paste(meta_unmixed$Species, meta_unmixed$Site_Date)

# Wavelengths
wl=as.numeric(gsub("nm","",names(spectra_unmixed))) # idem for spectra_raw

# Create spectral libraries
library_raw<-speclib(as.matrix(spectra_raw),wl)
library_unmixed<-speclib(as.matrix(spectra_unmixed),wl)
SI(library_raw)=meta_raw #'attribute' function has changed to 'SI' in version 0.6.0
SI(library_unmixed)=meta_unmixed

idSpeclib(library_raw)<-as.character(meta_raw[1,])  
idSpeclib(library_unmixed)<-as.character(meta_unmixed[1,])  

#--------------------------------------------------------------------------------

#####   2) Subset the libraries according to field/table and flower/green --> needed for visualisation in script "02. Spectral signature plots"

# Field measurements (F): without (G) or with flowers (F) 
library_FG<-subset(library_raw,Patch == "Field" & DvpPhase == "Green")
List_FG<- levels(factor((SI(library_FG))$Species)) #'attribute' function has changed to 'SI' in version 0.6.0
meta_FG<-SI(library_FG)
library_FF<-subset(library_raw,Patch == "Field" & DvpPhase == "Flower")
List_FF<- levels(factor((SI(library_FF))$Species))
meta_FF<-SI(library_FF)

# Raw table measurements (T): without (G) or with flowers (F)
library_TG_raw<-subset(library_raw,Patch == "Table" & DvpPhase == "Green")
List_TG_raw<- levels(factor((SI(library_TG_raw))$Species))
meta_TG_raw<-SI(library_TG_raw)
library_TF_raw<-subset(library_raw,Patch == "Table" & DvpPhase == "Flower")
List_TF_raw<- levels(factor((SI(library_TF_raw))$Species))
meta_TF_raw<-SI(library_TF_raw)

# Unmixed table measurements (T): without (G) or with flowers (F)
library_TG_unmixed<-subset(library_unmixed,Patch == "Table" & DvpPhase == "Green")
List_TG_unmixed<- levels(factor((SI(library_TG_unmixed))$Species))
meta_TG_unmixed<-SI(library_TG_unmixed)
library_TF_unmixed<-subset(library_unmixed,Patch == "Table" & DvpPhase == "Flower")
List_TF_unmixed<- levels(factor((SI(library_TF_unmixed))$Species))
meta_TF_unmixed<-SI(library_TF_unmixed)

#--------------------------------------------------------------------------------

#####   3) Create average libraries of field and table spectra 
#       --> needed for (dis)similarity index calculations in script "03 Inter- and intraspecific similarities"
#       --> needed for procrust analysis in "04. Optical trait comparisons"

# A) We have to make sure that only those spectra which appear in both the field and table dataset are included --> Determine intersect to make sure we have a paired dataset
meta_F<-meta_raw[which(meta_raw$Patch== "Field"),]
List_F<-levels(droplevels(meta_F$Sp_SiteDate_Phase))
List_T_unmixed<-levels(meta_unmixed$Sp_SiteDate_Phase) #raw might have more attributes (in case we were not able to unmix). We choose smallest dataset
Intersection<-intersect(List_F,List_T_unmixed)  # 30 correspondances

# B) AVERAGE FIELD and TABLE spectrum per SPECIES, SITE, DATE, AND CONDITION (Patch - Dvpmt phase)
#     This approach accounts for intraspecific variability between sites and dates
#     Needed for evaluations WITHIN species

#     Field
libr_F.intersect <- subset(library_raw, Patch == "Field")
idx_int <- SI(libr_F.intersect)$Sp_SiteDate_Phase %in% Intersection
libr_F.intersect<-libr_F.intersect[idx_int,]
# libr_F.intersect <- subset(library_raw, Patch == "Field" & Sp_SiteDate_Phase %in% Intersection)
meta_F.intersect<-SI(libr_F.intersect)
meta_F.intersect$Sp_SiteDate_Phase<-droplevels(meta_F.intersect$Sp_SiteDate_Phase)
meta_F.intersect$Species_Cond<-droplevels(meta_F.intersect$Species_Cond)
meta_F.intersect$Species_Site_Date<-droplevels(meta_F.intersect$Species_Site_Date)
meta_F.intersect$Species_Nr_Site_Date<-droplevels(meta_F.intersect$Species_Nr_Site_Date)

SI(libr_F.intersect)<-meta_F.intersect

libr_F.intersect.mean <- apply(libr_F.intersect, FUN = mean, bySI = "Sp_SiteDate_Phase") #SI is not transfered, so we add it in the following lines
list_F = levels(meta_F.intersect$Sp_SiteDate_Phase) #equals 'Intersection'
index_F<-which(meta_F.intersect$Sp_SiteDate_Phase == list_F[1])
drop<-c("Nr","Species_Nr","Height") #These colums should be ignored because they are not relevant for the mean dataset
meta_F.intersect.mean<-meta_F.intersect[index_F[1], !names(meta_F.intersect) %in% drop] 
for (i in 2:length(list_F)){
  index_F<-which(meta_F.intersect$Sp_SiteDate_Phase == list_F[i])
  meta_F.intersect.mean[i,]<-meta_F.intersect[index_F[1],!names(meta_F.intersect) %in% drop]
}
SI(libr_F.intersect.mean)=meta_F.intersect.mean

#     Table

###   Another preprocessing round for table data
# After comparison I found out that raw table intersection library contained 1 more record than the unmixed library
# I think there is one raw spectrum which has not been unmixed (because no fraction cover calculated?)
# Remove this point in order to be able to correctly compare raw vs. unmixed
#     Table raw
libr_T_raw.intersect1 <- subset(library_raw, Patch == "Table")
idx_int <- SI(libr_T_raw.intersect1)$Sp_SiteDate_Phase %in% Intersection
libr_T_raw.intersect1<-libr_T_raw.intersect1[idx_int,]
# libr_T_raw.intersect1<-subset(library_raw, Patch == "Table" & Sp_SiteDate_Phase %in% Intersection)
meta_T_raw.intersect1<-SI(libr_T_raw.intersect1)
#     Table unmixed
libr_T_unmixed.intersect1 <- subset(library_unmixed, Patch == "Table")
idx_int <- SI(libr_T_unmixed.intersect1)$Sp_SiteDate_Phase %in% Intersection
libr_T_unmixed.intersect1<-libr_T_unmixed.intersect1[idx_int,]
# libr_T_unmixed.intersect1<-subset(library_unmixed, Patch == "Table" & Sp_SiteDate_Phase %in% Intersection)
meta_T_unmixed.intersect1<-SI(libr_T_unmixed.intersect1)
Intersection2<-intersect(levels(meta_T_raw.intersect1$Code),levels(meta_T_unmixed.intersect1$First_measurement))

###   Further processing
#     Table raw
idx_int <- SI(libr_T_raw.intersect1)$Code %in% Intersection2
libr_T_raw.intersect<-libr_T_raw.intersect1[idx_int,]
# libr_T_raw.intersect<-subset(libr_T_raw.intersect1, Code %in% Intersection2)
meta_T_raw.intersect<-SI(libr_T_raw.intersect)
meta_T_raw.intersect$Species_Cond<-droplevels(meta_T_raw.intersect$Species_Cond)
meta_T_raw.intersect$Species_Site_Date<-droplevels(meta_T_raw.intersect$Species_Site_Date)
meta_T_raw.intersect$Sp_SiteDate_Phase<-droplevels(meta_T_raw.intersect$Sp_SiteDate_Phase)
meta_T_raw.intersect$Species_Nr_Site_Date<-droplevels(meta_T_raw.intersect$Species_Nr_Site_Date)

SI(libr_T_raw.intersect)<-meta_T_raw.intersect

libr_T_raw.intersect.mean <- apply(libr_T_raw.intersect, FUN = mean, bySI = "Sp_SiteDate_Phase") #SI is not transfered, so we add it in the following lines
list_T_raw = levels(meta_T_raw.intersect$Sp_SiteDate_Phase) #equals 'Intersection'
index_T_raw<-which(meta_T_raw.intersect$Sp_SiteDate_Phase == list_T_raw[1])
drop<-c("Nr","Species_Nr","Height") #These colums should be ignored because they are not relevant for the mean dataset
meta_T_raw.intersect.mean<-meta_T_raw.intersect[index_T_raw[1], !names(meta_T_raw.intersect) %in% drop] 
for (i in 2:length(list_T_raw)){
  index_T_raw<-which(meta_T_raw.intersect$Sp_SiteDate_Phase == list_T_raw[i])
  meta_T_raw.intersect.mean[i,]<-meta_T_raw.intersect[index_T_raw[1],!names(meta_T_raw.intersect) %in% drop]
}
SI(libr_T_raw.intersect.mean)=meta_T_raw.intersect.mean

#     Table unmixed
idx_int <- SI(libr_T_unmixed.intersect1)$First_measurement %in% Intersection2
libr_T_unmixed.intersect<-libr_T_unmixed.intersect1[idx_int,]
# libr_T_unmixed.intersect<-subset(libr_T_unmixed.intersect1, First_measurement %in% Intersection2)
meta_T_unmixed.intersect<-SI(libr_T_unmixed.intersect)
# libr_T_unmixed.intersect <- subset(library_unmixed, Patch == "Table")
# idx_int <- SI(libr_T_unmixed.intersect)$Sp_SiteDate_Phase %in% Intersection
# libr_T_unmixed.intersect<-libr_T_unmixed.intersect[idx_int,]
# libr_T_unmixed.intersect<-subset(library_unmixed, Patch == "Table" & Sp_SiteDate_Phase %in% Intersection)
meta_T_unmixed.intersect<-SI(libr_T_unmixed.intersect)
meta_T_unmixed.intersect$Species_Cond<-droplevels(meta_T_unmixed.intersect$Species_Cond)
meta_T_unmixed.intersect$Species_Site_Date<-droplevels(meta_T_unmixed.intersect$Species_Site_Date)
meta_T_unmixed.intersect$Sp_SiteDate_Phase<-droplevels(meta_T_unmixed.intersect$Sp_SiteDate_Phase)
meta_T_unmixed.intersect$Species_Nr_Site_Date<-droplevels(meta_T_unmixed.intersect$Species_Nr_Site_Date)

SI(libr_T_unmixed.intersect)<-meta_T_unmixed.intersect

libr_T_unmixed.intersect.mean <- apply(libr_T_unmixed.intersect, FUN = mean, bySI = "Sp_SiteDate_Phase") #SI is not transfered, so we add it in the following lines
list_T_unmixed = levels(meta_T_unmixed.intersect$Sp_SiteDate_Phase) #equals 'Intersection'
index_T_unmixed<-which(meta_T_unmixed.intersect$Sp_SiteDate_Phase == list_T_unmixed[1])
drop<-c("Nr","Species_Nr","Height") #These colums should be ignored because they are not relevant for the mean dataset
meta_T_unmixed.intersect.mean<-meta_T_unmixed.intersect[index_T_unmixed[1], !names(meta_T_unmixed.intersect) %in% drop] 
for (i in 2:length(list_T_unmixed)){
  index_T_unmixed<-which(meta_T_unmixed.intersect$Sp_SiteDate_Phase == list_T_unmixed[i])
  meta_T_unmixed.intersect.mean[i,]<-meta_T_unmixed.intersect[index_T_unmixed[1],!names(meta_T_unmixed.intersect) %in% drop]
}
SI(libr_T_unmixed.intersect.mean)=meta_T_unmixed.intersect.mean

#  C) AVERAGE FIELD and TABLE spectrum per SPECIES AND CONDITION (Patch - Dvpmt phase)
#     This approach results in 1 overall mean spectrum per species and condition (= no intraspecific variation)
#     Needed for evaluations BETWEEN species
#     2 possible approaches for averaging: 
#         - average libr_F.intersect based on Species_Cond or Sp_Phase
#         - average libr_F.intersect.mean based on Species_Cond or Sp_Phase

#     Field
libr_F.intersect.mean.overall<-apply(libr_F.intersect.mean, FUN = mean, bySI = "Species_Cond") # mean = average of mean spectra at each site
libr_F.intersect.mean.overall.ind<-apply(libr_F.intersect, FUN = mean, bySI = "Species_Cond") # mean = average of all individual spectra

list_F.overall = levels(meta_F.intersect.mean$Species_Cond)
index_F.overall<-which(meta_F.intersect.mean$Species_Cond == list_F.overall[1])
drop<-c("Nr","Species_Nr","Height") #These colums should be ignored because they are not relevant for the mean dataset
meta_F.intersect.mean.overall<-meta_F.intersect.mean[index_F.overall[1], !names(meta_F.intersect.mean) %in% drop] 
meta_F.intersect.mean.overall.ind <- meta_F.intersect.mean.overall
for (i in 2:length(list_F.overall)){
  index_F.overall<-which(meta_F.intersect.mean$Species_Cond == list_F.overall[i])
  meta_F.intersect.mean.overall[i,]<-meta_F.intersect.mean[index_F.overall[1],!names(meta_F.intersect.mean) %in% drop]
  meta_F.intersect.mean.overall.ind[i,]<-meta_F.intersect.mean[index_F.overall[1],!names(meta_F.intersect.mean) %in% drop]
}
SI(libr_F.intersect.mean.overall)=meta_F.intersect.mean.overall
SI(libr_F.intersect.mean.overall.ind)=meta_F.intersect.mean.overall.ind

#     Table raw
libr_T_raw.intersect.mean.overall<-apply(libr_T_raw.intersect.mean, FUN = mean, bySI = "Species_Cond")
libr_T_raw.intersect.mean.overall.ind<-apply(libr_T_raw.intersect, FUN = mean, bySI = "Species_Cond")
list_T_raw.overall = levels(meta_T_raw.intersect.mean$Species_Cond)
index_T_raw.overall<-which(meta_T_raw.intersect.mean$Species_Cond == list_T_raw.overall[1])
drop<-c("Nr","Species_Nr","Height") #These colums should be ignored because they are not relevant for the mean dataset
meta_T_raw.intersect.mean.overall<-meta_T_raw.intersect.mean[index_T_raw.overall[1], !names(meta_T_raw.intersect.mean) %in% drop] 
meta_T_raw.intersect.mean.overall.ind <- meta_T_raw.intersect.mean.overall
for (i in 2:length(list_T_raw.overall)){
  index_T_raw.overall<-which(meta_T_raw.intersect.mean$Species_Cond == list_T_raw.overall[i])
  meta_T_raw.intersect.mean.overall[i,]<-meta_T_raw.intersect.mean[index_T_raw.overall[1],!names(meta_T_raw.intersect.mean) %in% drop]
  meta_T_raw.intersect.mean.overall.ind[i,]<-meta_T_raw.intersect.mean[index_T_raw.overall[1],!names(meta_T_raw.intersect.mean) %in% drop]
}
SI(libr_T_raw.intersect.mean.overall)=meta_T_raw.intersect.mean.overall
SI(libr_T_raw.intersect.mean.overall.ind)=meta_T_raw.intersect.mean.overall.ind

#     Table unmixed
libr_T_unmixed.intersect.mean.overall<-apply(libr_T_unmixed.intersect.mean, FUN = mean, bySI = "Species_Cond")
libr_T_unmixed.intersect.mean.overall.ind<-apply(libr_T_unmixed.intersect, FUN = mean, bySI = "Species_Cond")
list_T_unmixed.overall = levels(meta_T_unmixed.intersect.mean$Species_Cond)
index_T_unmixed.overall<-which(meta_T_unmixed.intersect.mean$Species_Cond == list_T_unmixed.overall[1])
drop<-c("Nr","Species_Nr","Height") #These colums should be ignored because they are not relevant for the mean dataset
meta_T_unmixed.intersect.mean.overall<-meta_T_unmixed.intersect.mean[index_T_unmixed.overall[1], !names(meta_T_unmixed.intersect.mean) %in% drop] 
meta_T_unmixed.intersect.mean.overall.ind <- meta_T_unmixed.intersect.mean.overall
for (i in 2:length(list_T_unmixed.overall)){
  index_T_unmixed.overall<-which(meta_T_unmixed.intersect.mean$Species_Cond == list_T_unmixed.overall[i])
  meta_T_unmixed.intersect.mean.overall[i,]<-meta_T_unmixed.intersect.mean[index_T_unmixed.overall[1],!names(meta_T_unmixed.intersect.mean) %in% drop]
  meta_T_unmixed.intersect.mean.overall.ind[i,]<-meta_T_unmixed.intersect.mean[index_T_unmixed.overall[1],!names(meta_T_unmixed.intersect.mean) %in% drop]
}
SI(libr_T_unmixed.intersect.mean.overall)=meta_T_unmixed.intersect.mean.overall
SI(libr_T_unmixed.intersect.mean.overall.ind)=meta_T_unmixed.intersect.mean.overall.ind

#  D) Prepocessing data

# Also store original spectra (= before preprocessing)
libr_F.intersect_orig<-libr_F.intersect
libr_T_raw.intersect_orig<-libr_T_raw.intersect
libr_T_unmixed.intersect_orig<-libr_T_unmixed.intersect

libr_F.intersect.mean_orig<-libr_F.intersect.mean
libr_T_raw.intersect.mean_orig<-libr_T_raw.intersect.mean
libr_T_unmixed.intersect.mean_orig<-libr_T_unmixed.intersect.mean

# Smooth: Apply Savitsky-Golay smoothing to all libraries 
  # libr_F.intersect<- smoothSpeclib(libr_F.intersect,method="sgolay", n=51)
  # libr_F.intersect.mean<- smoothSpeclib(libr_F.intersect.mean,method="sgolay", n=51)
  # libr_T_raw.intersect<- smoothSpeclib(libr_T_raw.intersect,method="sgolay", n=51)
  # libr_T_unmixed.intersect<- smoothSpeclib(libr_T_unmixed.intersect,method="sgolay", n=51)
  # libr_T_raw.intersect.mean<- smoothSpeclib(libr_T_raw.intersect.mean,method="sgolay", n=51)
  # libr_T_unmixed.intersect.mean<- smoothSpeclib(libr_T_unmixed.intersect.mean,method="sgolay", n=51)
  # 
  # libr_F.intersect.mean.overall<- smoothSpeclib(libr_F.intersect.mean.overall,method="sgolay", n=51)
  # libr_T_raw.intersect.mean.overall<- smoothSpeclib(libr_T_raw.intersect.mean.overall,method="sgolay", n=51)
  # libr_T_unmixed.intersect.mean.overall<- smoothSpeclib(libr_T_unmixed.intersect.mean.overall,method="sgolay", n=51)
  # libr_F.intersect.mean.overall.ind<- smoothSpeclib(libr_F.intersect.mean.overall.ind,method="sgolay", n=51)
  # libr_T_raw.intersect.mean.overall.ind<- smoothSpeclib(libr_T_raw.intersect.mean.overall.ind,method="sgolay", n=51)
  # libr_T_unmixed.intersect.mean.overall.ind<- smoothSpeclib(libr_T_unmixed.intersect.mean.overall.ind,method="sgolay", n=51)
  
  library(stringr) #str_sub
  smooth_asd_sv <- function(speclib){
    SI(speclib)$Instrument <- str_sub(SI(speclib)[,1],-3,-1)
    libr_sv<-subset(speclib,Instrument == ".sv")
    libr_asd<-subset(speclib,Instrument == "asd")
    libr_sv<- smoothSpeclib(libr_sv,method="sgolay", n=101)
    libr_asd<- smoothSpeclib(libr_asd,method="sgolay", n=51)
    spectra_sv <- spectra(libr_sv)
    meta_sv <- SI(libr_sv)
    spectra_asd <- spectra(libr_asd)
    meta_asd <- SI(libr_asd)
    spectra_total <- rbind(spectra_asd,spectra_sv)
    wl <- seq(350,2500,1)
    colnames(spectra_total)<-wl
    meta_total <- rbind(meta_asd,meta_sv)
    Refl <- cbind(meta_total, spectra_total)
    Refl <- Refl[order(Refl$ID),] # order the observations according to ID
    meta_ordered <- Refl[,c(1:which(names(Refl)=="Instrument"))]
    speclib_smooth <- speclib(as.matrix(Refl[,-c(1:which(names(Refl)=="Instrument"))]), wl)
    SI(speclib_smooth)<-meta_ordered
    return(speclib_smooth)
  }
  libr_F.intersect<- smooth_asd_sv(libr_F.intersect)
  libr_F.intersect.mean<- smooth_asd_sv(libr_F.intersect.mean)
  libr_T_raw.intersect<- smooth_asd_sv(libr_T_raw.intersect)
  libr_T_unmixed.intersect<- smooth_asd_sv(libr_T_unmixed.intersect)
  libr_T_raw.intersect.mean<- smooth_asd_sv(libr_T_raw.intersect.mean)
  libr_T_unmixed.intersect.mean<- smooth_asd_sv(libr_T_unmixed.intersect.mean)
  
  libr_F.intersect.mean.overall<- smooth_asd_sv(libr_F.intersect.mean.overall)
  libr_T_raw.intersect.mean.overall<- smooth_asd_sv(libr_T_raw.intersect.mean.overall)
  libr_T_unmixed.intersect.mean.overall<- smooth_asd_sv(libr_T_unmixed.intersect.mean.overall)
  libr_F.intersect.mean.overall.ind<- smooth_asd_sv(libr_F.intersect.mean.overall.ind)
  libr_T_raw.intersect.mean.overall.ind<- smooth_asd_sv(libr_T_raw.intersect.mean.overall.ind)
  libr_T_unmixed.intersect.mean.overall.ind<- smooth_asd_sv(libr_T_unmixed.intersect.mean.overall.ind)
  
# 1st and 2nd derivatives
  libr_F_1deriv.intersect <- derivative.speclib(libr_F.intersect, m = 1)
  libr_F_1deriv.intersect.mean <- derivative.speclib(libr_F.intersect.mean, m = 1)
  libr_F_2deriv.intersect <- derivative.speclib(libr_F.intersect, m = 2)
  libr_F_2deriv.intersect.mean <- derivative.speclib(libr_F.intersect.mean, m = 2)
  libr_T_raw_1deriv.intersect <- derivative.speclib(libr_T_raw.intersect, m = 1)
  libr_T_raw_1deriv.intersect.mean <- derivative.speclib(libr_T_raw.intersect.mean, m = 1)
  libr_T_raw_2deriv.intersect <- derivative.speclib(libr_T_raw.intersect, m = 2)
  libr_T_raw_2deriv.intersect.mean <- derivative.speclib(libr_T_raw.intersect.mean, m = 2)
  libr_T_unmixed_1deriv.intersect <- derivative.speclib(libr_T_unmixed.intersect, m = 1)
  libr_T_unmixed_1deriv.intersect.mean <- derivative.speclib(libr_T_unmixed.intersect.mean, m = 1)
  libr_T_unmixed_2deriv.intersect <- derivative.speclib(libr_T_unmixed.intersect, m = 2)
  libr_T_unmixed_2deriv.intersect.mean <- derivative.speclib(libr_T_unmixed.intersect.mean, m = 2)

  
  # Mask atmospheric noise bands
# End of the mask = 2501, because the smoothing added an extra wavelength.
# This 'extended masking' is important when transforming a speclib back into a matrix
  mask(libr_F.intersect) <- c(349,400,1340,1460,1780,1970,2400,2501)
  mask(libr_F.intersect.mean) <- c(349,400,1340,1460,1780,1970,2400,2501)
  mask(libr_T_raw.intersect) <-  c(349,400,1340,1460,1780,1970,2400,2501)
  mask(libr_T_unmixed.intersect) <-  c(349,400,1340,1460,1780,1970,2400,2501)
  mask(libr_T_raw.intersect.mean)<- c(349,400,1340,1460,1780,1970,2400,2501)
  mask(libr_T_unmixed.intersect.mean)<- c(349,400,1340,1460,1780,1970,2400,2501)
  
  mask(libr_F.intersect.mean.overall) <- c(349,400,1340,1460,1780,1970,2400,2501)
  mask(libr_T_raw.intersect.mean.overall) <-  c(349,400,1340,1460,1780,1970,2400,2501)
  mask(libr_T_unmixed.intersect.mean.overall) <-  c(349,400,1340,1460,1780,1970,2400,2501)
  mask(libr_F.intersect.mean.overall.ind) <- c(349,400,1340,1460,1780,1970,2400,2501)
  mask(libr_T_raw.intersect.mean.overall.ind) <-  c(349,400,1340,1460,1780,1970,2400,2501)
  mask(libr_T_unmixed.intersect.mean.overall.ind) <-  c(349,400,1340,1460,1780,1970,2400,2501)
  
  mask(libr_F_1deriv.intersect) <- c(349,400,1340,1460,1780,1970,2400,2501)
  mask(libr_F_1deriv.intersect.mean) <- c(349,400,1340,1460,1780,1970,2400,2501)
  mask(libr_F_2deriv.intersect) <- c(349,400,1340,1460,1780,1970,2400,2501)
  mask(libr_F_2deriv.intersect.mean) <- c(349,400,1340,1460,1780,1970,2400,2501)
  mask(libr_T_raw_1deriv.intersect) <- c(349,400,1340,1460,1780,1970,2400,2501)
  mask(libr_T_raw_1deriv.intersect.mean) <- c(349,400,1340,1460,1780,1970,2400,2501)
  mask(libr_T_raw_2deriv.intersect) <- c(349,400,1340,1460,1780,1970,2400,2501)
  mask(libr_T_raw_2deriv.intersect.mean) <- c(349,400,1340,1460,1780,1970,2400,2501)
  mask(libr_T_unmixed_1deriv.intersect) <- c(349,400,1340,1460,1780,1970,2400,2501)
  mask(libr_T_unmixed_1deriv.intersect.mean) <- c(349,400,1340,1460,1780,1970,2400,2501)
  mask(libr_T_unmixed_2deriv.intersect) <- c(349,400,1340,1460,1780,1970,2400,2501)
  mask(libr_T_unmixed_2deriv.intersect.mean) <- c(349,400,1340,1460,1780,1970,2400,2501)
  
#--------------------------------------------------------------------------------

#####   4) Create 'green' libraries of field and table spectra --> needed for optical traits and ANOVA in "04. Optical trait comparisons"
  
# A) We have to make sure that only those spectra which appear in both the field and table dataset are included --> Determine intersect to make sure we have a paired dataset
 
# Option 1: no preprocessed data 
  List_FieldG<-levels(droplevels(meta_FG$Sp_SiteDate_Phase))
  List_TableG_raw<-levels(droplevels(meta_TG_raw$Sp_SiteDate_Phase))
  List_TableG_unmixed<-levels(droplevels(meta_TG_unmixed$Sp_SiteDate_Phase))
  Intersection_green<-intersect(List_FieldG,List_TableG_unmixed) #raw might have more attributes (in case we were not able to unmix). We choose smallest dataset
  
  libr_FG.intersect_orig2 <- subset(library_raw, Patch == "Field" & DvpPhase == "Green" )
  idx_int <- SI(libr_FG.intersect_orig2)$Sp_SiteDate_Phase %in% Intersection_green
  libr_FG.intersect_orig2<-libr_FG.intersect_orig2[idx_int,]
  # libr_FG.intersect_orig2<-subset(library_raw,Patch == "Field" & DvpPhase == "Green" & Sp_SiteDate_Phase %in% Intersection_green)
  libr_TG_raw.intersect_orig2 <- subset(library_raw, Patch == "Table" & DvpPhase == "Green" )
  idx_int <- SI(libr_TG_raw.intersect_orig2)$Sp_SiteDate_Phase %in% Intersection_green
  libr_TG_raw.intersect_orig2<-libr_TG_raw.intersect_orig2[idx_int,]
  # libr_TG_raw.intersect_orig2<-subset(library_raw,Patch == "Table" & DvpPhase == "Green" & Sp_SiteDate_Phase %in% Intersection_green)
  libr_TG_unmixed.intersect_orig2 <- subset(library_unmixed, Patch == "Table" & DvpPhase == "Green" )
  idx_int <- SI(libr_TG_unmixed.intersect_orig2)$Sp_SiteDate_Phase %in% Intersection_green
  libr_TG_unmixed.intersect_orig2<-libr_TG_unmixed.intersect_orig2[idx_int,]
  # libr_TG_unmixed.intersect_orig2<-subset(library_unmixed,Patch == "Table" & DvpPhase == "Green" & Sp_SiteDate_Phase %in% Intersection_green)

  libr_FG.intersect_orig<-subset(libr_F.intersect_orig, DvpPhase == "Green")
  libr_TG_raw.intersect_orig<-subset(libr_T_raw.intersect_orig, DvpPhase == "Green")
  libr_TG_unmixed.intersect_orig<-subset(libr_T_unmixed.intersect_orig, DvpPhase == "Green")
  meta_FG.intersect_orig<-SI(libr_FG.intersect_orig)
  meta_TG_raw.intersect_orig<-SI(libr_TG_raw.intersect_orig)
  meta_TG_unmixed.intersect_orig<-SI(libr_TG_unmixed.intersect_orig)

  libr_FG.intersect.mean_orig<-subset(libr_F.intersect.mean_orig, DvpPhase == "Green")
  libr_TG_raw.intersect.mean_orig<-subset(libr_T_raw.intersect.mean_orig, DvpPhase == "Green")
  libr_TG_unmixed.intersect.mean_orig<-subset(libr_T_unmixed.intersect.mean_orig, DvpPhase == "Green")
  meta_FG.intersect.mean_orig<-SI(libr_FG.intersect.mean_orig)
  meta_TG_raw.intersect.mean_orig<-SI(libr_TG_raw.intersect.mean_orig)
  meta_TG_unmixed.intersect.mean_orig<-SI(libr_TG_unmixed.intersect.mean_orig)
  
  
# Option 2: preprocess data = smoothed and without noise bands 
# All individual meausurements
  libr_FG.intersect<-subset(libr_F.intersect, DvpPhase == "Green")
  libr_TG_raw.intersect<-subset(libr_T_raw.intersect, DvpPhase == "Green")
  libr_TG_unmixed.intersect<-subset(libr_T_unmixed.intersect, DvpPhase == "Green")
  meta_FG.intersect<-SI(libr_FG.intersect)
  meta_TG_raw.intersect<-SI(libr_TG_raw.intersect)
  meta_TG_unmixed.intersect<-SI(libr_TG_unmixed.intersect)
  
  libr_FG_1deriv.intersect<-subset(libr_F_1deriv.intersect, DvpPhase == "Green")
  libr_FG_2deriv.intersect<-subset(libr_F_2deriv.intersect, DvpPhase == "Green")
  libr_TG_raw_1deriv.intersect<-subset(libr_T_raw_1deriv.intersect, DvpPhase == "Green")
  libr_TG_raw_2deriv.intersect<-subset(libr_T_raw_2deriv.intersect, DvpPhase == "Green")
  libr_TG_unmixed_1deriv.intersect<-subset(libr_T_unmixed_1deriv.intersect, DvpPhase == "Green")
  libr_TG_unmixed_2deriv.intersect<-subset(libr_T_unmixed_2deriv.intersect, DvpPhase == "Green")
  meta_FG_1deriv.intersect<-SI(libr_F_1deriv.intersect)
  meta_FG_2deriv.intersect<-SI(libr_F_2deriv.intersect)
  meta_TG_raw_1deriv.intersect<-SI(libr_TG_raw_1deriv.intersect)
  meta_TG_raw_2deriv.intersect<-SI(libr_TG_raw_2deriv.intersect)
  meta_TG_unmixed_1deriv.intersect<-SI(libr_TG_unmixed_1deriv.intersect)
  meta_TG_unmixed_2deriv.intersect<-SI(libr_TG_unmixed_2deriv.intersect)
  
  # Mean spectrum per species for each Site_Date
  libr_FG.intersect.mean<-subset(libr_F.intersect.mean, DvpPhase == "Green")
  libr_TG_raw.intersect.mean<-subset(libr_T_raw.intersect.mean, DvpPhase == "Green")
  libr_TG_unmixed.intersect.mean<-subset(libr_T_unmixed.intersect.mean, DvpPhase == "Green")
  meta_FG.intersect.mean<-SI(libr_FG.intersect.mean)
  meta_TG_raw.intersect.mean<-SI(libr_TG_raw.intersect.mean)
  meta_TG_unmixed.intersect.mean<-SI(libr_TG_unmixed.intersect.mean)
  
  libr_FG_1deriv.intersect.mean<-subset(libr_F_1deriv.intersect.mean, DvpPhase == "Green")
  libr_FG_2deriv.intersect.mean<-subset(libr_F_2deriv.intersect.mean, DvpPhase == "Green")
  libr_TG_raw_1deriv.intersect.mean<-subset(libr_T_raw_1deriv.intersect.mean, DvpPhase == "Green")
  libr_TG_raw_2deriv.intersect.mean<-subset(libr_T_raw_2deriv.intersect.mean, DvpPhase == "Green")
  libr_TG_unmixed_1deriv.intersect.mean<-subset(libr_T_unmixed_1deriv.intersect.mean, DvpPhase == "Green")
  libr_TG_unmixed_2deriv.intersect.mean<-subset(libr_T_unmixed_2deriv.intersect.mean, DvpPhase == "Green")
  meta_FG_1deriv.intersect.mean<-SI(libr_F_1deriv.intersect.mean)
  meta_FG_2deriv.intersect.mean<-SI(libr_F_2deriv.intersect.mean)
  meta_TG_raw_1deriv.intersect.mean<-SI(libr_TG_raw_1deriv.intersect.mean)
  meta_TG_raw_2deriv.intersect.mean<-SI(libr_TG_raw_2deriv.intersect.mean)
  meta_TG_unmixed_1deriv.intersect.mean<-SI(libr_TG_unmixed_1deriv.intersect.mean)
  meta_TG_unmixed_2deriv.intersect.mean<-SI(libr_TG_unmixed_2deriv.intersect.mean)
  
  
# merge spectra and metadata for saving 
# Original data, without preprocessing
FG.intersect_TableField_orig<-cbind(meta_FG.intersect_orig,spectra(libr_FG.intersect_orig)) 
names(FG.intersect_TableField_orig)[which(names(FG.intersect_TableField_orig) == '1')[1]:ncol(FG.intersect_TableField_orig)]<-as.character(paste("nm",wavelength(libr_FG.intersect_orig),sep=''))
TG_raw.intersect_TableField_orig<-cbind(meta_TG_raw.intersect_orig,spectra(libr_TG_raw.intersect_orig)) 
names(TG_raw.intersect_TableField_orig)[which(names(TG_raw.intersect_TableField_orig) == '1')[1]:ncol(TG_raw.intersect_TableField_orig)]<-as.character(paste("nm",wavelength(libr_TG_raw.intersect_orig),sep=''))
TG_unmixed.intersect_TableField_orig<-cbind(meta_TG_unmixed.intersect_orig,spectra(libr_TG_unmixed.intersect_orig)) 
names(TG_unmixed.intersect_TableField_orig)[which(names(TG_unmixed.intersect_TableField_orig) == '1')[1]:ncol(TG_unmixed.intersect_TableField_orig)]<-as.character(paste("nm",wavelength(libr_TG_unmixed.intersect_orig),sep=''))
  
FG.intersect.mean_TableField_orig<-cbind(meta_FG.intersect.mean_orig,spectra(libr_FG.intersect.mean_orig)) 
names(FG.intersect.mean_TableField_orig)[which(names(FG.intersect.mean_TableField_orig) == '1')[1]:ncol(FG.intersect.mean_TableField_orig)]<-as.character(paste("nm",wavelength(libr_FG.intersect.mean_orig),sep=''))
TG_raw.intersect.mean_TableField_orig<-cbind(meta_TG_raw.intersect.mean_orig,spectra(libr_TG_raw.intersect.mean_orig)) 
names(TG_raw.intersect.mean_TableField_orig)[which(names(TG_raw.intersect.mean_TableField_orig) == '1')[1]:ncol(TG_raw.intersect.mean_TableField_orig)]<-as.character(paste("nm",wavelength(libr_TG_raw.intersect.mean_orig),sep=''))
TG_unmixed.intersect.mean_TableField_orig<-cbind(meta_TG_unmixed.intersect.mean_orig,spectra(libr_TG_unmixed.intersect.mean_orig)) 
names(TG_unmixed.intersect.mean_TableField_orig)[which(names(TG_unmixed.intersect.mean_TableField_orig) == '1')[1]:ncol(TG_unmixed.intersect.mean_TableField_orig)]<-as.character(paste("nm",wavelength(libr_TG_raw.intersect.mean_orig),sep=''))
  
# Preprocessed data: smoothing and removal of noise bands
FG.intersect_TableField<-cbind(meta_FG.intersect,spectra(libr_FG.intersect)) 
names(FG.intersect_TableField)[which(names(FG.intersect_TableField) == '1')[1]:ncol(FG.intersect_TableField)]<-as.character(paste("nm",wavelength(libr_FG.intersect),sep=''))
TG_raw.intersect_TableField<-cbind(meta_TG_raw.intersect,spectra(libr_TG_raw.intersect)) 
names(TG_raw.intersect_TableField)[which(names(TG_raw.intersect_TableField) == '1')[1]:ncol(TG_raw.intersect_TableField)]<-as.character(paste("nm",wavelength(libr_TG_raw.intersect),sep=''))
TG_unmixed.intersect_TableField<-cbind(meta_TG_unmixed.intersect,spectra(libr_TG_unmixed.intersect)) 
names(TG_unmixed.intersect_TableField)[which(names(TG_unmixed.intersect_TableField) == '1')[1]:ncol(TG_unmixed.intersect_TableField)]<-as.character(paste("nm",wavelength(libr_TG_unmixed.intersect),sep=''))
 
FG.intersect.mean_TableField<-cbind(meta_FG.intersect.mean,spectra(libr_FG.intersect.mean)) 
names(FG.intersect.mean_TableField)[which(names(FG.intersect.mean_TableField) == '1')[1]:ncol(FG.intersect.mean_TableField)]<-as.character(paste("nm",wavelength(libr_FG.intersect.mean),sep=''))
TG_raw.intersect.mean_TableField<-cbind(meta_TG_raw.intersect.mean,spectra(libr_TG_raw.intersect.mean)) 
names(TG_raw.intersect.mean_TableField)[which(names(TG_raw.intersect.mean_TableField) == '1')[1]:ncol(TG_raw.intersect.mean_TableField)]<-as.character(paste("nm",wavelength(libr_TG_raw.intersect.mean),sep=''))
TG_unmixed.intersect.mean_TableField<-cbind(meta_TG_unmixed.intersect.mean,spectra(libr_TG_unmixed.intersect.mean)) 
names(TG_unmixed.intersect.mean_TableField)[which(names(TG_unmixed.intersect.mean_TableField) == '1')[1]:ncol(TG_unmixed.intersect.mean_TableField)]<-as.character(paste("nm",wavelength(libr_TG_raw.intersect.mean),sep=''))
  

#--------------------------------------------------------------------------------
# Create list of all needed data and  store it 

Black_datasets <- list(
          # data needed in 01b. and 04.
          "libr_FG.intersect" = libr_FG.intersect, # intersection of field data with table data (based on Sp_SiteDate_Phase)
          "libr_TG_raw.intersect" = libr_TG_raw.intersect,
          "libr_TG_unmixed.intersect" = libr_TG_unmixed.intersect,
          
          "libr_FG.intersect.mean" = libr_FG.intersect.mean, # mean of intersection data bySI = "Sp_SiteDate_Phase"
          "libr_TG_raw.intersect.mean" = libr_TG_raw.intersect.mean,
          "libr_TG_unmixed.intersect.mean" = libr_TG_unmixed.intersect.mean,
          
          # data needed in 02.
          # (
          "meta_FG" = meta_FG,
          "meta_TG_raw" = meta_TG_raw,
          "meta_TG_unmixed" = meta_TG_unmixed,
          "library_FG" = library_FG,
          "library_TG_raw" = library_TG_raw,
          "library_TG_unmixed" = library_TG_unmixed,
                 
          "meta_FF" = meta_FF,
          "meta_TF_raw" = meta_TF_raw,
          "meta_TF_unmixed" = meta_TF_unmixed,
          "library_FF" = library_FF,
          "library_TF_raw" = library_TF_raw,
          "library_TF_unmixed" = library_TF_unmixed, # )
          
          "libr_F.intersect_orig" = libr_F.intersect_orig,
          "libr_T_raw.intersect_orig" = libr_T_raw.intersect_orig,
          "libr_T_unmixed.intersect_orig" = libr_T_unmixed.intersect_orig,
          
          "libr_F.intersect.mean_orig" = libr_F.intersect.mean_orig,
          "libr_T_raw.intersect.mean_orig" = libr_T_raw.intersect.mean_orig,
          "libr_T_unmixed.intersect.mean_orig" = libr_T_unmixed.intersect.mean_orig,
          
          # data needed in 03.
          "libr_F.intersect" = libr_F.intersect,
          "libr_T_raw.intersect" = libr_T_raw.intersect,
          "libr_T_unmixed.intersect" = libr_T_unmixed.intersect,
          
          "libr_F.intersect.mean" = libr_F.intersect.mean,
          "libr_T_raw.intersect.mean" = libr_T_raw.intersect.mean,
          "libr_T_unmixed.intersect.mean" = libr_T_unmixed.intersect.mean,
          
          "libr_F.intersect.mean.overall" = libr_F.intersect.mean.overall, # mean of intersect.mean data bySI = "Species_Cond"
          "libr_T_raw.intersect.mean.overall" = libr_T_raw.intersect.mean.overall,
          "libr_T_unmixed.intersect.mean.overall" = libr_T_unmixed.intersect.mean.overall,
          "libr_F.intersect.mean.overall.ind" = libr_F.intersect.mean.overall.ind, # mean of intersect data bySI = "Species_Cond"
          "libr_T_raw.intersect.mean.overall.ind" = libr_T_raw.intersect.mean.overall.ind,
          "libr_T_unmixed.intersect.mean.overall.ind" = libr_T_unmixed.intersect.mean.overall.ind,
          
          # datasets instead of libraries
          "FG.intersect_TableField_orig" = FG.intersect_TableField_orig,
          "TG_raw.intersect_TableField_orig" = TG_raw.intersect_TableField_orig,
          "TG_unmixed.intersect_TableField_orig" = TG_unmixed.intersect_TableField_orig,
          "FG.intersect.mean_TableField_orig"  = FG.intersect.mean_TableField_orig,
          "TG_raw.intersect.mean_TableField_orig" = TG_raw.intersect.mean_TableField_orig,
          "TG_unmixed.intersect.mean_TableField_orig" = TG_unmixed.intersect.mean_TableField_orig,
           
          "FG.intersect_TableField" = FG.intersect_TableField,
          "TG_raw.intersect_TableField" = TG_raw.intersect_TableField,
          "TG_unmixed.intersect_TableField" = TG_unmixed.intersect_TableField,
          "FG.intersect.mean_TableField"  = FG.intersect.mean_TableField,
          "TG_raw.intersect.mean_TableField" = TG_raw.intersect.mean_TableField,
          "TG_unmixed.intersect.mean_TableField" = TG_unmixed.intersect.mean_TableField)

setwd(wd3)
saveRDS(Black_datasets,file=paste(paste("datasets_",fracperc,sep=jumpcorr),".rds",sep=""))

#--------------------------------------------------------------------------------
# Save metadata to folder of Black table letter: afterwards add number of pictures manually 
setwd(wd4)
write.table(meta_F.intersect,file=paste(paste("meta_F.intersect_",fracperc,sep=jumpcorr),".csv",sep=""))
write.table(meta_T_raw.intersect,file=paste(paste("meta_T_raw.intersect_",fracperc,sep=jumpcorr),".csv",sep=""))
write.table(meta_T_unmixed.intersect,file=paste(paste("meta_T_unmixed.intersect_",fracperc,sep=jumpcorr),".csv",sep=""))

return(Black_datasets)

}

##################################################################################
# Apply the function to datasets with different jump corrections

wd1 =  "C:/Users/u0091812/Box Sync/DATA/Hyperspectral data/Data with metadata/3. Aggregated spectra/"
wd2 = "C:/Users/u0091812/Box Sync/DATA/Hyperspectral data/Data with metadata/4. Signal unmixed spectra/"
wd3 = "C:/Users/u0091812/Box Sync/03. Black table procedure/Data/Spectral data"
wd4 = "C:/Users/u0091812/Box Sync/03. Black table procedure/Data/Spectral metadata"

Datasets_Nocorr <-Black_data(wd1,wd2,wd3,wd4,"Nocorr", fracperc = "",remove="1")
# Datasets_jcorr_add_SWIR1<-Black_data(wd1,wd2,wd3,wd4,"jcorr_add_SWIR1", fracperc = "",remove="1")
# Datasets_jcorr_mult_SWIR1<-Black_data(wd1,wd2,wd3,wd4,"jcorr_mult_SWIR1", fracperc = "",remove="1")
# Datasets_jcorr_mult_SWIR12<-Black_data(wd1,wd2,wd3,wd4,"jcorr_mult_SWIR12", fracperc = "",remove="1")
# Datasets_jcorr_mult_VNIR<-Black_data(wd1,wd2,wd3,wd4,"jcorr_mult_VNIR", fracperc = "",remove="1")
# Datasets_Nocorr_0.9 <-Black_data(wd1,wd2,wd3,wd4,"Nocorr", fracperc = "_frac0.9",remove="1")
