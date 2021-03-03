###############################################################
###                                                        ###
###   01b. Paired dataset creation - field - table - lab   ###
###                                                        ###
##############################################################

### The objective of this script is to create datasets in which only those records which appear in both the field and table hyperspectral and trait dataset are kept.

# We call the obtained datasets 'paired dataset' because they are identical regarding Species_SiteDate. 

### The input datasets are the average spectra of the bouquets/patches of either
# - the raw spectra
# - the unmixed spectra
# -> these data were averaged per site and date in "O1. Black table"

# In addition we need the species x traits dataset

# Script by Elisa Van Cleemput, 2017-2018
#--------------------------------------------------------------------------------

# Load spectral data: use script "01. Black table data" (We only need the green measurements) or open .rds data
paired_metadata <- function (wd1,wd2,wd3,jumpcorr) {

setwd(wd1)
datasets<-readRDS(paste("datasets_",".rds",sep=jumpcorr))

meta_FG<-SI(datasets$libr_FG.intersect)
meta_FG.mean<-SI(datasets$libr_FG.intersect.mean)
meta_TG_raw<-SI(datasets$libr_TG_raw.intersect)
meta_TG_raw.mean<-SI(datasets$libr_TG_raw.intersect.mean)
meta_TG_unmixed<-SI(datasets$libr_TG_unmixed.intersect)
meta_TG_unmixed.mean<-SI(datasets$libr_TG_unmixed.intersect.mean)

#--------------------------------------------------------------------------------

library(xlsx)

#--------------------------------------------------------------------------------
# Choose dataset for analyses and load it from the directory (dependent on requirements of further analyses)

# Intraspecific variation in plant functional traits
setwd(wd2)
PFT_intra_all<-read.xlsx("SpeciesxTraits_v11.xlsx", 5)
rownames(PFT_intra_all)<-PFT_intra_all[,4]
PFT_intra_all[PFT_intra_all==0]<-NA
PFT_intra_all[PFT_intra_all=="#N/A"]<-NA
PFT_intra_all<-PFT_intra_all[,1:which(names(PFT_intra_all) == 'Ni')]
PFT_Prosail<-PFT_intra_all[,c("Site","Species.code","Site_Species_Date","SLA","EWT","Chltotal_area","Carot_area")]
cols.dont.want<-c("SSD","Seed_No","Seed_weight")
PFT_intra <- PFT_intra_all[, ! names(PFT_intra_all) %in% cols.dont.want, drop = F]
PFT_intra<-PFT_intra[complete.cases(PFT_intra),] # remove all columns with NAN


# Optiopnal ---
# Only keep the first measurements on each site (v11): 102 --> 76 records (tabblad 5) 
PFT_Prosail<-PFT_Prosail[complete.cases(PFT_Prosail),] # remove all columns with NAN

# PFT_meta<-PFT_Prosail[,1:3]
# PFT_traits<-PFT_Prosail[,4:7]

PFT_meta <- PFT_intra[,1:6]
PFT_traits <- PFT_intra[,7:ncol(PFT_intra)]

#--------------------------------------------------------------------------------

# Make sure all datasets have the same terminology.
# This will ensure us to make parallel figures, and to run a general script for both datasets

colnames(meta_FG)[which(names(meta_FG) == "Species")] <- "Species.code"
colnames(meta_FG.mean)[which(names(meta_FG.mean) == "Species")] <- "Species.code"
colnames(meta_TG_raw)[which(names(meta_TG_raw) == "Species")] <- "Species.code"
colnames(meta_TG_raw.mean)[which(names(meta_TG_raw.mean) == "Species")] <- "Species.code"
colnames(meta_TG_unmixed)[which(names(meta_TG_unmixed) == "Species")] <- "Species.code"
colnames(meta_TG_unmixed.mean)[which(names(meta_TG_unmixed.mean) == "Species")] <- "Species.code"

levels(PFT_meta$Site)[levels(PFT_meta$Site) == "BD"] <- "BattSD"
levels(PFT_meta$Site)[levels(PFT_meta$Site) == "BL"] <- "BattSP"
levels(PFT_meta$Site)[levels(PFT_meta$Site) == "BN"] <- "BattN"
levels(PFT_meta$Site)[levels(PFT_meta$Site) == "DI"] <- "DoBeI"
levels(PFT_meta$Site)[levels(PFT_meta$Site) == "PI"] <- "PapenIn"
levels(PFT_meta$Site)[levels(PFT_meta$Site) == "DS"] <- "DoBeS"
levels(PFT_meta$Site)[levels(PFT_meta$Site) == "GK"] <- "DeMatCem"
levels(PFT_meta$Site)[levels(PFT_meta$Site) == "GV"] <- "DeMatSpo"
levels(PFT_meta$Site)[levels(PFT_meta$Site) == "HB"] <- "HaachtBev"
levels(PFT_meta$Site)[levels(PFT_meta$Site) == "HH"] <- "HaachtHol"
levels(PFT_meta$Site)[levels(PFT_meta$Site) == "HK"] <- "HeverleeKapel"

library(stringr)
right = function(text, num_char) {
  substr(text, nchar(as.character(text)) - (num_char-1), nchar(as.character(text)))
}
PFT_meta$Date<-right(PFT_meta$Site_Species_Date,8)
PFT_meta$ID<-paste(PFT_meta$Species.code,PFT_meta$Site,as.factor(PFT_meta$Date))


# Combine the modified metadata with the trait values again for use in later analyses
PFT<-cbind(PFT_meta,PFT_traits)
# spec_FG<-cbind(meta_FG,spectra_FG)

# meta_TG_raw and meta_TG_unmixed contain exactly the same records.
# They also containt the same records as meta_FG (this was done in "01. Black table")
# This can be verified with following statements:
  # test<-intersect(meta_TG_raw$ID,meta_TG_unmixed$ID)
  # checktest<-PFT[!(meta_TG_raw$ID %in% test), ]
  # checktest2<-PFT[!(meta_TG_unmixed$ID %in% test), ]
  # test2<-intersect(meta_FG$ID,meta_TG_unmixed$ID)
  # checktest3<-PFT[!(meta_FG$ID %in% test2), ]
  # checktest4<-PFT[!(meta_TG_raw$ID %in% test2), ]

  # test.mean<-intersect(meta_TG_raw.mean$ID,meta_TG_unmixed.mean$ID)
  # checktest.mean<-PFT[!(meta_TG_raw.mean$ID %in% test.mean), ]
  # checktest2.mean<-PFT[!(meta_TG_unmixed.mean$ID %in% test.mean), ]
  # test2.mean<-intersect(meta_FG.mean$ID,meta_TG_unmixed.mean$ID)
  # checktest3.mean<-PFT[!(meta_FG.mean$ID %in% test2.mean), ]
  # checktest4.mean<-PFT[!(meta_TG_raw.mean$ID %in% test2.mean), ]

#--------------------------------------------------------------------------------

# Determine which records coincide in all datasets

PFT$ID<-factor(PFT$ID)
levels(PFT$ID)[levels(PFT$ID) == "EQUPAL BattSP 20160719"] <- "EQUPAL BattSD 20160719"
levels(PFT$ID)[levels(PFT$ID) == "AGRSP BattSP 20160719"] <- "AGRSP BattSD 20160719"
levels(PFT$ID)[levels(PFT$ID) == "HOLLAN BattSP 20160719"] <- "HOLLAN BattSD 20160719"
levels(PFT$ID)[levels(PFT$ID) == "RANREP BattSP 20160719"] <- "RANREP BattSD 20160719"
levels(PFT$ID)[levels(PFT$ID) == "SYMOFF HeverleeKapel 20160813"] <- "SYMOFF HeverleeKapel 20160816"
levels(PFT$ID)[levels(PFT$ID) == "HYPPER HeverleeKapel 20160813"] <-  "HYPPER HeverleeKapel 20160816"

Intersection<-intersect(PFT$ID,meta_FG$ID) # should be the same as intersection with TG_raw or TG_unmixed, all measurements or means

# The datasets have 14 rows in common 
# Check which functional measurements do not have a specral equivalent on the field and table
check<-PFT[!(PFT$ID %in% Intersection), ]
# Check which field spectra do not have a functional measurement equivalent
check2<-levels(as.factor(meta_FG[!(meta_FG$ID %in% Intersection),c("ID") ]))
# 4 spectral measurements without trait data --> these were taken in teh second measuring round (august)

PFT_paired<-PFT[PFT$ID %in% Intersection, ]
PFT_paired<-PFT_paired[order(PFT_paired$ID),]
library_FG_paired<-subset(datasets$libr_FG.intersect,ID %in% Intersection)

meta_FG.mean_paired<-meta_FG.mean[meta_FG.mean$ID %in% Intersection, ]
meta_FG.mean_paired<-meta_FG.mean_paired[order(meta_FG.mean_paired$ID),]
library_FG.mean_paired<-subset(datasets$libr_FG.intersect.mean,ID %in% Intersection)

meta_TG_raw_paired<-meta_TG_raw[meta_TG_raw$ID %in% Intersection, ]
meta_TG_raw_paired<-meta_TG_raw_paired[order(meta_TG_raw_paired$ID),]
library_TG_raw_paired<-subset(datasets$libr_TG_raw.intersect,ID %in% Intersection)

meta_TG_raw.mean_paired<-meta_TG_raw.mean[meta_TG_raw.mean$ID %in% Intersection, ]
meta_TG_raw.mean_paired<-meta_TG_raw.mean_paired[order(meta_TG_raw.mean_paired$ID),]
library_TG_raw.mean_paired<-subset(datasets$libr_TG_raw.intersect.mean,ID %in% Intersection)

meta_TG_unmixed_paired<-meta_TG_unmixed[meta_TG_unmixed$ID %in% Intersection, ]
meta_TG_unmixed_paired<-meta_TG_unmixed_paired[order(meta_TG_unmixed_paired$ID),]
library_TG_unmixed_paired<-subset(datasets$libr_TG_unmixed.intersect,ID %in% Intersection)

meta_TG_unmixed.mean_paired<-meta_TG_unmixed.mean[meta_TG_unmixed.mean$ID %in% Intersection, ]
meta_TG_unmixed.mean_paired<-meta_TG_unmixed.mean_paired[order(meta_TG_unmixed.mean_paired$ID),]
library_TG_unmixed.mean_paired<-subset(datasets$libr_TG_unmixed.intersect.mean,ID %in% Intersection)

Paired_libraries <- list ("library_FG_paired" = library_FG_paired, 
                          "library_FG.mean_paired" = library_FG.mean_paired,
                          "library_TG_raw_paired" = library_TG_raw_paired,
                          "library_TG_raw.mean_paired" = library_TG_raw.mean_paired,
                          "library_TG_unmixed_paired" = library_TG_unmixed_paired,
                          "library_TG_unmixed.mean_paired" = library_TG_unmixed.mean_paired)
  
#--------------------------------------------------------------------------------
# Save the twin datasets

setwd(wd3)

write.table(PFT_paired,file="PFT_paired.csv",sep=";",row.names = F)
write.table(Intersection,file="Intersection_Field_Table_Lab.csv",sep=";",row.names = F)

saveRDS(Paired_libraries,file=paste("Spectra_paired_",".rds",sep=jumpcorr))

}

wd1 = "C:/Users/u0091812/Box Sync/03. Black table procedure/Data/Spectral data"
wd2 = "C:/Users/u0091812/Box Sync/DATA/Species x traits/Species x traits/"
wd3 = "C:/Users/u0091812/Box Sync/03. Black table procedure/Data/Paired Field-Table-Lab list"

paired_metadata(wd1,wd2,wd3,"Nocorr")

