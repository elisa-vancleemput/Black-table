#######################################################
###                                                 ###
###          02. Spectral signature plot            ###
###                                                 ###
#######################################################

### The objective of this script is to plot spectral signatures of species meausured in the field and on the table.
# Hereby we aim to visually compare the various measurement methods.

# More specifically we want to compare following measurements:
#   a) unmixed table spectra vs. raw table spectra
#   b) field spectra vs.(unmixed) table spectra --> do the unmixed spectra represent the field better than raw table spectra?

### The input of these steps are the libraries created in "01. Black table data":
#   = combinations of 
#         - field OR table_raw OR table_unmixed; and
#         - Green OR Flower

# For each species we will create separate plots
# The script consists of 3 parts:
#   1) specify Species and Development phase (Green or Flower) for plotting
#       + generate plot specifications
#   2) Plot mean spectrum for each Site_Date: field and unmixed spectra - field and raw spectra (objective b)
#   3) Plot all individual spectra taken at each Site_Date

# Script by Elisa Van Cleemput, 2018
#--------------------------------------------------------------------------------

# Clean workspace

rm(list=ls())

#--------------------------------------------------------------------------------

# Load libraries
library(hsdar)
library(colorRamps)

#--------------------------------------------------------------------------------
# Figure per species: colour indicates Site_Date, line type indicates field vs. table
# Make this combination separately for raw and unmixed table spectra, or include the three options by three different line types
# Make this figure separately for Flower and flower spectra

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

spectral_plots <- function (datasets, wd1, wd2, wd3, jumpcorr, SpeciesName, DvpPhase, VIbands, VIbandpos) {

# 1) Specify species and datasets for visualisation (development phase: Green or Flower)
# if (DvpPhase == "Green"){
#   dataset_Field<-datasets$meta_FG
#   dataset_Table_unmixed<-datasets$meta_TG_unmixed
#   dataset_Table_raw<-datasets$meta_TG_raw
#   Libr1a<-datasets$library_FG
#   Libr2a<-datasets$library_TG_unmixed
#   Libr3a<-datasets$library_TG_raw
# } else if (DvpPhase == "Flower"){
#   dataset_Field<-datasets$meta_FF
#   dataset_Table_unmixed<-datasets$meta_TF_unmixed
#   dataset_Table_raw<-datasets$meta_TF_raw
#   Libr1a<-datasets$library_FF
#   Libr2a<-datasets$library_TF_unmixed
#   Libr3a<-datasets$library_TF_raw
# }
  
  Libr1a<-datasets$libr_F.intersect_orig
  Libr2a<-datasets$libr_T_unmixed.intersect_orig
  Libr3a<-datasets$libr_T_raw.intersect_orig
  if (DvpPhase == "Green"){
    Libr1a<-subset(Libr1a, DvpPhase == "Green")
    Libr2a<-subset(Libr2a, DvpPhase == "Green")
    Libr3a<-subset(Libr3a, DvpPhase == "Green")
  } else if (DvpPhase == "Flower"){
    Libr1a<-subset(Libr1a, DvpPhase == "Flower")
    Libr2a<-subset(Libr2a, DvpPhase == "Flower")
    Libr3a<-subset(Libr3a, DvpPhase == "Flower")
  }
  dataset_Field<-SI(Libr1a)
  dataset_Table_unmixed<-SI(Libr2a)
  dataset_Table_raw<-SI(Libr3a)
  
  # Mean per day:
  Libr1b<-datasets$libr_F.intersect.mean_orig
  Libr2b<-datasets$libr_T_unmixed.intersect.mean_orig
  Libr3b<-datasets$libr_T_raw.intersect.mean_orig
  if (DvpPhase == "Green"){
    Libr1mean<-subset(Libr1b, DvpPhase == "Green")
    Libr2mean<-subset(Libr2b, DvpPhase == "Green")
    Libr3mean<-subset(Libr3b, DvpPhase == "Green")
  } else if (DvpPhase == "Flower"){
    Libr1mean<-subset(Libr1b, DvpPhase == "Flower")
    Libr2mean<-subset(Libr2b, DvpPhase == "Flower")
    Libr3mean<-subset(Libr3b, DvpPhase == "Flower")
  }
  dataset_Field.mean<-SI(Libr1mean)
  dataset_Table_unmixed.mean<-SI(Libr2mean)
  dataset_Table_raw.mean<-SI(Libr3mean)
  
  
  # 2) Make a list with all Site_Date for which the specified species contains data in the field dataset
  # meta_Field_Species<-dataset_Field[which(dataset_Field$Species==SpeciesName),]
  meta_Field_Species<-eval(substitute(dataset_Field[which(dataset_Field$Species==SpeciesName),], list(s = SpeciesName)))
  List_Site_Date_Field<-levels(droplevels(meta_Field_Species$Site_Date))

# 3) Check if all these attribtutes have a table equivalent. Determine intersect to make sure we have a paired dataset
  # meta_Table_Species<-dataset_Table_unmixed[which(dataset_Table_unmixed$Species==SpeciesName),] #raw might have more attributes (in case we were not able to unmix). We choose smallest dataset
  meta_Table_Species<-eval(substitute(dataset_Table_unmixed[which(dataset_Table_unmixed$Species==SpeciesName),], list(s = SpeciesName)))
  List_Site_Date_Table<-levels(droplevels(meta_Table_Species$Site_Date))
  
  List_Site_Date<-intersect(List_Site_Date_Field,List_Site_Date_Table)

# 4) Determine the number of records per Site_Date attribute
  # dataset_Field_subset<-dataset_Field[dataset_Field$Species==SpeciesName & dataset_Field$Site_Date %in% List_Site_Date,]
  # dataset_Table_unmixed_subset<-dataset_Table_unmixed[dataset_Table_unmixed$Species==SpeciesName & dataset_Table_unmixed$Site_Date %in% List_Site_Date,]
  # dataset_Table_raw_subset<-dataset_Table_raw[dataset_Table_raw$Species==SpeciesName & dataset_Table_raw$Site_Date %in% List_Site_Date,]
  dataset_Field_subset<-eval(substitute(dataset_Field[which(dataset_Field$Species==SpeciesName),], list(s = SpeciesName)))
  idx_int<-dataset_Field_subset$Site_Date %in% List_Site_Date
  dataset_Field_subset<-dataset_Field_subset[idx_int,]
  dataset_Table_unmixed_subset<-eval(substitute(dataset_Table_unmixed[which(dataset_Table_unmixed$Species==SpeciesName),], list(s = SpeciesName)))
  idx_int<-dataset_Table_unmixed_subset$Site_Date %in% List_Site_Date
  dataset_Table_unmixed_subset<-dataset_Table_unmixed_subset[idx_int,]
  dataset_Table_raw_subset<-eval(substitute(dataset_Table_raw[which(dataset_Table_raw$Species==SpeciesName),], list(s = SpeciesName)))
  idx_int<-dataset_Table_raw_subset$Site_Date %in% List_Site_Date
  dataset_Table_raw_subset<-dataset_Table_raw_subset[idx_int,]
  
  library(dplyr)
  dataset_Field_subset.count<-count(dataset_Field_subset,Site_Date)
  dataset_Table_unmixed_subset.count<-count(dataset_Table_unmixed_subset,Site_Date)
  dataset_Table_raw_subset.count<-count(dataset_Table_raw_subset,Site_Date)
  
  dataset_Field_subset.totalcount<-count(dataset_Field_subset,Species)
  dataset_Table_unmixed_subset.totalcount<-count(dataset_Table_unmixed_subset,Species)
  dataset_Table_raw_subset.totalcount<-count(dataset_Table_raw_subset,Species)
  
  
  dataset_Field_subset.mean<-eval(substitute(dataset_Field.mean[which(dataset_Field.mean$Species==SpeciesName),], list(s = SpeciesName)))
  idx_int<-dataset_Field_subset.mean$Site_Date %in% List_Site_Date
  dataset_Field_subset.mean<-dataset_Field_subset.mean[idx_int,]
  dataset_Table_unmixed_subset.mean<-eval(substitute(dataset_Table_unmixed.mean[which(dataset_Table_unmixed.mean$Species==SpeciesName),], list(s = SpeciesName)))
  idx_int<-dataset_Table_unmixed_subset.mean$Site_Date %in% List_Site_Date
  dataset_Table_unmixed_subset.mean<-dataset_Table_unmixed_subset.mean[idx_int,]
  dataset_Table_raw_subset.mean<-eval(substitute(dataset_Table_raw.mean[which(dataset_Table_raw.mean$Species==SpeciesName),], list(s = SpeciesName)))
  idx_int<-dataset_Table_raw_subset.mean$Site_Date %in% List_Site_Date
  dataset_Table_raw_subset.mean<-dataset_Table_raw_subset.mean[idx_int,]
  
  dataset_Field_subset.totalcount.mean<-count(dataset_Field_subset.mean,Species)
  dataset_Table_unmixed_subset.totalcount.mean<-count(dataset_Table_unmixed_subset.mean,Species)
  dataset_Table_raw_subset.totalcount.mean<-count(dataset_Table_raw_subset.mean,Species)

# 5) Determine the index of the measurements fitting the above specified criteria 
  index_Field<-which(dataset_Field$Species==SpeciesName & dataset_Field$Site_Date %in% List_Site_Date)
  dataset_Field_subset$index<-index_Field
  index_Table_unmixed<-which(dataset_Table_unmixed$Species==SpeciesName & dataset_Table_unmixed$Site_Date %in% List_Site_Date)
  dataset_Table_unmixed_subset$index<-index_Table_unmixed
  index_Table_raw<-which(dataset_Table_raw$Species==SpeciesName & dataset_Table_raw$Site_Date %in% List_Site_Date)
  dataset_Table_raw_subset$index<-index_Table_raw

# 6) Preprocess the libraries to be plotted
  # Smooth: Apply Savitsky-Golay smoothing to all libraries 
  # Libr1<- smoothSpeclib(Libr1a,method="sgolay", n=51)
  # Libr2<- smoothSpeclib(Libr2a,method="sgolay", n=51)
  # Libr3<- smoothSpeclib(Libr3a,method="sgolay", n=51)
  # Libr1mean<- smoothSpeclib(Libr1mean,method="sgolay", n=51)
  # Libr2mean<- smoothSpeclib(Libr2mean,method="sgolay", n=51)
  # Libr3mean<- smoothSpeclib(Libr3mean,method="sgolay", n=51)
  
  Libr1<- smooth_asd_sv(Libr1a)
  Libr2<- smooth_asd_sv(Libr2a)
  Libr3<- smooth_asd_sv(Libr3a)
  Libr1mean<- smooth_asd_sv(Libr1mean)
  Libr2mean<- smooth_asd_sv(Libr2mean)
  Libr3mean<- smooth_asd_sv(Libr3mean)
  
  
  # Mask atmospheric noise bands
  mask(Libr1) <-  c(349,400,1340,1460,1780,1970,2400,2501)
  mask(Libr2) <-  c(349,400,1340,1460,1780,1970,2400,2501)
  mask(Libr3) <-  c(349,400,1340,1460,1780,1970,2400,2501)

  mask(Libr1mean) <-  c(349,400,1340,1460,1780,1970,2400,2501)
  mask(Libr2mean) <-  c(349,400,1340,1460,1780,1970,2400,2501)
  mask(Libr3mean) <-  c(349,400,1340,1460,1780,1970,2400,2501)
  
#--------------------------------------------------------------------------------

########### Plot mean spectrum for each Site_Date: field and unmixed spectra OR field and raw spectra
if (VIbands == "No") {
  
  # Define color scheme for plotting
  colors = c("orange","pink","red","blue","gray","chocolate","blueviolet","darkorange","wheat","black","Flower","yellow","darkorchid","gold")
  # library("RColorBrewer")
  # colors = brewer.pal(12, "Paired")
  # colors=brewer.pal(6,"Dark2")
  colors = colorRamps::primary.colors(20)
  
  setwd(wd1)
jpeg(paste(paste(paste(SpeciesName,DvpPhase,sep="_"),jumpcorr,sep="_Field+Unmixed_"),".jpeg",sep=""),res=300,width=14,height=7,units='in')
# jpeg(paste(paste(SpeciesName,DvpPhase,sep="_"),".jpeg",sep="_Field+Unmixed"),res=300,width=14,height=7,units='in')
# jpeg("SOLGIG_Flower_Mean_Field+Unmixed.jpeg",res=300,width=14,height=7,units='in')
# x11(width=1280,height=720)
plot(c(),c(),axes=FALSE,ylim=c(0,1.0),xlim=c(350,2500),ylab='',xlab='', xaxs="i", yaxs="i")
for (i in 1:length(List_Site_Date)){
  Libr1_sp<-Libr1[SI(Libr1)$Species == SpeciesName,]
  Libr1_plot<-Libr1_sp[SI(Libr1_sp)$Site_Date == List_Site_Date[i],]
  plot(Libr1_plot, FUN="mean", lwd=2,lty=1, ylim=c(0,1), col=c(colors[i],colors[i]), new=FALSE)
  Libr2_sp<-Libr2[SI(Libr2)$Species == SpeciesName,]
  Libr2_plot<-Libr2_sp[SI(Libr2_sp)$Site_Date == List_Site_Date[i],]
  plot(Libr2_plot, FUN="mean", lwd=2,lty=2, ylim=c(0,1), col=c(colors[i],colors[i]), new=FALSE)
  # plot(subset(Libr1, Species == SpeciesName & Site_Date == List_Site_Date[i]), FUN="mean", lwd=2,lty=1, ylim=c(0,1), col=c(colors[i],colors[i]), new=FALSE)
  # plot(subset(Libr2,  Species == SpeciesName & Site_Date == List_Site_Date[i]), FUN="mean", lwd=2,lty=2, ylim=c(0,1), col=c(colors[i],colors[i]), new=FALSE)
  # # plot(subset(Libr3, Species == SpeciesName & Site_Date == List_Site_Date[i]), FUN="mean", lwd=2,lty=3, ylim=c(0,1), col=c(colors[i],colors[i]), new=FALSE)
} 
legend(1300,0.99,List_Site_Date,lty=rep(1,length(List_Site_Date)),cex=1.2,lwd=2,bty="n",col=colors[1:length(List_Site_Date)])
legend(1750,0.99,paste("#",dataset_Field_subset.count[[2]]),cex=1.2,bty="n")
legend(1900,0.99,rep("Unmixed table bouquet spectrum",length(List_Site_Date)),lty=rep(2,length(List_Site_Date)),cex=1.2, lwd=2,bty="n",col=colors[1:length(List_Site_Date)])
legend(2350,0.99,paste("#",dataset_Table_unmixed_subset.count[[2]]),cex=1.2,bty="n")
# legend(1900,0.99,rep("Raw Table spectrum",length(List_Site_Date)),lty=rep(3,length(List_Site_Date)),cex=1.2, lwd=2,bty="n",col=colors[1:length(List_Site_Date)])
# legend(2280,0.99,paste("#",dataset_Table_raw_subset.count[[2]]),cex=1.2,bty="n")
axis(1,seq(350,2500,150),seq(350,2500,150),font=2)
axis(2,seq(0.0,1.05,0.1),seq(0.0,1.05,0.1),font=2)
mtext('Reflectance',side=2,font=2,line=3)
mtext('Wavelength (nm)',side=1,font=2,line=3)
title(SpeciesName)
dev.off()
# add confidence intervals???

jpeg(paste(paste(paste(SpeciesName,DvpPhase,sep="_"),jumpcorr,sep="_Field+Raw_"),".jpeg",sep=""),res=300,width=14,height=7,units='in')
# jpeg(paste(paste(SpeciesName,DvpPhase,sep="_"),".jpeg",sep="_Field+Raw"),res=300,width=14,height=7,units='in')
# jpeg("SOLGIG_Flower_Mean_Field+Raw.jpeg",res=300,width=14,height=7,units='in')
# x11(width=1280,height=720)
plot(c(),c(),axes=FALSE,ylim=c(0,1.0),xlim=c(350,2500),ylab='',xlab='', xaxs="i", yaxs="i")
for (i in 1:length(List_Site_Date)){
  Libr1_sp<-Libr1[SI(Libr1)$Species == SpeciesName,]
  Libr1_plot<-Libr1_sp[SI(Libr1_sp)$Site_Date == List_Site_Date[i],]
  plot(Libr1_plot, FUN="mean", lwd=2,lty=1, ylim=c(0,1), col=c(colors[i],colors[i]), new=FALSE)
  Libr3_sp<-Libr3[SI(Libr3)$Species == SpeciesName,]
  Libr3_plot<-Libr3_sp[SI(Libr3_sp)$Site_Date == List_Site_Date[i],]
  plot(Libr3_plot, FUN="mean", lwd=2,lty=3, ylim=c(0,1), col=c(colors[i],colors[i]), new=FALSE)
  # plot(subset(Libr1, Species == SpeciesName & Site_Date == List_Site_Date[i]), FUN="mean", lwd=2,lty=1, ylim=c(0,1), col=c(colors[i],colors[i]), new=FALSE)
  # plot(subset(Libr3, Species == SpeciesName & Site_Date == List_Site_Date[i]), FUN="mean", lwd=2,lty=3, ylim=c(0,1), col=c(colors[i],colors[i]), new=FALSE)
} 
legend(1300,0.99,List_Site_Date,lty=rep(1,length(List_Site_Date)),cex=1.2,lwd=2,bty="n",col=colors[1:length(List_Site_Date)])
legend(1750,0.99,paste("#",dataset_Field_subset.count[[2]]),cex=1.2,bty="n")
legend(1900,0.99,rep("Raw table bouquet spectrum",length(List_Site_Date)),lty=rep(3,length(List_Site_Date)),cex=1.2, lwd=2,bty="n",col=colors[1:length(List_Site_Date)])
legend(2280,0.99,paste("#",dataset_Table_raw_subset.count[[2]]),cex=1.2,bty="n")
axis(1,seq(350,2500,150),seq(350,2500,150),font=2)
axis(2,seq(0.0,1.05,0.1),seq(0.0,1.05,0.1),font=2)
mtext('Reflectance',side=2,font=2,line=3)
mtext('Wavelength (nm)',side=1,font=2,line=3)
title(SpeciesName)
dev.off()

#--------------------------------------------------------------------------------

########### Plot all individual spectra taken at each Site_Date
setwd(wd2)
# Field spectra
# colorstry=brewer.pal(6,"Dark2") #----> the plotted colors are different. However, the stored color codes are identic. How is this possible?
colorset_Field <- data.frame(levels(droplevels(dataset_Field_subset$Site_Date)),
                             color=colors[1:nlevels(droplevels(dataset_Field_subset$Site_Date))])
colnames(colorset_Field)<-c("Site_Date","color")
dataset_Field_subset_colorset <- merge(dataset_Field_subset, colorset_Field, by="Site_Date")
dataset_Field_subset_colorset<-dataset_Field_subset_colorset[order(dataset_Field_subset_colorset$Site_Date),]

jpeg(paste(paste(paste(SpeciesName,DvpPhase,sep="_"),jumpcorr,sep="_Field_"),".jpeg",sep=""),res=300,width=14,height=7,units='in')
# jpeg(paste(paste(SpeciesName,DvpPhase,sep="_"),".jpeg",sep="_Field"),res=300,width=14,height=7,units='in')
# jpeg("SOLGIG_Flower_Field.jpeg",res=300,width=14,height=7,units='in')
# x11(width=1280,height=720)
plot(c(),c(),axes=FALSE,ylim=c(0,1.05),xlim=c(350,2500),ylab='',xlab='', xaxs="i", yaxs="i")
for (i in 1:length(index_Field)){
  plot(Libr1, FUN=dataset_Field_subset_colorset$index[i],ylim=c(0,1), col=dataset_Field_subset_colorset[i,"color"],lty=1,lwd=2,new=FALSE)
}
legend(1400,0.99,dataset_Field_subset_colorset[,"Site_Date"],lty=rep(1,length(index_Field)),
       cex=0.8,lwd=2,bty="n", col=dataset_Field_subset_colorset[,"color"],ncol=2)
axis(1,seq(350,2500,150),seq(350,2500,150),font=2)
axis(2,seq(0.0,1.05,0.1),seq(0.0,1.05,0.1),font=2)
mtext('Reflectance',side=2,font=2,line=3)
mtext('Wavelength (nm)',side=1,font=2,line=3)
title(SpeciesName)
dev.off()

# Unmixed table spectra
colorset_Table_unmixed <- data.frame(levels(droplevels(dataset_Table_unmixed_subset$Site_Date)),
                                     color=colors[1:nlevels(droplevels(dataset_Table_unmixed_subset$Site_Date))])
colnames(colorset_Table_unmixed)<-c("Site_Date","color")
dataset_Table_unmixed_subset_colorset <- merge(dataset_Table_unmixed_subset, colorset_Table_unmixed, by="Site_Date")
dataset_Table_unmixed_subset_colorset<-dataset_Table_unmixed_subset_colorset[order(dataset_Table_unmixed_subset_colorset$Site_Date),]

jpeg(paste(paste(paste(SpeciesName,DvpPhase,sep="_"),jumpcorr,sep="_Unmixed_"),".jpeg",sep=""),res=300,width=14,height=7,units='in')
# jpeg(paste(paste(SpeciesName,DvpPhase,sep="_"),".jpeg",sep="_Unmixed"),res=300,width=14,height=7,units='in')
# jpeg("SOLGIG_Flower_Unmixed.jpeg",res=300,width=14,height=7,units='in')
# x11(width=1280,height=720)
plot(c(),c(),axes=FALSE,ylim=c(0,1.05),xlim=c(350,2500),ylab='',xlab='', xaxs="i", yaxs="i")
for (i in 1:length(index_Table_unmixed)){
  plot(Libr2, FUN=dataset_Table_unmixed_subset_colorset$index[i],ylim=c(0,1), col=dataset_Table_unmixed_subset_colorset[i,"color"],lty=2,lwd=2,new=FALSE)
}
legend(1400,0.99,paste(dataset_Table_unmixed_subset_colorset[,"Site_Date"],"Unmixed table bouquet spectrum"),lty=rep(2,length(index_Table_unmixed)),
       cex=0.8,lwd=2,bty="n", col=dataset_Table_unmixed_subset_colorset[,"color"],ncol=2)
# legend(2100,0.99,paste("fraction_veg = ",round(dataset_Table_unmixed_subset_colorset[,"Fraction_veg"],digits=2)),
#        cex=0.8,bty="n"
axis(1,seq(350,2500,150),seq(350,2500,150),font=2)
axis(2,seq(0.0,1.05,0.1),seq(0.0,1.05,0.1),font=2)
mtext('Reflectance',side=2,font=2,line=3)
mtext('Wavelength (nm)',side=1,font=2,line=3)
title(SpeciesName)
dev.off()

# Raw table spectra
colorset_Table_raw <- data.frame(levels(droplevels(dataset_Table_raw_subset$Site_Date)),
                                 color=colors[1:nlevels(droplevels(dataset_Table_raw_subset$Site_Date))])
colnames(colorset_Table_raw)<-c("Site_Date","color")
dataset_Table_raw_subset_colorset <- merge(dataset_Table_raw_subset, colorset_Table_raw, by="Site_Date")
dataset_Table_raw_subset_colorset<-dataset_Table_raw_subset_colorset[order(dataset_Table_raw_subset_colorset$Site_Date),]

jpeg(paste(paste(paste(SpeciesName,DvpPhase,sep="_"),jumpcorr,sep="_Raw_"),".jpeg",sep=""),res=300,width=14,height=7,units='in')
# jpeg(paste(paste(SpeciesName,DvpPhase,sep="_"),".jpeg",sep="Raw"),res=300,width=14,height=7,units='in')
# jpeg("SOLGIG_Flower_Raw.jpeg",res=300,width=14,height=7,units='in')
# x11(width=1280,height=720)
plot(c(),c(),axes=FALSE,ylim=c(0,1.05),xlim=c(350,2500),ylab='',xlab='', xaxs="i", yaxs="i")
for (i in 1:length(index_Table_raw)){
  plot(Libr3, FUN=dataset_Table_raw_subset_colorset$index[i],ylim=c(0,1), col=dataset_Table_raw_subset_colorset[i,"color"],lty=3,lwd=2,new=FALSE)
}
legend(1400,0.99,paste(dataset_Table_unmixed_subset_colorset[,"Site_Date"],"Raw table bouquet spectrum"),lty=rep(3,length(index_Table_raw)),
       cex=0.8,lwd=2,bty="n", col=dataset_Table_unmixed_subset_colorset[,"color"],ncol=2)
axis(1,seq(350,2500,150),seq(350,2500,150),font=2)
axis(2,seq(0.0,1.05,0.1),seq(0.0,1.05,0.1),font=2)
mtext('Reflectance',side=2,font=2,line=3)
mtext('Wavelength (nm)',side=1,font=2,line=3)
title(SpeciesName)
dev.off()

#--------------------------------------------------------------------------------

########### Plot mean and sd of all individual spectra of a species
library("RColorBrewer")
colors = brewer.pal(3, "Dark2") # Color codes are: #1b9e77, #d95f02, #7570b3
colors_sd<-c("#b3e2cd","#fdcdac","#cbd5e8") #a softer combination of the same colors defined above (http://colorbrewer2.org/#type=qualitative&scheme=Pastel2&n=3)
transp<-0.3

# ATTENTION: Better take average of mean spectra per day than average of all bouquets because a day with lots of measurements will then control the overall average
# --> _MeanOverDays instead of _MeanOverBouqets
# Drawback: IN case of only 1 day we won't have a sd interval...
setwd(wd3)
jpeg(paste(paste(paste(SpeciesName,DvpPhase,sep="_"),jumpcorr,sep="_"),"_MeanOverBouqets.jpeg",sep=""),res=300,width=12.8,height=7.2,units='in')
# jpeg(paste(paste(SpeciesName,DvpPhase,sep="_"),".jpeg",sep=""),res=300,width=12.8,height=7.2,units='in')
# jpeg("SOLGIG_Flower.jpeg",res=300,width=12.8,height=7.2,units='in')
# pdf("SOLGIG_Flower.pdf",width=12.8,height=7.2)
# x11(width=1280,height=720)
plot(c(),c(),axes=FALSE,ylim=c(0,1.0),xlim=c(350,2500),ylab='',xlab='', xaxs="i", yaxs="i")
  # 1) Colour the area between the sd-borders
  wl_plot <- c(seq(400, 1340, 1), seq(1460,1780,1), seq (1970, 2400, 1))
  xx<-c(wl_plot,(rev(wl_plot)))
  Libr3_sp<-Libr3[SI(Libr3)$Species == SpeciesName,]
  idx_int <- SI(Libr3_sp)$Site_Date %in% List_Site_Date
  Libr3_plot<-Libr3_sp[idx_int,]
  mean_spec_Libr3 <- apply(Libr3_plot, FUN = mean, na.rm = TRUE)
  sd_spec_Libr3 <- apply(Libr3_plot, FUN = sd, na.rm = TRUE)
  # mean_spec_Libr3 <- apply(subset(Libr3, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = mean, na.rm = TRUE)
  # sd_spec_Libr3   <- apply(subset(Libr3, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = sd, na.rm = TRUE)
  spectra2plot_Libr3 <- rbind(spectra(mean_spec_Libr3) + spectra(sd_spec_Libr3),
                              spectra(mean_spec_Libr3),
                              spectra(mean_spec_Libr3) - spectra(sd_spec_Libr3))
  yy_Libr3<-c(spectra2plot_Libr3[1,],rev(spectra2plot_Libr3[3,]))
  # yy<-c(spectra2plot_Libr3[1,c(56:1001,1101:1446,1646:2051)],rev(spectra2plot_Libr3[3,c(56:1001,1101:1446,1646:2051)])) # the same result
  polygon(xx,yy_Libr3,col=alpha(c(colors_sd[3]),transp),border=NA,new=F)

  Libr2_sp<-Libr2[SI(Libr2)$Species == SpeciesName,]
  idx_int <- SI(Libr2_sp)$Site_Date %in% List_Site_Date
  Libr2_plot<-Libr2_sp[idx_int,]
  mean_spec_Libr2 <- apply(Libr2_plot, FUN = mean, na.rm = TRUE)
  sd_spec_Libr2 <- apply(Libr2_plot, FUN = sd, na.rm = TRUE)
  # mean_spec_Libr2 <- apply(subset(Libr2, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = mean, na.rm = TRUE)
  # sd_spec_Libr2   <- apply(subset(Libr2, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = sd, na.rm = TRUE)
  spectra2plot_Libr2 <- rbind(spectra(mean_spec_Libr2) + spectra(sd_spec_Libr2),
                              spectra(mean_spec_Libr2),
                              spectra(mean_spec_Libr2) - spectra(sd_spec_Libr2))
  yy_Libr2<-c(spectra2plot_Libr2[1,],rev(spectra2plot_Libr2[3,]))
  polygon(xx,yy_Libr2,col=alpha(c(colors_sd[2]),transp),border=NA,new=F)

  Libr1_sp<-Libr1[SI(Libr1)$Species == SpeciesName,]
  idx_int <- SI(Libr1_sp)$Site_Date %in% List_Site_Date
  Libr1_plot<-Libr1_sp[idx_int,]
  mean_spec_Libr1 <- apply(Libr1_plot, FUN = mean, na.rm = TRUE)
  sd_spec_Libr1 <- apply(Libr1_plot, FUN = sd, na.rm = TRUE)
  # mean_spec_Libr1 <- apply(subset(Libr1, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = mean, na.rm = TRUE)
  # sd_spec_Libr1   <- apply(subset(Libr1, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = sd, na.rm = TRUE)
  spectra2plot_Libr1 <- rbind(spectra(mean_spec_Libr1) + spectra(sd_spec_Libr1),
                              spectra(mean_spec_Libr1),
                              spectra(mean_spec_Libr1) - spectra(sd_spec_Libr1))
  yy_Libr1<-c(spectra2plot_Libr1[1,],rev(spectra2plot_Libr1[3,]))
  polygon(xx,yy_Libr1,col=alpha(c(colors_sd[1]),transp),border=NA,new=F)
  xx_mask1<-c(seq(1340,1460,1),rev(seq(1340,1460,1)))
  xx_mask2<-c(seq(1780,1970,1),rev(seq(1780,1970,1)))
  yy_mask1<-c(seq(0,0,length.out=length(xx_mask1)/2),seq(1,1,length.out=length(xx_mask1)/2))
  yy_mask2<-c(seq(0,0,length.out=length(xx_mask2)/2),seq(1,1,length.out=length(xx_mask2)/2))
  polygon(xx_mask1,yy_mask1,col="white",border=NA,new=F)
  polygon(xx_mask2,yy_mask2,col="white",border=NA,new=F)

  # 2) Plot mean spectrum and sd on top of the areas
  # plot(subset(Libr1, Species == SpeciesName & Site_Date %in% List_Site_Date), lwd=2, ylim=c(0,1), col=c(colors[1]), new=FALSE)
  # plot(subset(Libr3, Species == SpeciesName & Site_Date %in% List_Site_Date), lwd=2, ylim=c(0,1), col=c(colors[3]), new=FALSE)
  # plot(subset(Libr2, Species == SpeciesName & Site_Date %in% List_Site_Date), lwd=2, ylim=c(0,1), col=c(colors[2]), new=FALSE)
  plot(Libr1_plot, lwd=2, ylim=c(0,1), col=c(colors[1]), new=FALSE)
  plot(Libr2_plot, lwd=2, ylim=c(0,1), col=c(colors[2]), new=FALSE)
  plot(Libr3_plot, lwd=2, ylim=c(0,1), col=c(colors[3]), new=FALSE)
  
  # 3) Plot graphical details: axes, legend, etc.
  legend(1850,0.99,"Field patch",cex=1.2,lwd=2,bty="n",col=colors[1])
  legend(1850,0.89,"Unmixed table bouquet",cex=1.2,lwd=2,bty="n",col=colors[2])
  legend(1850,0.79,"Raw table bouquet",cex=1.2,lwd=2,bty="n",col=colors[3])
  legend(2280,0.99,paste("#",dataset_Field_subset.totalcount[[2]]),cex=1.2,bty="n")
  legend(2280,0.89,paste("#",dataset_Table_unmixed_subset.totalcount[[2]]),cex=1.2,bty="n")
  legend(2280,0.79,paste("#",dataset_Table_raw_subset.totalcount[[2]]),cex=1.2,bty="n")
  axis(1,seq(350,2500,150),seq(350,2500,150),font=2)
  axis(2,seq(0.0,1.05,0.1),seq(0.0,1.05,0.1),font=2,las=2)
  mtext('Reflectance',side=2,font=2,line=3)
  mtext('Wavelength (nm)',side=1,font=2,line=3)
  title(paste(SpeciesName, DvpPhase))
dev.off()

setwd(wd4)
jpeg(paste(paste(paste(SpeciesName,DvpPhase,sep="_"),jumpcorr,sep="_"),"_MeanOverDays.jpeg",sep=""),res=300,width=12.8,height=7.2,units='in')
# x11(width=1280,height=720)
plot(c(),c(),axes=FALSE,ylim=c(0,1.0),xlim=c(350,2500),ylab='',xlab='', xaxs="i", yaxs="i")
  # 1) Colour the area between the sd-borders
  wl_plot <- c(seq(400, 1340, 1), seq(1460,1780,1), seq (1970, 2400, 1))
  xx<-c(wl_plot,(rev(wl_plot)))
  Libr3_sp<-Libr3mean[SI(Libr3mean)$Species == SpeciesName,]
  idx_int <- SI(Libr3_sp)$Site_Date %in% List_Site_Date
  Libr3_plot<-Libr3_sp[idx_int,]
  mean_spec_Libr3 <- apply(Libr3_plot, FUN = mean, na.rm = TRUE)
  sd_spec_Libr3 <- apply(Libr3_plot, FUN = sd, na.rm = TRUE)
  # mean_spec_Libr3 <- apply(subset(Libr3, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = mean, na.rm = TRUE)
  # sd_spec_Libr3   <- apply(subset(Libr3, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = sd, na.rm = TRUE)
  spectra2plot_Libr3 <- rbind(spectra(mean_spec_Libr3) + spectra(sd_spec_Libr3),
                              spectra(mean_spec_Libr3),
                              spectra(mean_spec_Libr3) - spectra(sd_spec_Libr3))
  yy_Libr3<-c(spectra2plot_Libr3[1,],rev(spectra2plot_Libr3[3,]))
  # yy<-c(spectra2plot_Libr3[1,c(56:1001,1101:1446,1646:2051)],rev(spectra2plot_Libr3[3,c(56:1001,1101:1446,1646:2051)])) # the same result
  polygon(xx,yy_Libr3,col=alpha(c(colors_sd[3]),transp),border=NA,new=F)
  
  Libr2_sp<-Libr2mean[SI(Libr2mean)$Species == SpeciesName,]
  idx_int <- SI(Libr2_sp)$Site_Date %in% List_Site_Date
  Libr2_plot<-Libr2_sp[idx_int,]
  mean_spec_Libr2 <- apply(Libr2_plot, FUN = mean, na.rm = TRUE)
  sd_spec_Libr2 <- apply(Libr2_plot, FUN = sd, na.rm = TRUE)
  # mean_spec_Libr2 <- apply(subset(Libr2, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = mean, na.rm = TRUE)
  # sd_spec_Libr2   <- apply(subset(Libr2, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = sd, na.rm = TRUE)
  spectra2plot_Libr2 <- rbind(spectra(mean_spec_Libr2) + spectra(sd_spec_Libr2),
                              spectra(mean_spec_Libr2),
                              spectra(mean_spec_Libr2) - spectra(sd_spec_Libr2))
  yy_Libr2<-c(spectra2plot_Libr2[1,],rev(spectra2plot_Libr2[3,]))
  polygon(xx,yy_Libr2,col=alpha(c(colors_sd[2]),transp),border=NA,new=F)
  
  Libr1_sp<-Libr1mean[SI(Libr1mean)$Species == SpeciesName,]
  idx_int <- SI(Libr1_sp)$Site_Date %in% List_Site_Date
  Libr1_plot<-Libr1_sp[idx_int,]
  mean_spec_Libr1 <- apply(Libr1_plot, FUN = mean, na.rm = TRUE)
  sd_spec_Libr1 <- apply(Libr1_plot, FUN = sd, na.rm = TRUE)
  # mean_spec_Libr1 <- apply(subset(Libr1, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = mean, na.rm = TRUE)
  # sd_spec_Libr1   <- apply(subset(Libr1, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = sd, na.rm = TRUE)
  spectra2plot_Libr1 <- rbind(spectra(mean_spec_Libr1) + spectra(sd_spec_Libr1),
                              spectra(mean_spec_Libr1),
                              spectra(mean_spec_Libr1) - spectra(sd_spec_Libr1))
  yy_Libr1<-c(spectra2plot_Libr1[1,],rev(spectra2plot_Libr1[3,]))
  polygon(xx,yy_Libr1,col=alpha(c(colors_sd[1]),transp),border=NA,new=F)
  xx_mask1<-c(seq(1340,1460,1),rev(seq(1340,1460,1)))
  xx_mask2<-c(seq(1780,1970,1),rev(seq(1780,1970,1)))
  yy_mask1<-c(seq(0,0,length.out=length(xx_mask1)/2),seq(1,1,length.out=length(xx_mask1)/2))
  yy_mask2<-c(seq(0,0,length.out=length(xx_mask2)/2),seq(1,1,length.out=length(xx_mask2)/2))
  polygon(xx_mask1,yy_mask1,col="white",border=NA,new=F)
  polygon(xx_mask2,yy_mask2,col="white",border=NA,new=F)

  # 2) Plot mean spectrum and sd on top of the areas
  # plot(subset(Libr1, Species == SpeciesName & Site_Date %in% List_Site_Date), lwd=2, ylim=c(0,1), col=c(colors[1]), new=FALSE)
  # plot(subset(Libr3, Species == SpeciesName & Site_Date %in% List_Site_Date), lwd=2, ylim=c(0,1), col=c(colors[3]), new=FALSE)
  # plot(subset(Libr2, Species == SpeciesName & Site_Date %in% List_Site_Date), lwd=2, ylim=c(0,1), col=c(colors[2]), new=FALSE)
  plot(Libr1_plot, lwd=2, ylim=c(0,1), col=c(colors[1]), new=FALSE)
  plot(Libr2_plot, lwd=2, ylim=c(0,1), col=c(colors[2]), new=FALSE)
  plot(Libr3_plot, lwd=2, ylim=c(0,1), col=c(colors[3]), new=FALSE)
  
  # 3) Plot graphical details: axes, legend, etc.
  legend(1850,0.99,"Field patch",cex=1.2,lwd=2,bty="n",col=colors[1])
  legend(1850,0.89,"Unmixed table bouquet",cex=1.2,lwd=2,bty="n",col=colors[2])
  legend(1850,0.79,"Raw table bouquet",cex=1.2,lwd=2,bty="n",col=colors[3])
  legend(2280,0.99,paste("#",dataset_Field_subset.totalcount.mean[[2]]),cex=1.2,bty="n")
  legend(2280,0.89,paste("#",dataset_Table_unmixed_subset.totalcount.mean[[2]]),cex=1.2,bty="n")
  legend(2280,0.79,paste("#",dataset_Table_raw_subset.totalcount.mean[[2]]),cex=1.2,bty="n")
  axis(1,seq(350,2500,150),seq(350,2500,150),font=2)
  axis(2,seq(0.0,1.05,0.1),seq(0.0,1.05,0.1),font=2,las=2)
  mtext('Reflectance',side=2,font=2,line=3)
  mtext('Wavelength (nm)',side=1,font=2,line=3)
  title(paste(SpeciesName, DvpPhase))
dev.off()

}
else if (VIbands == "Yes") {
  ########### Plot mean and sd of all individual spectra of a species
  library("RColorBrewer")
  colors = brewer.pal(3, "Dark2") # Color codes are: #1b9e77, #d95f02, #7570b3
  colors_sd<-c("#b3e2cd","#fdcdac","#cbd5e8") #a softer combination of the same colors defined above (http://colorbrewer2.org/#type=qualitative&scheme=Pastel2&n=3)
  transp<-0.3
  
  # ATTENTION: Better take average of mean spectra per day than average of all bouquets because a day with lots of measurements will then control the overall average
  # --> _MeanOverDays instead of _MeanOverBouqets
  # Drawback: IN case of only 1 day we won't have a sd interval...
  
  setwd(wd3)
  # jpeg(paste(paste(paste(SpeciesName,DvpPhase,sep="_"),jumpcorr,sep="_"),"_MeanOverBouqets.jpeg",sep=""),res=300,width=12.8,height=7.2,units='in')
  x11(width=1280,height=720)
  plot(c(),c(),axes=FALSE,ylim=c(0,1.0),xlim=c(350,2500),ylab='',xlab='', xaxs="i", yaxs="i")
  # 1) Colour the area between the sd-borders
  wl_plot <- c(seq(400, 1340, 1), seq(1460,1780,1), seq (1970, 2400, 1))
  xx<-c(wl_plot,(rev(wl_plot)))
  Libr3_sp<-Libr3[SI(Libr3)$Species == SpeciesName,]
  idx_int <- SI(Libr3_sp)$Site_Date %in% List_Site_Date
  Libr3_plot<-Libr3_sp[idx_int,]
  mean_spec_Libr3 <- apply(Libr3_plot, FUN = mean, na.rm = TRUE)
  sd_spec_Libr3 <- apply(Libr3_plot, FUN = sd, na.rm = TRUE)
  # mean_spec_Libr3 <- apply(subset(Libr3, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = mean, na.rm = TRUE)
  # sd_spec_Libr3   <- apply(subset(Libr3, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = sd, na.rm = TRUE)
  spectra2plot_Libr3 <- rbind(spectra(mean_spec_Libr3) + spectra(sd_spec_Libr3),
                              spectra(mean_spec_Libr3),
                              spectra(mean_spec_Libr3) - spectra(sd_spec_Libr3))
  yy_Libr3<-c(spectra2plot_Libr3[1,],rev(spectra2plot_Libr3[3,]))
  # yy<-c(spectra2plot_Libr3[1,c(56:1001,1101:1446,1646:2051)],rev(spectra2plot_Libr3[3,c(56:1001,1101:1446,1646:2051)])) # the same result
  polygon(xx,yy_Libr3,col=alpha(c(colors_sd[3]),transp),border=NA,new=F)
  
  Libr2_sp<-Libr2[SI(Libr2)$Species == SpeciesName,]
  idx_int <- SI(Libr2_sp)$Site_Date %in% List_Site_Date
  Libr2_plot<-Libr2_sp[idx_int,]
  mean_spec_Libr2 <- apply(Libr2_plot, FUN = mean, na.rm = TRUE)
  sd_spec_Libr2 <- apply(Libr2_plot, FUN = sd, na.rm = TRUE)
  # mean_spec_Libr2 <- apply(subset(Libr2, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = mean, na.rm = TRUE)
  # sd_spec_Libr2   <- apply(subset(Libr2, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = sd, na.rm = TRUE)
  spectra2plot_Libr2 <- rbind(spectra(mean_spec_Libr2) + spectra(sd_spec_Libr2),
                              spectra(mean_spec_Libr2),
                              spectra(mean_spec_Libr2) - spectra(sd_spec_Libr2))
  yy_Libr2<-c(spectra2plot_Libr2[1,],rev(spectra2plot_Libr2[3,]))
  polygon(xx,yy_Libr2,col=alpha(c(colors_sd[2]),transp),border=NA,new=F)
  
  Libr1_sp<-Libr1[SI(Libr1)$Species == SpeciesName,]
  idx_int <- SI(Libr1_sp)$Site_Date %in% List_Site_Date
  Libr1_plot<-Libr1_sp[idx_int,]
  mean_spec_Libr1 <- apply(Libr1_plot, FUN = mean, na.rm = TRUE)
  sd_spec_Libr1 <- apply(Libr1_plot, FUN = sd, na.rm = TRUE)
  # mean_spec_Libr1 <- apply(subset(Libr1, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = mean, na.rm = TRUE)
  # sd_spec_Libr1   <- apply(subset(Libr1, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = sd, na.rm = TRUE)
  spectra2plot_Libr1 <- rbind(spectra(mean_spec_Libr1) + spectra(sd_spec_Libr1),
                              spectra(mean_spec_Libr1),
                              spectra(mean_spec_Libr1) - spectra(sd_spec_Libr1))
  yy_Libr1<-c(spectra2plot_Libr1[1,],rev(spectra2plot_Libr1[3,]))
  polygon(xx,yy_Libr1,col=alpha(c(colors_sd[1]),transp),border=NA,new=F)
  xx_mask1<-c(seq(1340,1460,1),rev(seq(1340,1460,1)))
  xx_mask2<-c(seq(1780,1970,1),rev(seq(1780,1970,1)))
  yy_mask1<-c(seq(0,0,length.out=length(xx_mask1)/2),seq(1,1,length.out=length(xx_mask1)/2))
  yy_mask2<-c(seq(0,0,length.out=length(xx_mask2)/2),seq(1,1,length.out=length(xx_mask2)/2))
  polygon(xx_mask1,yy_mask1,col="white",border=NA,new=F)
  polygon(xx_mask2,yy_mask2,col="white",border=NA,new=F)
  
  # 2) Plot mean spectrum and sd on top of the areas
  # plot(subset(Libr1, Species == SpeciesName & Site_Date %in% List_Site_Date), lwd=2, ylim=c(0,1), col=c(colors[1]), new=FALSE)
  # plot(subset(Libr3, Species == SpeciesName & Site_Date %in% List_Site_Date), lwd=2, ylim=c(0,1), col=c(colors[3]), new=FALSE)
  # plot(subset(Libr2, Species == SpeciesName & Site_Date %in% List_Site_Date), lwd=2, ylim=c(0,1), col=c(colors[2]), new=FALSE)
  plot(Libr1_plot, lwd=2, ylim=c(0,1), col=c(colors[1]), new=FALSE)
  plot(Libr2_plot, lwd=2, ylim=c(0,1), col=c(colors[2]), new=FALSE)
  plot(Libr3_plot, lwd=2, ylim=c(0,1), col=c(colors[3]), new=FALSE)
  
  # 3) Plot graphical details: axes, legend, etc.
  legend(1850,0.99,"Field patch",cex=1.2,lwd=2,bty="n",col=colors[1])
  legend(1850,0.89,"Unmixed table bouquet",cex=1.2,lwd=2,bty="n",col=colors[2])
  legend(1850,0.79,"Raw table bouquet",cex=1.2,lwd=2,bty="n",col=colors[3])
  legend(2280,0.99,paste("#",dataset_Field_subset.totalcount[[2]]),cex=1.2,bty="n")
  legend(2280,0.89,paste("#",dataset_Table_unmixed_subset.totalcount[[2]]),cex=1.2,bty="n")
  legend(2280,0.79,paste("#",dataset_Table_raw_subset.totalcount[[2]]),cex=1.2,bty="n")
  axis(1,seq(350,2500,150),seq(350,2500,150),font=2)
  axis(2,seq(0.0,1.05,0.1),seq(0.0,1.05,0.1),font=2,las=2)
  mtext('Reflectance',side=2,font=2,line=3)
  mtext('Wavelength (nm)',side=1,font=2,line=3)
  title(paste(SpeciesName, DvpPhase))
  
  # 4) Plot VI wavelength positions
  abline(v=VIbandpos)
  # dev.off()
  
   setwd(wd3)
  # jpeg(paste(paste(paste(SpeciesName,DvpPhase,sep="_"),jumpcorr,sep="_"),"_MeanOverDays.jpeg",sep=""),res=300,width=12.8,height=7.2,units='in')
  x11(width=1280,height=720)
  plot(c(),c(),axes=FALSE,ylim=c(0,1.0),xlim=c(350,2500),ylab='',xlab='', xaxs="i", yaxs="i")
  # 1) Colour the area between the sd-borders
  wl_plot <- c(seq(400, 1340, 1), seq(1460,1780,1), seq (1970, 2400, 1))
  xx<-c(wl_plot,(rev(wl_plot)))
  Libr3_sp<-Libr3mean[SI(Libr3mean)$Species == SpeciesName,]
  idx_int <- SI(Libr3_sp)$Site_Date %in% List_Site_Date
  Libr3_plot<-Libr3_sp[idx_int,]
  mean_spec_Libr3 <- apply(Libr3_plot, FUN = mean, na.rm = TRUE)
  sd_spec_Libr3 <- apply(Libr3_plot, FUN = sd, na.rm = TRUE)
  # mean_spec_Libr3 <- apply(subset(Libr3, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = mean, na.rm = TRUE)
  # sd_spec_Libr3   <- apply(subset(Libr3, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = sd, na.rm = TRUE)
  spectra2plot_Libr3 <- rbind(spectra(mean_spec_Libr3) + spectra(sd_spec_Libr3),
                              spectra(mean_spec_Libr3),
                              spectra(mean_spec_Libr3) - spectra(sd_spec_Libr3))
  yy_Libr3<-c(spectra2plot_Libr3[1,],rev(spectra2plot_Libr3[3,]))
  # yy<-c(spectra2plot_Libr3[1,c(56:1001,1101:1446,1646:2051)],rev(spectra2plot_Libr3[3,c(56:1001,1101:1446,1646:2051)])) # the same result
  polygon(xx,yy_Libr3,col=alpha(c(colors_sd[3]),transp),border=NA,new=F)
  
  Libr2_sp<-Libr2mean[SI(Libr2mean)$Species == SpeciesName,]
  idx_int <- SI(Libr2_sp)$Site_Date %in% List_Site_Date
  Libr2_plot<-Libr2_sp[idx_int,]
  mean_spec_Libr2 <- apply(Libr2_plot, FUN = mean, na.rm = TRUE)
  sd_spec_Libr2 <- apply(Libr2_plot, FUN = sd, na.rm = TRUE)
  # mean_spec_Libr2 <- apply(subset(Libr2, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = mean, na.rm = TRUE)
  # sd_spec_Libr2   <- apply(subset(Libr2, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = sd, na.rm = TRUE)
  spectra2plot_Libr2 <- rbind(spectra(mean_spec_Libr2) + spectra(sd_spec_Libr2),
                              spectra(mean_spec_Libr2),
                              spectra(mean_spec_Libr2) - spectra(sd_spec_Libr2))
  yy_Libr2<-c(spectra2plot_Libr2[1,],rev(spectra2plot_Libr2[3,]))
  polygon(xx,yy_Libr2,col=alpha(c(colors_sd[2]),transp),border=NA,new=F)
  
  Libr1_sp<-Libr1mean[SI(Libr1mean)$Species == SpeciesName,]
  idx_int <- SI(Libr1_sp)$Site_Date %in% List_Site_Date
  Libr1_plot<-Libr1_sp[idx_int,]
  mean_spec_Libr1 <- apply(Libr1_plot, FUN = mean, na.rm = TRUE)
  sd_spec_Libr1 <- apply(Libr1_plot, FUN = sd, na.rm = TRUE)
  # mean_spec_Libr1 <- apply(subset(Libr1, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = mean, na.rm = TRUE)
  # sd_spec_Libr1   <- apply(subset(Libr1, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = sd, na.rm = TRUE)
  spectra2plot_Libr1 <- rbind(spectra(mean_spec_Libr1) + spectra(sd_spec_Libr1),
                              spectra(mean_spec_Libr1),
                              spectra(mean_spec_Libr1) - spectra(sd_spec_Libr1))
  yy_Libr1<-c(spectra2plot_Libr1[1,],rev(spectra2plot_Libr1[3,]))
  polygon(xx,yy_Libr1,col=alpha(c(colors_sd[1]),transp),border=NA,new=F)
  xx_mask1<-c(seq(1340,1460,1),rev(seq(1340,1460,1)))
  xx_mask2<-c(seq(1780,1970,1),rev(seq(1780,1970,1)))
  yy_mask1<-c(seq(0,0,length.out=length(xx_mask1)/2),seq(1,1,length.out=length(xx_mask1)/2))
  yy_mask2<-c(seq(0,0,length.out=length(xx_mask2)/2),seq(1,1,length.out=length(xx_mask2)/2))
  polygon(xx_mask1,yy_mask1,col="white",border=NA,new=F)
  polygon(xx_mask2,yy_mask2,col="white",border=NA,new=F)
  
  # 2) Plot mean spectrum and sd on top of the areas
  # plot(subset(Libr1, Species == SpeciesName & Site_Date %in% List_Site_Date), lwd=2, ylim=c(0,1), col=c(colors[1]), new=FALSE)
  # plot(subset(Libr3, Species == SpeciesName & Site_Date %in% List_Site_Date), lwd=2, ylim=c(0,1), col=c(colors[3]), new=FALSE)
  # plot(subset(Libr2, Species == SpeciesName & Site_Date %in% List_Site_Date), lwd=2, ylim=c(0,1), col=c(colors[2]), new=FALSE)
  plot(Libr1_plot, lwd=2, ylim=c(0,1), col=c(colors[1]), new=FALSE)
  plot(Libr2_plot, lwd=2, ylim=c(0,1), col=c(colors[2]), new=FALSE)
  plot(Libr3_plot, lwd=2, ylim=c(0,1), col=c(colors[3]), new=FALSE)
  
  # 3) Plot graphical details: axes, legend, etc.
  legend(1850,0.99,"Field patch",cex=1.2,lwd=2,bty="n",col=colors[1])
  legend(1850,0.89,"Unmixed table bouquet",cex=1.2,lwd=2,bty="n",col=colors[2])
  legend(1850,0.79,"Raw table bouquet",cex=1.2,lwd=2,bty="n",col=colors[3])
  legend(2280,0.99,paste("#",dataset_Field_subset.totalcount.mean[[2]]),cex=1.2,bty="n")
  legend(2280,0.89,paste("#",dataset_Table_unmixed_subset.totalcount.mean[[2]]),cex=1.2,bty="n")
  legend(2280,0.79,paste("#",dataset_Table_raw_subset.totalcount.mean[[2]]),cex=1.2,bty="n")
  axis(1,seq(350,2500,150),seq(350,2500,150),font=2)
  axis(2,seq(0.0,1.05,0.1),seq(0.0,1.05,0.1),font=2,las=2)
  mtext('Reflectance',side=2,font=2,line=3)
  mtext('Wavelength (nm)',side=1,font=2,line=3)
  title(paste(SpeciesName, DvpPhase))
  
  # 4) Plot VI wavelength positions
  abline(v=c(925,710,754,709,681,750,705,550),col="green",lty=2)
  abline(v=c(539,490,807),col="orange",lty=2)
  abline(v=c(790,705,860,720,523,583,1645,1715,531,570),col="yellow",lty=2)
  abline(v=c(857,1241),col="blue",lty=2)
  abline(v=c(1662,1732,1550,1750,2260,1490),col="brown",lty=2)
  # abline(v=VIbandpos)
  # text(VIbandpos,1,label=VIbandpos,pos=2,srt=90)
  # dev.off()
  }
}
# Apply spectral_plots to all species measured
spectral_plots_apply <- function(wd1, wd2, wd3, wd4, jumpcorr, VIbands) {
  setwd(wd4)
  datasets<-readRDS(paste("datasets_",".rds",sep=jumpcorr))  
  
  levels <- levels(droplevels(SI(datasets$libr_T_unmixed.intersect)$Sp_Phase))
  SpeciesNames <- substr(levels,1,6)
  if (VIbands == "No"){
    for (i in 1:length(levels)){
      if (grepl("Flower",levels[i],fixed=T)=="TRUE"){
        spectral_plots(datasets=datasets, wd1, wd2, wd3, jumpcorr=jumpcorr, SpeciesNames[i], DvpPhase="Flower", VIbands="No",VIbandpos=c(0))
      } else if (grepl("Green",levels[i],fixed=T)=="TRUE"){
        spectral_plots(datasets=datasets, wd1, wd2, wd3, jumpcorr=jumpcorr, SpeciesNames[i], DvpPhase="Green", VIbands="No",VIbandpos=c(0))
      }}
  }
  else if (VIbands == "Yes"){
    VIbandpositions = c(925,710,754,709,681,750,705,550,539,490,807,490,730,705,860,720,523,583,1645,1715,857,1241,1662,1732,1550,1750,2260,1490,531,570)
    spectral_plots(datasets=datasets, wd1, wd2, wd3, jumpcorr=jumpcorr, "URTDIO", DvpPhase="Green", VIbands="Yes", VIbandpos=VIbandpositions)
  }
  
}

# Plot field + raw data and field + unmixed data
spectral_plots_separate_raw_unmix <- function (datasets, wd5, jumpcorr, SpeciesName, DvpPhase) {
  
  # 1) Specify species and datasets for visualisation (development phase: Green or Flower)
  # if (DvpPhase == "Green"){
  #   dataset_Field<-datasets$meta_FG
  #   dataset_Table_unmixed<-datasets$meta_TG_unmixed
  #   dataset_Table_raw<-datasets$meta_TG_raw
  #   Libr1a<-datasets$library_FG
  #   Libr2a<-datasets$library_TG_unmixed
  #   Libr3a<-datasets$library_TG_raw
  # } else if (DvpPhase == "Flower"){
  #   dataset_Field<-datasets$meta_FF
  #   dataset_Table_unmixed<-datasets$meta_TF_unmixed
  #   dataset_Table_raw<-datasets$meta_TF_raw
  #   Libr1a<-datasets$library_FF
  #   Libr2a<-datasets$library_TF_unmixed
  #   Libr3a<-datasets$library_TF_raw
  # }
  
  Libr1a<-datasets$libr_F.intersect_orig
  Libr2a<-datasets$libr_T_unmixed.intersect_orig
  Libr3a<-datasets$libr_T_raw.intersect_orig
  if (DvpPhase == "Green"){
    Libr1a<-subset(Libr1a, DvpPhase == "Green")
    Libr2a<-subset(Libr2a, DvpPhase == "Green")
    Libr3a<-subset(Libr3a, DvpPhase == "Green")
  } else if (DvpPhase == "Flower"){
    Libr1a<-subset(Libr1a, DvpPhase == "Flower")
    Libr2a<-subset(Libr2a, DvpPhase == "Flower")
    Libr3a<-subset(Libr3a, DvpPhase == "Flower")
  }
  dataset_Field<-SI(Libr1a)
  dataset_Table_unmixed<-SI(Libr2a)
  dataset_Table_raw<-SI(Libr3a)
  
  # Mean per day:
  Libr1b<-datasets$libr_F.intersect.mean_orig
  Libr2b<-datasets$libr_T_unmixed.intersect.mean_orig
  Libr3b<-datasets$libr_T_raw.intersect.mean_orig
  if (DvpPhase == "Green"){
    Libr1mean<-subset(Libr1b, DvpPhase == "Green")
    Libr2mean<-subset(Libr2b, DvpPhase == "Green")
    Libr3mean<-subset(Libr3b, DvpPhase == "Green")
  } else if (DvpPhase == "Flower"){
    Libr1mean<-subset(Libr1b, DvpPhase == "Flower")
    Libr2mean<-subset(Libr2b, DvpPhase == "Flower")
    Libr3mean<-subset(Libr3b, DvpPhase == "Flower")
  }
  dataset_Field.mean<-SI(Libr1mean)
  dataset_Table_unmixed.mean<-SI(Libr2mean)
  dataset_Table_raw.mean<-SI(Libr3mean)
  
  
  # 2) Make a list with all Site_Date for which the specified species contains data in the field dataset
  # meta_Field_Species<-dataset_Field[which(dataset_Field$Species==SpeciesName),]
  meta_Field_Species<-eval(substitute(dataset_Field[which(dataset_Field$Species==SpeciesName),], list(s = SpeciesName)))
  List_Site_Date_Field<-levels(droplevels(meta_Field_Species$Site_Date))
  
  # 3) Check if all these attribtutes have a table equivalent. Determine intersect to make sure we have a paired dataset
  # meta_Table_Species<-dataset_Table_unmixed[which(dataset_Table_unmixed$Species==SpeciesName),] #raw might have more attributes (in case we were not able to unmix). We choose smallest dataset
  meta_Table_Species<-eval(substitute(dataset_Table_unmixed[which(dataset_Table_unmixed$Species==SpeciesName),], list(s = SpeciesName)))
  List_Site_Date_Table<-levels(droplevels(meta_Table_Species$Site_Date))
  
  List_Site_Date<-intersect(List_Site_Date_Field,List_Site_Date_Table)
  
  # 4) Determine the number of records per Site_Date attribute
  # dataset_Field_subset<-dataset_Field[dataset_Field$Species==SpeciesName & dataset_Field$Site_Date %in% List_Site_Date,]
  # dataset_Table_unmixed_subset<-dataset_Table_unmixed[dataset_Table_unmixed$Species==SpeciesName & dataset_Table_unmixed$Site_Date %in% List_Site_Date,]
  # dataset_Table_raw_subset<-dataset_Table_raw[dataset_Table_raw$Species==SpeciesName & dataset_Table_raw$Site_Date %in% List_Site_Date,]
  dataset_Field_subset<-eval(substitute(dataset_Field[which(dataset_Field$Species==SpeciesName),], list(s = SpeciesName)))
  idx_int<-dataset_Field_subset$Site_Date %in% List_Site_Date
  dataset_Field_subset<-dataset_Field_subset[idx_int,]
  dataset_Table_unmixed_subset<-eval(substitute(dataset_Table_unmixed[which(dataset_Table_unmixed$Species==SpeciesName),], list(s = SpeciesName)))
  idx_int<-dataset_Table_unmixed_subset$Site_Date %in% List_Site_Date
  dataset_Table_unmixed_subset<-dataset_Table_unmixed_subset[idx_int,]
  dataset_Table_raw_subset<-eval(substitute(dataset_Table_raw[which(dataset_Table_raw$Species==SpeciesName),], list(s = SpeciesName)))
  idx_int<-dataset_Table_raw_subset$Site_Date %in% List_Site_Date
  dataset_Table_raw_subset<-dataset_Table_raw_subset[idx_int,]
  
  library(dplyr)
  dataset_Field_subset.count<-count(dataset_Field_subset,Site_Date)
  dataset_Table_unmixed_subset.count<-count(dataset_Table_unmixed_subset,Site_Date)
  dataset_Table_raw_subset.count<-count(dataset_Table_raw_subset,Site_Date)
  
  dataset_Field_subset.totalcount<-count(dataset_Field_subset,Species)
  dataset_Table_unmixed_subset.totalcount<-count(dataset_Table_unmixed_subset,Species)
  dataset_Table_raw_subset.totalcount<-count(dataset_Table_raw_subset,Species)
  
  
  dataset_Field_subset.mean<-eval(substitute(dataset_Field.mean[which(dataset_Field.mean$Species==SpeciesName),], list(s = SpeciesName)))
  idx_int<-dataset_Field_subset.mean$Site_Date %in% List_Site_Date
  dataset_Field_subset.mean<-dataset_Field_subset.mean[idx_int,]
  dataset_Table_unmixed_subset.mean<-eval(substitute(dataset_Table_unmixed.mean[which(dataset_Table_unmixed.mean$Species==SpeciesName),], list(s = SpeciesName)))
  idx_int<-dataset_Table_unmixed_subset.mean$Site_Date %in% List_Site_Date
  dataset_Table_unmixed_subset.mean<-dataset_Table_unmixed_subset.mean[idx_int,]
  dataset_Table_raw_subset.mean<-eval(substitute(dataset_Table_raw.mean[which(dataset_Table_raw.mean$Species==SpeciesName),], list(s = SpeciesName)))
  idx_int<-dataset_Table_raw_subset.mean$Site_Date %in% List_Site_Date
  dataset_Table_raw_subset.mean<-dataset_Table_raw_subset.mean[idx_int,]
  
  dataset_Field_subset.totalcount.mean<-count(dataset_Field_subset.mean,Species)
  dataset_Table_unmixed_subset.totalcount.mean<-count(dataset_Table_unmixed_subset.mean,Species)
  dataset_Table_raw_subset.totalcount.mean<-count(dataset_Table_raw_subset.mean,Species)
  
  # 5) Determine the index of the meausurements fitting the above specified criteria 
  index_Field<-which(dataset_Field$Species==SpeciesName & dataset_Field$Site_Date %in% List_Site_Date)
  dataset_Field_subset$index<-index_Field
  index_Table_unmixed<-which(dataset_Table_unmixed$Species==SpeciesName & dataset_Table_unmixed$Site_Date %in% List_Site_Date)
  dataset_Table_unmixed_subset$index<-index_Table_unmixed
  index_Table_raw<-which(dataset_Table_raw$Species==SpeciesName & dataset_Table_raw$Site_Date %in% List_Site_Date)
  dataset_Table_raw_subset$index<-index_Table_raw
  
  # 6) Preprocess the libraries to be plotted
  # Smooth: Apply Savitsky-Golay smoothing to all libraries 
  # Libr1<- smoothSpeclib(Libr1a,method="sgolay", n=51)
  # Libr2<- smoothSpeclib(Libr2a,method="sgolay", n=51)
  # Libr3<- smoothSpeclib(Libr3a,method="sgolay", n=51)
  # Libr1mean<- smoothSpeclib(Libr1mean,method="sgolay", n=51)
  # Libr2mean<- smoothSpeclib(Libr2mean,method="sgolay", n=51)
  # Libr3mean<- smoothSpeclib(Libr3mean,method="sgolay", n=51)
  
  Libr1 <- smooth_asd_sv(Libr1a)
  Libr2 <- smooth_asd_sv(Libr2a)
  Libr3 <- smooth_asd_sv(Libr3a)
  Libr1mean <- smooth_asd_sv(Libr1mean)
  Libr2mean <- smooth_asd_sv(Libr2mean)
  Libr3mean <- smooth_asd_sv(Libr3mean)
  
  # Mask atmospheric noise bands
  mask(Libr1) <-  c(349,400,1340,1460,1780,1970,2400,2501)
  mask(Libr2) <-  c(349,400,1340,1460,1780,1970,2400,2501)
  mask(Libr3) <-  c(349,400,1340,1460,1780,1970,2400,2501)
  
  mask(Libr1mean) <-  c(349,400,1340,1460,1780,1970,2400,2501)
  mask(Libr2mean) <-  c(349,400,1340,1460,1780,1970,2400,2501)
  mask(Libr3mean) <-  c(349,400,1340,1460,1780,1970,2400,2501)
  
  #--------------------------------------------------------------------------------
  
    ########### Plot mean and sd of all individual spectra of a species
    library("RColorBrewer")
    colors = brewer.pal(3, "Dark2") # Color codes are: #1b9e77, #d95f02, #7570b3
    colors_sd<-c("#b3e2cd","#fdcdac","#cbd5e8") #a softer combination of the same colors defined above (http://colorbrewer2.org/#type=qualitative&scheme=Pastel2&n=3)
    transp<-0.3
    
    # ATTENTION: Better take average of mean spectra per day than average of all bouquets because a day with lots of measurements will then control the overall average
    # --> _MeanOverDays instead of _MeanOverBouqets
    # Drawback: IN case of only 1 day we won't have a sd interval...
     
    setwd(wd5)
    # jpeg(paste(paste(paste(SpeciesName,DvpPhase,sep="_"),jumpcorr,sep="_"),"_MeanOverBouqets_raw.jpeg",sep=""),res=300,width=7.5,height=6,units='in')
    jpeg(paste(paste(paste(SpeciesName,DvpPhase,sep="_"),jumpcorr,sep="_"),"_MeanOverBouqets_raw.jpeg",sep=""),res=600,width=7.5,height=6,units='in')
    # jpeg(paste(paste(SpeciesName,DvpPhase,sep="_"),".jpeg",sep=""),res=300,width=12.8,height=7.2,units='in')
    # jpeg("SOLGIG_Flower.jpeg",res=300,width=12.8,height=7.2,units='in')
    # pdf("SOLGIG_Flower.pdf",width=12.8,height=7.2)
    # x11(width=1280,height=720)
    plot(c(),c(),axes=FALSE,ylim=c(0,1.0),xlim=c(350,2500),ylab='',xlab='', xaxs="i", yaxs="i")
    # 1) Colour the area between the sd-borders
    wl_plot <- c(seq(400, 1340, 1), seq(1460,1780,1), seq (1970, 2400, 1))
    xx<-c(wl_plot,(rev(wl_plot)))
    Libr3_sp<-Libr3[SI(Libr3)$Species == SpeciesName,]
    idx_int <- SI(Libr3_sp)$Site_Date %in% List_Site_Date
    Libr3_plot<-Libr3_sp[idx_int,]
    mean_spec_Libr3 <- apply(Libr3_plot, FUN = mean, na.rm = TRUE)
    sd_spec_Libr3 <- apply(Libr3_plot, FUN = sd, na.rm = TRUE)
    # mean_spec_Libr3 <- apply(subset(Libr3, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = mean, na.rm = TRUE)
    # sd_spec_Libr3   <- apply(subset(Libr3, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = sd, na.rm = TRUE)
    spectra2plot_Libr3 <- rbind(spectra(mean_spec_Libr3) + spectra(sd_spec_Libr3),
                                spectra(mean_spec_Libr3),
                                spectra(mean_spec_Libr3) - spectra(sd_spec_Libr3))
    yy_Libr3<-c(spectra2plot_Libr3[1,],rev(spectra2plot_Libr3[3,]))
    # yy<-c(spectra2plot_Libr3[1,c(56:1001,1101:1446,1646:2051)],rev(spectra2plot_Libr3[3,c(56:1001,1101:1446,1646:2051)])) # the same result
    polygon(xx,yy_Libr3,col=alpha(c(colors_sd[3]),transp),border=NA,new=F)
    
    Libr1_sp<-Libr1[SI(Libr1)$Species == SpeciesName,]
    idx_int <- SI(Libr1_sp)$Site_Date %in% List_Site_Date
    Libr1_plot<-Libr1_sp[idx_int,]
    mean_spec_Libr1 <- apply(Libr1_plot, FUN = mean, na.rm = TRUE)
    sd_spec_Libr1 <- apply(Libr1_plot, FUN = sd, na.rm = TRUE)
    # mean_spec_Libr1 <- apply(subset(Libr1, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = mean, na.rm = TRUE)
    # sd_spec_Libr1   <- apply(subset(Libr1, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = sd, na.rm = TRUE)
    spectra2plot_Libr1 <- rbind(spectra(mean_spec_Libr1) + spectra(sd_spec_Libr1),
                                spectra(mean_spec_Libr1),
                                spectra(mean_spec_Libr1) - spectra(sd_spec_Libr1))
    yy_Libr1<-c(spectra2plot_Libr1[1,],rev(spectra2plot_Libr1[3,]))
    polygon(xx,yy_Libr1,col=alpha(c(colors_sd[1]),transp),border=NA,new=F)
    xx_mask1<-c(seq(1340,1460,1),rev(seq(1340,1460,1)))
    xx_mask2<-c(seq(1780,1970,1),rev(seq(1780,1970,1)))
    yy_mask1<-c(seq(0,0,length.out=length(xx_mask1)/2),seq(1,1,length.out=length(xx_mask1)/2))
    yy_mask2<-c(seq(0,0,length.out=length(xx_mask2)/2),seq(1,1,length.out=length(xx_mask2)/2))
    polygon(xx_mask1,yy_mask1,col="white",border=NA,new=F)
    polygon(xx_mask2,yy_mask2,col="white",border=NA,new=F)
    
    # 2) Plot mean spectrum and sd on top of the areas
    # plot(subset(Libr1, Species == SpeciesName & Site_Date %in% List_Site_Date), lwd=2, ylim=c(0,1), col=c(colors[1]), new=FALSE)
    # plot(subset(Libr3, Species == SpeciesName & Site_Date %in% List_Site_Date), lwd=2, ylim=c(0,1), col=c(colors[3]), new=FALSE)
    plot(Libr1_plot, lwd=2, ylim=c(0,1), col=c(colors[1]), new=FALSE)
    plot(Libr3_plot, lwd=2, ylim=c(0,1), col=c(colors[3]), new=FALSE)
    
    # 3) Plot graphical details: axes, legend, etc.

    legend(800,0.99,"Field patch reflectance",cex=1.4,lwd=2,bty="n",col=colors[1])
    legend(800,0.89,"Original table bouqet reflectance",cex=1.4,lwd=2,bty="n",col=colors[3])
    legend(2200,0.99,paste("#",dataset_Field_subset.totalcount[[2]]),cex=1.4,bty="n")
    legend(2200,0.89,paste("#",dataset_Table_raw_subset.totalcount[[2]]),cex=1.4,bty="n")
    axis(1,seq(350,2500,150),seq(350,2500,150),font=2,cex.axis=1.2)
    axis(2,seq(0.0,1.05,0.1),seq(0.0,1.05,0.1),font=2,las=2,cex.axis=1.2)
    mtext('Reflectance',side=2,font=2,line=3,cex=1.3)
    mtext('Wavelength (nm)',side=1,font=2,line=3,cex=1.3)
    # title(paste(SpeciesName, DvpPhase))
    # title("(a)",adj=0.02,cex.main=2)
    dev.off()
    
       
    setwd(wd5)    
    jpeg(paste(paste(paste(SpeciesName,DvpPhase,sep="_"),jumpcorr,sep="_"),"_MeanOverBouqets_unmixed.jpeg",sep=""),res=600,width=7.5,height=6,units='in')
    # jpeg(paste(paste(SpeciesName,DvpPhase,sep="_"),".jpeg",sep=""),res=300,width=12.8,height=7.2,units='in')
    # jpeg("SOLGIG_Flower.jpeg",res=300,width=12.8,height=7.2,units='in')
    # pdf("SOLGIG_Flower.pdf",width=12.8,height=7.2)
    # x11(width=1280,height=720)
    plot(c(),c(),axes=FALSE,ylim=c(0,1.0),xlim=c(350,2500),ylab='',xlab='', xaxs="i", yaxs="i")
    # 1) Colour the area between the sd-borders
    wl_plot <- c(seq(400, 1340, 1), seq(1460,1780,1), seq (1970, 2400, 1))
    xx<-c(wl_plot,(rev(wl_plot)))
    Libr2_sp<-Libr2[SI(Libr2)$Species == SpeciesName,]
    idx_int <- SI(Libr2_sp)$Site_Date %in% List_Site_Date
    Libr2_plot<-Libr2_sp[idx_int,]
    mean_spec_Libr2 <- apply(Libr2_plot, FUN = mean, na.rm = TRUE)
    sd_spec_Libr2 <- apply(Libr2_plot, FUN = sd, na.rm = TRUE)
    # mean_spec_Libr2 <- apply(subset(Libr2, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = mean, na.rm = TRUE)
    # sd_spec_Libr2   <- apply(subset(Libr2, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = sd, na.rm = TRUE)
    spectra2plot_Libr2 <- rbind(spectra(mean_spec_Libr2) + spectra(sd_spec_Libr2),
                                spectra(mean_spec_Libr2),
                                spectra(mean_spec_Libr2) - spectra(sd_spec_Libr2))
    yy_Libr2<-c(spectra2plot_Libr2[1,],rev(spectra2plot_Libr2[3,]))
    polygon(xx,yy_Libr2,col=alpha(c(colors_sd[2]),transp),border=NA,new=F)
    
    Libr1_sp<-Libr1[SI(Libr1)$Species == SpeciesName,]
    idx_int <- SI(Libr1_sp)$Site_Date %in% List_Site_Date
    Libr1_plot<-Libr1_sp[idx_int,]
    mean_spec_Libr1 <- apply(Libr1_plot, FUN = mean, na.rm = TRUE)
    sd_spec_Libr1 <- apply(Libr1_plot, FUN = sd, na.rm = TRUE)
    # mean_spec_Libr1 <- apply(subset(Libr1, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = mean, na.rm = TRUE)
    # sd_spec_Libr1   <- apply(subset(Libr1, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = sd, na.rm = TRUE)
    spectra2plot_Libr1 <- rbind(spectra(mean_spec_Libr1) + spectra(sd_spec_Libr1),
                                spectra(mean_spec_Libr1),
                                spectra(mean_spec_Libr1) - spectra(sd_spec_Libr1))
    yy_Libr1<-c(spectra2plot_Libr1[1,],rev(spectra2plot_Libr1[3,]))
    polygon(xx,yy_Libr1,col=alpha(c(colors_sd[1]),transp),border=NA,new=F)
    xx_mask1<-c(seq(1340,1460,1),rev(seq(1340,1460,1)))
    xx_mask2<-c(seq(1780,1970,1),rev(seq(1780,1970,1)))
    yy_mask1<-c(seq(0,0,length.out=length(xx_mask1)/2),seq(1,1,length.out=length(xx_mask1)/2))
    yy_mask2<-c(seq(0,0,length.out=length(xx_mask2)/2),seq(1,1,length.out=length(xx_mask2)/2))
    polygon(xx_mask1,yy_mask1,col="white",border=NA,new=F)
    polygon(xx_mask2,yy_mask2,col="white",border=NA,new=F)
    
    # 2) Plot mean spectrum and sd on top of the areas
    # plot(subset(Libr1, Species == SpeciesName & Site_Date %in% List_Site_Date), lwd=2, ylim=c(0,1), col=c(colors[1]), new=FALSE)
    # plot(subset(Libr2, Species == SpeciesName & Site_Date %in% List_Site_Date), lwd=2, ylim=c(0,1), col=c(colors[2]), new=FALSE)
    plot(Libr1_plot, lwd=2, ylim=c(0,1), col=c(colors[1]), new=FALSE)
    plot(Libr2_plot, lwd=2, ylim=c(0,1), col=c(colors[2]), new=FALSE)

    # 3) Plot graphical details: axes, legend, etc.
    legend(800,0.99,"Field patch reflectance",cex=1.4,lwd=2,bty="n",col=colors[1])
    legend(800,0.89,"Unmixed table bouqet reflectance",cex=1.4,lwd=2,bty="n",col=colors[2])
    legend(2200,0.99,paste("#",dataset_Field_subset.totalcount[[2]]),cex=1.4,bty="n")
    legend(2200,0.89,paste("#",dataset_Table_unmixed_subset.totalcount[[2]]),cex=1.4,bty="n")
    axis(1,seq(350,2500,150),seq(350,2500,150),font=2,cex.axis=1.2)
    axis(2,seq(0.0,1.05,0.1),seq(0.0,1.05,0.1),font=2,las=2,cex.axis=1.2)
    mtext('Reflectance',side=2,font=2,line=3,cex=1.3)
    mtext('Wavelength (nm)',side=1,font=2,line=3,cex=1.3)
    # title(paste(SpeciesName, DvpPhase))
    # title("(b)",adj=0.02,cex.main=2)
    dev.off()
    
    
    setwd(wd5)
    jpeg(paste(paste(paste(SpeciesName,DvpPhase,sep="_"),jumpcorr,sep="_"),"_MeanOverDays_raw.jpeg",sep=""),res=600,width=7.5,height=6,units='in')
    # x11(width=1280,height=720)
    plot(c(),c(),axes=FALSE,ylim=c(0,1.0),xlim=c(350,2500),ylab='',xlab='', xaxs="i", yaxs="i")
    # 1) Colour the area between the sd-borders
    wl_plot <- c(seq(400, 1340, 1), seq(1460,1780,1), seq (1970, 2400, 1))
    xx<-c(wl_plot,(rev(wl_plot)))
    Libr3_sp<-Libr3mean[SI(Libr3mean)$Species == SpeciesName,]
    idx_int <- SI(Libr3_sp)$Site_Date %in% List_Site_Date
    Libr3_plot<-Libr3_sp[idx_int,]
    mean_spec_Libr3 <- apply(Libr3_plot, FUN = mean, na.rm = TRUE)
    sd_spec_Libr3 <- apply(Libr3_plot, FUN = sd, na.rm = TRUE)
    # mean_spec_Libr3 <- apply(subset(Libr3, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = mean, na.rm = TRUE)
    # sd_spec_Libr3   <- apply(subset(Libr3, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = sd, na.rm = TRUE)
    spectra2plot_Libr3 <- rbind(spectra(mean_spec_Libr3) + spectra(sd_spec_Libr3),
                                spectra(mean_spec_Libr3),
                                spectra(mean_spec_Libr3) - spectra(sd_spec_Libr3))
    yy_Libr3<-c(spectra2plot_Libr3[1,],rev(spectra2plot_Libr3[3,]))
    # yy<-c(spectra2plot_Libr3[1,c(56:1001,1101:1446,1646:2051)],rev(spectra2plot_Libr3[3,c(56:1001,1101:1446,1646:2051)])) # the same result
    polygon(xx,yy_Libr3,col=alpha(c(colors_sd[3]),transp),border=NA,new=F)
    
    Libr1_sp<-Libr1mean[SI(Libr1mean)$Species == SpeciesName,]
    idx_int <- SI(Libr1_sp)$Site_Date %in% List_Site_Date
    Libr1_plot<-Libr1_sp[idx_int,]
    mean_spec_Libr1 <- apply(Libr1_plot, FUN = mean, na.rm = TRUE)
    sd_spec_Libr1 <- apply(Libr1_plot, FUN = sd, na.rm = TRUE)
    # mean_spec_Libr1 <- apply(subset(Libr1, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = mean, na.rm = TRUE)
    # sd_spec_Libr1   <- apply(subset(Libr1, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = sd, na.rm = TRUE)
    spectra2plot_Libr1 <- rbind(spectra(mean_spec_Libr1) + spectra(sd_spec_Libr1),
                                spectra(mean_spec_Libr1),
                                spectra(mean_spec_Libr1) - spectra(sd_spec_Libr1))
    yy_Libr1<-c(spectra2plot_Libr1[1,],rev(spectra2plot_Libr1[3,]))
    polygon(xx,yy_Libr1,col=alpha(c(colors_sd[1]),transp),border=NA,new=F)
    xx_mask1<-c(seq(1340,1460,1),rev(seq(1340,1460,1)))
    xx_mask2<-c(seq(1780,1970,1),rev(seq(1780,1970,1)))
    yy_mask1<-c(seq(0,0,length.out=length(xx_mask1)/2),seq(1,1,length.out=length(xx_mask1)/2))
    yy_mask2<-c(seq(0,0,length.out=length(xx_mask2)/2),seq(1,1,length.out=length(xx_mask2)/2))
    polygon(xx_mask1,yy_mask1,col="white",border=NA,new=F)
    polygon(xx_mask2,yy_mask2,col="white",border=NA,new=F)
    
    # 2) Plot mean spectrum and sd on top of the areas
    # plot(subset(Libr1, Species == SpeciesName & Site_Date %in% List_Site_Date), lwd=2, ylim=c(0,1), col=c(colors[1]), new=FALSE)
    # plot(subset(Libr3, Species == SpeciesName & Site_Date %in% List_Site_Date), lwd=2, ylim=c(0,1), col=c(colors[3]), new=FALSE)
    plot(Libr1_plot, lwd=2, ylim=c(0,1), col=c(colors[1]), new=FALSE)
    plot(Libr3_plot, lwd=2, ylim=c(0,1), col=c(colors[3]), new=FALSE)
    
    # 3) Plot graphical details: axes, legend, etc.
    legend(800,0.99,"Field patch reflectance",cex=1.4,lwd=2,bty="n",col=colors[1])
    legend(800,0.89,"Original table bouquet reflectance",cex=1.4,lwd=2,bty="n",col=colors[3])
    legend(2200,0.99,paste("#",dataset_Field_subset.totalcount.mean[[2]]),cex=1.4,bty="n")
    legend(2200,0.89,paste("#",dataset_Table_raw_subset.totalcount.mean[[2]]),cex=1.4,bty="n")
    axis(1,seq(350,2500,150),seq(350,2500,150),font=2,cex.axis=1.2)
    axis(2,seq(0.0,1.05,0.1),seq(0.0,1.05,0.1),font=2,las=2,cex.axis=1.2)
    mtext('Reflectance',side=2,font=2,line=3,cex=1.3)
    mtext('Wavelength (nm)',side=1,font=2,line=3,cex=1.3)
    # title(paste(SpeciesName, DvpPhase))
    # title("(a)",adj=0.02,cex.main=2)
    dev.off()
    
    
    setwd(wd5)
    jpeg(paste(paste(paste(SpeciesName,DvpPhase,sep="_"),jumpcorr,sep="_"),"_MeanOverDays_unmixed.jpeg",sep=""),res=600,width=7.5,height=6,units='in')
    # x11(width=1280,height=720)
    plot(c(),c(),axes=FALSE,ylim=c(0,1.0),xlim=c(350,2500),ylab='',xlab='', xaxs="i", yaxs="i")
    # 1) Colour the area between the sd-borders
    wl_plot <- c(seq(400, 1340, 1), seq(1460,1780,1), seq (1970, 2400, 1))
    xx<-c(wl_plot,(rev(wl_plot)))
    Libr2_sp<-Libr2mean[SI(Libr2mean)$Species == SpeciesName,]
    idx_int <- SI(Libr2_sp)$Site_Date %in% List_Site_Date
    Libr2_plot<-Libr2_sp[idx_int,]
    mean_spec_Libr2 <- apply(Libr2_plot, FUN = mean, na.rm = TRUE)
    sd_spec_Libr2 <- apply(Libr2_plot, FUN = sd, na.rm = TRUE)
    # mean_spec_Libr2 <- apply(subset(Libr2, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = mean, na.rm = TRUE)
    # sd_spec_Libr2   <- apply(subset(Libr2, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = sd, na.rm = TRUE)
    spectra2plot_Libr2 <- rbind(spectra(mean_spec_Libr2) + spectra(sd_spec_Libr2),
                                spectra(mean_spec_Libr2),
                                spectra(mean_spec_Libr2) - spectra(sd_spec_Libr2))
    yy_Libr2<-c(spectra2plot_Libr2[1,],rev(spectra2plot_Libr2[3,]))
    polygon(xx,yy_Libr2,col=alpha(c(colors_sd[2]),transp),border=NA,new=F)
    
    Libr1_sp<-Libr1mean[SI(Libr1mean)$Species == SpeciesName,]
    idx_int <- SI(Libr1_sp)$Site_Date %in% List_Site_Date
    Libr1_plot<-Libr1_sp[idx_int,]
    mean_spec_Libr1 <- apply(Libr1_plot, FUN = mean, na.rm = TRUE)
    sd_spec_Libr1 <- apply(Libr1_plot, FUN = sd, na.rm = TRUE)
    # mean_spec_Libr1 <- apply(subset(Libr1, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = mean, na.rm = TRUE)
    # sd_spec_Libr1   <- apply(subset(Libr1, Species == SpeciesName & Site_Date %in% List_Site_Date), FUN = sd, na.rm = TRUE)
    spectra2plot_Libr1 <- rbind(spectra(mean_spec_Libr1) + spectra(sd_spec_Libr1),
                                spectra(mean_spec_Libr1),
                                spectra(mean_spec_Libr1) - spectra(sd_spec_Libr1))
    yy_Libr1<-c(spectra2plot_Libr1[1,],rev(spectra2plot_Libr1[3,]))
    polygon(xx,yy_Libr1,col=alpha(c(colors_sd[1]),transp),border=NA,new=F)
    xx_mask1<-c(seq(1340,1460,1),rev(seq(1340,1460,1)))
    xx_mask2<-c(seq(1780,1970,1),rev(seq(1780,1970,1)))
    yy_mask1<-c(seq(0,0,length.out=length(xx_mask1)/2),seq(1,1,length.out=length(xx_mask1)/2))
    yy_mask2<-c(seq(0,0,length.out=length(xx_mask2)/2),seq(1,1,length.out=length(xx_mask2)/2))
    polygon(xx_mask1,yy_mask1,col="white",border=NA,new=F)
    polygon(xx_mask2,yy_mask2,col="white",border=NA,new=F)
    
    # 2) Plot mean spectrum and sd on top of the areas
    # plot(subset(Libr1, Species == SpeciesName & Site_Date %in% List_Site_Date), lwd=2, ylim=c(0,1), col=c(colors[1]), new=FALSE)
    # plot(subset(Libr2, Species == SpeciesName & Site_Date %in% List_Site_Date), lwd=2, ylim=c(0,1), col=c(colors[2]), new=FALSE)
    plot(Libr1_plot, lwd=2, ylim=c(0,1), col=c(colors[1]), new=FALSE)
    plot(Libr2_plot, lwd=2, ylim=c(0,1), col=c(colors[2]), new=FALSE)

    # 3) Plot graphical details: axes, legend, etc.
    legend(800,0.99,"Field patch reflectance",cex=1.4,lwd=2,bty="n",col=colors[1])
    legend(800,0.89,"Unmixed table bouquet reflectance",cex=1.4,lwd=2,bty="n",col=colors[2])
    legend(2200,0.99,paste("#",dataset_Field_subset.totalcount.mean[[2]]),cex=1.4,bty="n")
    legend(2200,0.89,paste("#",dataset_Table_unmixed_subset.totalcount.mean[[2]]),cex=1.4,bty="n")
    axis(1,seq(350,2500,150),seq(350,2500,150),font=2,cex.axis=1.2)
    axis(2,seq(0.0,1.05,0.1),seq(0.0,1.05,0.1),font=2,las=2,cex.axis=1.2)
    mtext('Reflectance',side=2,font=2,line=3,cex=1.3)
    mtext('Wavelength (nm)',side=1,font=2,line=3,cex=1.3)
    # title(paste(SpeciesName, DvpPhase))
    # title("(b)",adj=0.02,cex.main=2)
    dev.off()
}
# Apply spectral_plots_separate_raw_unmix to all species measured
spectral_plots_separate_apply <- function(wd4, wd5, jumpcorr) {
  setwd(wd4)
  datasets<-readRDS(paste("datasets_",".rds",sep=jumpcorr))  
  
  levels <- levels(droplevels(SI(datasets$libr_T_unmixed.intersect)$Sp_Phase))
  SpeciesNames <- substr(levels,1,6)
  for (i in 1:length(levels)){
      if (grepl("Flower",levels[i],fixed=T)=="TRUE"){
        spectral_plots_separate_raw_unmix(datasets=datasets, wd5, jumpcorr=jumpcorr, SpeciesNames[i], DvpPhase="Flower")
      } else if (grepl("Green",levels[i],fixed=T)=="TRUE"){
        spectral_plots_separate_raw_unmix(datasets=datasets, wd5, jumpcorr=jumpcorr, SpeciesNames[i], DvpPhase="Green")
      }
  }
}

# Create plot with mean of the 42 paired field-table spectra and indicate bands used for VI's
spectral_plots_VIbands <- function (wd4, wd6, jumpcorr) {
  
  # 1) Specify datasets for visualisation
  setwd(wd4)
  datasets<-readRDS(paste("datasets_",".rds",sep=jumpcorr))  
  # Mean per day:
  Libr1mean<-datasets$libr_F.intersect.mean_orig
  Libr2mean<-datasets$libr_T_unmixed.intersect.mean_orig
  Libr3mean<-datasets$libr_T_raw.intersect.mean_orig
  dataset_Field.mean<-SI(Libr1mean)
  dataset_Table_unmixed.mean<-SI(Libr2mean)
  dataset_Table_raw.mean<-SI(Libr3mean)
  
  # 2) Make a list with all Site_Date contained in the field dataset
  List_Site_Date_Field<-levels(droplevels(dataset_Field.mean$Site_Date))
  
  # 3) Check if all these attribtutes have a table equivalent. Determine intersect to make sure we have a paired dataset
  List_Site_Date_Table<-levels(droplevels(dataset_Table_unmixed.mean$Site_Date))
  
  List_Site_Date<-intersect(List_Site_Date_Field,List_Site_Date_Table)
  
  # 4) Determine the number of records per Site_Date attribute
  library(dplyr)
  dataset_Field.mean.count<-count(dataset_Field.mean,Sp_Phase)
  dataset_Table_unmixed.mean.count<-count(dataset_Table_unmixed.mean,Sp_Phase)
  dataset_Table_raw_subset.count<-count(dataset_Table_raw.mean,Sp_Phase)
  
  dataset_Field.mean.totalcount<-count(dataset_Field.mean)
  dataset_Table_unmixed.mean.totalcount<-count(dataset_Table_unmixed.mean)
  dataset_Table_raw.mean.totalcount<-count(dataset_Table_raw.mean)
  
  # 5) Determine the index of the meausurements fitting the above specified criteria 
  # nvt
  
  # 6) Preprocess the libraries to be plotted
  # Smooth: Apply Savitsky-Golay smoothing to all libraries 
  # Libr1mean<- smoothSpeclib(Libr1mean,method="sgolay", n=51)
  # Libr2mean<- smoothSpeclib(Libr2mean,method="sgolay", n=51)
  # Libr3mean<- smoothSpeclib(Libr3mean,method="sgolay", n=51)
  
  Libr1mean <- smooth_asd_sv(Libr1mean)
  Libr2mean <- smooth_asd_sv(Libr2mean)
  Libr3mean <- smooth_asd_sv(Libr3mean)
  
  # Mask atmospheric noise bands
  mask(Libr1mean) <-  c(349,400,1340,1460,1780,1970,2400,2501)
  mask(Libr2mean) <-  c(349,400,1340,1460,1780,1970,2400,2501)
  mask(Libr3mean) <-  c(349,400,1340,1460,1780,1970,2400,2501)
  
  #--------------------------------------------------------------------------------
  
  ########### Plot mean and sd of all mean spectra of species
  library("RColorBrewer")
  colors = brewer.pal(3, "Dark2") # Color codes are: #1b9e77, #d95f02, #7570b3
  colors_sd<-c("#b3e2cd","#fdcdac","#cbd5e8") #a softer combination of the same colors defined above (http://colorbrewer2.org/#type=qualitative&scheme=Pastel2&n=3)
  colors_VI_pal<-brewer.pal(9,"Pastel1")
  colors_VI<-colors_VI_pal[c(3,1,2,7)] # chl - carotenoids - water - LMA
  # colors_VI_pal<-brewer.pal(8,"Accent")
  # colors_VI<-colors_VI_pal[c(1,2,5,8)]
  transp<-0.3
  
  # ATTENTION: Better take average of mean spectra per day than average of all bouquets because a day with lots of measurements will then control the overall average
  # --> _MeanOverDays instead of _MeanOverBouqets
  # Drawback: IN case of only 1 day we won't have a sd interval...
  # Here, we take the mean of MeanOverDays
  
  setwd(wd6)
  jpeg("MeanSpectra_VIbands.jpeg",res=300,width=12.8,height=7.2,units='in')
  # x11(width=1280,height=720)
  plot(c(),c(),axes=FALSE,ylim=c(0,1.0),xlim=c(350,2500),ylab='',xlab='', xaxs="i", yaxs="i")
  
  # 4) Plot VI wavelength positions
  abline(v=c(925,710,750,705,550,754,709,681,531,570,860,720),col=colors_VI[1],lty=3,lwd=2) #chlorophyll, nitrogen and phosphorus
  abline(v=c(539,490,807),col=colors_VI[2],lty=3,lwd=2) # carotenoids
  abline(v=c(860,1240,858,2130,531,570,700,670,800),col=colors_VI[3],lty=3,lwd=2) # water
  abline(v=c(1662,1732,1550,1750,1649,1722),col=colors_VI[4],lty=3,lwd=2) # LMA
  
  
  # 1) Colour the area between the sd-borders
  wl_plot <- c(seq(400, 1340, 1), seq(1460,1780,1), seq (1970, 2400, 1))
  xx<-c(wl_plot,(rev(wl_plot)))
  mean_spec_Libr3 <- apply(Libr3mean, FUN = mean, na.rm = TRUE)
  sd_spec_Libr3 <- apply(Libr3mean, FUN = sd, na.rm = TRUE)
  spectra2plot_Libr3 <- rbind(spectra(mean_spec_Libr3) + spectra(sd_spec_Libr3),
                              spectra(mean_spec_Libr3),
                              spectra(mean_spec_Libr3) - spectra(sd_spec_Libr3))
  yy_Libr3<-c(spectra2plot_Libr3[1,],rev(spectra2plot_Libr3[3,]))
  mean_spec_Libr2 <- apply(Libr2mean, FUN = mean, na.rm = TRUE)
  sd_spec_Libr2 <- apply(Libr2mean, FUN = sd, na.rm = TRUE)
  spectra2plot_Libr2 <- rbind(spectra(mean_spec_Libr2) + spectra(sd_spec_Libr2),
                              spectra(mean_spec_Libr2),
                              spectra(mean_spec_Libr2) - spectra(sd_spec_Libr2))
  yy_Libr2<-c(spectra2plot_Libr2[1,],rev(spectra2plot_Libr2[3,]))
  mean_spec_Libr1 <- apply(Libr1mean, FUN = mean, na.rm = TRUE)
  sd_spec_Libr1 <- apply(Libr1mean, FUN = sd, na.rm = TRUE)
  spectra2plot_Libr1 <- rbind(spectra(mean_spec_Libr1) + spectra(sd_spec_Libr1),
                              spectra(mean_spec_Libr1),
                              spectra(mean_spec_Libr1) - spectra(sd_spec_Libr1))
  yy_Libr1<-c(spectra2plot_Libr1[1,],rev(spectra2plot_Libr1[3,]))
  # First plot white polygons on top of the vertical lines
  polygon(xx,yy_Libr3,col="white",border=NA,new=F)
  polygon(xx,yy_Libr2,col="white",border=NA,new=F)
  polygon(xx,yy_Libr1,col="white",border=NA,new=F)
  # Then plot coloured polygons on top
  polygon(xx,yy_Libr3,col=alpha(c(colors_sd[3]),transp),border=NA,new=F)
  polygon(xx,yy_Libr2,col=alpha(c(colors_sd[2]),transp),border=NA,new=F)
  polygon(xx,yy_Libr1,col=alpha(c(colors_sd[1]),transp),border=NA,new=F)
  # Mask the atmospheric noise regions
  xx_mask1<-c(seq(1340,1460,1),rev(seq(1340,1460,1)))
  xx_mask2<-c(seq(1780,1970,1),rev(seq(1780,1970,1)))
  yy_mask1<-c(seq(0,0,length.out=length(xx_mask1)/2),seq(1,1,length.out=length(xx_mask1)/2))
  yy_mask2<-c(seq(0,0,length.out=length(xx_mask2)/2),seq(1,1,length.out=length(xx_mask2)/2))
  polygon(xx_mask1,yy_mask1,col="white",border=NA,new=F)
  polygon(xx_mask2,yy_mask2,col="white",border=NA,new=F)
  
  # 2) Plot mean spectrum and sd on top of the areas
  plot(Libr1mean, lwd=2, ylim=c(0,1), col=c(colors[1]), new=FALSE)
  plot(Libr2mean, lwd=2, ylim=c(0,1), col=c(colors[2]), new=FALSE)
  plot(Libr3mean, lwd=2, ylim=c(0,1), col=c(colors[3]), new=FALSE)
  
  # 3) Plot graphical details: axes, legend, etc.
  legend(1850,0.99,"Field spectra",cex=1.2,lwd=2,col=colors[1],bg="white",box.col="white") #,bty="n"
  legend(1850,0.92,"Unmixed table spectra",cex=1.2,lwd=2,col=colors[2],bg="white",box.col="white")
  legend(1850,0.85,"Original table spectra",cex=1.2,lwd=2,col=colors[3],bg="white",box.col="white")
  legend(2350,0.99,paste("#",dataset_Field.mean.totalcount),cex=1.2,bg="white",box.col="white")
  legend(2350,0.92,paste("#",dataset_Table_unmixed.mean.totalcount),cex=1.2,bg="white",box.col="white")
  legend(2350,0.85,paste("#",dataset_Table_raw.mean.totalcount),cex=1.2,bg="white",box.col="white")
  
  legend(1850,0.70,"Chlorophyll, N and P",cex=1.2,lwd=2.5,lty=3,col=colors_VI[1],bg="white",box.col="white")
  legend(1850,0.62,"Carotenoids",cex=1.2,lwd=2.5,lty=3,col=colors_VI[2],bg="white",box.col="white")
  legend(1850,0.55,"Water",cex=1.2,lwd=2.5,lty=3,col=colors_VI[3],bg="white",box.col="white")
  legend(1850,0.48,"LMA",cex=1.2,lwd=2.5,lty=3,col=colors_VI[4],bg="white",box.col="white")
  
  axis(1,seq(350,2500,150),seq(350,2500,150),font=2,cex.axis=1.2)
  axis(2,seq(0.0,1.05,0.1),seq(0.0,1.05,0.1),font=2,las=2,cex.axis=1.2)
  mtext('Reflectance',side=2,font=2,line=3,cex=1.3)
  mtext('Wavelength (nm)',side=1,font=2,line=3,cex=1.3)

  dev.off()
}

#--------------------------------------------------------------------------------
wd1 = "C:/Users/u0091812/Box Sync/03. Black table procedure/Figures/Field vs. Table/2016+2018/Mean spectra per day"
wd2 = "C:/Users/u0091812/Box Sync/03. Black table procedure/Figures/Field vs. Table/2016+2018/Individual spectra"
wd3 = "C:/Users/u0091812/Box Sync/03. Black table procedure/Figures/Field vs. Table/2016+2018/Mean spectra overall"
wd4 = "C:/Users/u0091812/Box Sync/03. Black table procedure/Data/Spectral data"
wd5 = "C:/Users/u0091812/Box Sync/03. Black table procedure/Figures/Field vs. Table/2016+2018/Mean spectra overall - raw and unmix separate"
wd6 = "C:/Users/u0091812/Box Sync/03. Black table procedure/Figures/Field vs. Table/2016+2018/Mean and VI bands"
wd6 = "C:/Users/u0091812/Box Sync/03. Black table procedure/Figures/Field vs. Table/2016+2018/test"


spectral_plots_apply(wd1, wd2, wd3, wd4, "Nocorr", VIbands="No")
# spectral_plots_apply(wd1, wd2, wd3, wd4, "jcorr_add_SWIR1", VIbands="No")
# spectral_plots_apply(wd1, wd2, wd3, wd4, "jcorr_mult_SWIR1", VIbands="No")
# spectral_plots_apply(wd1, wd2, wd3, wd4, "jcorr_mult_SWIR12", VIbands="NoNo")
# spectral_plots_apply(wd1, wd2, wd3, wd4, "jcorr_mult_VNIR", VIbands="No")
# spectral_plots_apply("wd1, wd2, wd3, wd4, Nocorr_frac0.9", VIbands="No")

spectral_plots_separate_apply(wd4, wd5, "Nocorr")

#  Plot figure on which bands are indicated that are used in calculation of VI's
# spectral_plots_apply("wd1, wd2, wd3, wd4, "Nocorr", VIbands="Yes")
spectral_plots_VIbands(wd4, wd6, "Nocorr")

