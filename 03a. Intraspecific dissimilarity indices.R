###################################################################################
###                                                                             ###
###       03a. Intraspecific similarities based on (dis)similarity indices      ###
###                                                                             ###
###################################################################################

### The objective of this script is to QUANTITATIVELY evaluate the performance of the black table
# WITHIN SPECIES (intraspecific) we investigate
#   1) to what extent (unmixed) table spectra resemble field spectra: shape and distance (dis)similarity indices on entire spectra -> does not provide info if no comparison
#   2) the effect of vegetation fraction on the unmixing quality: linear model (dis)similarity indices ~ fraction_veg

### The input of these steps are the libraries created in "01. Black table data":
#   For assessing (dis)similarities WIHTIN species we need:
#       - AVERAGE FIELD spectrum per SPECIES, SITE, DATE, AND CONDITION (Patch - Dvpmt phase):
#         libr_F.intersect.mean
#       - INDIVIDUAL table (average bouquet) spectra:
#         libr_T_unmixed.intersect
#         libr_T_raw.intersect

# ! These signatures have been preprocessed with following specifications:
#   - smoothing: savitsky-golay filter with n = 51 and 101 for asd and sv respectively
#   - removal of atmospheric noise windows:  c(349,400,1340,1460,1780,1970,2400,2501)

# Script by Elisa Van Cleemput, 2018
#--------------------------------------------------------------------------------

# Clean workspace

rm(list=ls())

#--------------------------------------------------------------------------------
######## WITHIN species: mean of field + individual table spectra
# calculate indices for both raw and unmixed spectra

library(resemble) # SAM calculation with fDiss function
library(tempR) # City-block distances
library(hsdar)
library(stats)
library(cowplot) # get_legend
# library(lme4) # lmer function
library(lmerTest) # lmer function from lme4 package + p-values
# library(nlme) # lme function
library(MuMIn) # r.squaredGLMM
library(car) #Anova, qqp
library(randomcoloR) # distinctColorPalette

#--------------------------------------------------------------------------------
# FUNCTIONS

# Self-defined distance functions. Do they coincide with results from package?
# https://stats.stackexchange.com/questions/270492/coding-mahalanobis-and-manhattan-distance-with-r
euclideanDist <- function(a, b){
  d<-(a-b)^2
  d<-sum(d)
  d<-sqrt(d)
  # return(d)
}
euclideanDist2 <- function(a, b){
  d<-(a-b)^2
  d<-sum(d)
  d<-d/length(a)
  d<-sqrt(d)
  # return(d)
}
manhattanDist_tempR<-function(a,b){ # algorithm in tempR function. No theoretical equivalent found
  d.mat<-abs(a-b)
  d<-sum(d.mat)
  return(list(distance=d, size= prod(dim(d.mat))))
}
manhattanDist<-function(a,b){ # 
  # If d is divided by sum(a,b), output = Bray-Curtis or Sorensen dissimilarity 
  # If not: output : Manhattan distance
  d.mat<-abs(a-b)
  d<-sum(d.mat)
  # d<-d/sum(a,b)
  return(d)
}

# Functions to calculate intraspecific (dis)similarities --> output = meta_T_raw.intersect and meta_T_unmixed.intersect
intraDissim_proc <- function(datasets,spectra_T_proc.intersect,meta_T_proc.intersect,proc){
  meta_F.intersect.mean<-SI(datasets$libr_F.intersect.mean_orig)
  libr_F.intersect.mean<-datasets$libr_F.intersect.mean
  if (proc == "unmixed") {
    libr_T_proc.intersect<-datasets$libr_T_unmixed.intersect
  } else if (proc == "raw") {
    libr_T_proc.intersect<-datasets$libr_T_raw.intersect
  }
  
  # Created empty vectors to be filled with (dis)similarity values
  error_proc<-spectra_T_proc.intersect
  scm_proc<-rep(0,nrow(error_proc))
  sam_proc<-rep(0,nrow(error_proc))
  sam2_proc<-rep(0,nrow(error_proc))
  sid_proc<-rep(0,nrow(error_proc))
  euc_proc<-rep(0,nrow(error_proc))
  euc2_proc<-rep(0,nrow(error_proc))
  euc3_proc<-rep(0,nrow(error_proc))
  manabs_proc<-rep(0,nrow(error_proc))
  manabs2_proc<-rep(0,nrow(error_proc))
  manabs3_proc<-rep(0,nrow(error_proc))
  manrel_proc<-rep(0,nrow(error_proc))
  manrel2_proc<-rep(0,nrow(error_proc))
  bcdabs_proc<-rep(0,nrow(error_proc))
  bcdrel_proc<-rep(0,nrow(error_proc))
  test<-rep(0,nrow(error_proc))
  
  Fraction_veg_proc <- rep(0,nrow(error_proc))

  # ------ proc table spectra
  for (i in 1:nrow(meta_T_proc.intersect)){
    # Determine field measurement corresponding to table measurement
    meta_Field_spec =meta_F.intersect.mean[meta_F.intersect.mean$Sp_SiteDate_Phase %in% meta_T_proc.intersect[i,"Sp_SiteDate_Phase"],]
    idx <- SI(libr_F.intersect.mean)$Sp_SiteDate_Phase == meta_Field_spec$Sp_SiteDate_Phase
    Field_spec <- libr_F.intersect.mean[idx,]
    # Field_spec<-subset(libr_F.intersect.mean,Sp_SiteDate_Phase == meta_Field_spec$Sp_SiteDate_Phase)
    spectra_F<-spectra(Field_spec)
    
    # Error: absolute amplitude deviation of table spectrum from field spectrum at each wavelength
    # perform smoothing first!!!!
    error_proc[i,] <- spectra_T_proc.intersect[i,] - spectra_F
    
    # SCM: Spectral Correlation Measure (shape measure) + t-test
    # package: resemble (corDiss function) --> Is this the correct function????????
    scm_proc[i] <- corDiss(spectra_F,t(as.matrix(spectra_T_proc.intersect[i,])),center=F,scaled=F) # resemble  package: input is matrix or data.frame
    
    # SAM = Spectral Angle Meausure: between table spectrum and corresponding field spectrum 
    # 2 packages available with the same results: resemble (fDiss function) and hsdar (sam function)
    # sam_proc[i] <- fDiss(t(as.matrix(spectra_T_proc.intersect[i,])),spectra_F,method="cosine",center=F,scaled=F) # resemble  package: input is matrix or data.frame
    sam_proc[i] <- fDiss(spectra_F,t(as.matrix(spectra_T_proc.intersect[i,])),method="cosine",center=F,scaled=F) # resemble  package: input is matrix or data.frame
    # Or, using hsdar package:
    
# Table_spec<-subset(libr_T_proc.intersect,Species_Nr_Site_Date == meta_T_proc.intersect$Species_Nr_Site_Date[i])
    # interest<- meta_T_proc.intersect[i,"Species_Nr_Site_Date"]
    # Table_spec<-hsdar::subset(libr_T_proc.intersect,Species_Nr_Site_Date == interest)
    # sam2_proc[i] <- sam(Field_spec,Table_spec) #hsdar package = sam_rawi 
    # sam2_proc[i] <- sam(Table_spec,Field_spec) #hsdar package = sam_raw
    
    
    # SID =  Spectral Information Divergence
    # package: resemble (sid function)
    sid_proc[i] <- sid(spectra_F,t(as.matrix(spectra_T_proc.intersect[i,])),mode="density",center=F,scaled=F)$sid # resemble  package: input is matrix or data.frame
    
    #EUCLIDEAN DISTANCE
    euc2_proc[i] <- dist(rbind(spectra_F,t(as.matrix(spectra_T_proc.intersect[i,]))),method="euclidean")
    
    
    
    # EUCLIDEAN DISTANCE: different methods deliver different results!
    euc_proc[i] <- euclideanDist(spectra_F,t(as.matrix(spectra_T_proc.intersect[i,])))
    euc2_proc[i] <- dist(rbind(spectra_F,t(as.matrix(spectra_T_proc.intersect[i,]))),method="euclidean")
    euc3_proc[i] <- stats::dist(rbind(spectra_F,t(as.matrix(spectra_T_proc.intersect[i,]))), method = "euclidean")
    # # The fDiss function does not involve correct formula: no sum under the square root
    # euc4_proc[i] <- fDiss(spectra_F,t(as.matrix(spectra_T_proc.intersect[i,])),method="euclid",center=F,scaled=F) # resemble  package: input is matrix or data.frame
    # # Alternative formulas according to Van der Meer et al. 2006 ---> different results!!!
    # euc5_proc[i] <- 2*sqrt(1-cos(sam_proc[i]))
    # euc6_proc[i] <- 2*sin(sam_proc[i]/2)
    # # # Wang et al. 2018 average squared difference under the sqrt; they hence calculate RMS difference
    # euc7_proc[i] <- euclideanDist2(spectra_F,t(as.matrix(spectra_T_proc.intersect[i,])))
    
    
    # Sorensen similarity or Curtis-Bray dissimilarity AND Manhattan or city-block DISTANCE
    manabs_proc[i] <- dist(rbind(spectra_F,t(as.matrix(spectra_T_proc.intersect[i,]))),method="manhattan") # = absolute Manhattan distance 
    manabs2_proc[i]<-manhattanDist(spectra_F,t(as.matrix(spectra_T_proc.intersect[i,]))) 
    
    manrel_proc[i] <- dist(rbind(spectra_F/sum(spectra_F),t(as.matrix(spectra_T_proc.intersect[i,]))/sum(t(as.matrix(spectra_T_proc.intersect[i,])))),method="manhattan") # = relative Manhattan distance 
    manrel2_proc[i]<-manhattanDist(spectra_F/sum(spectra_F),t(as.matrix(spectra_T_proc.intersect[i,]))/sum(t(as.matrix(spectra_T_proc.intersect[i,])))) 
    
    bcdabs_proc[i]<-manabs_proc[i]/(sum(spectra_F,t(as.matrix(spectra_T_proc.intersect[i,])))) # absolute Bray-Curtis or Sorenson dissimilarity
    bcdrel_proc[i]<-manrel_proc[i]/2 #relative  Bray-Curtis or Sorenson dissimilarity
    
    
    manabs3_proc[i] <- dist.city.block(spectra_F,t(as.matrix(spectra_T_proc.intersect[i,]))) #tempR package --> different result
    Md<-manhattanDist_tempR(spectra_F,t(as.matrix(spectra_T_proc.intersect[i,]))) 
    test[i]<-Md$distance/Md$size # = man_proc (= tempR algorithm)
    
    Fraction_veg_proc[i]<-meta_T_proc.intersect[i,"Fraction_veg"]
  }
  
  meta_T_proc.intersect$sam<-sam_proc
  meta_T_proc.intersect$scm<-scm_proc
  meta_T_proc.intersect$sid<-sid_proc
  meta_T_proc.intersect$euc<-euc2_proc
  meta_T_proc.intersect$manabs<-manabs_proc
  meta_T_proc.intersect$manrel<-manrel_proc
  meta_T_proc.intersect$bcdabs<-bcdabs_proc
  meta_T_proc.intersect$bcdrel<-bcdrel_proc
  
  return(meta_T_proc.intersect)
}

# Function that plots (dis)similarities values in function of cover fraction + save fig
Graph_FractionDis <- function(wd2,jumpcorr,dataset,dataset_twin,dist_proc,dist,title){
  ### ordinary linear model
  regr<-lm(dataset[[dist]] ~ dataset[["Fraction_veg"]])
  # regr<-lm(datasetdist~datasetFracveg)  # ,data=dataset
  setwd(wd2)
  jpeg(paste(paste("Intraspecific_","_",sep=dist_proc),"_LinAssum.jpeg",sep=jumpcorr),res=600,width=12,height=12,units='in')
  par(mfrow=c(2,2))
  plot(regr) # same as: plot(gvlma::gvlma.lm(regr,alphalevel=0.05))
  dev.off()
  # str(gvlma::gvlma.lm(regr,alphalevel=0.05)$GlobalTest)
      # GlobalStat4:       global test about appropriateness of linear model assumptions
      # DirectionalStat1:  detect Skewness = measure of the lack of symmetry
      # DirectionalStat2:  detect Kurtosis = measure of whether the data are heavy-tailed or light-tailed relative to normal distribution
      # DirectionalStat3:  detect nonlinear Link function
      # DirectionalStat4:  detect Heteroscedasticity
      #  Pena, EA and Slate, EH (2006). "Global validation of linear model assumptions," J.\ Amer.\ Statist.\ Assoc., 101(473):341-354.
  # hist(regr$residuals)
  p_slope_regr<-summary(regr)$coefficients[2,4] # same as: summary(gvlma::gvlma.lm(regr,alphalevel=0.05))
  R2 <- summary(regr)$r.squared
  if (p_slope_regr < 0.001) {
    p_slope_regr <- "< 0.001"
  } else {
    p_slope_regr <- paste("=",round(p_slope_regr,digits=3),sep=" ")
  }
  
  ### linear mixed model
  # regr_mixed <- lmerTest::lmer(dataset[[dist]] ~ dataset[["Fraction_veg"]] + (1|dataset[["Sp_SiteDate_Phase"]]), data = dataset)
  dataset2<-dataset
  dataset2$dist<-dataset[[dist]]
  regr_mixed <- lmerTest::lmer(dist ~ Fraction_veg + (1|Sp_SiteDate_Phase), data=dataset2)
  regr_mixed <- lmerTest::lmer(dist ~ Fraction_veg + (1|Species), data=dataset2)
  # regr_mixed <- lmerTest::lmer(dataset[[dist]] ~ dataset$Fraction_veg + (1|dataset$Sp_SiteDate_Phase))
  # regr_mixed <- lmerTest::lmer(datasetdist ~ datasetFracveg + (1|dataset[["Sp_SiteDate_Phase"]]))
  p_slope_regr_mixed <-summary(regr_mixed)$coefficients[2,5]
  R2m<-r.squaredGLMM(regr_mixed)[,"R2m"]  # Evaluate model goodness-of-fit
  R2c<-r.squaredGLMM(regr_mixed)[,"R2c"]
  # AIC(regr_mixed)
  # anova(regr_mixed) #same as: anova(regr_mixed,type="marginal") and drop1(regr_mixed, test = "Chisq")
  # Anova(regr_mixed,test.statistic=c("Chisq"))
  # Anova(regr_mixed,test.statistic=c("F"),type="3")
  # Anova(regr_mixed,test.statistic=c("F"),type="2") # same as: Anova(regr_mixed,test.statistic=c("F")) 

  ymin <- min(dataset[[dist]],dataset_twin[[dist]])
  ymax <- max(dataset[[dist]],dataset_twin[[dist]])

  # palette_color <- distinctColorPalette(16) # find a good color palette. E.g.:
  palette_color<-c("#CFE1D4","#72BCD0","#7C91D1","#85DD92","#DC5E7E","#B236E8","#7961DC","#6EE253",
                   "#D5DD95","#D18AD3","#D3E152","#E09B49","#78E3D3","#D959CC","#DBB9D9","#CDA08A")
  palette_shape <- c(10,16)
  
  # x11(width=720,height=720)
  g <- ggplot(dataset,aes(dataset[["Fraction_veg"]] ,dataset[[dist]]))+
        # geom_point(aes(color=Sp_Phase),size=3)+
        geom_point(aes(colour=Species, shape=DvpPhase), size = 3)+
        scale_colour_manual(name='Species',values=palette_color) +
        scale_shape_manual(name='Phenology',values=palette_shape)+
        xlab("Vegetation fraction in FOV")+
        ylab(title)+
        annotate("text", x = max(dataset[["Fraction_veg"]]), y = Inf ,hjust=1, vjust = 2,size=7, #dataset$Fraction_veg
                 label = paste("R²",round(R2,digits=2),sep=" = "), parse=F)+
        annotate("text", x = max(dataset[["Fraction_veg"]]),y= Inf,hjust=1, vjust=5,size=7,
                 label = paste("p",p_slope_regr), parse=F)+
        # annotate("text", x = min(dataset[["Fraction_veg"]]), y = Inf ,hjust=0, vjust = 2,size=7, #dataset$Fraction_veg
        #          label = paste("R²m",round(R2m,digits=2),sep=" = "), parse=F)+
        # annotate("text", x = min(dataset[["Fraction_veg"]]), y = Inf ,hjust=0, vjust = 5,size=7, #dataset$Fraction_veg
        #          label = paste("R²c",round(R2c,digits=2),sep=" = "), parse=F)+
        # annotate("text", x = min(dataset[["Fraction_veg"]]),y= Inf,hjust=0, vjust=8,size=7,
        #          label = paste("p slope = ",round(p_slope_regr_mixed,digits=3)), parse=F)+
        # annotate("text", x = 0.3, y = 0.15,hjust=0,
        # label = paste("nRMSE",round(regr_sam_sum$sigma/mean(meta_T_unmixed.intersect$sam),digits=3),sep=" = "), parse=F)+
        # ! sigma represents pearson RMSE (restricted max likelihood estimator of residual variance), >< statistical RMSE, which is biased downward
        # https://stackoverflow.com/questions/43123462/how-to-obtain-rmse-out-of-lm-result
        # geom_abline(slope=coef(regr_sam)[2],intercept=coef(regr_sam)[1])+
        scale_y_continuous(limits=c(ymin,ymax))+
        theme_bw()+
        # theme(axis.text=element_text(size=18),axis.title=element_text(size=20),
        #       legend.text=element_text(size=20),legend.title=element_text(size=18))
        theme(axis.text=element_text(size=16),axis.title=element_text(size=18))+
        theme(legend.position="none")
  if (p_slope_regr < 0.05) {
    g <- g +geom_smooth(method=lm,se=T,color="black") # se = 0.95 confidence interval
  }
  print(g) 
  setwd(wd2)
  ggsave(paste(paste("Intraspecific_","_",sep=dist_proc),".jpeg",sep=jumpcorr),width = 15,height = 15,units=("cm"),dpi=200)
}

# Function to plot legend:
Graph_FractionDis_legend <- function (wd2,dataset,datasetFracveg,datasetdist){
  regr<-lm(datasetdist~datasetFracveg)  # ,data=dataset
  regr_sum<-summary(regr) # same as: summary(gvlma::gvlma.lm(regr,alphalevel=0.05))
  
  palette_color<-c("#CFE1D4","#72BCD0","#7C91D1","#85DD92","#DC5E7E","#B236E8","#7961DC","#6EE253",
                   "#D5DD95","#D18AD3","#D3E152","#E09B49","#78E3D3","#D959CC","#DBB9D9","#CDA08A")
  palette_shape <- c(10,16)
  
    g<-ggplot(dataset,aes(datasetFracveg,datasetdist))+
    # g<-ggplot(meta_T_unmixed.intersect,aes(Fraction_veg,euc))+
    # geom_point(aes(color=Sp_Phase),size=3)+
      geom_point(aes(colour=Species, shape=DvpPhase), size = 3)+
      scale_colour_manual(name='Species',values=palette_color) +
      scale_shape_manual(name='Phenology',values=palette_shape)+
      
    xlab("Vegetation fracion in FOV")+
    ylab("Euclidean distance")+
    annotate("text", x = min(datasetFracveg), y = Inf ,hjust=0, vjust = 2,size=7,
             label = paste("R²",round(regr_sum$r.squared,digits=2),sep=" = "), parse=F)+
    annotate("text", x = min(datasetFracveg),y= Inf,hjust=0, vjust=5,size=7,
             label = paste("p slope = ",round(regr_sum$coefficients[2,4],digits=2)), parse=F)+
    geom_smooth(method=lm,se=T,color="black")+
    theme_bw()+
    theme(axis.text=element_text(size=12),axis.title=element_text(size=16),
          legend.text=element_text(size=18),legend.title=element_text(size=20))+
    guides(col=guide_legend(ncol=2)) +
    theme(legend.key.size=unit(2,'lines'))
  legend<-get_legend(g)
  ggdraw(plot_grid(legend,ncol=1))
  setwd(wd2)
  ggsave("Intraspecific_legend.jpeg",width = 20,height = 15,units=("cm"),dpi=600)
}

# Apply intraDissim_proc on raw and unmixed spectra and save graphs
intra_dissim <- function (wd1, wd2, jumpcorr,legend){
  setwd(wd1)
  datasets<-readRDS(paste("datasets_",".rds",sep=jumpcorr))

  # Create output vectors and matrices
  spectra_T_raw.intersect<-spectra(datasets$libr_T_raw.intersect)
  spectra_T_unmixed.intersect<-spectra(datasets$libr_T_unmixed.intersect)
  meta_T_raw.intersect<-SI(datasets$libr_T_raw.intersect)
  meta_T_unmixed.intersect<-SI(datasets$libr_T_unmixed.intersect)
  meta_T_raw.intersect$Fraction_veg<-meta_T_unmixed.intersect$Fraction_veg # Add the variable 'vegetation fraction' also to to the raw table dataset

  # Calculate (dis)similarities
  meta_T_raw.intersect <- intraDissim_proc(datasets=datasets,
                                           spectra_T_proc.intersect=spectra_T_raw.intersect,
                                           meta_T_proc.intersect=meta_T_raw.intersect,
                                           proc="raw")
  meta_T_unmixed.intersect <- intraDissim_proc(datasets=datasets,
                                               spectra_T_proc.intersect=spectra_T_unmixed.intersect,
                                               meta_T_proc.intersect=meta_T_unmixed.intersect,
                                               proc="unmixed")
  
  # Visualise
  Graph_FractionDis(wd2,jumpcorr,meta_T_unmixed.intersect,meta_T_raw.intersect,"sam_unmixed","sam","SAM between field patch and \nunmixed table bouquet measurements")
  Graph_FractionDis(wd2,jumpcorr,meta_T_unmixed.intersect,meta_T_raw.intersect,"sid_unmixed","sid","SID between field patch and \nunmixed table bouquet measurements")
  Graph_FractionDis(wd2,jumpcorr,meta_T_unmixed.intersect,meta_T_raw.intersect,"euc_unmixed","euc","Euclidean distance between field patch \nand unmixed table bouquet measurements")
  Graph_FractionDis(wd2,jumpcorr,meta_T_unmixed.intersect,meta_T_raw.intersect,"manabs_unmixed","manabs","Manhattan distance between field patch \nand unmixed table bouquet measurements")
  Graph_FractionDis(wd2,jumpcorr,meta_T_unmixed.intersect,meta_T_raw.intersect,"manrel_unmixed","manrel","Relative Manhattan distance between field patch \nand unmixed table bouquet measurements")
  Graph_FractionDis(wd2,jumpcorr,meta_T_unmixed.intersect,meta_T_raw.intersect,"bcdabs_unmixed","bcdabs","Bray-Curtis distance between field patch \nand unmixed table bouquet measurements")
  Graph_FractionDis(wd2,jumpcorr,meta_T_unmixed.intersect,meta_T_raw.intersect,"bcdrel_unmixed","bcdrel","Relative Bray-Curtis distance between field patch \nand unmixed table bouquet measurements")
  Graph_FractionDis(wd2,jumpcorr,meta_T_unmixed.intersect,meta_T_raw.intersect,"scm_unmixed","scm","Spectral correlation measure between field patch \nand unmixed table bouquet measurements")
  
  Graph_FractionDis(wd2,jumpcorr,meta_T_raw.intersect,meta_T_unmixed.intersect,"sam_raw","sam","SAM between field patch and \noriginal table bouquet measurements")
  Graph_FractionDis(wd2,jumpcorr,meta_T_raw.intersect,meta_T_unmixed.intersect,"sid_raw","sid","SID between field patch and \noriginal table bouquet measurements")
  Graph_FractionDis(wd2,jumpcorr,meta_T_raw.intersect,meta_T_unmixed.intersect,"euc_raw","euc","Euclidean distance between field patch \nand original table bouquet measurements")
  Graph_FractionDis(wd2,jumpcorr,meta_T_raw.intersect,meta_T_unmixed.intersect,"manabs_raw","manabs","Manhattan distance between field patch \nand original table bouquet measurements")
  Graph_FractionDis(wd2,jumpcorr,meta_T_raw.intersect,meta_T_unmixed.intersect,"manrel_raw","manrel","Relative Manhattan distance between field patch \nand original table bouquet measurements")
  Graph_FractionDis(wd2,jumpcorr,meta_T_raw.intersect,meta_T_unmixed.intersect,"bcdabs_raw","bcdabs","Bray-Curtis distance between field patch \nand original table bouquet measurements")
  Graph_FractionDis(wd2,jumpcorr,meta_T_raw.intersect,meta_T_unmixed.intersect,"bcdrel_raw","bcdrel","Relatice Bray-Curtis distance between field patch \nand original table bouquet measurements")
  Graph_FractionDis(wd2,jumpcorr,meta_T_raw.intersect,meta_T_unmixed.intersect,"scm_raw","scm","Spectral correlation measure between field patch \nand original table bouquet measurements")
  
  if (legend == 1){
    Graph_FractionDis_legend(wd2,meta_T_unmixed.intersect,meta_T_unmixed.intersect$Fraction_veg,meta_T_unmixed.intersect$euc)
  }
  
  metalist <- list ("meta_T_raw.intersect" = meta_T_raw.intersect,
                    "meta_T_unmixed.intersect" = meta_T_unmixed.intersect,
                    "jumpcorr"=jumpcorr)
  
  return(metalist)
} # Also add argument for unmixing choice -> focal length ???

#--------------
# Pearson r
Pearson_FractionDis <- function(datasetdist_raw, datasetdist_unmixed){
  Dissim_names<-c("sam","scm","sid","euc","manabs","manrel","bcdabs","bcdrel")
  Pearson_Frac <- data.frame(matrix(ncol = 5, nrow = length(Dissim_names)))
  colnames(Pearson_Frac) <- c("Dissim","Frac-Raw-r","Frac-Raw-sign","Frac-Unmixed-r","Frac-Unmixed-sign")
  Pearson_Frac[,1] <- Dissim_names
  
  for (i in 1:length(Dissim_names)) {
    Pearson_Frac[i,"Frac-Raw-r"]<-cor.test(datasetdist_raw[[Dissim_names[i]]],datasetdist_raw[["Fraction_veg"]],method='pearson')$estimate
    Pearson_Frac[i,"Frac-Raw-sign"]<-cor.test(datasetdist_raw[[Dissim_names[i]]],datasetdist_raw[["Fraction_veg"]],method='pearson')$p.value
    Pearson_Frac[i,"Frac-Unmixed-r"]<-cor.test(datasetdist_unmixed[[Dissim_names[i]]],datasetdist_unmixed[["Fraction_veg"]],method='pearson')$estimate
    Pearson_Frac[i,"Frac-Unmixed-sign"]<-cor.test(datasetdist_unmixed[[Dissim_names[i]]],datasetdist_unmixed[["Fraction_veg"]],method='pearson')$p.value
  }
  
  return(Pearson_Frac)
}
Pearson_FractionDis_apply <- function(wd1,jumpcorr) {
  setwd(wd1)
  datasets<-readRDS(paste("datasets_",".rds",sep=jumpcorr))
  
  # Create output vectors and matrices
  spectra_T_raw.intersect<-spectra(datasets$libr_T_raw.intersect)
  spectra_T_unmixed.intersect<-spectra(datasets$libr_T_unmixed.intersect)
  meta_T_raw.intersect<-SI(datasets$libr_T_raw.intersect)
  meta_T_unmixed.intersect<-SI(datasets$libr_T_unmixed.intersect)
  meta_T_raw.intersect$Fraction_veg<-meta_T_unmixed.intersect$Fraction_veg # Add the variable 'vegetation fraction' also to to the raw table dataset
  
  # Calculate (dis)similarities
  meta_T_raw.intersect <- intraDissim_proc(datasets,spectra_T_raw.intersect,meta_T_raw.intersect,"raw")
  meta_T_unmixed.intersect <- intraDissim_proc(datasets,spectra_T_unmixed.intersect,meta_T_unmixed.intersect,"unmixed")
  
  output <- Pearson_FractionDis(meta_T_raw.intersect, meta_T_unmixed.intersect)
}

# Spearman r
Spearman_FractionDis <- function(datasetdist_raw, datasetdist_unmixed){
  Dissim_names<-c("sam","scm","sid","euc","manabs","manrel","bcdabs","bcdrel")
  Spearman_Frac <- data.frame(matrix(ncol = 5, nrow = length(Dissim_names)))
  colnames(Spearman_Frac) <- c("Dissim","Frac-Raw-r","Frac-Raw-sign","Frac-Unmixed-r","Frac-Unmixed-sign")
  Spearman_Frac[,1] <- Dissim_names
  
  for (i in 1:length(Dissim_names)) {
    Spearman_Frac[i,"Frac-Raw-r"]<-cor.test(datasetdist_raw[[Dissim_names[i]]],datasetdist_raw[["Fraction_veg"]],method='spearman')$estimate
    Spearman_Frac[i,"Frac-Raw-sign"]<-cor.test(datasetdist_raw[[Dissim_names[i]]],datasetdist_raw[["Fraction_veg"]],method='spearman')$p.value
    Spearman_Frac[i,"Frac-Unmixed-r"]<-cor.test(datasetdist_unmixed[[Dissim_names[i]]],datasetdist_unmixed[["Fraction_veg"]],method='spearman')$estimate
    Spearman_Frac[i,"Frac-Unmixed-sign"]<-cor.test(datasetdist_unmixed[[Dissim_names[i]]],datasetdist_unmixed[["Fraction_veg"]],method='spearman')$p.value
  }
  
  return(Spearman_Frac)
}
Spearman_FractionDis_apply <- function(wd1,jumpcorr) {
  setwd(wd1)
  datasets<-readRDS(paste("datasets_",".rds",sep=jumpcorr))
  
  # Create output vectors and matrices
  spectra_T_raw.intersect<-spectra(datasets$libr_T_raw.intersect)
  spectra_T_unmixed.intersect<-spectra(datasets$libr_T_unmixed.intersect)
  meta_T_raw.intersect<-SI(datasets$libr_T_raw.intersect)
  meta_T_unmixed.intersect<-SI(datasets$libr_T_unmixed.intersect)
  meta_T_raw.intersect$Fraction_veg<-meta_T_unmixed.intersect$Fraction_veg # Add the variable 'vegetation fraction' also to to the raw table dataset
  
  # Calculate (dis)similarities
  meta_T_raw.intersect <- intraDissim_proc(datasets,spectra_T_raw.intersect,meta_T_raw.intersect,"raw")
  meta_T_unmixed.intersect <- intraDissim_proc(datasets,spectra_T_unmixed.intersect,meta_T_unmixed.intersect,"unmixed")
  
  output <- Spearman_FractionDis(meta_T_raw.intersect, meta_T_unmixed.intersect)
}

# Kendall r
Kendall_FractionDis <- function(datasetdist_raw, datasetdist_unmixed){
  Dissim_names<-c("sam","scm","sid","euc","manabs","manrel","bcdabs","bcdrel")
  Kendall_Frac <- data.frame(matrix(ncol = 5, nrow = length(Dissim_names)))
  colnames(Kendall_Frac) <- c("Dissim","Frac-Raw-r","Frac-Raw-sign","Frac-Unmixed-r","Frac-Unmixed-sign")
  Kendall_Frac[,1] <- Dissim_names
  
  for (i in 1:length(Dissim_names)) {
    Kendall_Frac[i,"Frac-Raw-r"]<-cor.test(datasetdist_raw[[Dissim_names[i]]],datasetdist_raw[["Fraction_veg"]],method='kendall')$estimate
    Kendall_Frac[i,"Frac-Raw-sign"]<-cor.test(datasetdist_raw[[Dissim_names[i]]],datasetdist_raw[["Fraction_veg"]],method='kendall')$p.value
    Kendall_Frac[i,"Frac-Unmixed-r"]<-cor.test(datasetdist_unmixed[[Dissim_names[i]]],datasetdist_unmixed[["Fraction_veg"]],method='kendall')$estimate
    Kendall_Frac[i,"Frac-Unmixed-sign"]<-cor.test(datasetdist_unmixed[[Dissim_names[i]]],datasetdist_unmixed[["Fraction_veg"]],method='kendall')$p.value
  }
  
  return(Kendall_Frac)
}
Kendall_FractionDis_apply <- function(wd1,jumpcorr) {
  setwd(wd1)
  datasets<-readRDS(paste("datasets_",".rds",sep=jumpcorr))
  
  # Create output vectors and matrices
  spectra_T_raw.intersect<-spectra(datasets$libr_T_raw.intersect)
  spectra_T_unmixed.intersect<-spectra(datasets$libr_T_unmixed.intersect)
  meta_T_raw.intersect<-SI(datasets$libr_T_raw.intersect)
  meta_T_unmixed.intersect<-SI(datasets$libr_T_unmixed.intersect)
  meta_T_raw.intersect$Fraction_veg<-meta_T_unmixed.intersect$Fraction_veg # Add the variable 'vegetation fraction' also to to the raw table dataset
  
  # Calculate (dis)similarities
  meta_T_raw.intersect <- intraDissim_proc(datasets,spectra_T_raw.intersect,meta_T_raw.intersect,"raw")
  meta_T_unmixed.intersect <- intraDissim_proc(datasets,spectra_T_unmixed.intersect,meta_T_unmixed.intersect,"unmixed")
  
  output <- Kendall_FractionDis(meta_T_raw.intersect, meta_T_unmixed.intersect)
}

#--------------
# Functions to create histograms and apply them to every (dis)similarity metric and species (cover fractions)
Hist_Dissim <- function(wd3,dataset,datasetdist,dist_proc,jumpcorr) {
  setwd(wd3)
  ggplot(dataset, aes(x=datasetdist))+
    geom_histogram(binwidth = 0.008,color="black", fill="gray")+
    # geom_histogram(aes(y=..density..), colour="black", fill="gray")+
    # geom_density(alpha=.2, fill="#FF6666") +
    xlim(0,0.17)+
    ylim(0,10)+
    # ggtitle("SAM between umixed table spectrum and field spectrum")+
    theme_classic()
  ggsave(paste(paste("Hist_","_",sep=dist_proc),".jpeg",sep=jumpcorr),width = 15,height = 10,units=("cm"),dpi=600)
}
Hist_Dissim_apply <- function (wd3,metalist) {
  meta_unmixed <- metalist$meta_T_unmixed.intersect
  meta_raw <- metalist$meta_T_raw.intersect
  jumpcorr <- metalist$jumpcorr
  
  Hist_Dissim(wd3,meta_unmixed, meta_unmixed$sam, "sam_unmixed",jumpcorr)
  Hist_Dissim(wd3,meta_unmixed, meta_unmixed$sam, "sid_unmixed",jumpcorr)
  Hist_Dissim(wd3,meta_unmixed, meta_unmixed$sam, "euc_unmixed",jumpcorr)
  Hist_Dissim(wd3,meta_unmixed, meta_unmixed$sam, "manabs_unmixed",jumpcorr)
  Hist_Dissim(wd3,meta_unmixed, meta_unmixed$sam, "manrel_unmixed",jumpcorr)
  Hist_Dissim(wd3,meta_unmixed, meta_unmixed$sam, "sorabs_unmixed",jumpcorr)
  Hist_Dissim(wd3,meta_unmixed, meta_unmixed$sam, "scm_unmixed",jumpcorr)
  
  Hist_Dissim(wd3,meta_raw, meta_raw$sam, "sam_raw",jumpcorr)
  Hist_Dissim(wd3,meta_raw, meta_raw$sam, "sid_raw",jumpcorr)
  Hist_Dissim(wd3,meta_raw, meta_raw$sam, "euc_raw",jumpcorr)
  Hist_Dissim(wd3,meta_raw, meta_raw$sam, "manabs_raw",jumpcorr)
  Hist_Dissim(wd3,meta_raw, meta_raw$sam, "manrel_raw",jumpcorr)
  Hist_Dissim(wd3,meta_raw, meta_raw$sam, "sorabs_raw",jumpcorr)
  Hist_Dissim(wd3,meta_raw, meta_raw$sam, "scm_raw",jumpcorr)
}
Hist_FractionVeg <- function(wd4,dataset,jumpcorr,Species_Dvp) {
  Hist_SpDvp =dataset[dataset$Species_Cond == Species_Dvp ,]
  
  setwd(wd4)
  ggplot(Hist_SpDvp, aes(x=Fraction_veg))+
    geom_histogram(binwidth = 0.05,color="black", fill="gray")+
    xlim(0,1.05)+
    xlab("Vegation fraction on the black table")+
    ylab("Number of bouquets")+
    theme_classic()
  ggsave(paste("Hist_FractionVeg_",".jpeg",sep=Species_Dvp),width = 15,height = 10,units=("cm"),dpi=600)
  
  ggsave(paste(paste("Hist_FractionVeg_","_",sep=Species_Dvp),".jpeg",sep=jumpcorr),width = 15,height = 10,units=("cm"),dpi=600)
  
} # might need input for every argument
Hist_FractionVeg_apply <- function(wd4,metalist) {
  meta_unmixed <- metalist$meta_T_unmixed.intersect
  jumpcorr <- metalist$jumpcorr
  
  levels<-levels(droplevels(meta_unmixed$Species_Cond))
  for (i in 1:length(levels)) {
    Hist_FractionVeg(wd4,meta_unmixed,jumpcorr,levels[i])
  }
}

#-------------------------------------------------------------------------------- 
# APPLY all of the previous functions on the spectral data, with jump correction method  
wd1 = "C:/Users/u0091812/Box Sync/03. Black table procedure/Data/Spectral data"
wd2 = "C:/Users/u0091812/Box Sync/03. Black table procedure/Figures/Intraspecific/Regression cover fraction vs similarity"
wd3 = "C:/Users/u0091812/Box Sync/03. Black table procedure/Figures/Intraspecific/Histograms of similarity measures"
wd4 = "C:/Users/u0091812/Box Sync/Black table Letter/Figures/Intraspecific/Histograms of vegetation fractions"

intra_dissim_Nocorr<-intra_dissim(wd1, wd2, "Nocorr",legend="1")
Hist_Dissim_apply(wd3,intra_dissim_Nocorr)
Hist_FractionVeg_apply(wd4,intra_dissim_Nocorr)

Pearson_Nocorr<-Pearson_FractionDis_apply(wd1,"Nocorr")
# https://stats.stackexchange.com/questions/156995/correlation-between-a-percentage-and-a-real-number
Spearman_Nocorr<-Spearman_FractionDis_apply(wd1,"Nocorr")
Kendall_Nocorr<-Kendall_FractionDis_apply(wd1,"Nocorr")

# intra_dissim_jcorr_add_SWIR1<-intra_dissim(wd1, wd2, "jcorr_add_SWIR1",legend="0")
# Hist_Dissim_apply(wd3,intra_dissim_jcorr_add_SWIR1)
# Hist_FractionVeg_apply(wd4,intra_dissim_jcorr_add_SWIR1)

# intra_dissim_jcorr_mult_SWIR1<-intra_dissim(wd1, wd2, "jcorr_mult_SWIR1",legend="0")
# Hist_Dissim_apply(wd3,intra_dissim_jcorr_mult_SWIR1)
# Hist_FractionVeg_apply(wd4,intra_dissim_jcorr_mult_SWIR1)

# intra_dissim_jcorr_mult_SWIR12<-intra_dissim(wd1, wd2, "jcorr_mult_SWIR12",legend="0")
# Hist_Dissim_apply(wd3,intra_dissim_jcorr_mult_SWIR12)
# Hist_FractionVeg_apply(wd4,intra_dissim_jcorr_mult_SWIR12)

# intra_dissim_jcorr_mult_VNIR<-intra_dissim(wd1, wd2, "jcorr_mult_VNIR",legend="0")
# Hist_Dissim_apply(wd3,intra_dissim_jcorr_mult_VNIR)
# Hist_FractionVeg_apply(wd4,intra_dissim_jcorr_mult_VNIR)

# intra_dissim_Nocorr<-intra_dissim(wd1, wd2, "Nocorr_frac0.9",legend="0")
# Hist_Dissim_apply(wd3,intra_dissim_Nocorr)
# Hist_FractionVeg_apply(wd4,intra_dissim_Nocorr)
