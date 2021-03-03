###################################################################################
###                                                                             ###
###        03b. Interspecific similarities based on (dis)similarity indices     ###
###                                                                             ###
###################################################################################


### The objective of this script is to QUANTITATIVELY evaluate the performance of the black table
# BETWEEN SPECIES (interspecific) we investigate
#   1) the magnitude of interspecific spectral differences calculated on field vs. table spectra: matrix of spectral similarity indices
#   2) if interspecific spectral differences are preserved on the table: mantel test between the similarity matrices of field and table spectra

### The input of these steps are the libraries created in "01. Black table data":
#   For assessing (dis)similarities BETWEEN species we need:
#       - AVERAGE FIELD and TABLE spectrum per SPECIES AND CONDITION (Patch - Dvpmt phase):
#         libr_F.intersect.mean.overall
#         libr_T_raw.intersect.mean.overall
#         libr_T_unmixed.intersect.mean.overall

# ! These signatures have been preprocessed with following specifications:
#   - smoothing: savitsky-golay filter with n = 51 and 101 for asd and sv respectively
#   - removal of atmospheric noise windows:  c(349,400,1340,1460,1780,1970,2400,2501)

# Script by Elisa Van Cleemput, 2018
#--------------------------------------------------------------------------------
# Clean workspace

rm(list=ls())

#--------------------------------------------------------------------------------

######## BETWEEN species: overall mean of field + overall mean of table spectra
# calculate indices for both raw and unmixed spectra
# Imitate Figure 6. from Jiménez & Díaz-Delgado (2015) + Figure 1 and SI-3 from Schweiger et al. (2018)

library(hsdar) # dist.speclib, sam_distance
library(resemble) # fDiss, sid
library(vegan) # mantel
library(reshape2) # 'melt' function
library(ggplot2) #create heatmap
library(lemon) # facet_rep_wrap

#--------------------------------------------------------------------------------
# FUNCTIONS

# (dis)similarity matrix calculation
interDissim_calc <- function(em,libr){
  sam<-fDiss(em,method="cosine",center=F,scaled=F)
  sam2<-dist.speclib(libr,method="sam") # the same values as sam, but does not include the 1:1 diagonal in the results 
  sam3<-sam_distance(libr) # the same values as sam:-)
  sid<-sid(em,mode="density",center=F,scaled=F)$sid
  scm<- corDiss(em,center=F,scaled=F)
  euc<-as.matrix(dist(em,method="euclidean",diag=T,upper=T))
  manabs<-as.matrix(dist(em,method="manhattan",diag=T,upper=T))
  bcdabs<-manabs
  for (i in 1:nrow(bcdabs)) {
    for (j in 1:nrow(bcdabs)){
      bcdabs[i,j] <- dist(rbind(em[i,],em[j,]),method="manhattan")/sum(em[i,],em[j,])
    }
    
  }
  InterDissim <- list("SAM" = sam,
                      # "SAM2" = sam2,
                      "SID" = sid,
                      "SCM" = scm,
                      "EUC" = euc,
                      "MAN_ABS" = manabs,
                      "BCD_ABS" = bcdabs)

  labels<-SI(libr)$Sp_Phase
  InterDissim<-lapply(InterDissim, function(x) {rownames(x) <- labels; x})
  InterDissim<-lapply(InterDissim, function(x) {colnames(x) <- labels; x})
  
  Dissim_names<-c("SAM","SID","SCM","EUC","MAN_ABS","BCD_ABS")

  output<- list("InterDissim" = InterDissim,
                "Dissim_names" = Dissim_names)
  return(output)
}

# Mantel test between field and table (dis)similarity matrix --> also applicable on similarity matrices?
mantel_FieldTable <-function(interDissim_field,interDissim_table_raw,interDissim_table_unmixed,Dissim_names) {
  mantel_Dissim <- data.frame(matrix(ncol = 5, nrow = length(interDissim_field)))
  colnames(mantel_Dissim) <- c("Dissim","Field-Raw-r","Field-Raw-sign","Field-Unmixed-r","Field-Unmixed-sign")
  mantel_Dissim[,1] <- Dissim_names
  
  for (i in 1:length(Dissim_names)) {
    mantel_Dissim[i,"Field-Raw-r"] <- mantel(interDissim_field[[Dissim_names[i]]],interDissim_table_raw[[Dissim_names[i]]])$statistic
    mantel_Dissim[i,"Field-Raw-sign"] <- mantel(interDissim_field[[Dissim_names[i]]],interDissim_table_raw[[Dissim_names[i]]])$signif
    mantel_Dissim[i,"Field-Unmixed-r"] <- mantel(interDissim_field[[Dissim_names[i]]],interDissim_table_unmixed[[Dissim_names[i]]])$statistic
    mantel_Dissim[i,"Field-Unmixed-sign"] <- mantel(interDissim_field[[Dissim_names[i]]],interDissim_table_unmixed[[Dissim_names[i]]])$signif
  }
  # sam vs. sam2: r is the same, but significance and quantiles differ. Which input is prefered?: dist objects or matrices?
  
  return(mantel_Dissim)
}

### Visualise (dis)similarity matrices: similar to Figure 6. from Jiménez & Díaz-Delgado (2015)
#  http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
get_lower_tri<-function(matrix){
  matrix[upper.tri(matrix)] <- NA
  return(matrix)
} # Get lower triangle of the matrix
get_upper_tri <- function(matrix){
  matrix[lower.tri(matrix)]<- NA
  # rownames(matrix)<-axislabels
  # colnames(matrix)<-axislabels
  return(matrix)
} # Get upper triangle of the matrix
DisMatrix_vis <- function(wd3,jumpcorr,em_field,dist,distproc,MatrixField,MatrixTable){
  mincolor<-min(min(MatrixField),min(MatrixField))
  maxcolor<-max(max(MatrixField),max(MatrixField))
  midcolor<-(mincolor+maxcolor)/2
  
  axislabels<-rownames(em_field)
  
  theme<- theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    axis.line=element_blank(),
    legend.title=element_blank())

  # x11(width=720,height=720)
  m<-ggplot(melt(get_upper_tri(MatrixField), na.rm=T)) +
    geom_tile(aes(x=Var1, y=Var2, fill=value),color="white") +
    # geom_text(aes(Var1, Var2, label = round(value,digits=1)), color = "black", size = 4) +
    geom_tile(data=melt(get_lower_tri(MatrixTable), na.rm=T), aes(x=Var1, y=Var2, fill=value),color="white") +
    geom_abline(intercept = 0, slope = 1, color="white",size=1) +
    # geom_text(data=melt(get_lower_tri(sam_table_raw), na.rm=T),aes(X1, X2, label = round(value,digits=1)), color = "black", size = 4) +
    # scale_fill_gradient2(low = "black", high = "lightgrey", mid = "grey",
    #                      midpoint = 0.17, limit = c(0,0.35), space = "Lab",
    #                      name=dist) +
    scale_fill_gradient2(low = "black", high = "lightgrey", mid = "grey", 
                         midpoint = midcolor,limit = c(mincolor,maxcolor),
                         space = "Lab") +
    guides(fill = guide_colorbar(barwidth = 2, barheight = 10))+
    theme(axis.text.x = element_text(angle = 55, vjust = 1,size = 16, hjust = 1),
          axis.text.y = element_text(size=16),
          legend.text=element_text(size=16),
          legend.title=element_text(size=20),
          legend.spacing.x = unit(0.3, 'cm'))+
    scale_x_discrete(labels=axislabels)+
    scale_y_discrete(labels=axislabels)+
    coord_fixed() +
    # labs(title="field (upper triangle) and raw table SAM (lower triangle)",
    # subtitle=paste("Mantel r",round(m_FRSAM$statistic,digits=3),sep=" = "))+
    # ggtitle("field (upper triangle) and raw table SAM (lower triangle)") +
    # theme(plot.title = element_text(size=10))+
    theme
  # m
  setwd(wd3)
  ggsave(paste(paste("Matrix_","_",sep=distproc),".jpeg",sep=jumpcorr),width = 20,height = 20,units=("cm"),dpi=600)
  
  # upper triangle = field data + lower triangle = raw table data
}

### Visualise as regression (similar to Fig. 1 from Schweiger et al. 2018): pairwise similarities 
DisRegr_vis <- function(wd4,jumpcorr,em_field,dist,MatrixField,MatrixTable_unmixed,MatrixTable_raw,title, intercept) {
  axislabels<-rownames(em_field)

  rownames(MatrixField)<-axislabels
  colnames(MatrixField)<-axislabels
  data_full<-cbind(melt(MatrixField),
                   data.frame(melt(MatrixTable_unmixed)$value),
                   data.frame(melt(MatrixTable_raw)$value))
  colnames(data_full)<-c("sp1","sp2","field","unmixed","raw")

    # Remove zeros (= values on diagonal) = dissimilarity of species with itself
  data_full<-subset(data_full,data_full[,'field']>0.00005)

  # Linear regression 
  if (intercept == "_with_intercept") {
    regr_unmixed<-lm(unmixed ~ field, data = data_full) # with intercept
    p_slope_regr_unmixed<-summary(regr_unmixed)$coefficients[2,4] # same as: summary(gvlma::gvlma.lm(regr_unmixed,alphalevel=0.05))
    t_val_unmixed=coef(summary(regr_unmixed))[2,3]
    # str(gvlma::gvlma.lm(regr_unmixed,alphalevel=0.05)$GlobalTest)
    # GlobalStat4:       global test about appropriateness of linear model assumptions
    # DirectionalStat1:  detect Skewness = measure of the lack of symmetry
    # DirectionalStat2:  detect Kurtosis = measure of whether the data are heavy-tailed or light-tailed relative to normal distribution
    # DirectionalStat3:  detect nonlinear Link function
    # DirectionalStat4:  detect Heteroscedasticity
    #  Pena, EA and Slate, EH (2006). "Global validation of linear model assumptions," J.\ Amer.\ Statist.\ Assoc., 101(473):341-354.
  } else if (intercept == "_without_intercept"){
    regr_unmixed<-lm(unmixed ~ field -1, data = data_full) # without intercept
    p_slope_regr_unmixed<-summary(regr_unmixed)$coefficients[1,4] # same as: summary(gvlma::gvlma.lm(regr_unmixed,alphalevel=0.05))
    t_val_unmixed=coef(summary(regr_unmixed))[1,3]
  }
  setwd(wd4)
    jpeg(paste(paste(paste("Interspecific_","_field_unmixed_",sep=dist),intercept,sep=jumpcorr),"_LinAssum.jpeg",sep=""),res=300,width=12,height=12,units='in')
    par(mfrow=c(2,2))
    plot(regr_unmixed) # same as: plot(gvlma::gvlma.lm(regr,alphalevel=0.05))
    dev.off()
    # hist(regr$residuals)
    R2_unmixed <- summary(regr_unmixed)$r.squared
    RMSE_unmixed <- sqrt(mean(regr_unmixed$residuals^2))
    
    if (intercept == "_with_intercept") {
      regr_raw<-lm(raw ~ field, data = data_full)
      p_slope_regr_raw<-summary(regr_raw)$coefficients[2,4] # same as: summary(gvlma::gvlma.lm(regr_unmixed,alphalevel=0.05))
      t_val_raw <- coef(summary(regr_raw))[2,3]
      # str(gvlma::gvlma.lm(regr_raw,alphalevel=0.05)$GlobalTest)
    } else if (intercept == "_without_intercept"){
      regr_raw<-lm(raw ~ field -1, data = data_full)
      p_slope_regr_raw<-summary(regr_raw)$coefficients[1,4] # same as: summary(gvlma::gvlma.lm(regr_unmixed,alphalevel=0.05))
      t_val_raw <- coef(summary(regr_raw))[1,3]
    }
    
    setwd(wd4)
    jpeg(paste(paste(paste("Interspecific_","_field_raw_",sep=dist),intercept,sep=jumpcorr),"_LinAssum.jpeg",sep=""),res=300,width=12,height=12,units='in')
    par(mfrow=c(2,2))
    plot(regr_raw) # same as: plot(gvlma::gvlma.lm(regr,alphalevel=0.05))
    dev.off()
    R2_raw <- summary(regr_raw)$r.squared
    RMSE_raw <- sqrt(mean(regr_raw$residuals^2))
    
    if (p_slope_regr_unmixed < 0.001) {
      p_slope_regr_unmixed <- "< 0.001"
    } else {
      p_slope_regr_unmixed <- paste("=",round(p_slope_regr_unmixed,digits=3),sep=" ")
    }
    if (p_slope_regr_raw < 0.001) {
      p_slope_regr_raw <- "< 0.001"
    } else {
      p_slope_regr_raw <- paste("=",round(p_slope_regr_raw,digits=3),sep=" ")
    }
    
    stat_raw <- data.frame(R2_raw=R2_raw,
                           p.val_raw=p_slope_regr_raw,
                           t.val_raw=t_val_raw,
                           df.regression_raw=summary(regr_raw)$df[1],
                           df.residual_raw=summary(regr_raw)$df[2])
    stat_unmixed <- data.frame(R2_unmixed=R2_unmixed,
                               p.val_unmixed=p_slope_regr_unmixed,
                               t.val_unmixed=t_val_unmixed,
                               df.regression_unmixed=summary(regr_unmixed)$df[1],
                               df.residual_unmixed=summary(regr_unmixed)$df[2])
    overview <- cbind(stat_raw,stat_unmixed)
    overview_stat <- list ("regr_raw"=regr_raw,
                           "regr_unmixed"=regr_unmixed,
                           overview)
    
  # palette_color <- distinctColorPalette(16) # find a good color palette. E.g.:
  palette_color<-c("#CFE1D4","#72BCD0","#7C91D1","#85DD92","#DC5E7E","#B236E8","#7961DC","#6EE253",
                   "#D5DD95","#D18AD3","#D3E152","#E09B49","#78E3D3","#D959CC","#DBB9D9","#CDA08A")
  palette_shape <- c(10,16)
    
  data_full$sp1_name<-as.factor(substr(data_full$sp1,1,6))
  data_full$sp2_name<-as.factor(substr(data_full$sp2,1,6))
  data_full$sp2_dvp<-as.factor(substr(data_full$sp2,8,13))
  
  axismax<-max(data_full$field,data_full$raw,data_full$unmixed)
  # axismin<-min(data_full$field,data_full$raw,data_full$unmixed)
  axismin<-0
  
  # Unmixed vs. field
 g <- ggplot(data_full,aes(field,unmixed))+
    geom_point()+
    #   geom_point(aes(colour=sp2_name, shape=sp2_dvp), size = 3)+
    #   scale_colour_manual(name='Species',values=palette_color) +
    #   scale_shape_manual(name='Phenology',values=palette_shape)+
    # scale_colour_discrete(name = "Species and \nphenological stage")+
    xlab(paste(title,"between field \n patch measurements"))+
    ylab(paste(title,"between unmixed \n table bouquet meausurements"))+
    # geom_smooth(aes(group=sp2,colour=sp2_name,linetype=sp2_dvp), method="lm",se=F)+
    scale_colour_manual(name='Species and \nphenological stage',values=palette_color) +
    scale_linetype_manual(values=c("twodash","longdash"))+
    scale_x_continuous(limits=c(axismin,axismax))+
    scale_y_continuous(limits=c(axismin,axismax))+
    annotate("text", x = axismin, y = Inf ,hjust=0, vjust = 2,size=7,
             label = paste("R²",round(R2_unmixed,digits=2),sep=" = "), parse=F)+
    annotate("text", x = axismin, y = Inf ,hjust=-1, vjust = 2,size=7,
              label = paste("RMSE",round(RMSE_unmixed,digits=2),sep=" = "), parse=F)+
    annotate("text", x = axismin, y = Inf ,hjust=0, vjust = 5,size=7,
             label = paste("p",p_slope_regr_unmixed,sep=" "), parse=F)+
    # annotate("text", x = axismin, y = Inf ,hjust=0, vjust = 5,size=7,
    #          label = paste("p slope",round(p_slope_regr_unmixed,digits=3),sep=" = "), parse=F)+
    theme(legend.key.size=unit(2,'lines'))+
    theme_bw()+
    theme(plot.margin=unit(c(10,10,6,6),"points"))+
    theme(axis.text=element_text(size=14),axis.title=element_text(size=18))+
    # legend.text=element_text(size=16),legend.title=element_text(size=18))+
    theme(legend.position="none")
  
  if (intercept == "_without_intercept") {
    g <- g + geom_smooth(aes(group=sp2,colour=sp2_name,linetype=sp2_dvp), method="lm",se=F, formula=y~x-1) +
             geom_smooth(method="lm", colour="black", formula=y~x-1)
  } else if (intercept =="_with_intercept"){
    g <- g + geom_smooth(aes(group=sp2,colour=sp2_name,linetype=sp2_dvp), method="lm",se=F) +
               geom_smooth(method="lm", colour="black")
  }
 print(g)
  setwd(wd4)
  ggsave(paste(paste(paste("Regression_overall_","_field_unmixed_",sep=dist),intercept,sep=jumpcorr),".jpeg",sep=""),width = 16,height = 15,units=("cm"),dpi=600)

  # Raw vs. field
  q <- ggplot(data_full,aes(field,raw))+
    geom_point()+
    xlab(paste(title,"between field \n patch measurements"))+
    ylab(paste(title,"between original \n table bouquet meausurements"))+
    # geom_smooth(aes(group=sp2,colour=sp2_name,linetype=sp2_dvp), method="lm",se=F)+
    scale_colour_manual(name = "Species and \nphenological stage",values=palette_color)+
    scale_linetype_manual(values=c("twodash","longdash"))+
    # geom_smooth(method="lm", colour="black")+
    scale_x_continuous(limits=c(axismin,axismax))+
    scale_y_continuous(limits=c(axismin,axismax))+
    annotate("text", x = axismin, y = Inf ,hjust=0, vjust = 2,size=7,
             label = paste("R²",round(R2_raw,digits=2),sep=" = "), parse=F)+
    annotate("text", x = axismin, y = Inf ,hjust=-1, vjust = 2,size=7,
             label = paste("RMSE",round(RMSE_raw,digits=2),sep=" = "), parse=F)+
    annotate("text", x = axismin, y = Inf ,hjust=0, vjust = 5,size=7,
             label = paste("p",p_slope_regr_raw,sep=" "), parse=F)+
    theme(legend.key.size=unit(2,'lines'))+
    theme_bw()+
    theme(plot.margin=unit(c(10,10,6,6),"points"))+
    theme(axis.text=element_text(size=14),axis.title=element_text(size=18))+
    theme(legend.position="none")
  if (intercept == "_without_intercept") {
    q <- q + geom_smooth(aes(group=sp2,colour=sp2_name,linetype=sp2_dvp), method="lm",se=F, formula=y~x-1) +
      geom_smooth(method="lm", colour="black", formula=y~x-1)
  } else if (intercept =="_with_intercept"){
    q <- q + geom_smooth(aes(group=sp2,colour=sp2_name,linetype=sp2_dvp), method="lm",se=F) +
      geom_smooth(method="lm", colour="black")
  }
  print(q)
  setwd(wd4)
  ggsave(paste(paste(paste("Regression_overall_","_field_raw_",sep=dist),intercept,sep=jumpcorr),".jpeg",sep=""),width = 16,height = 15,units=("cm"),dpi=600)
  
  palette_color<-c("#CFE1D4","#72BCD0","#7C91D1","#85DD92","#DC5E7E","#B236E8","#7961DC","#6EE253",
                   "#D5DD95","#D18AD3","#D3E152","#E09B49","#78E3D3","#D959CC","#DBB9D9","#CDA08A",
                   "#FFCC00")
  
  # Save legend separately
  # x11(width=1280,height=720)
  r<- ggplot(data_full,aes(field,unmixed))+
      geom_point()+
      xlab(paste(title,"between field bouquet measurements"))+
      ylab(paste(title,"between table meausurements"))+
      geom_smooth(aes(group=sp2,colour=sp2_name,linetype=sp2_dvp), method="lm",se=F)+
      scale_linetype_manual(name = "Phenology", values=c("twodash","longdash"), 
                            guide=guide_legend(override.aes = list(colour=c("black", "black"))))+
      scale_colour_manual(name = "Species",values=palette_color)+
      geom_smooth(method="lm", colour="black")+
      scale_x_continuous(limits=c(axismin,axismax))+
      scale_y_continuous(limits=c(axismin,axismax))+
      theme_bw()+
      theme(axis.text=element_text(size=16),axis.title=element_text(size=18),
            legend.text=element_text(size=22),legend.title=element_text(size=24))+
      guides(col=guide_legend(ncol=2)) +
      theme(legend.key.size=unit(3,'lines'))
  library(cowplot)
  legend<-get_legend(r)
  ggdraw(plot_grid(legend,ncol=1))
  setwd(wd4)
  ggsave("Regression_legend.jpeg",width = 15,height = 20,units=("cm"),dpi=600)

  return(overview_stat)
}

### Visualise as regression for each species separately (similar to Fig. SI-3 from Schweiger et al. 2018): pairwise similarities 
DisRegr_apart_vis <- function(wd4,jumpcorr,em_field,dist,MatrixField,MatrixTable_unmixed,MatrixTable_raw,title, intercept) {
  axislabels<-rownames(em_field)
  
  rownames(MatrixField)<-axislabels
  colnames(MatrixField)<-axislabels
  data_full<-cbind(melt(MatrixField),
                   data.frame(melt(MatrixTable_unmixed)$value),
                   data.frame(melt(MatrixTable_raw)$value))
  colnames(data_full)<-c("sp1","sp2","field","unmixed","raw")
  data_full$sp2_name<-as.factor(substr(data_full$sp2,1,6))
  data_full$sp2_dvp<-as.factor(substr(data_full$sp2,8,13))
  # Remove zeros (= values on diagonal) = dissimilarity of species with itself
  data_full<-subset(data_full,data_full[,'field']>0.00005)
  
  ## Create a list containing the p-values and R² values for regressions on each combination of species and DvpPhase (sp2)
  # https://stackoverflow.com/questions/30792685/draw-geom-smooth-only-for-fits-that-are-significant
  stat_raw = lapply(levels(data_full$sp2), function(i) {
    if(nrow(data_full[data_full[,"sp2"]==i,]) > 1) {
      if (intercept == "_with_intercept") {
        data.frame(sp2=i, 
                   R2=summary(lm(raw ~ field, data = data_full[data_full[,"sp2"]==i,]))$r.squared,
                   RMSE = sqrt(mean(lm(raw ~ field, data = data_full[data_full[,"sp2"]==i,])$residuals^2)),
                   p.val=coef(summary(lm(raw ~ field, data = data_full[data_full[,"sp2"]==i,])))[2,4],
                   t.val=coef(summary(lm(raw ~ field, data = data_full[data_full[,"sp2"]==i,])))[2,3],
                   df.regression=summary(lm(raw ~ field, data = data_full[data_full[,"sp2"]==i,]))$df[1],
                   df.residual=summary(lm(raw ~ field, data = data_full[data_full[,"sp2"]==i,]))$df[2])
      } else if (intercept == "_without_intercept"){
        data.frame(sp2=i, 
                   R2=summary(lm(raw ~ field -1, data = data_full[data_full[,"sp2"]==i,]))$r.squared,
                   RMSE = sqrt(mean(lm(raw ~ field -1, data = data_full[data_full[,"sp2"]==i,])$residuals^2)),
                   p.val=coef(summary(lm(raw ~ field -1, data = data_full[data_full[,"sp2"]==i,])))[1,4],
                   t.val=coef(summary(lm(raw ~ field -1, data = data_full[data_full[,"sp2"]==i,])))[1,3],
                   df.regression=summary(lm(raw ~ field -1, data = data_full[data_full[,"sp2"]==i,]))$df[1],
                   df.residual=summary(lm(raw ~ field -1, data = data_full[data_full[,"sp2"]==i,]))$df[2])
      }
    }
  })
  stat_raw = do.call(rbind, stat_raw)   # Turn list into a data frame
  stat_raw2 = stat_raw[stat_raw$p.val < 0.05, ]   # Keep only rows with p.val < 0.05
  stat_raw$RMSE_print <- stat_raw$RMSE
  stat_raw[which(stat_raw$RMSE < 0.01),"RMSE_print"] <- paste("RMSE < 0.01")
  stat_raw[which(stat_raw$RMSE > 0.01),"RMSE_print"] <- paste("RMSE =",round(stat_raw$RMSE,digits=2),sep=" ")
  data_full.subset_raw = data_full[data_full[,"sp2"] %in% stat_raw2$sp2, ]
  
  stat_unmixed = lapply(levels(data_full$sp2), function(i) {
    if(nrow(data_full[data_full[,"sp2"]==i,]) > 1) {
      if (intercept == "_with_intercept") {
        data.frame(sp2=i, 
                   R2=summary(lm(unmixed ~ field, data = data_full[data_full[,"sp2"]==i,]))$r.squared,
                   RMSE = sqrt(mean(lm(unmixed ~ field, data = data_full[data_full[,"sp2"]==i,])$residuals^2)),
                   p.val=coef(summary(lm(unmixed ~ field, data = data_full[data_full[,"sp2"]==i,])))[2,4],
                   t.val=coef(summary(lm(unmixed ~ field, data = data_full[data_full[,"sp2"]==i,])))[2,3],
                   df.regression=summary(lm(unmixed ~ field, data = data_full[data_full[,"sp2"]==i,]))$df[1],
                   df.residual=summary(lm(unmixed ~ field, data = data_full[data_full[,"sp2"]==i,]))$df[2])  
      } else if (intercept == "_without_intercept") {
        data.frame(sp2=i, 
                   R2=summary(lm(unmixed ~ field -1, data = data_full[data_full[,"sp2"]==i,]))$r.squared,
                   RMSE = sqrt(mean(lm(unmixed ~ field -1, data = data_full[data_full[,"sp2"]==i,])$residuals^2)),
                   p.val=coef(summary(lm(unmixed ~ field -1, data = data_full[data_full[,"sp2"]==i,])))[1,4],
                   t.val=coef(summary(lm(unmixed ~ field -1, data = data_full[data_full[,"sp2"]==i,])))[1,3],
                   df.regression=summary(lm(unmixed ~ field -1, data = data_full[data_full[,"sp2"]==i,]))$df[1],
                   df.residual=summary(lm(unmixed ~ field -1, data = data_full[data_full[,"sp2"]==i,]))$df[2])
      }
    }
  })
  stat_unmixed = do.call(rbind, stat_unmixed)   # Turn list into a data frame
  stat_unmixed2 = stat_unmixed[stat_unmixed$p.val < 0.05, ]   # Keep only rows with p.val < 0.05
  stat_unmixed$RMSE_print <- stat_unmixed$RMSE
  stat_unmixed[which(stat_unmixed$RMSE < 0.01),"RMSE_print"] <- paste("RMSE < 0.01")
  stat_unmixed[which(stat_unmixed$RMSE > 0.01),"RMSE_print"] <- paste("RMSE =",round(stat_unmixed$RMSE,digits=2),sep=" ")
  data_full.subset_unmixed = data_full[data_full[,"sp2"] %in% stat_unmixed2$sp2, ]
  
  palette_shape <- c(10,16)
  # palette_color <- distinctColorPalette(16) # find a good color palette. E.g.:
  palette_color<-c("#CFE1D4","#72BCD0","#7C91D1","#85DD92","#DC5E7E","#B236E8","#7961DC","#6EE253",
                   "#D5DD95","#D18AD3","#D3E152","#E09B49","#78E3D3","#D959CC","#DBB9D9","#CDA08A")
  palette_color<-cbind(palette_color,levels(data_full$sp2_name))
  colnames(palette_color)<-c("palette_color","sp2_name")
  palette_color_raw<-palette_color[palette_color[,"sp2_name"] %in% data_full.subset_raw$sp2_name,]
  palette_color_unmixed<-palette_color[palette_color[,"sp2_name"] %in% data_full.subset_unmixed$sp2_name,]
  
  axismax<-max(data_full$field,data_full$raw,data_full$unmixed)
  # axismin<-min(data_full$field,data_full$raw,data_full$unmixed)
  axismin<-0
  
  # Raw vs. field
  if (is.matrix(palette_color_raw) == T){
    col<-palette_color_raw[,"palette_color"]
  } else {
    col <- palette_color_raw[["palette_color"]]
  }
  g <- ggplot(data_full,aes(field,raw))+
        geom_point(aes(shape=sp2_dvp))+
        scale_shape_manual(name='Phenology',values=palette_shape)+
        xlab(paste(title,"between field patch measurements"))+
        ylab(paste(title,"between original table bouquet meausurements"))+
        facet_rep_wrap(~sp2,repeat.tick.labels = T)+
        # facet_wrap(~sp2,scales="free")+
        # geom_smooth(data=data_full.subset_raw,aes(color=sp2_name),method="lm",se=T) + 
        # geom_smooth(aes(color=sp2_name),method="lm")+
        scale_colour_manual(name='Species',values=col) +
        scale_x_continuous(limits=c(axismin,axismax))+
        scale_y_continuous(limits=c(axismin-10,axismax+10))+
        coord_cartesian(xlim=c(axismin,axismax), ylim=c(axismin,axismax)) +
        theme_classic()+
        theme(strip.background = element_rect(colour="white", fill="white"))+
        theme(legend.position="none")+
        geom_text(data = stat_raw,
                  mapping = aes(x= min(data_full$field), y = Inf, label = paste("R² = ",round(R2,digits=1))),
                  hjust = 0, vjust = 1, size = 3) +
        geom_text(data = stat_raw,
                  mapping = aes(x= min(data_full$field), y = Inf, label = RMSE_print),
                  hjust = 0, vjust = 2.5, size = 3) 
  
  if (intercept == "_without_intercept"){
    g <- g + geom_smooth(data=data_full.subset_raw,aes(color=sp2_name),method="lm",se=T, formula = y~x-1)
  } else if (intercept == "_with_intercept"){
    g <- g + geom_smooth(data=data_full.subset_raw,aes(color=sp2_name),method="lm",se=T)
  }
  print(g)
  setwd(wd4)
  ggsave(paste(paste(paste("Regression_pairwise_","_field_raw_",sep=dist),intercept,sep=jumpcorr),".jpeg",sep=""),width = 20,height = 15,units=("cm"),dpi=600)
  
  # Unmixed vs.field 
  if (is.matrix(palette_color_unmixed) == T){
    col<-palette_color_unmixed[,"palette_color"]
  } else {
    col <- palette_color_unmixed[["palette_color"]]
  }
  q <- ggplot(data_full,aes(field,unmixed))+
        geom_point(aes(shape=sp2_dvp))+
        scale_shape_manual(name='Phenology',values=palette_shape)+
        xlab(paste(title,"between field patch measurements"))+
        ylab(paste(title,"between unmixed table bouquet meausurements"))+
        facet_rep_wrap(~sp2,repeat.tick.labels = T)+
        # facet_wrap(~sp2,scales="free")+
        # geom_smooth(data=data_full.subset_unmixed,aes(color=sp2_name),method="lm",se=T) + 
        # geom_smooth(aes(color=sp2_name),method="lm")+
        scale_colour_manual(name='Species',values=col) +
        scale_x_continuous(limits=c(axismin,axismax))+
        scale_y_continuous(limits=c(axismin-10,axismax+10))+
        coord_cartesian(xlim=c(axismin,axismax), ylim=c(axismin,axismax)) +
        theme_classic()+
        theme(strip.background = element_rect(colour="white", fill="white"))+
        theme(legend.position="none")+
        geom_text(data = stat_unmixed,
                  mapping = aes(x= min(data_full$field), y = Inf, label = paste("R² = ",round(R2,digits=1))),
                  hjust = 0, vjust = 1, size = 3) +
        geom_text(data = stat_unmixed,
                  mapping = aes(x= min(data_full$field), y = Inf, label = RMSE_print),
                  hjust = 0, vjust = 2.5, size = 3) 
      
  if (intercept == "_without_intercept"){
    q <- q + geom_smooth(data=data_full.subset_unmixed,aes(color=sp2_name),method="lm",se=T, formula = y~x-1)
  } else if (intercept == "_with_intercept"){
    q <- q + geom_smooth(data=data_full.subset_unmixed,aes(color=sp2_name),method="lm",se=T)
  }
  print(q)
  setwd(wd4)
  ggsave(paste(paste(paste("Regression_pairwise_","_field_unmixed_",sep=dist),intercept,sep=jumpcorr),".jpeg",sep=""),width = 20,height = 15,units=("cm"),dpi=600)
  
  # Raw vs.field: old, simple code
  # ylim<-max(data_full$field)
  # ggplot(data_full,aes(field,raw))+
  #   geom_point()+
  #   xlab(paste(dist," between field measurements"))+
  #   ylab(paste(dist," between table meausurement"))+
  #   # scale_x_continuous(limits=c(0,ylim+0.02)) +
  #   # scale_y_continuous(limits=c(0,ylim+0.02)) +
  #   # scale_y_continuous(limits=c(0,max(samdata_full$field)))+
  #   geom_smooth(aes(color=sp2),method="lm")+
  #   facet_rep_wrap(~sp2,repeat.tick.labels = T)+
  #   # facet_wrap(~sp2,scales="free")+
  #   theme_classic()+
  #   theme(strip.background = element_rect(colour="white", fill="white"))+
  #   theme(legend.position="none")
  # setwd("C:/Users/u0091812/Box Sync/03. Black table procedure/Figures/Interspecific/Regression between dissimilarity metrics")
  # ggsave(paste(paste("Regression_pairwise_","_field_raw_",sep=dist),".jpeg",sep=jumpcorr),width = 20,height = 15,units=("cm"),dpi=600)
  
  overview_stat <- cbind(stat_raw,stat_unmixed)
  return(overview_stat)
}


#--------------------------------------------------------------------------------
# APPLY functions

# Calculate dissimiliraty matrices and mantel r's + save statistics
interDissim_mantel_apply <- function (wd1,wd2,jumpcorr) {
  setwd(wd1)
  datasets<-readRDS(paste("datasets_",".rds",sep=jumpcorr))
  
  # spectra data in matrix format --> needed for certain (dis)similarity calculations
  em_field<-spectra(datasets$libr_F.intersect.mean.overall)
  meta_F.intersect.mean.overall<-SI(datasets$libr_F.intersect.mean.overall)
  rownames(em_field)<-meta_F.intersect.mean.overall$Sp_Phase
  em_table_raw<-spectra(datasets$libr_T_raw.intersect.mean.overall)
  em_table_unmixed<-spectra(datasets$libr_T_unmixed.intersect.mean.overall)
  libr_field<-datasets$libr_F.intersect.mean.overall
  libr_table_raw<-datasets$libr_T_raw.intersect.mean.overall
  libr_table_unmixed<-datasets$libr_T_unmixed.intersect.mean.overall
  
  em_field.ind<-spectra(datasets$libr_F.intersect.mean.overall.ind)
  meta_F.intersect.mean.overall.ind<-SI(datasets$libr_F.intersect.mean.overall.ind)
  rownames(em_field.ind)<-meta_F.intersect.mean.overall.ind$Sp_Phase
  em_table_raw.ind<-spectra(datasets$libr_T_raw.intersect.mean.overall.ind)
  em_table_unmixed.ind<-spectra(datasets$libr_T_unmixed.intersect.mean.overall.ind)
  libr_field.ind<-datasets$libr_F.intersect.mean.overall.ind
  libr_table_raw.ind<-datasets$libr_T_raw.intersect.mean.overall.ind
  libr_table_unmixed.ind<-datasets$libr_T_unmixed.intersect.mean.overall.ind
  
  interDissim_field<-interDissim_calc(em_field,libr_field)$InterDissim
  interDissim_table_raw<-interDissim_calc(em_table_raw,libr_table_raw)$InterDissim
  interDissim_table_unmixed<-interDissim_calc(em_table_unmixed,libr_table_unmixed)$InterDissim
  Dissim_names<-interDissim_calc(em_field,libr_field)$Dissim_names
  
  interDissim_field.ind<-interDissim_calc(em_field.ind,libr_field.ind)$InterDissim
  interDissim_table_raw.ind<-interDissim_calc(em_table_raw.ind,libr_table_raw.ind)$InterDissim
  interDissim_table_unmixed.ind<-interDissim_calc(em_table_unmixed.ind,libr_table_unmixed.ind)$InterDissim
  Dissim_names.ind<-interDissim_calc(em_field.ind,libr_field.ind)$Dissim_names
  
  # Calculate mantel r between dissimilarity matrices and store
  mantel_Dissim <- mantel_FieldTable(interDissim_field,interDissim_table_raw,interDissim_table_unmixed,Dissim_names) 
  mantel_Dissim.ind <- mantel_FieldTable(interDissim_field.ind,interDissim_table_raw.ind,interDissim_table_unmixed.ind,Dissim_names.ind) 
  
  setwd(wd2)
  write.table(mantel_Dissim,paste("Mantel_Dissim_",".csv",sep=jumpcorr),row.names=F)
  write.table(mantel_Dissim.ind,paste("Mantel_Dissim.ind_",".csv",sep=jumpcorr),row.names=F)
  
  interDissim_mantel_apply_ouput <- list ("em_field" = em_field,
                                          "Dissim_names" = Dissim_names,
                                          "interDissim_field" = interDissim_field,
                                          "interDissim_table_raw" = interDissim_table_raw,
                                          "interDissim_table_unmixed" = interDissim_table_unmixed,
                                          "mantel_Dissim" = mantel_Dissim,
                                          
                                          "em_field.ind" = em_field.ind,
                                          "Dissim_names.ind" = Dissim_names.ind,
                                          "interDissim_field.ind" = interDissim_field.ind,
                                          "interDissim_table_raw.ind" = interDissim_table_raw.ind,
                                          "interDissim_table_unmixed.ind" = interDissim_table_unmixed.ind,
                                          "mantel_Dissim.ind" = mantel_Dissim.ind)
  
  return(interDissim_mantel_apply_ouput)
}

# Visualise Dissim matrices and save them
DisMatrix_vis_apply <- function(wd3,jumpcorr,interDissim_mantel_apply_ouput) {
  em_field <- interDissim_mantel_apply_ouput$em_field
  Dissim_names <- interDissim_mantel_apply_ouput$Dissim_names
  # Dissim_names<-Dissim_names[!Dissim_names == "SAM2"]
  interDissim_field <- interDissim_mantel_apply_ouput$interDissim_field
  interDissim_table_raw <- interDissim_mantel_apply_ouput$interDissim_table_raw
  interDissim_table_unmixed <- interDissim_mantel_apply_ouput$interDissim_table_unmixed
  
  for (i in 1:length(Dissim_names)) {
    DisMatrix_vis(wd3,jumpcorr,em_field,Dissim_names[i],paste(Dissim_names[i],"_field_raw"),interDissim_field[[Dissim_names[i]]],interDissim_table_raw[[Dissim_names[i]]])
    DisMatrix_vis(wd3,jumpcorr,em_field,Dissim_names[i],paste(Dissim_names[i],"_field_unmixed"),interDissim_field[[Dissim_names[i]]],interDissim_table_unmixed[[Dissim_names[i]]])
  }  
  # upper triangle = field data + lower triangle = raw table data
}
DisRegr_vis_apply <- function(wd4,jumpcorr,interDissim_mantel_apply_ouput, intercept) {
  em_field <- interDissim_mantel_apply_ouput$em_field
  Dissim_names <- interDissim_mantel_apply_ouput$Dissim_names
  # Dissim_names<-Dissim_names[!Dissim_names == "SAM2"]
  interDissim_field <- interDissim_mantel_apply_ouput$interDissim_field
  interDissim_table_raw <- interDissim_mantel_apply_ouput$interDissim_table_raw
  interDissim_table_unmixed <- interDissim_mantel_apply_ouput$interDissim_table_unmixed
  
  titles <- c("SAM "," \nSID "," \nSCM ","Euclidean distance \n", "Manhattan distance", "Bray-Curtis distance \n")
  listofstats <- list()
  for (i in 1:length(Dissim_names)) {
    stat_dissim<- DisRegr_vis(wd4,jumpcorr,em_field,Dissim_names[i],interDissim_field[[Dissim_names[i]]],
                             interDissim_table_unmixed[[Dissim_names[i]]],interDissim_table_raw[[Dissim_names[i]]],titles[i],intercept=intercept)
    listofstats[[i]] <- stat_dissim
    }
  names(listofstats)<-Dissim_names
  
  setwd(wd4)
  saveRDS(listofstats,file=paste("Statistics_overall_",".rds",sep=jumpcorr))
  
  return(listofstats)
}
DisRegr_apart_vis_apply <- function(wd4,jumpcorr,interDissim_mantel_apply_ouput, intercept) {
  em_field <- interDissim_mantel_apply_ouput$em_field
  Dissim_names <- interDissim_mantel_apply_ouput$Dissim_names
  # Dissim_names<-Dissim_names[!Dissim_names == "SAM2"]
  interDissim_field <- interDissim_mantel_apply_ouput$interDissim_field
  interDissim_table_raw <- interDissim_mantel_apply_ouput$interDissim_table_raw
  interDissim_table_unmixed <- interDissim_mantel_apply_ouput$interDissim_table_unmixed
  
  titles <- c("SAM ","SID ","SCM ","Euclidean distance ", "Manhattan distance ", "Bray-Curtis distance ")
  listofstats <- list()
  for (i in 1:length(Dissim_names)) {
      stat_dissim <- DisRegr_apart_vis(wd4,jumpcorr,em_field,
                                      Dissim_names[i],
                                      interDissim_field[[Dissim_names[i]]],
                                      interDissim_table_unmixed[[Dissim_names[i]]],
                                      interDissim_table_raw[[Dissim_names[i]]],
                                      titles[i],
                                      intercept)
      listofstats[[i]] <- stat_dissim
  }
  names(listofstats)<-Dissim_names
  
  setwd(wd4)
  saveRDS(listofstats,file=paste("Statistics_pairwise_",".rds",sep=jumpcorr))

  return(listofstats)
}

#--------------------------------------------------------------------------------
# APPLY to the data

wd1 = "C:/Users/u0091812/Box Sync/03. Black table procedure/Data/Spectral data"
wd2 = "C:/Users/u0091812/Box Sync/03. Black table procedure/Figures/Interspecific/Mantel between dissimilarity matrices/"
wd3 = "C:/Users/u0091812/Box Sync/03. Black table procedure/Figures/Interspecific/Matrices"
wd4 = "C:/Users/u0091812/Box Sync/03. Black table procedure/Figures/Interspecific/Regression between dissimilarity metrics"


Dissim_Nocorr <- interDissim_mantel_apply(wd1,wd2,"Nocorr")
DisMatrix_vis_apply(wd3,"Nocorr",Dissim_Nocorr)
DisRegr_Nocorr_stats_df1 <- DisRegr_vis_apply(wd4,"Nocorr",Dissim_Nocorr, intercept="_without_intercept")
DisRegr_apart_Nocorr_stats_df1 <- DisRegr_apart_vis_apply(wd4,"Nocorr",Dissim_Nocorr, intercept="_without_intercept")
DisRegr_Nocorr_stats_df2 <- DisRegr_vis_apply(wd4,"Nocorr",Dissim_Nocorr, intercept="_with_intercept")
DisRegr_apart_Nocorr_stats_df2 <- DisRegr_apart_vis_apply(wd4,"Nocorr",Dissim_Nocorr, intercept="_with_intercept")

#--------------------------------------------------------------------------------
#  Read created rds files 

jumpcorr <- "Nocorr"

setwd(wd4)
stat_overall <- readRDS(paste("Statistics_overall_",".rds",sep=jumpcorr))
stat_pairwise <- readRDS(paste("Statistics_pairwise_",".rds",sep=jumpcorr))
