#######################################################
###                                                 ###
###        04. Interspecific VI comparisons         ###
###                                                 ###
#######################################################


# In this script we have a closer look at predicted optical traits (vegetation indices; VIs) from field and table spectra.
# We ask:
#   Is the relative position of species in the optical trait space preserved on the table?: pca on optical field and table traits + procrustes analysis
#     (procrustes is only possible when # observations field = # observation table, which is not the case,
#     so we have to make sure to have only one species_SiteDate point for Procrustes analysis)


# There are several methods to quantify optical traits: we focus here on
#   - vegetation indices
# For these analysis we will only use green table data, so no flower spectra, because we believe it only makes sense to predict optical leaf traits from green spectra.

### The input of these steps are the following libraries created in "01. Black table data":
#   For visualising full optical trait space (PCA) 
#       - all INDIVIDUAL bouquets and patch measurements:
#         libr_FG.intersect
#         libr_TG_raw.intersect
#         libr_TG_unmixed.intersect
#   For visualising and comparing optical trait space (PCA + procrustes) --> or use individual meausurements and calculate mean and sd
#       - AVERAGE FIELD and TABLE spectrum per SPECIES AND SITE_DATE 
#         libr_FG.intersect.mean
#         libr_TG_raw.intersect
#         libr_TG_unmixed.intersect

# ! These signatures have been preprocessed with following specifications:
#   - smoothing: savitsky-golay filter with n = 51 and 101 for asd and sv respectively
#   - removal of atmospheric noise windows:  c(349,400,1340,1460,1780,1970,2400,2501)

# Script by Elisa Van Cleemput, 2018
#--------------------------------------------------------------------------------
# Clean workspace

rm(list=ls())

library(hsdar)
# library(devtools)
# install_github("ggbiplot", "vqv")
# library(ggbiplot)
library(randomcoloR) # distinctColorPalette
library(vegan) #envfit and procrustes
library(gginnards) # move_layers # https://cran.r-project.org/web/packages/gginnards/vignettes/user-guide-2.html
library(ggplot2)
library(grid)

#--------------------------------------------------------------------------------

# 1) Function that selects and calculates set of vegetation indices: consult excel sheet 'Vegetation indices_Per trait'
### The following indices seem to be good candidates according to the literature:
    # Chlorophyll total: NDVI(925,710) - MTCI - MCARI2/OSAVI2 - TCARI2/OSAVI2
    # Carotenoids: SR1 - SR2
    # Nitrogen: SR(730, 705) - NDVI(860, 720) - MCARI2/OSAVI2 - MTCI - REP_LE
    # Phosphorus: NDVI(523, 583) - NDVI(1645,1715) - NDVI(860, 720) - SR(730, 705)
    # Potassium: NDVI(523, 583) - NDVI(1645,1715)
    # LMA: NDVI(1662,1732) - NDVI(1550,1750) - NDVI(2260,1490)
    # Water:
    # Photosynthetic efficiency: PRI
# Check references of the indices implemented in the hsdar package
# hsdardocs("References.pdf")
VI_calc <- function (wd1,jumpcorr) {
    setwd(wd1)
    datasets<-readRDS(paste("datasets_",".rds",sep=jumpcorr))
  
  indices<-c("(R925-R710)/(R925+R710)","TCARI2/OSAVI2","MCARI2/OSAVI2","MTCI",
             "REP_LE","PRI","(R860-R720)/(R860+R720)",
             "R539/R490","R807/R490",
             "(R860-R1240)/(R860+R1240)","(R858-R2130)/(R858+R2130)","PRI_norm",
             "(R1662-R1732)/(R1662+R1732)","(R1550-R1750)/(R1550+R1750)","(R1649-R1722)/(R1649+R1722)")
  
  # All individual measurements
  FG_indices <- vegindex(datasets$libr_FG.intersect,indices)
  TG_raw_indices <- vegindex(datasets$libr_TG_raw.intersect,indices)
  TG_unmixed_indices <- vegindex(datasets$libr_TG_unmixed.intersect,indices)
  # Mean spectrum per species for each Site_Date
  FG_indices.mean <- vegindex(datasets$libr_FG.intersect.mean,indices)
  TG_raw_indices.mean <- vegindex(datasets$libr_TG_raw.intersect.mean,indices)
  TG_unmixed_indices.mean <- vegindex(datasets$libr_TG_unmixed.intersect.mean,indices)
  
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


# 2) Function: PCA on indices: all indices or 1 index per trait???
# myggbiplot (based on ggbiplot function) creates PCA figures
myggbiplot <- function (pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE, 
                        obs.scale = 1 - scale, var.scale = scale, groups = NULL, 
                        ellipse = FALSE, ellipse.prob = 0.68, labels = NULL, labels.size = 3, 
                        alpha = 1, var.axes = TRUE, circle = FALSE, circle.prob = 0.69, 
                        varname.size = 3, varname.adjust = 1.5, varname.abbrev = FALSE, 
                        col_arrows,var.labels,
                        ...) {
  library(ggplot2)
  library(plyr)
  library(scales)
  library(grid)
  library(ggrepel) #geom_text_repel
  
  stopifnot(length(choices) == 2)
  if (inherits(pcobj, "prcomp")) {
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$rotation
  }
  else if (inherits(pcobj, "princomp")) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(pcobj$scores, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$loadings
  }
  else if (inherits(pcobj, "PCA")) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(pcobj$ind$coord, 2, 1/(d * nobs.factor), FUN = "*")
    v <- sweep(pcobj$var$coord, 2, sqrt(pcobj$eig[1:ncol(pcobj$var$coord), 
                                                  1]), FUN = "/")
  }
  else if (inherits(pcobj, "lda")) {
    nobs.factor <- sqrt(pcobj$N)
    d <- pcobj$svd
    u <- predict(pcobj)$x/nobs.factor
    v <- pcobj$scaling
    d.total <- sum(d^2)
  }
  else {
    stop("Expected a object of class prcomp, princomp, PCA, or lda")
  }
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale, 
                              FUN = "*"))
  v <- sweep(v, 2, d^var.scale, FUN = "*")
  df.v <- as.data.frame(v[, choices])
  names(df.u) <- c("xvar", "yvar")
  names(df.v) <- names(df.u)
  if (pc.biplot) {
    df.u <- df.u * nobs.factor
  }
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  v.scale <- rowSums(v^2)
  df.v <- r * df.v/sqrt(max(v.scale))
  if (obs.scale == 0) {
    u.axis.labs <- paste("standardized PC", choices, sep = "")
  }
  else {
    u.axis.labs <- paste("PC", choices, sep = "")
  }
  u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%% explained var.)", 
                                            100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
  if (!is.null(labels)) {
    df.u$labels <- labels
  }
  if (!is.null(groups)) {
    df.u$groups <- groups
  }
  
  # Variable names
  if (varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  }
  else {
    df.v$varname <- rownames(v)
  }
  
  # Variables for text label placement
  df.v$angle <- with(df.v, (180/pi) * atan(yvar/xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar))/2)
  
  # Base plot
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + xlab(u.axis.labs[1]) + 
    ylab(u.axis.labs[2]) + coord_equal()

  if (var.axes) {
    if (circle) {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, 
                                                length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * 
                             sin(theta))
      g <- g + geom_path(data = circle, color = muted("white"), 
                         size = 1/2, alpha = 1/3)
    }
    # g <- g + geom_segment(data = df.v, aes(x = 0, y = 0, 
    #                                        xend = xvar, yend = yvar), arrow = arrow(length = unit(1/2, 
    #                                                                                               "picas")), color = muted("red"))
    g <- g + geom_segment(data = df.v, aes(x = 0, y = 0,
                                           xend = xvar, yend = yvar), arrow = arrow(length = unit(1/2, 
                                                                                                  "picas")), color = col_arrows)
  }
  if (!is.null(df.u$labels)) {
    if (!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups), 
                         size = labels.size)
    }
    else {
      g <- g + geom_text(aes(label = labels), size = labels.size)
    }
  }
  else {
    if (!is.null(df.u$groups)) {
      g <- g + geom_point(aes(color = groups), alpha = alpha)
    }
    else {
      g <- g + geom_point(alpha = alpha)
    }
  }
  if (!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    ell <- ddply(df.u, "groups", function(x) {
      if (nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, 
                       mu, FUN = "+"), groups = x$groups[1])
    })
    names(ell)[1:2] <- c("xvar", "yvar")
    g <- g + geom_path(data = ell, aes(color = groups, group = groups))
  }
  if (var.labels) {
    # g <- g + geom_text(data = df.v, aes(label = varname, 
    #                                     x = xvar, y = yvar, angle = angle, hjust = hjust), 
    #                    color = "darkred", size = varname.size)
    # g <- g + geom_text(data = df.v, aes(label = varname,
    #                                     x = xvar, y = yvar, angle = angle, hjust = hjust),
    #                    color = col_arrows, size = varname.size)
    g <- g + geom_text_repel(data = df.v, aes(label = varname,
                                              x = xvar, y = yvar, angle = angle),
                             color = col_arrows, size = varname.size)
    
  }
  return(g)
}
# there are several options for PCA: princomp function or prcomp function (results and figures are the same :-))
# Save figure with and without variable labels, and add labels to the latter figure in inkscape
VI_pca <- function (wd1,wd2,wd3,jumpcorr,datasets_VI,VI_data,libr_data,labels) {
  graphics.off()
  
  setwd(wd1)
  datasets<-readRDS(paste("datasets_",".rds",sep=jumpcorr))
  libr<-datasets[[libr_data]]
  meta<-SI(libr)
  groups<-SI(libr)$Species
  VI_datavalues<-datasets_VI[[VI_data]]
  
  # 2 ways of calculating pca
  prin<-princomp((VI_datavalues), cor = TRUE, scores = TRUE)
  prc<-prcomp(VI_datavalues,scale=TRUE,center=TRUE)
  
  summaries <- list("var_prc" = summary(prc),
                    "var_prin" = summary(prin),
                    "loadings_prc" = prc$rotation,# variable loadings = correlation between variable and PC
                    "loadings_prin" = prin$loadings,
                    "scores_prc" = prc$x, # coordinates of the observations on the PC's
                    "scores_prin" = prin$scores)
  
  # --------- plot results using ggbiplot  --> based on prc
  # x11(width=1500,height=720)
  setwd(wd2)
  # ggsave(paste(paste("Screeplot_","_",sep=VI_data),".jpeg",sep=jumpcorr),width = 20,height = 15,units=("cm"),dpi=600)
  jpeg(paste(paste("Screeplot_","_",sep=VI_data),".jpeg",sep=jumpcorr),res=300,width=14,height=7,units='in')
  par(mfrow=c(1,2))
  screeplot(prc) # It seems reasonable to withhold 3 PC's (also when looking at proportion of variance they explain)
  abline(h=1, col='blue')
  plot(prc,type='l')
  abline(h=1, col='blue')
  dev.off()
  
  if (labels == "species") {
    vis<-meta[,'Species']  # Indicate what should be visualised on the plot:
    setwd(wd2)
  } else if (labels == "sensors") {
    library(stringr) #str_sub
    vis<-str_sub(meta[,1] ,-3,-1) # Indicate what should be visualised on the plot:
    setwd(wd3)
  }
  
  # n <- nlevels(droplevels(vis))
  # palette_color <- distinctColorPalette(n) # find a good color palette. E.g.:
  palette_color<-c("#73E0B0","#E06762","#CCBAD8","#E0D64D","#804BDD","#D3DFCD","#8AE659","#7AC8D7","#DC8DC8","#C3DB8D","#7A88D4","#DA50CB","#D4A67E")
  # palette_shape <- c(0,1,2,3,4,5,6,7,8,15,16,17,18)
  palette_shape <- c(15,16,17,18,15,16,17,18,15,16,17,18,15)
  # x11(width=1500,height=720)
  g <- myggbiplot(prc, pc.biplot=T, obs.scale = 1, var.scale = 1, #choices = 2:3,
                  ellipse = F, circle = F, var.axes=T,alpha=0, # specify groups = vis and remove alpha = 0 when no shape differences
                  col_arrows = "dimgray",var.labels=T) + 
    geom_point(aes(colour=vis, shape=vis), size = 3) +
    scale_color_manual(name='Species',values=palette_color) +
    scale_shape_manual(name='Species',values=palette_shape)+
    theme(legend.direction = 'vertical', 
          legend.position = 'right') +
    theme(legend.position="none")+
    scale_x_continuous(breaks=seq(floor(min(prc$x[,1]))+floor(min(prc$x[,1]))%%2,ceiling(max(prc$x[,1]))+ceiling(max(prc$x[,1]))%%2,2))+
    scale_y_continuous(breaks=seq(floor(min(prc$x[,2]))+floor(min(prc$x[,2]))%%2,ceiling(max(prc$x[,2]))+ceiling(max(prc$x[,2]))%%2,2))+
    # scale_x_continuous(sec.axis = dup_axis()) + scale_y_continuous(sec.axis = dup_axis()) +
    theme_bw()+
    # ggtitle("Mean Optical traits of unmixed Table Green measurements")+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  g_reorder <- move_layers(g,idx=4L,position="bottom")
  print(g_reorder)

  ggsave(paste(paste("pca12_","_",sep=VI_data),".jpeg",sep=jumpcorr),width = 20,height = 15,units=("cm"),dpi=600)
  
  p <- myggbiplot(prc, obs.scale = 1, var.scale = 1, #choices = 2:3,
                  ellipse = F, circle = F, var.axes=T,alpha=0, # specify groups = vis and remove alpha = 0 when no shape differences
                  col_arrows = "dimgray",var.labels=F) + 
    geom_point(aes(colour=vis, shape=vis), size = 3) +
    scale_color_manual(name='Species',values=palette_color) +
    scale_shape_manual(name='Species',values=palette_shape)+
    theme(legend.direction = 'vertical', 
          legend.position = 'right') +
    theme(legend.position="none")+
    scale_x_continuous(breaks=seq(floor(min(prc$x[,1]))+floor(min(prc$x[,1]))%%2,ceiling(max(prc$x[,1]))+ceiling(max(prc$x[,1]))%%2,2))+
    scale_y_continuous(breaks=seq(floor(min(prc$x[,2]))+floor(min(prc$x[,2]))%%2,ceiling(max(prc$x[,2]))+ceiling(max(prc$x[,2]))%%2,2))+
    theme_bw()+
    # ggtitle("Mean Optical traits of unmixed Table Green measurements")+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  p_reorder <- move_layers(p,idx=3L,position="bottom")
  print(p_reorder)
  # setwd("C:/Users/u0091812/Box Sync/03. Black table procedure/Figures/Interspecific/VI_PCA")
  ggsave(paste(paste("pca12_","_",sep=VI_data),"_NoVarLabels.jpeg",sep=jumpcorr),width = 20,height = 15,units=("cm"),dpi=600)
  
  # x11(width=1500,height=720)
  # library(rgl)
  # jpeg("test.jpg")
  # # x11(width=1500,height=720)
  # plot3d(prin$scores[,1:3], col=vis)
  # dev.off()
  
  # Save legend separately
  m <- myggbiplot(prc, obs.scale = 1, var.scale = 1, #choices = 2:3,
                  ellipse = F, circle = F, var.axes=T,alpha=0, # specify groups = vis and remove alpha = 0 when no shape differences
                  col_arrows = "black",var.labels=F) + 
    geom_point(aes(colour=vis, shape=vis), size = 3) +
    scale_color_manual(name='Species',values=palette_color) +
    scale_shape_manual(name='Species',values=palette_shape)+
    theme(legend.direction = 'vertical', 
          legend.position = 'right') +
    theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  library(cowplot)
  legend<-get_legend(m)
  ggdraw(plot_grid(legend,ncol=1))
  # setwd("C:/Users/u0091812/Box Sync/03. Black table procedure/Figures/Interspecific/VI_PCA")
  ggsave("PCA_legend.jpeg",width = 8,height = 16,units=("cm"),dpi=600)
  
  # --------- plot results with parts of the code adopted from Díaz et al.2015
  pc12<-prin$scores[,1:2]
  # pc12[,1]<-pc12[,1]*-1
  # pc12[,2]<-pc12[,2]*-1
  
  ll<-prin$loadings
  # x11(width=1500,height=720)
  plot(prin,type='l',col='red',main="screeplot")
  abline(h=1, col='blue')
  fit<-envfit(pc12, VI_datavalues) # use envfit(vegan package) for drawing arrows, can be also done using trait loadings
  fit2<-fit$vectors$arrows*-2.5 # drawing line segments in arrow opposites direction for pretty layout
  plot( pc12[,], pch=16, cex=0.5, col="black")
  plot(fit, cex=0.90, col=1, labels=list(vectors = colnames(VI_datavalues)))
  segments(0,0, fit2[,1], fit2[,2], col=1, lty=2, lwd=1)
  mtext("PC1", cex=0.75, side=1, line=0.5, adj=1)
  mtext("PC2", cex=0.75, side=2, line=0.5, at=4.7) #, las=2)
  
  return(summaries)
}
VI_pca_apply <- function (wd1,wd2,wd3,jumpcorr,VI_output,labels) {
  sum_VI_FG <- VI_pca(wd1,wd2,wd3,jumpcorr,VI_output,"VI_FG","libr_FG.intersect",labels)
  sum_VI_TG_raw <- VI_pca(wd1,wd2,wd3,jumpcorr,VI_output,"VI_TG_raw","libr_TG_raw.intersect",labels)
  sum_VI_TG_unmixed <- VI_pca(wd1,wd2,wd3,jumpcorr,VI_output,"VI_TG_unmixed","libr_TG_unmixed.intersect",labels)
  
  sum_VI_FG.mean <- VI_pca(wd1,wd2,wd3,jumpcorr,VI_output,"VI_FG.mean","libr_FG.intersect.mean",labels)
  sum_VI_TG_raw.mean <- VI_pca(wd1,wd2,wd3,jumpcorr,VI_output,"VI_TG_raw.mean","libr_TG_raw.intersect.mean",labels)
  sum_VI_TG_unmixed.mean <- VI_pca(wd1,wd2,wd3,jumpcorr,VI_output,"VI_TG_unmixed.mean","libr_TG_unmixed.intersect.mean",labels)
  
  pca_summaries <- list("sum_VI_FG" = sum_VI_FG,
                        "sum_VI_TG_raw" = sum_VI_TG_raw,
                        "sum_VI_TG_unmixed" = sum_VI_TG_unmixed,
                        "sum_VI_FG.mean" = sum_VI_FG.mean,
                        "sum_VI_TG_raw.mean" = sum_VI_TG_raw.mean,
                        "sum_VI_TG_unmixed.mean" = sum_VI_TG_unmixed.mean)
}


# 3) Function that performs and visualises procrustes analysis between two PCA results
#  procrustes on mean of each Species_Site_Date
VI_pca_procrustes <- function (wd1, wd4, jumpcorr, pca_spec) {
  # Two ordinations can be very similar but this similarity may be masked as a result of the two ordinations
  # having different scalings, orientations and signs
  # the main difference between basic Procrustes and PROTEST is that protest is always symmetric whereas procrustes defaults to non-symmetric
  
  # table measurements are rotated to the field measurements. Arrows from field to table
  proc_F_Traw<-procrustes(pca_spec$sum_VI_FG.mean$scores_prc,pca_spec$sum_VI_TG_raw.mean$scores_prc,symmetric=F)
  proc_F_Tunmixed<-procrustes(pca_spec$sum_VI_FG.mean$scores_prc,pca_spec$sum_VI_TG_unmixed.mean$scores_prc,symmetric=F)
  prot_F_Traw<-protest(pca_spec$sum_VI_FG.mean$scores_prc,pca_spec$sum_VI_TG_raw.mean$scores_prc)
  prot_F_Tunmixed<-protest(pca_spec$sum_VI_FG.mean$scores_prc,pca_spec$sum_VI_TG_unmixed.mean$scores_prc)
  
  # Simple plots
  par(mfrow=c(1,2))
  plot(proc_F_Traw,kind="1")
  plot(proc_F_Traw,kind="2") # residual plot
  
  plot(proc_F_Tunmixed)

  # Advanced plots: separate colour for each species (https://stackoverflow.com/questions/30325739/ggplot2-for-procrustes-rotation-in-vegan)
  # Procrustes residual plots
  setwd(wd1)
  datasets<-readRDS(paste("datasets_",".rds",sep=jumpcorr))
  groups<-SI(datasets$libr_FG.intersect.mean)$Species
  groups2<-SI(datasets$libr_TG_raw.intersect.mean)$Species
  groups3<-SI(datasets$libr_TG_unmixed.intersect.mean)$Species
plotpro_F_Traw.sym <- data.frame(pca2x=prot_F_Traw$Yrot[,1],
                                     pca2y=prot_F_Traw$Yrot[,2],
                                     pca1x=prot_F_Traw$X[,1],
                                     pca1y=prot_F_Traw$X[,2],
                                     species=groups, species2=groups2)
plotpro_F_Traw.asym <- data.frame(pca2x=proc_F_Traw$Yrot[,1],
                                     pca2y=proc_F_Traw$Yrot[,2],
                                     pca1x=proc_F_Traw$X[,1],
                                     pca1y=proc_F_Traw$X[,2],
                                     species=groups, species2=groups2)
plotpro_F_Tunmixed.sym <- data.frame(pca2x=prot_F_Tunmixed$Yrot[,1],
                                     pca2y=prot_F_Tunmixed$Yrot[,2],
                                     pca1x=prot_F_Tunmixed$X[,1],
                                     pca1y=prot_F_Tunmixed$X[,2],
                                     species=groups, species3=groups3)
plotpro_F_Tunmixed.asym <- data.frame(pca2x=proc_F_Tunmixed$Yrot[,1],
                                      pca2y=proc_F_Tunmixed$Yrot[,2],
                                      pca1x=proc_F_Tunmixed$X[,1],
                                      pca1y=proc_F_Tunmixed$X[,2],
                                      species=groups, species3=groups3)

setwd(wd4)
palette_color<-c("#73E0B0","#E06762","#CCBAD8","#E0D64D","#804BDD","#D3DFCD","#8AE659","#7AC8D7","#DC8DC8","#C3DB8D","#7A88D4","#DA50CB","#D4A67E")
palette_shape <- c(15,16,17,18,15,16,17,18,15,16,17,18,15)

proc_plot <- function(plotpro_data,protest,sym){
  # x11(width=720,height=720)
  u <- ggplot(plotpro_data) +
    geom_point(aes(x=pca1x, y=pca1y, colour=groups, shape=groups),size=3) +
    # geom_point(aes(x=pca2x, y=pca2y, colour=groups, shape=groups),size=3) +
    geom_segment(aes(x=pca1x,y=pca1y,xend=pca2x,yend=pca2y,colour=groups),arrow=arrow(length=unit(0.2,"cm")),size=1)+
    scale_color_manual(name='Species',values=palette_color) +
    scale_shape_manual(name='Species',values=palette_shape)+
    annotate("text", x = Inf, y = -Inf,hjust=1.1,vjust=-3,size=5,
             label = paste("correlation",round(protest$t0,digits=2),sep=" = "), parse=F)+ # Add r and signifcance: find out how they are stored
    annotate("text", x = Inf, y = -Inf,hjust=1.1,vjust=-1,size=5,
             label = paste("p value",round(protest$signif,digits=3),sep=" = "), parse=F)+ # Add r and signifcance: find out how they are stored
    labs(x="PC1", y="PC2") +
    theme_bw()+ theme(text = element_text(size=16))+
    theme(legend.key.size=unit(2,'lines')) +
    theme(legend.position="none")+
    coord_equal()
  
  if (sym == "asym") {
     u <- u + scale_x_continuous(breaks=seq(floor(min(plotpro_data$pca1x,plotpro_data$pca2x))+floor(min(plotpro_data$pca1x,plotpro_data$pca2x))%%2,ceiling(max(plotpro_data$pca1x,plotpro_data$pca2x))+ceiling(max(plotpro_data$pca1x,plotpro_data$pca2x))%%2,2))+
              scale_y_continuous(breaks=seq(floor(min(plotpro_data$pca1y,plotpro_data$pca2y))+floor(min(plotpro_data$pca1y,plotpro_data$pca2y))%%2,ceiling(max(plotpro_data$pca1y,plotpro_data$pca2y))+ceiling(max(plotpro_data$pca1y,plotpro_data$pca2y))%%2,2))
    } 
  else if (sym == "sym") {
     u <- u + scale_x_continuous(breaks=seq(round(min(plotpro_data$pca1x,plotpro_data$pca2x),digits=1),max(round(plotpro_data$pca1x,plotpro_data$pca2x),digits=1),0.1))+
              scale_y_continuous(breaks=seq(round(min(plotpro_data$pca1y,plotpro_data$pca2y),digits=1),round(max(plotpro_data$pca1y,plotpro_data$pca2y),digits=1),0.1))
    }
    # print(u)
}

u <- ggplot(plotpro_F_Traw.sym) +
  geom_point(aes(x=pca1x, y=pca1y, colour=groups, shape=groups),size=3) +
  geom_segment(aes(x=pca1x,y=pca1y,xend=pca2x,yend=pca2y,colour=groups),arrow=arrow(length=unit(0.2,"cm")),size=1)+
  scale_color_manual(name='Species',values=palette_color) +
  scale_shape_manual(name='Species',values=palette_shape)+
  theme_bw()+ theme(text = element_text(size=16))+
  theme(legend.key.size=unit(2,'lines')) +
  coord_equal()
library(cowplot)
legend<-get_legend(u)
ggdraw(plot_grid(legend,ncol=1))
setwd(wd4)
ggsave("Procrustes_legend.jpeg",width = 8,height = 16,units=("cm"),dpi=600)


proc_plot(plotpro_F_Traw.sym,prot_F_Traw,sym="sym")
ggsave(paste("proc_F_Traw_sym_",".jpeg",sep=jumpcorr),width = 20,height = 15,units=("cm"),dpi=600)
proc_plot(plotpro_F_Traw.asym,prot_F_Traw,sym="asym")
ggsave(paste("proc_F_Traw_asym_",".jpeg",sep=jumpcorr),width = 20,height = 15,units=("cm"),dpi=600)
proc_plot(plotpro_F_Tunmixed.sym,prot_F_Tunmixed,sym="sym")
ggsave(paste("proc_F_Tunmixed_sym_",".jpeg",sep=jumpcorr),width = 20,height = 15,units=("cm"),dpi=600)
proc_plot(plotpro_F_Tunmixed.asym,prot_F_Tunmixed,sym="asym")
ggsave(paste("proc_F_Tunmixed_asym_",".jpeg",sep=jumpcorr),width = 20,height = 15,units=("cm"),dpi=600)

}

#--------------------------------------------------------------------------------

# APPLY the functions
wd1 = "C:/Users/u0091812/Box Sync/03. Black table procedure/Data/Spectral data"
wd2 = "C:/Users/u0091812/Box Sync/03. Black table procedure/Figures/Interspecific/VI_PCA"
wd3 = "C:/Users/u0091812/Box Sync/03. Black table procedure/Figures/Interspecific/VI_PCA_instrument"
wd4 = "C:/Users/u0091812/Box Sync/03. Black table procedure/Figures/Interspecific/VI_PCA_procrustes"

VI_Nocorr <- VI_calc(wd1,"Nocorr")
pca_Nocorr <- VI_pca_apply(wd1,wd2,wd3,"Nocorr",VI_Nocorr, labels="species")
VI_pca_procrustes(wd1, wd4,"Nocorr",pca_Nocorr)

# The same pca as pca_Nocorrbut with instruments coloured differently: figure included in response to reviewers
pca_Nocorr_sensors <- VI_pca_apply(wd1,wd2,wd3,"Nocorr",VI_Nocorr, labels="sensors")


# Store PCa information as csv
setwd(wd2)
write.table(pca_Nocorr$sum_VI_FG$loadings_prc,"pca_Nocorr_FG.csv",col.names=T,row.names=T)
write.table(pca_Nocorr$sum_VI_TG_raw$loadings_prc,"pca_Nocorr_TG_raw.csv",col.names=T,row.names=T)
write.table(pca_Nocorr$sum_VI_TG_unmixed$loadings_prc,"pca_Nocorr_TG_unmixed.csv",col.names=T,row.names=T)
write.table(pca_Nocorr$sum_VI_FG.mean$loadings_prc,"pca_Nocorr_FG.mean.csv",col.names=T,row.names=T)
write.table(pca_Nocorr$sum_VI_TG_raw.mean$loadings_prc,"pca_Nocorr_TG_raw.mean.csv",col.names=T,row.names=T)
write.table(pca_Nocorr$sum_VI_TG_unmixed.mean$loadings_prc,"pca_Nocorr_TG_unmixed.mean.csv",col.names=T,row.names=T)

