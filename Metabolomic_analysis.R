# Clean global environment
rm(list = ls())

##Packages
library(vroom)
library(MetaboAnalystR)
library(plotly)
library(iheatmapr)
library(pacman)
library(dplyr)

##Set your home folder/working folder
#work computer
working<-("users/your_file_location")
setwd(working)
old_wd<-(working)


##Load data
#work computer

data <- read.csv("metabolomics.csv")



run_stat<- function (feature){
  setwd(old_wd)
  file_name<- paste(feature,"_2025_01_15_")
  new_folder<-dir.create(file_name)
  newfile<- cbind(data[,1],data[,feature], data[,105:ncol(data)])##adjust based on the metabolites in the file you want to analyze
  setwd(file_name)
  write.csv(newfile, paste(feature, ".csv", sep=""), row.names = FALSE)
  
  
  # Load MetaboAnalystR
  mSet<-InitDataObjects("pktable", "stat", FALSE)
  
  #make sure that your file is clean and that there are no duplicates
  mSet<-Read.TextData(mSet,paste(feature,".csv",sep = ""),"rowu", "disc")
  mSet<-SanityCheckData(mSet);
  
  
  # normalization parameters
  mSet<-ReplaceMin(mSet);
  mSet<-PreparePrenormData(mSet);
  mSet<-(Normalization(mSet, "MedianNorm", "NULL", "AutoNorm", ratio=FALSE, ratioNum=20))
  
  #normalization plots optional
  mSet<-PlotNormSummary(mSet, "norm_0_", format ="png", dpi=72, width=NA);
  # mSet<-PlotSampleNormSummary(mSet, "snorm_0_", format = "png", dpi=72, width=NA);
  
  
  ##Perform if only 2 groups else multi group
  # Perform fold-change analysis on uploaded data, unpaired
  if(mSet$dataSet$cls.num==2){
    # Plot fold-change analysis
    mSet<-FC.Anal(mSet, 2.0, 0, FALSE); 
    mSet<-PlotFC(mSet, "fc_0_", "png", 72, width=NA);
    
    # Perform T-test (parametric)
    mSet<-Ttests.Anal(mSet, nonpar=F, threshp=0.05, paired=FALSE, equal.var=TRUE, "fdr", FALSE);
    
    ##For all groups Orthoganl Partial Least Squares 
    # Perform oPLS-DA analysis
    mSet<-OPLSR.Anal(mSet, reg=TRUE);
    
    # Create a 2D oPLS-DA score plot
    mSet<-PlotOPLS2DScore(mSet, "opls_score2d_0_", format = "png", dpi=72, width=NA, 1,2,0.95,1,0);
    
    # Create a significant features plot
    mSet<-PlotOPLS.Splot(mSet, "opls_splot_0_", "all", "png", 72, width=NA);
    
    # Create a plot of features ranked by importMCVe
    mSet<-PlotOPLS.Imp(mSet, "opls_imp_0_", "png", 72, width=NA, "vip", "tscore", 15,FALSE);
    
    # Create a plot of the model overview
    mSet<-PlotOPLS.MDL(mSet, "opls_mdl_0_", format = "png", dpi=72, width=NA);
    
    # Perform and plot oPLS-DA permutation 
    mSet<-OPLSDA.Permut(mSet, 100);
    
    mSet<-PlotOPLS.Permutation(mSet, "opls_perm_2_", format = "png", dpi=72, width=NA);
    
    # Perform hierarchical clustering and plot heat map
    #mSet<-PlotHeatMap(mSet, "heatmap_0_", "png", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", 8, "overview", T, T, NULL, T, F, T, T, T);
    
    # Perform PCA analysis
    mSet<-PCA.Anal(mSet)
    
    # Create PCA overview
    mSet<-PlotPCAPairSummary(mSet, "pca_pair_0_", format = "png", dpi = 72, width=NA, 5)
    
    # Create PCA scree plot
    mSet<-PlotPCAScree(mSet, "pca_scree_0_", "png", dpi = 72, width=NA, 5)
    
    # Create a 2D PCA score plot
    mSet<-PlotPCA2DScore(mSet, "pca_score2d_0_", format = "png", dpi=72, width=NA, 1, 2, 0.95, 1, 0)
    
    # Create a 3D PCA score plot
    mSet<-PlotPCA3DScoreImg(mSet, "pca_score3d_0_", "png", 72, width=NA, 1,2,3, 40)
    
    # Create a PCA loadings Plots
    mSet<-PlotPCALoading(mSet, "pca_loading_0_", "png", 72, width=NA, 1,2);
    
    # Create a PCA Biplot
    mSet<-PlotPCABiplot(mSet, "pca_biplot_0_", format = "png", dpi = 72, width=NA, 1, 2)
    
    
    setwd(old_wd)
  } else{
    
    ##For all groups Orthoganl Partial Least Squares
    # Perform oPLS-DA analysis
    mSet<-OPLSR.Anal(mSet, reg=TRUE);
    
    # Create a 2D oPLS-DA score plot
    mSet<-PlotOPLS2DScore(mSet, "opls_score2d_0_", format = "png", dpi=72, width=NA, 1,2,0.95,1,0);
    
    # Create a significant features plot
    mSet<-PlotOPLS.Splot(mSet, "opls_splot_0_", "all", "png", 72, width=NA);
    
    # Create a plot of features ranked by importMCVe
    mSet<-PlotOPLS.Imp(mSet, "opls_imp_0_", "png", 72, width=NA, "vip", "tscore", 15,FALSE);
    
    # Create a plot of the model overview
    mSet<-PlotOPLS.MDL(mSet, "opls_mdl_0_", format = "png", dpi=72, width=NA);
    
    # Perform and plot oPLS-DA permutation
    mSet<-OPLSDA.Permut(mSet, 100);
    
    mSet<-PlotOPLS.Permutation(mSet, "opls_perm_2_", format = "png", dpi=72, width=NA);
    
    # Perform PCA analysis
    mSet<-PCA.Anal(mSet)
    
    # Create PCA overview
    mSet<-PlotPCAPairSummary(mSet, "pca_pair_0_", format = "png", dpi = 72, width=NA, 5)
    
    # Create PCA scree plot
    mSet<-PlotPCAScree(mSet, "pca_scree_0_", "png", dpi = 72, width=NA, 5)
    
    # Create a 2D PCA score plot
    mSet<-PlotPCA2DScore(mSet, "pca_score2d_0_", format = "png", dpi=72, width=NA, 1, 2, 0.95, 1, 0)
    
    # Create a 3D PCA score plot
    mSet<-PlotPCA3DScoreImg(mSet, "pca_score3d_0_", "png", 72, width=NA, 1,2,3, 40)
    
    # Create a PCA loadings Plots
    mSet<-PlotPCALoading(mSet, "pca_loading_0_", "png", 72, width=NA, 1,2);
    
    # Create a PCA Biplot
    mSet<-PlotPCABiplot(mSet, "pca_biplot_0_", format = "png", dpi = 72, width=NA, 1, 2)
    
    # Perform hierarchical clustering and plot heat map
    mSet<-PlotHeatMap(mSet, "heatmap_0_", "png", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", 8, "overview", T, T, NULL, T, F, T, T, T);
    
    # Perform ANOVA
    mSet <- ANOVA.Anal(mSet, F, 0.05, "fisher");
    # Plot ANOVA
    mSet <- PlotANOVA(mSet, "aov_0_", "png", 72, width=NA);
    
    # mSet<-Calculate.ANOVA.posthoc(mSet, "fisher", 0.05)
    setwd(old_wd)
  }
  
}

run_stat("ELN_risk")
