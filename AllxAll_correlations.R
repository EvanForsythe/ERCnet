#! /usr/local/bin/Rscript --vanilla --default-packages=utils

#USAGE
#

#Load packages
package_list<-c("igraph")

#Loop to check if package is installed and libraried
for(p in 1:length(package_list)){
  if (!require(package_list[p], character.only = TRUE)) {
    install.packages(package_list[p], dependencies = TRUE)
    library(package_list[p], character.only=TRUE)
  }
}

getScriptPath <- function(){
  cmd.args <- commandArgs()
  m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
  script.dir <- dirname(regmatches(cmd.args, m))
  if(length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
  if(length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
  return(script.dir)
}
# Set the script path:
working_dir<- paste0(getScriptPath(),"/") #added the "/" at the end so paste0 commands below work.
#working_dir<-"/Users/esforsythe/Documents/Work/Bioinformatics/ERC_networks/Analysis/ERCnet_dev/"
setwd(working_dir)

#Because the above function returns "." as the working dir, using this command to set the full path (to avoid issues below)
working_dir<-paste0(getwd(),"/")

#Read in arguments
args = commandArgs(trailingOnly=TRUE)

#Get the job name (used to identify the proper output folder)
jobname<-args[1]
#jobname<-"TEST"
out_dir<-paste0("OUT_", jobname, "/")

### ERC
#Read results in (store as separate variable)
bxb_measure_df_res<-read.table(file = paste0(working_dir, out_dir, "BL_results/bxb_BLs_normalized.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
r2t_measure_df_res<-read.table(file = paste0(working_dir, out_dir, "BL_results/r2t_BLs_normalized.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#Input_N1_df<-read.table(file = paste0(working_dir, "BL_results/N1_input_file.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#Get all pairwise comparisons (without repeats) (this helps avoid nested for loops)
ERC_df<-as.data.frame(t(combn(unique(bxb_measure_df_res$HOG_ID), 2)))
names(ERC_df)<-c("GeneA_HOG", "GeneB_HOG")

#Add columns for results
ERC_results_df<-cbind(ERC_df, 
                      data.frame(Overlapping_branches_BXB=NA, R2_BXB=NA, Slope_BXB=NA, 
                                 Pearson_P_BXB=NA, Spearman_P_BXB=NA, 
                                 Overlapping_branches_R2T=NA, R2_R2T=NA, Slope_R2T=NA, 
                                 Pearson_P_R2T=NA, Spearman_P_R2T=NA))

#Print a messages
print("Beginning all-by-all pairwise ERC correlation analysis...")

#Use cat because paste+print doesn't recognize \n
cat("\nNumber of genes: ", nrow(bxb_measure_df_res), 
    "\nNumber of pair-wise comparisons: ", nrow(ERC_df))

#Loop through the rows
for(e in 1:nrow(ERC_results_df)){
  ## BXB
  #Make a datafrome of the test-set
  BXB_df_temp<-data.frame(
    geneA=paste(bxb_measure_df_res[which(bxb_measure_df_res$HOG_ID==paste(ERC_results_df$GeneA_HOG[e])),]),
    geneB=paste(bxb_measure_df_res[which(bxb_measure_df_res$HOG_ID==paste(ERC_results_df$GeneB_HOG[e])),])
  )
  
  #Clean up dataframe
  BXB_df_temp[-1,]
  BXB_df_temp$geneA<-as.numeric(paste(BXB_df_temp$geneA))
  BXB_df_temp$geneB<-as.numeric(paste(BXB_df_temp$geneB))
  BXB_df_temp_complete<-BXB_df_temp[complete.cases(BXB_df_temp), ]
  
  #Number of overlapping branches
  ERC_results_df$Overlapping_branches_BXB[e]<-nrow(BXB_df_temp_complete)
  
  #Check if there are at least three points (3 are needed for a linear correlation)
  if(ERC_results_df$Overlapping_branches_BXB[e]>2){
    
    #Get lm
    BXB_lm_temp<-lm(BXB_df_temp_complete$geneA~BXB_df_temp_complete$geneB)
    #R2
    ERC_results_df$R2_BXB[e]<-summary(BXB_lm_temp)$r.squared
    #Slope
    ERC_results_df$Slope_BXB[e]<-BXB_lm_temp$coefficients[2]
    #pearson
    ERC_results_df$Pearson_P_BXB[e]<-cor.test(x = BXB_df_temp_complete$geneA, y = BXB_df_temp_complete$geneB, method = "pearson")$p.value
    #spearman
    ERC_results_df$Spearman_P_BXB[e]<-cor.test(x = BXB_df_temp_complete$geneA, y = BXB_df_temp_complete$geneB, method = "spearman")$p.value
  }#End points check if
  
  ## R2T
  #Make a test datafrome
  R2T_df_temp<-data.frame(
    geneA=paste(r2t_measure_df_res[which(r2t_measure_df_res$HOG_ID==paste(ERC_results_df$GeneA_HOG[e])),]),
    geneB=paste(r2t_measure_df_res[which(r2t_measure_df_res$HOG_ID==paste(ERC_results_df$GeneB_HOG[e])),])
  )
  
  #Clean up dataframe
  R2T_df_temp[-1,]
  R2T_df_temp$geneA<-as.numeric(paste(R2T_df_temp$geneA))
  R2T_df_temp$geneB<-as.numeric(paste(R2T_df_temp$geneB))
  R2T_df_temp_complete<-R2T_df_temp[complete.cases(R2T_df_temp), ]
  
  #Number of overlapping branches
  ERC_results_df$Overlapping_branches_R2T[e]<-nrow(R2T_df_temp_complete)
  
  #Check if there are at least three points (3 are needed for a linear correlation)
  if(ERC_results_df$Overlapping_branches_R2T[e]>2){
    
    #Get lm
    R2T_lm_temp<-lm(R2T_df_temp_complete$geneA~R2T_df_temp_complete$geneB)
    #R2
    ERC_results_df$R2_R2T[e]<-summary(R2T_lm_temp)$r.squared
    #Slope
    ERC_results_df$Slope_R2T[e]<-R2T_lm_temp$coefficients[2]
    #pearson
    ERC_results_df$Pearson_P_R2T[e]<-cor.test(x = R2T_df_temp_complete$geneA, y = R2T_df_temp_complete$geneB, method = "pearson")$p.value
    #spearman
    ERC_results_df$Spearman_P_R2T[e]<-cor.test(x = R2T_df_temp_complete$geneA, y = R2T_df_temp_complete$geneB, method = "spearman")$p.value
  }#End points check if
  
}#End pairwise combo loop (variable = e)

#Print figures describing the common branches between pairs of genes
#bxb
pdf(file = paste0(working_dir, out_dir, "Results_ERC/BxB_overlap_hist.pdf"), width=5, height = 5)

hist(ERC_results_df$Overlapping_branches_BXB, 
     breaks = (max(ERC_results_df$Overlapping_branches_BXB)-min(ERC_results_df$Overlapping_branches_BXB)+1), 
     main = "Branch-by-branch", xlab = "Branches in common between \nthe two genes being analyzed")

dev.off()

#r2t
pdf(file = paste0(working_dir, out_dir, "Results_ERC/R2T_overlap_hist.pdf"), width=5, height = 5)

hist(ERC_results_df$Overlapping_branches_R2T, 
     breaks = (max(ERC_results_df$Overlapping_branches_R2T)-min(ERC_results_df$Overlapping_branches_R2T)+1),
     main = "Root-to-tip", xlab = "Branches in common between \nthe two genes being analyzed")

dev.off()

##NETWORKS
#Branch by branch
#Trim down to only look at significant correlation
ERC_results_df_hits_bxb<-subset(ERC_results_df, Pearson_P_BXB<0.01)[,c(-1,-2)]

#make a network graph
network_graph_bxb<-graph.data.frame(ERC_results_df_hits_bxb, directed = FALSE)

#save pdf
pdf(file = paste0(working_dir, out_dir, "Results_ERC/BxB_network.pdf"), width=5, height = 5)

#plot the graph
plot.igraph(network_graph_bxb, vertex.size=10, vertex.label=NA)

dev.off()

#Root to tip
#Trim down to only look at significant correlation
ERC_results_df_hits_r2t<-subset(ERC_results_df, Pearson_P_R2T<0.01)[,c(-1,-2)]

#make a network graph
network_graph_r2t<-graph.data.frame(ERC_results_df_hits_r2t, directed = FALSE)

#save pdf
pdf(file = paste0(working_dir, out_dir, "Results_ERC/R2T_network.pdf"), width=5, height = 5)

#plot the graph
plot.igraph(network_graph_r2t, vertex.size=10, vertex.label=NA)

dev.off()




