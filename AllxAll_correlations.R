#! /usr/local/bin/Rscript --vanilla --default-packages=utils

#USAGE
#

##Load packages
#package_list<-c()
#
##Loop to check if package is installed and libraried
#for(p in 1:length(package_list)){
#  if (!require(package_list[p], character.only = TRUE)) {
#    install.packages(package_list[p], dependencies = TRUE)
#    library(package_list[p], character.only=TRUE)
#  }
#}

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
#jobname<-"Clptest"
#jobname<-"FULL"
out_dir<-paste0("OUT_", jobname, "/")

#Get the focal species
foc_sp<-args[2]
#foc_sp<-"A_thaliana_prot"

#read in ERC results file
ERC_results_df<-read.table(paste0(working_dir, out_dir, "ERC_results/ERC_results.tsv"), sep = "\t", header = TRUE)
#ERC_results_df<-read.table(paste0(working_dir, out_dir, "ERC_results/ERC_results.tsv"), sep = "\t", header = TRUE, nrows=1000000)

#Clean data to write tsv file of the ERC results
#Reorder by R2T pearson R2
ERC_results_df_order<-ERC_results_df[order(ERC_results_df$Pearson_P_R2T),]

#Read in file that contains sequence IDs
All_HOGs_df<-read.table(file = paste0(working_dir, out_dir, "Filtered_genefam_dataset.csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)

#Remove the "N1." (or N2. N3. etc...) from the string
All_HOGs_df$HOG<-sapply(strsplit(as.character(All_HOGs_df$HOG), "\\."), `[`, 2)

All_HOGs_df_trimmed<-cbind(All_HOGs_df[,1:3], data.frame(
  Focal_sp_ID=sapply(strsplit(All_HOGs_df[,which(names(All_HOGs_df) == foc_sp)], ","), `[`, 1)
))

#Add space to include gene IDs
ERC_results_df_plus<-cbind(ERC_results_df_order, data.frame(GeneA_ID=NA, GeneB_ID=NA))

#Loop through rows and add gene ID in
for(n in 1:nrow(ERC_results_df_plus)){
  ERC_results_df_plus$GeneA_ID[n]<-All_HOGs_df_trimmed$Focal_sp_ID[which(All_HOGs_df_trimmed$HOG == ERC_results_df_plus$GeneA_HOG[n])]
  ERC_results_df_plus$GeneB_ID[n]<-All_HOGs_df_trimmed$Focal_sp_ID[which(All_HOGs_df_trimmed$HOG == ERC_results_df_plus$GeneB_HOG[n])]
  #if(n%%1000 == 0){
  #  print(n)
  #  }
  }

#Write the table
write.table(ERC_results_df_plus, file = paste0(working_dir, out_dir, "ERC_results/ERC_results_ordered_full.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

#Print figures describing the common branches between pairs of genes
#bxb
pdf(file = paste0(working_dir, out_dir, "ERC_results/BxB_overlap_hist.pdf"), width=5, height = 5)

hist(ERC_results_df$Overlapping_branches_BXB, 
     breaks = (max(na.omit(as.numeric(paste(ERC_results_df$Overlapping_branches_BXB))))-min(na.omit(as.numeric(paste(ERC_results_df$Overlapping_branches_BXB)))+1)), 
     main = "Branch-by-branch", xlab = "Branches in common between \nthe two genes being analyzed")

dev.off()

#r2t
pdf(file = paste0(working_dir, out_dir, "ERC_results/R2T_overlap_hist.pdf"), width=5, height = 5)

hist(ERC_results_df$Overlapping_branches_R2T, 
     breaks = (max(na.omit(as.numeric(paste(ERC_results_df$Overlapping_branches_R2T))))-min(na.omit(as.numeric(paste(ERC_results_df$Overlapping_branches_R2T)))+1)),
     main = "Root-to-tip", xlab = "Branches in common between \nthe two genes being analyzed")

dev.off()


#Print figures describing coorelation between different strategies

#Linear model
mod0<-lm(Pearson_R2_R2T~Pearson_R2_BXB, data=ERC_results_df_plus)
modsum0 = summary(mod0)

#r2t vs bxb (Pearson)
jpeg(file = paste0(working_dir, out_dir, "ERC_results/Branch_length_strategies_pearson.jpg"), quality = 75)

plot(ERC_results_df_plus$Pearson_R2_BXB, ERC_results_df_plus$Pearson_R2_R2T,
     xlab ="Branch-by-branch R-squared",
     ylab ="Root-to-tip R-squared",
     col = rgb(red = 0, green = 0, blue = 0, alpha = 0.2),
     main = paste0("Pearson R-squared values calculated from \nR2T vs BXB branch lengths", "\nR2=", format(modsum0$adj.r.squared, digits = 3))
     )

#add trend line
abline(mod0)

dev.off()


#r2t vs bxb (Spearman)
#Linear model
mod1<-lm(Spearman_R2_R2T~Spearman_R2_BXB, data=ERC_results_df_plus)
modsum1 = summary(mod1)

jpeg(file = paste0(working_dir, out_dir, "ERC_results/Branch_length_strategies_spearman.jpg"), quality = 75)

plot(ERC_results_df_plus$Spearman_R2_BXB, ERC_results_df_plus$Spearman_R2_R2T,
     xlab ="Branch-by-branch R-squared",
     ylab ="Root-to-tip R-squared",
     col = rgb(red = 0, green = 0, blue = 0, alpha = 0.2),
     main = paste0("Spearman R-squared values calculated from \nR2T vs BXB branch lengths", "\nR2=", format(modsum1$adj.r.squared, digits = 3))
     )

#add trend line
abline(mod1)

dev.off()


#Pearson vs Spearman (R2T)
#Linear model
mod2<-lm(Spearman_R2_R2T~Pearson_R2_R2T, data=ERC_results_df_plus)
modsum2 = summary(mod2)

jpeg(file = paste0(working_dir, out_dir, "ERC_results/Correlation_strategies_R2T.jpg"), quality = 75)

plot(ERC_results_df_plus$Pearson_R2_R2T, ERC_results_df_plus$Spearman_R2_R2T,
     xlab ="Pearson R-squared",
     ylab ="Spearman R-squared",
     col = rgb(red = 0, green = 0, blue = 0, alpha = 0.2),
     main = paste0("R2T R-squared values calculated from \nPearson vs Spearman correlation", "\nR2=", format(modsum2$adj.r.squared, digits = 3))
     )
     
#add trend line
abline(mod2)

dev.off()

#Pearson vs Spearman (BXB)
#Linear model
mod3<-lm(Spearman_R2_BXB~Pearson_R2_BXB, data=ERC_results_df_plus)
modsum3 = summary(mod3)

jpeg(file = paste0(working_dir, out_dir, "ERC_results/Correlation_strategies_BXB.jpg"), quality = 75)

plot(ERC_results_df_plus$Pearson_R2_BXB, ERC_results_df_plus$Spearman_R2_BXB,
     xlab ="Pearson R-squared",
     ylab ="Spearnam R-squared",
     col = rgb(red = 0, green = 0, blue = 0, alpha = 0.2),
     main = paste0("BXB R-squared values calculated from \nPearson vs Spearman correlation", "\nR2=", format(modsum3$adj.r.squared, digits = 3))
     )

#add trend line
abline(mod3)

dev.off()


#Print figures describing the distribution of R-squared values

# plot(density(na.omit(as.numeric(paste(ERC_results_df$Pearson_R2_BXB)))), col = "dark red")
# lines(density(na.omit(as.numeric(paste(ERC_results_df$Spearman_R2_BXB)))), col = "light red")
# lines(density(na.omit(as.numeric(paste(ERC_results_df$Pearson_R2_R2T)))), col = "dark blue")
# lines(density(na.omit(as.numeric(paste(ERC_results_df$Spearman_R2_R2T)))), col = "light blue")

#Note: these are not very compelling figures. Skipping for now. 
