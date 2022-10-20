#! /usr/local/bin/Rscript --vanilla --default-packages=utils

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

#Which file to run stats on
stats_file<-args[2]
#stats_file<-"hits"

#read in designated ERC results file
if(stats_file=="full"){
  ERC_results_df<-read.table(paste0(working_dir, out_dir, "ERC_results/ERC_results.tsv"), sep = "\t", header = TRUE)
}else if(stats_file=="hits")
  ERC_results_df<-read.table(paste0(working_dir, out_dir, "ERC_results/ERC_results_potential_hits.tsv"), sep = "\t", header = TRUE)

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
mod0<-lm(Pearson_R2_R2T~Pearson_R2_BXB, data=ERC_results_df)
modsum0 = summary(mod0)

#r2t vs bxb (Pearson)
jpeg(file = paste0(working_dir, out_dir, "ERC_results/Branch_length_strategies_pearson.jpg"), quality = 75)

plot(ERC_results_df$Pearson_R2_BXB, ERC_results_df$Pearson_R2_R2T,
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
mod1<-lm(Spearman_R2_R2T~Spearman_R2_BXB, data=ERC_results_df)
modsum1 = summary(mod1)

jpeg(file = paste0(working_dir, out_dir, "ERC_results/Branch_length_strategies_spearman.jpg"), quality = 75)

plot(ERC_results_df$Spearman_R2_BXB, ERC_results_df$Spearman_R2_R2T,
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
mod2<-lm(Spearman_R2_R2T~Pearson_R2_R2T, data=ERC_results_df)
modsum2 = summary(mod2)

jpeg(file = paste0(working_dir, out_dir, "ERC_results/Correlation_strategies_R2T.jpg"), quality = 75)

plot(ERC_results_dfs$Pearson_R2_R2T, ERC_results_df$Spearman_R2_R2T,
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
mod3<-lm(Spearman_R2_BXB~Pearson_R2_BXB, data=ERC_results_df)
modsum3 = summary(mod3)

jpeg(file = paste0(working_dir, out_dir, "ERC_results/Correlation_strategies_BXB.jpg"), quality = 75)

plot(ERC_results_df$Pearson_R2_BXB, ERC_results_df$Spearman_R2_BXB,
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
