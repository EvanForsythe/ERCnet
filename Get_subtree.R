#! /usr/local/bin/Rscript --vanilla --default-packages=utils

#Load packages
library("ape")
library("stringr")

#Set the working directory with the path to this script
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

#Get OG tree path
OG_file_path<-paste(args[1])
#OG_file_path<-"/Users/esforsythe/Documents/Work/Bioinformatics/ERC_networks/Analysis/Orthofinder/Results_Oct15/Resolved_Gene_Trees/"
#OG_file_path<-"/Users/esforsythe/Documents/Work/Bioinformatics/ERC_networks/Analysis/Orthofinder/Plant_cell/Results_Feb15/Resolved_Gene_Trees/"

#Read in the csv file of filtered gene gene families
out_dir<-paste(args[2])
#out_dir<-"OUT_Clptest/"
out_dir_full<-paste0(working_dir, out_dir)

#Read in file
HOGs_df<-read.csv(paste0(out_dir_full, "Filtered_genefam_dataset.csv"))

#Loop through all the rows in the file
for(h in 1:nrow(HOGs_df)){
#Get the name of the outfile
HOG_temp<-paste0(HOGs_df$HOG[h])
  
#Remove the N1 (or N2, N3, etc...) string
HOG_temp<-unlist(str_split(HOG_temp, "\\."))[2]
  
#Read in tree of OG
full_tree_temp<-read.tree(file = paste0(OG_file_path, HOGs_df$OG[h], "_tree.txt"))

#Prune the tree to contain only keeper tips
#Get argument and split into a vector

#Get list of tips
raw_strings<-paste(HOGs_df[h,4:ncol(HOGs_df)])
if(length(which(raw_strings==""))>0){
  raw_strings<-raw_strings[-which(raw_strings=="")]
}

if(length(which(raw_strings=="NA"))>0){
  raw_strings<-raw_strings[-which(raw_strings=="NA")]
}

keeper_tips<-unlist(strsplit(str_replace_all(toString(raw_strings), " ", ""), split = ","))

#Check if this was a pruned alignment
if(length(list.files(path = paste0(out_dir_full, "Aln_pruning/"), pattern = paste0(HOG_temp, "_ALN_RETAINED.txt")))==1){
  prune_seqs<-read.table(file = paste0(out_dir_full, "Aln_pruning/", list.files(path = paste0(out_dir_full, "Aln_pruning/"), pattern = paste0(HOG_temp, "_ALN_RETAINED.txt"))))$V1
  
  #loop through seqs that need to be pruned from subtree
  for(p in 1:length(prune_seqs)){
    keeper_tips<-keeper_tips[-which(keeper_tips==prune_seqs[p])]
    }#End for loop
  }#End prune if statement

#Get the subtree
subtree_temp<-keep.tip(phy = full_tree_temp, tip =c(keeper_tips))

#Add suffix
out_file_name<-paste0(HOG_temp, "_tree.txt")

#Check if subtrees has polytomies
if(is.binary(subtree_temp)){
  #write the subtree
  write.tree(subtree_temp, file = paste0(out_dir_full,"HOG_subtrees/", out_file_name))
}else{
  cat(paste0(out_file_name, '\n'),
    file = paste0(out_dir_full, "Non-binary_subtrees.txt"), append = TRUE)
  }

}#end for loop
