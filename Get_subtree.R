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

#Read in the csv file of filtered gene gene families
out_dir<-paste(args[2])
#out_dir<-"OUT_TPC_TEST/"
out_dir_full<-paste0(working_dir, out_dir)

#Read in file
HOGs_df<-read.csv(paste0(out_dir_full, "Filtered_genefam_dataset.csv"))

#Loop through all the rows in the file
for(h in 1:nrow(HOGs_df)){

#Read in tree of OG
full_tree_temp<-read.tree(file = paste0(OG_file_path, HOGs_df$OG[h], "_tree.txt"))

# #Development
# keeper_tips<-unlist(strsplit(
#   "Atha_AT5G57750,Bole_Bol009773,Brap_I04829,Bstr_26833s0480,Cgra_2848s0056,Crub_0008s1826,Esal_10014578m,Spar_Sp6g19430",
#   ","))

#Prune the tree to contain only keeper tips
#Get argument and split into a vector
#keeper_tips<-unlist(strsplit(args[3], ','))

#Get list of tips
raw_strings<-paste(HOGs_df[h,4:ncol(HOGs_df)])
if(length(which(raw_strings==""))>0){
  raw_strings<-raw_strings[-which(raw_strings=="")]
}

if(length(which(raw_strings=="NA"))>0){
  raw_strings<-raw_strings[-which(raw_strings=="NA")]
}

keeper_tips<-unlist(strsplit(str_replace_all(toString(raw_strings), " ", ""), split = ","))

#Get the subtree
subtree_temp<-keep.tip(phy = full_tree_temp, tip =c(keeper_tips))

#Get the name of the outfile
out_file_name<-paste0(HOGs_df$HOG[h])

#Remove the N1 (or N2, N3, etc...) string
out_file_name<-unlist(str_split(out_file_name, "\\."))[2]

#Add suffix
out_file_name<-paste0(out_file_name, "_tree.txt")

#Check if subtrees has polytomies
if(is.binary(subtree_temp)){
  #write the subtree
  write.tree(subtree_temp, file = paste0(out_dir_full,"HOG_subtrees/", out_file_name))
}else{
  cat(paste0(out_file_name, '\n'),
    file = paste0(out_dir_full, "Non-binary_subtrees.txt"), append = TRUE)
  }

}#end for loop
