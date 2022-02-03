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
script_path<- getScriptPath()
setwd(script_path)

#Load packages
package_list<-c("ape")

#Loop to check if package is installed and libraried
for(p in 1:length(package_list)){
  if (!require(package_list[p], character.only = TRUE)) {
    install.packages(package_list[p], dependencies = TRUE)
    library(package_list[p], character.only=TRUE)
  }
}

#Read in arguments
args = commandArgs(trailingOnly=TRUE)

#Get OG tree path
OG_file_path<-paste(args[1])

#Read in tree of OG
full_tree_temp<-read.tree(file = OG_file_path)

# #Development
# keeper_tips<-unlist(strsplit(
#   "Atha_AT5G57750,Bole_Bol009773,Brap_I04829,Bstr_26833s0480,Cgra_2848s0056,Crub_0008s1826,Esal_10014578m,Spar_Sp6g19430",
#   ","))

#Prune the tree to contain only keeper tips
#Get argument and split into a vector
keeper_tips<-unlist(strsplit(args[3], ','))
#Get the subtree
subtree_temp<-keep.tip(phy = full_tree_temp, tip = keeper_tips)

#Get the name of the outfile
out_file<-args[2]

#Check if subtrees has polytomies
if(is.binary(subtree_temp)){
  #write the subtree
  write.tree(subtree_temp, file = out_file)
}else{
  cat(paste0(out_file, '\n'),
    file = "Non-binary_subtrees.txt", append = TRUE)
  }

