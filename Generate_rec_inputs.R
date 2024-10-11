#! /usr/local/bin/Rscript --vanilla --default-packages=utils

#USAGE
#Rscript Generate_rec_inputs.R test 

#Load packages
library("ape")
library("stringr")
library("phytools")

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

#Get the job name (used to identify the proper output folder)
jobname<-args[1]
#jobname<-"test"
#jobname<-"TPC_test"
out_dir<-paste0("OUT_", jobname, "/")

#Get list of trees to be reconciled
tree_input_list<-list.files(path = paste0(working_dir, out_dir, "BL_trees/"), pattern = ".treefile")
#tree_input_list<-str_replace(str_replace(tree_input_list, "RAxML_result.", ""), "_BL.txt", "")
tree_input_list<-str_replace(tree_input_list, "_BL.treefile", "")

print(length(tree_input_list))

#Loop through input files
for(d in 1:length(tree_input_list)){
  #d<-1
  #Get the BL tree
  BL_tree_temp<-read.tree(file = paste0(working_dir, out_dir, "BL_trees/",tree_input_list[d], "_BL.treefile")) 

  #print("Printing BL_tree")
  #newick_string <- write.tree(BL_tree_temp)
  #system(paste("echo", shQuote(newick_string), "| nw_display -"))


  #Get the reconciled/rooted tree
  subtree_temp2<-read.tree(file = paste0(working_dir, out_dir, "Rearranged_trees/", tree_input_list[d], "_BS.treefile_recs.nwk"))
  
  #Make new version of subtree with all branches=1
  subtree_temp2$edge.length<-rep(1,length(subtree_temp2$edge.length))

  #print("Printing Rearranged_trees/")
  #newick_string2 <- write.tree(subtree_temp2)
  #system(paste("echo", shQuote(newick_string2), "| nw_display -"))
  
  #Split into rooting
  subtrees_4rooting<-treeSlice(subtree_temp2, 0.01, trivial=FALSE, prompt=FALSE)

  if(length(subtrees_4rooting)==1){
    outtaxa<-subtree_temp2$tip.label[!subtree_temp2$tip.label %in% subtrees_4rooting[[1]]$tip.label]
  }else if(length(subtrees_4rooting)==2){

    t_one_count<-subtrees_4rooting[[1]]$tip.label
    t_two_count<-subtrees_4rooting[[2]]$tip.label

    if(length(t_one_count)<length(t_two_count)){
      outtaxa<-t_one_count
    }else{
      outtaxa<-t_two_count
    }
  } #End subtrees if statement

  #print("\nOutgroup: ")
  #print(outtaxa)

  #Root the BL tree
  BL_tree_root<-root(BL_tree_temp, outgroup = outtaxa, resolve.root = TRUE)

  #Add node labels to the tree
  BL_tree_root$node.label<-paste0("n", 1:BL_tree_root$Nnode)
  
  #Write the BL_tree (this is what is used for DLCpar reconciliation)
  write.tree(phy = BL_tree_root, file = paste0(working_dir, out_dir, "DLCpar/", tree_input_list[d], "_NODES_BL.txt"))
  
} #End DLCpar input generator file loop (variable = d)

