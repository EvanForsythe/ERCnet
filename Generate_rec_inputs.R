#! /usr/local/bin/Rscript --vanilla --default-packages=utils

#USAGE
#Rscript Generate_rec_inputs.R /Users/esforsythe/Documents/Work/Bioinformatics/ERC_networks/Analysis/Orthofinder/Results_Oct15/


#Load packages
package_list<-c("ape", "stringr", "phytools")

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
#jobname<-"TPC_test"
out_dir<-paste0("OUT_", jobname, "/")

#Get dir where orthofinder results live
OFpath<-args[2]
#OFpath<-"/Users/esforsythe/Documents/Work/Bioinformatics/ERC_networks/Analysis/Orthofinder/Results_Oct15/"
#OFpath<-"/Users/esforsythe/Documents/Work/Bioinformatics/ERC_networks/Analysis/Orthofinder/Plant_cell/Results_Feb15/"


## Reconciliation
#Check if reconciliation dir exists (create one if not)
if(!dir.exists(paste0(working_dir, out_dir, "DLCpar/"))){
  system(paste("mkdir ", paste0(working_dir, out_dir, "DLCpar/")))
}

#get a copy of the species tree from Orthofinder folder
file.copy(from = paste0(OFpath,"Species_Tree/SpeciesTree_rooted_node_labels.txt"),
          to = paste0(working_dir, out_dir, "DLCpar/SpeciesTree_rooted_node_labels.txt"))

#Read species mapping file
mapping_df<-read.table(paste0(out_dir, "Species_mapping.csv"), sep = ",", header = TRUE)

#Create the DCLpar files needed
#Make species map table
write.table(data.frame(gt=paste0(mapping_df$Prefix, "*"), st=paste0(mapping_df$SpeciesID)),
            file = paste0(working_dir, out_dir, "DLCpar/speciesIDs.smap"), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

#Get list of trees to be reconciled
tree_input_list<-list.files(path = paste0(working_dir, out_dir, "BL_trees/"), pattern = "RAxML_result.")

#Remove extra text in file names
tree_input_list<-str_replace(str_replace(tree_input_list, "RAxML_result.", ""), "_BL.txt", "")

#Loop through input files
for(d in 1:length(tree_input_list)){
  #d<-1
  #Get the BL tree
  BL_tree_temp<-read.tree(file = paste0(working_dir, out_dir, "BL_trees/RAxML_result.", tree_input_list[d], "_BL.txt")) 
  #Get the original subtree
  subtree_temp2<-read.tree(file = paste0(working_dir, out_dir, "HOG_subtrees/", tree_input_list[d], "_tree.txt")) 
  
  #Make new version of subtree with all branches=1
  subtree_temp2$edge.length<-rep(1,length(subtree_temp2$edge.length))
  
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
  
  #Root the BL tree
  BL_tree_root<-root(BL_tree_temp, outgroup = outtaxa, resolve.root = TRUE)
  
  #Add node labels to the tree
  BL_tree_root$node.label<-paste0("n", 1:BL_tree_root$Nnode)
  
  #Write the BL_tree (this is what is used for DLCpar reconciliation)
  write.tree(phy = BL_tree_root, file = paste0(working_dir, out_dir, "DLCpar/", tree_input_list[d], "_NODES_BL.txt"))
  
  # ##Run DLCpar
  # DLCpar_cmd<-paste0("dlcpar_search -s SpeciesTree_rooted_node_labels.txt -S speciesIDs.smap ", tree_input_list[d], "_NODES_BL.txt")
  # system(DLCpar_cmd)
  # 
  # if((d %% 100) == 0){
  #   print(paste0(d, " trees reconciled"))
  # }
  
}#End DLCpar input generator file loop (variable = d)





