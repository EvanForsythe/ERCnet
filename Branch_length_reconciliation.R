#! /usr/local/bin/Rscript --vanilla --default-packages=utils

#USAGE

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
#jobname<-"TEST"
#jobname<-"Clptest"
out_dir<-paste0("OUT_", jobname, "/")

###Use DLCpar output to pull out relevant branch lengths
#Check if BL results dir exists (create one if not)
if(!dir.exists(paste0(working_dir, out_dir, "BL_results/"))){
  system(paste("mkdir ", paste0(working_dir, out_dir, "BL_results/")))
}

#Read species tree
sp_tr<-read.tree(paste0(working_dir, out_dir, "DLCpar/SpeciesTree_rooted_node_labels.txt"))

#Table of all edges on the species trees
all_nodes<-c(sp_tr$tip.label, sp_tr$node.label)
sp_branches_df<-as.data.frame(apply(X = sp_tr$edge, MARGIN = c(1,2), FUN = function(x){return(all_nodes[x])}))
names(sp_branches_df)<-c("ancestor", "decendant")
#Add a blank column to collect the BL measurement from the gene tree
sp_branches_df['gene_tree_BL_measure']<-NA

#Read species mapping file
mapping_df<-read.table(paste0(out_dir, "Species_mapping.csv"), sep = ",", header = TRUE)

#Get the list of reconciliation files
rec_files_list<-list.files(path = paste0(working_dir, out_dir, "DLCpar/"), pattern = "_NODES_BL.txt.dlcpar.locus.recon")

#Remove extra text
rec_files_list<-str_replace(rec_files_list, ".dlcpar.locus.recon", "")

#Make a table to store BL results for all trees
BL_measure_df<-as.data.frame(matrix(,ncol = (nrow(sp_branches_df)+1), nrow = length(rec_files_list)))
names(BL_measure_df)<-c("HOG_ID", paste0(sp_branches_df$ancestor, "_to_", sp_branches_df$decendant))

#Make a table to store root-2-tip (r2t) BL results
r2t_measure_df<-as.data.frame(matrix(,ncol=length(sp_tr$tip.label), nrow = length(rec_files_list)))
names(r2t_measure_df)<-sp_tr$tip.label
#Change the name of the outgroup (it should be the first) to be a column for storing HOG ID
names(r2t_measure_df)[1]<-"HOG_ID"

#Start loop
for(m in 1:length(rec_files_list)){
  #read in reconciliation file
  rec_df<-read.table(file = paste0(working_dir, out_dir, "DLCpar/", rec_files_list[m], ".dlcpar.locus.recon"), header = FALSE, stringsAsFactors = FALSE)
  names(rec_df)<-c("Gene_tree_node", "Sp_tree_location", "Event")

  #Note about rec file: the nodes in here refer to the nodes from the 'locus tree', which doesn't have branch lengths
  #There is a table that can be used to replace locus tree nodes with the relevant node from the gene tree

  #Read conversion table
  con_df<-read.table(file = paste0(working_dir, out_dir, "DLCpar/",  rec_files_list[m], ".dlcpar.coal.recon"), header = FALSE, stringsAsFactors = FALSE)
  names(con_df)<-c("Gene_tree_node_names", "Locus_tree_node_names", "Not_sure")
  #Replace the node names
  rec_df_replace<-rec_df
  rec_df_replace$Gene_tree_node<-NA

  #Note:need to figure out why some gene tree nodes are missing
  for(r in 1:nrow(rec_df)){
    if(length(which(con_df$Locus_tree_node_names==rec_df$Gene_tree_node[r]))>0){
      rec_df_replace$Gene_tree_node[r]<-con_df$Gene_tree_node_names[which(con_df$Locus_tree_node_names==rec_df$Gene_tree_node[r])]
    }}

  rec_df_replace<-rec_df_replace[complete.cases(rec_df_replace),]

  #Make a place to store the BL measurements for a given tree
  sp_branches_df_temp<-sp_branches_df

  #Read in gene tree with BLs
  temp_BL_tree<-read.tree(file = paste0(working_dir, out_dir, "DLCpar/", rec_files_list[m]))

  #make a matrix of branch lengths on gene tree
  BL_distance_mat<-dist.nodes(temp_BL_tree)

  #Make a list of all gene tree nodes (internal and external)
  all_gt_nodes<-c(temp_BL_tree$tip.label, temp_BL_tree$node.label)

  #Add the names of all nodes (internal and external)
  rownames(BL_distance_mat)<-all_gt_nodes
  colnames(BL_distance_mat)<-all_gt_nodes

  #Get the non-duplication nodes in the rec table
  spec_gene_rec<-subset(rec_df_replace, Event!="dup")

  #loop through all species branches and extract relevant distances
  #b<-16
  for(b in 1:nrow(sp_branches_df_temp)){
    #Check if both nodes exists in reconciliation file
    if(length(which(spec_gene_rec$Sp_tree_location==sp_branches_df_temp$ancestor[b]))>0
       &&
       length(which(spec_gene_rec$Sp_tree_location==sp_branches_df_temp$decendant[b]))>0){

      gt_anc_node<-subset(spec_gene_rec, Sp_tree_location==sp_branches_df_temp$ancestor[b])$Gene_tree_node
      gt_decend_node<-subset(spec_gene_rec, Sp_tree_location==sp_branches_df_temp$decendant[b])$Gene_tree_node

      store_wt_av<-list()
      #Get the average branch distance
      for(a in 1:length(gt_anc_node)){
        #Only retain decend nodes that decend from the particular anc node
        #Get the average from all the dec nodes that are actual dec
        store_wt_av[length(store_wt_av)+1]<-mean(BL_distance_mat[gt_anc_node[a],
                                                                 intersect(gt_decend_node, all_gt_nodes[getDescendants(temp_BL_tree, which(all_gt_nodes==gt_anc_node[a]))])
        ])
        sp_branches_df_temp$gene_tree_BL_measure[b]<-mean(na.omit(unlist(store_wt_av)))
      }#end anc node loop (varaible = a)

    }#End if statement
  }#End loop through species tree branches (variable b)

  #Change any BL=0 branches to NA
  sp_branches_df_temp$gene_tree_BL_measure[which(sp_branches_df_temp$gene_tree_BL_measure==0)]<-NA

  BL_measure_df[m,]<-c(str_replace(rec_files_list[m], "_NODES_BL.txt", ""),sp_branches_df_temp$gene_tree_BL_measure)
  
  #Check if there is an outlier branch
  #Get the branch values
  bxb_vals_temp<-as.numeric(paste(na.omit(as.numeric(BL_measure_df[m,-1]))))
  #Ask if the longest branch is more than 10x the second longest branch
  
  if(length(bxb_vals_temp)>0){
    if(bxb_vals_temp[order(bxb_vals_temp,decreasing = TRUE)][1]>(10*bxb_vals_temp[order(bxb_vals_temp,decreasing = TRUE)][2])){
      #If there is an outlier branch, set all values to NA for this tree
      BL_measure_df[m,-1]<-NA
  }}
  
  ### pull out root to tip distances
  #Find the node on the gene tree that maps to the node of interest (selected with --Node/-n in the Phylogenomics.py step) on the species tree
  #To figure out with node is used on this run, read in one of the earlier ouput files and pull out the N1/N2/etc.. string
  
  #Read second line of seq_counts file to pull of the node
  con <- file(paste0(out_dir, "Seq_counts_per_species.csv"),"r")
  second_line <- readLines(con, n = 2)[2]
  close(con)
  nodeID<-unlist(str_split(second_line, "\\."))[1]
  
  gt_root_node<-rec_df_replace$Gene_tree_node[which(rec_df_replace$Sp_tree_location==nodeID)]

  #Loop through all the branches that need r2t measures
  for(x in 2:ncol(r2t_measure_df)){
    #Record the r2t value for each species
    r2t_measure_df[m,x]<-mean(BL_distance_mat[gt_root_node, grep(mapping_df$Prefix[which(mapping_df$SpeciesID==names(r2t_measure_df)[x])], colnames(BL_distance_mat))])
  }#End r2t loop (variable = x)

  r2t_measure_df[m,1]<-str_replace(rec_files_list[m], "_NODES_BL.txt", "")

  #Check if there is an outlier branch
  #Get the branch values
  r2t_vals_temp<-as.numeric(paste(na.omit(as.numeric(r2t_measure_df[m,-1]))))
  #Ask if the longest branch is more than 10x the second longest branch
  if(length(r2t_vals_temp)>0){
    if(r2t_vals_temp[order(r2t_vals_temp,decreasing = TRUE)][1]>(10*r2t_vals_temp[order(r2t_vals_temp,decreasing = TRUE)][2])){
      #If there is an outlier branch, set all values to NA for this tree
      r2t_measure_df[m,-1]<-NA
  }}#end double if
  
}#End subtree loop (variable m)

### Write the results
#Write the branch by branch (bxb) results
write.table(BL_measure_df, file = paste0(working_dir, out_dir, "BL_results/bxb_BLs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

#Write the root to tip (r2t) results
write.table(r2t_measure_df, file = paste0(working_dir, out_dir, "BL_results/r2t_BLs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

#Normalize branch lengths by genome-wide average
#Branch by branch
bxb_BL_norm<-BL_measure_df

i<-2:ncol(bxb_BL_norm)
bxb_BL_norm[,i] <- apply(bxb_BL_norm[,i], 2,function(x) as.numeric(as.character(x)))

for(s in 2:ncol(bxb_BL_norm)){
  bxb_BL_norm[,s]<-bxb_BL_norm[,s]/(mean(na.omit(bxb_BL_norm[,s])))
}

#Branch by branch
r2t_BL_norm<-r2t_measure_df

i<-2:ncol(r2t_BL_norm)
r2t_BL_norm[,i] <- apply(r2t_BL_norm[,i], 2,function(x) as.numeric(as.character(x)))

for(s in 2:ncol(r2t_BL_norm)){
  r2t_BL_norm[,s]<-r2t_BL_norm[,s]/(mean(na.omit(r2t_BL_norm[,s])))
}

### Write the results
#Write the branch by branch (bxb) results
write.table(bxb_BL_norm, file = paste0(working_dir, out_dir, "BL_results/bxb_BLs_normalized.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

#Write the root to tip (r2t) results
write.table(r2t_BL_norm, file = paste0(working_dir, out_dir, "BL_results/r2t_BLs_normalized.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

