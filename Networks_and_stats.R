#! /usr/local/bin/Rscript --vanilla --default-packages=utils
#Load packages
library("igraph")

getScriptPath <- function(){
  cmd.args <- commandArgs()
  m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
  script.dir <- dirname(regmatches(cmd.args, m))
  if(length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
  if(length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
  return(script.dir)
}


### Read in all user-defined variables here

#Read in arguments
args = commandArgs(trailingOnly=TRUE)

#Get the job name (used to identify the proper output folder)
jobname<-args[1]

#Get parameters for filtering data to be used in networks
BL_type<-paste(args[2])
filter_stat<-paste(args[3])
filter_stat_cutoff<-as.numeric(paste(args[4]))

#Get the method to be used
clust_method<-args[5]

#User defined trim cutoff (i.e. what minimum community size is displayed)
trim_cutoff <- as.numeric(paste(args[6]))

#Get the ID of the focal species
foc_sp<-paste(args[7])

###



# Set the script path:
working_dir<- paste0(getScriptPath(),"/") #added the "/" at the end so paste0 commands below work.

setwd(working_dir)

#Because the above function returns "." as the working dir, using this command to set the full path (to avoid issues below)
working_dir<-paste0(getwd(),"/")

out_dir<-paste0("OUT_", jobname, "/")

#Read in ERC correlation results
#Read the table
ERC_results_df<-read.table(file = paste0(working_dir, out_dir, "ERC_results/ERC_results.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)


#Filter data to be included in ERC network
if(BL_type == "bxb"){
  if(filter_stat=="pval"){
    ERC_hits_df<-subset(ERC_results_df, Pearson_P_BXB<filter_stat_cutoff)[,c(1,2)]
  }else if(filter_stat=="R2"){
    ERC_hits_df<-subset(ERC_results_df, R2_BXB>=filter_stat_cutoff)[,c(1,2)]
  }
}else if(BL_type == "r2t"){
  if(filter_stat=="pval"){
    ERC_hits_df<-subset(ERC_results_df, Pearson_P_R2T<filter_stat_cutoff)[,c(1,2)]
  }else if(filter_stat=="R2"){
    ERC_hits_df<-subset(ERC_results_df, R2_R2T>=filter_stat_cutoff)[,c(1,2)]
  }
}

#Print message
#Use cat because paste+print doesn't recognize \n
cat("\nNumber of significant correlations (according to filter parameters): ", nrow(ERC_hits_df), 
    "\nNumber of genes in network: ", length(unique(c(ERC_hits_df$GeneA_HOG, ERC_hits_df$GeneB_HOG))),
    "\n\n")

print("Printing networks...")

#Build network
network_graph<-graph.data.frame(ERC_hits_df, directed = FALSE)



#CLUSTER COMMUNITIES

#Cluster
if(clust_method == "fg"){
  comms<-cluster_fast_greedy(network_graph)
  algo_name <- "fast and greedy"
}else if(clust_method == "eb"){
  comms<-cluster_edge_betweenness(network_graph)
  algo_name <- "edge betweenness"
}else if(clust_method == "op"){
  comms<-cluster_optimal(network_graph)
  algo_name <- "optimal"
}else if(clust_method == "wt"){
  comms<-cluster_walktrap(network_graph)
  algo_name <- "walktrap"
}else{
  comms<-cluster_fast_greedy(network_graph)
  algo_name <- "fast and greedy"
}


#TRIM OR NO TRIM

if (trim_cutoff > 0){
 
  comms_keep_ids <- as.numeric(names(sizes(comms)[sizes(comms) >= trim_cutoff]))
  comms_keep_v_idxs <- which(comms$membership %in% comms_keep_ids)
  
  network_graph_final <- induced_subgraph(network_graph, V(network_graph)[comms_keep_v_idxs])
  
  # subset community objects
  comms_final <- comms
  comms_final$names <- comms$names[comms_keep_v_idxs]
  comms_final$membership <- comms$membership[comms_keep_v_idxs]
  comms_final$vcount <- length(comms_final$names)
  comms_final$modularity <- modularity(network_graph_final, comms_final$membership, E(network_graph_final)$weight)
  
  #Network layout and colors
  LO <- layout_nicely(network_graph)
  LO_final <- LO[comms_keep_v_idxs, ]
  comms_plot_col <- rainbow(length(communities(comms_final)), alpha = 0.3)[comms_keep_ids]
  comms_plot_border <- rainbow(length(communities(comms_final)), alpha = 1)[comms_keep_ids]
  
  #legend colors
   legend_color <- rainbow(length(comms_final), alpha = 0.3)[comms_keep_ids]
  
}else{
    network_graph_final <- network_graph
    comms_final <- comms
    comms_plot_col <- rainbow(length(comms_final), alpha = 0.3)
    comms_plot_border <- rainbow(length(comms_final), alpha = 1)
    #Network layout
    LO_final <- layout_nicely(network_graph_final)
    #legend colors
    legend_color <- rainbow(length(comms_final), alpha = 0.3)
}



#save pdf
pdf(file = paste0(working_dir, out_dir, "Network_analyses/ERC_network_",BL_type, "_", filter_stat, "_", filter_stat_cutoff, "_", clust_method,".pdf"), width=8, height = 8)

#Plot the graph with all communities
plot(comms_final, network_graph_final,
            vertex.label=NA,
            vertex.size=5,
            edge.color="black",
            V(network_graph)$color<-"black",
            mark.col = comms_plot_col,
            mark.border = comms_plot_border,
            layout=LO_final,
            main=paste0(length(comms_final), " communities clustered using ", algo_name, " algorithm"))

mysubtitle<-paste0("BL method: ",BL_type, "  |  ", "Filter stat: ", filter_stat, "  |  ", "Cutoff: ", filter_stat_cutoff)
mtext(side = 3, line = 0, at = 1, adj = 1, mysubtitle)

legend("topleft",
       legend = as.factor(1:length(comms_final)),
        fill = legend_color
)


dev.off()

#

## Get the genes in each community
#Create df
comms_df<-data.frame(HOG=comms_final$names, Community=comms_final$membership)

#Reorder
comms_df<-comms_df[order(comms_df$Community),]

#Read in the HOG file
All_HOGs_df<-read.table(file = paste0(working_dir, out_dir, "Filtered_genefam_dataset.csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)

#Remove the "N1." (or N2. N3. etc...) from the string
All_HOGs_df$HOG<-sapply(strsplit(as.character(All_HOGs_df$HOG), "\\."), `[`, 2)

#Join the dataframes
comms_w_IDs<-merge(x = comms_df, y = All_HOGs_df, by="HOG", all.x=TRUE)

#Get a dataframe with the focal sp IDs
foc_sp_df<-cbind(comms_w_IDs[,c(1,2,3)],
data.frame(
Focal_sp_ID=sapply(strsplit(comms_w_IDs[,which(names(comms_w_IDs) == foc_sp)], ","), `[`, 1)
))

#Remove prefix
foc_sp_df$Focal_sp_ID<-gsub(foc_sp_df$Focal_sp_ID, pattern = paste0(foc_sp, "_"), replacement = "")

#Write txt files
#All genes in the network (as a background for enrichment analyses)
write.table(paste(na.omit(unique(foc_sp_df$Focal_sp_ID))), 
            file = paste0(working_dir, out_dir, "Network_analyses/Communities/Comm_","BACKGROUND_", BL_type, "_", filter_stat, "_", filter_stat_cutoff, "_", clust_method,".txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#Loop through all communities
for(n in min(foc_sp_df$Community):max(foc_sp_df$Community)){

  #Write tsv files
  write.table(paste(na.omit(unique(subset(foc_sp_df$Focal_sp_ID, foc_sp_df$Community==n)))), 
              file = paste0(working_dir, out_dir, "Network_analyses/Communities/Comm_",sprintf("%04d", n), "_", BL_type, "_", filter_stat, "_", filter_stat_cutoff, "_", clust_method,".txt"), 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

} 
  
  


