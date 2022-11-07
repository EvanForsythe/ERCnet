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
RSquared<-paste(args[3])
PValue<-paste(args[4])

#Get the method to be used
clust_method<-args[5]

#User defined trim cutoff (i.e. what minimum community size is displayed)
trim_cutoff <- as.numeric(paste(args[6]))

#Get the ID of the focal species
foc_sp<-paste(args[7])

#####Testing
#jobname<-"quick_test"
#BL_type<-"r2t"
#RSquared<-0.5
#PValue<-0.05
#clust_method<-"fg"
#trim_cutoff<-0
#foc_sp<-"A_thaliana_prot"
#working_dir<-"/Users/esforsythe/Documents/Work/Bioinformatics/ERC_networks/Analysis/ERCnet_dev/"
###

# Set the script path:
working_dir<- paste0(getScriptPath(),"/") #added the "/" at the end so paste0 commands below work.

setwd(working_dir)

#Because the above function returns "." as the working dir, using this command to set the full path (to avoid issues below)
working_dir<-paste0(getwd(),"/")

out_dir<-paste0("OUT_", jobname, "/")

#Read in ERC correlation results
#Read the table
ERC_hits_df<-read.table(file = paste0(working_dir, out_dir, "ERC_results/Filtered_ERC_Results.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#Print message
#Use cat because paste+print doesn't recognize \n
cat("\nNumber of significant correlations (according to filter parameters): ", nrow(ERC_hits_df), 
    "\nNumber of genes in network: ", length(unique(c(ERC_hits_df$GeneA_HOG, ERC_hits_df$GeneB_HOG))),
    "\n\n")


print("Printing networks...")

#make a version of the data frame that contains only the HOGs
ERC_hits_df_4network<-ERC_hits_df[,grep("HOG", names(ERC_hits_df))]

#Build network
network_graph<-graph.data.frame(ERC_hits_df_4network, directed = FALSE)

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


#TRIM OR NO TRIM (if statement)

validCutoff = FALSE

if (trim_cutoff > 0){
  
  if (trim_cutoff < max(sizes(comms))) {
	validCutoff = TRUE
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
	} else {
		print("Trim Cutoff would filter all communities. Retaining all communities instead.")
	}
  
} 
if (!validCutoff | trim_cutoff <= 0){
    network_graph_final <- network_graph
    comms_final <- comms
    comms_plot_col <- rainbow(length(comms_final), alpha = 0.3)
    comms_plot_border <- rainbow(length(comms_final), alpha = 1)
    #Network layout
    LO_final <- layout_nicely(network_graph_final)
    #legend colors
    legend_color <- rainbow(length(comms_final), alpha = 0.3)
    trim_cutoff = 0
}


#save pdf
pdf(file = paste0(working_dir, out_dir, "Network_analyses/ERC_network_",BL_type, "_", "R2_", RSquared ,"_Pv_", PValue, "_", clust_method,"_trimcutoff_", trim_cutoff,".pdf"), width=8, height = 8)

#Plot the graph with all communities
plot(comms_final, network_graph_final,
            vertex.label=NA,
            vertex.size=2,
            edge.color=rgb(0,0,0,0.3),
            adjustcolor(V(network_graph_final)$color<-"black", alpha = .5),
            mark.col = comms_plot_col,
            mark.border = comms_plot_border,
            layout=LO_final,
            main=paste0(length(comms_final), " communities clustered using ", algo_name, " algorithm"))

mysubtitle<-paste0("BL method: ",BL_type, "  |  ", "Filter stats: ","R2 " ,RSquared ," P ", PValue , "  |  ")
mtext(side = 3, line = 0, at = 1, adj = 1, mysubtitle)

legend("topleft",
       legend = as.factor(1:length(comms_final)),
       fill = legend_color
)


dev.off()


## Get the genes in each community
#Create df
comms_df<-data.frame(HOG=comms_final$names, Community=comms_final$membership)

#Reorder
comms_df<-comms_df[order(comms_df$Community),]

#reformat to create table of HOGs and IDs
A_genes<-ERC_hits_df[,1:2]
B_genes<-ERC_hits_df[,3:4]
names(A_genes)<-c("HOG", "ID")
names(B_genes)<-c("HOG", "ID")
ID_df<-rbind(A_genes, B_genes)


#Join the dataframes
comms_w_IDs<-merge(x = comms_df, y = ID_df, by="HOG", all.x=TRUE)

#Take the first sequence ID when mulitple are present (this is akin to choosing a random representative)
comms_w_IDs$ID<-sapply(strsplit(comms_w_IDs$ID, ","),`[`, 1)


#Write txt files
#All genes in the network (as a background for enrichment analyses)
write.table(paste(na.omit(unique(comms_w_IDs$ID))),
           file = paste0(working_dir, out_dir, "Network_analyses/Communities/Comm_","BACKGROUND_", BL_type, "_", "R2_", RSquared ,"_Pv_", PValue, "_", clust_method,".txt"),
           sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#Loop through all communities
for(n in min(comms_w_IDs$Community):max(comms_w_IDs$Community)){

  #Write tsv files
  write.table(paste(na.omit(unique(subset(comms_w_IDs$ID, comms_w_IDs$Community==n)))), 
              file = paste0(working_dir, out_dir, "Network_analyses/Communities/Comm_",sprintf("%04d", n), "_", BL_type, "_", "R2_", RSquared ,"_Pv_", PValue, "_", clust_method,".txt"), 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

}
  
###Create centrality metrics spreadsheet###

#Get all HOG names pre-filter
all_HOG_names_df <- data.frame(HOG_ID=unique(c(ERC_hits_df$GeneA_HOG,ERC_hits_df$GeneB_HOG)))

#make HOG column name match in foc_sp_df
names(comms_w_IDs)[names(comms_w_IDs) == 'HOG'] <- 'HOG_ID'

#only keep relevant columns
subset_focal_sp_df <- subset(comms_w_IDs, select = c("HOG_ID", "ID"))

#Get degree info (Authority and hub scores omitted because they are the same as eigenvector in undirected graphs)
degree_cen<-degree(network_graph_final, normalized=TRUE)

#Make dataframe
network_stats_df <- data.frame(HOG_ID=names(degree_cen),
          Eigenvector_centrality=as.numeric(paste(evcent(network_graph_final)$vector)),
          Degree_centrality=as.numeric(paste(degree_cen)),
          Eccentricity_centrality=as.numeric(paste(eccentricity(network_graph_final))),
          Betweenness_centrality=as.numeric(paste(betweenness(network_graph_final))),
          Closeness_centrality=as.numeric(paste(closeness(network_graph_final))),
          Community_num=as.numeric(paste(membership(comms_final)))
)

#merge HOGS not present in graph
all_HOGs_and_focalsp_df <- merge(all_HOG_names_df,subset_focal_sp_df,by="HOG_ID",all=TRUE, sort = FALSE)

#merge
final_network_stats_df <- merge(all_HOGs_and_focalsp_df,network_stats_df,by="HOG_ID",all=TRUE, sort = FALSE)

write.csv(final_network_stats_df, file = paste0(working_dir, out_dir, "Network_analyses/Network_stats_metrics_",BL_type, "_", "R2_", RSquared ,"_Pv_", PValue, "_", clust_method,"_trimcutoff_", trim_cutoff,".csv"))


