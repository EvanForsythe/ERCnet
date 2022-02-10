#! /usr/local/bin/Rscript --vanilla --default-packages=utils
#Load packages
package_list<-c("igraph")

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
#jobname<-"BIGTEST"
out_dir<-paste0("OUT_", jobname, "/")

#Read in ERC correlation results
#Read the table
ERC_results_df<-read.table(file = paste0(working_dir, out_dir, "ERC_results/ERC_results.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#Get parameters for filtering data to be used in networks
BL_type<-paste(args[2])
filter_stat<-paste(args[3])
filter_stat_cutoff<-as.numeric(paste(args[4]))

###Development***
##BLtype
#BL_type<-"bxb"
#BL_type<-"r2t"
#
##Stats
#filter_stat<-"pval"
#filter_stat_cutoff<-0.001
#OR
#filter_stat<-"R2"
#filter_stat_cutoff<-0.8

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

#degree(network_graph)
#Cluster communities
#Get the method to be used
clust_method<-args[5]
#clust_method<-"fg"
#clust_method<-"eb"

#Cluster
if(clust_method == "fg"){
  comms<-cluster_fast_greedy(network_graph)
}else if(clust_method == "eb"){
  comms<-cluster_edge_betweenness(network_graph)
}

#dendPlot(comms, mode = "hclust")
#membership(comms)
#length(comms)

#save pdf
pdf(file = paste0(working_dir, out_dir, "Network_analyses/ERC_network_",BL_type, "_", filter_stat, "_", filter_stat_cutoff, "_", clust_method,".pdf"), width=8, height = 8)

#plot the graph
plot.igraph(network_graph, 
            vertex.label=NA,
            vertex.size=5,
            edge.color="black",
            vertex.color=membership(comms), 
            layout=layout.fruchterman.reingold,
            vertex.label.color="black",
            main=paste0("BL method: ",BL_type,
                        "\nFilter stat: ", filter_stat,
                        "\nCutoff: ", filter_stat_cutoff,
                        "\n", length(comms), " communities clustered by: ", clust_method))

dev.off()

## Get the genes in each community
#Create df
comms_df<-data.frame(HOG=comms$names, Community=comms$membership)

#Reorder
comms_df<-comms_df[order(comms_df$Community),]

#Read in the HOG file
All_HOGs_df<-read.table(file = paste0(working_dir, out_dir, "Filtered_genefam_dataset.csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE)

#Remove the "N1."
All_HOGs_df$HOG<-gsub(All_HOGs_df$HOG, pattern = "N1.", replacement = "")

#Join the dataframes
comms_w_IDs<-merge(x = comms_df, y = All_HOGs_df, by="HOG", all.x=TRUE)

#Get the ID of the focal species
foc_sp<-paste(args[6])
#foc_sp<-"Atha"

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
  
  
