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

#Get the filename from Network_analyses.py
fileName<-paste(args[8])
#Removes '.tsv' from fileName
fileName = substr(fileName, 1, nchar(fileName) - 4)
 

#Run functional category analysis
func_cat_bool<-paste(args[9])

#Label some nodes on the network?
lab_bool<-paste(args[10])


#####Testing
#jobname<-"ERC_Final"
#BL_type<-"R2T"
#RSquared<-0.5
#PValue<-0.0001
#clust_method<-"fg"
#trim_cutoff <- 0
#foc_sp<-"Atha"
#fileName<-"Filtered_ERC_results_R2T_3_0.0001_0.5"
#func_cat_bool<-"True"
#lab_bool<-"True"

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
ERC_hits_df<-read.table(file = paste0(working_dir, out_dir, "ERC_results/Filtered_results/", fileName,".tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

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
  	#LO <- layout_on_grid(network_graph)
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
    #LO_final <- layout_on_grid(network_graph_final)
    LO_final <- layout_nicely(network_graph_final)
    #legend colors
    legend_color <- rainbow(length(comms_final), alpha = 0.3)
    trim_cutoff = 0
}

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

#Filter out duplicated rows
final_network_stats_df_unique<-final_network_stats_df[!duplicated(final_network_stats_df$HOG_ID), ]

#reorder by community number and then by Degree Centrality and then by Eigenvector_centrality
final_network_stats_df_unique_reorder<-final_network_stats_df_unique[order(final_network_stats_df_unique$Community_num, -final_network_stats_df_unique$Degree_centrality, -final_network_stats_df_unique$Eigenvector_centrality), ]

#Sort the columns
final_network_stats_df_unique_reorder_rows<-final_network_stats_df_unique_reorder[, c("Community_num", names(final_network_stats_df_unique_reorder)[-length(names(final_network_stats_df_unique_reorder))])]

#Write the file
write.csv(final_network_stats_df_unique_reorder_rows, file = paste0(working_dir, out_dir, "Network_analyses/Network_stats_metrics_",fileName, "_", clust_method,"_trimcutoff_", trim_cutoff,".csv"), row.names = FALSE, quote = FALSE)

###END Create centrality metrics spreadsheet END###

### Add node label information
if(lab_bool=="True"){
  #Read in node label tsv
  node_labs_df<-read.table(file = paste0(working_dir, "Node_labels.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  #Merge the HOG info
  node_labs_HOGs<-merge(x = final_network_stats_df_unique_reorder_rows, y = node_labs_df, by = "ID", all.x = TRUE, all.y = FALSE)
  
  #Store as categorical data
  lab_values <- node_labs_HOGs$Label
  names(lab_values) <- node_labs_HOGs$HOG_ID
  
  #Add the attribute
  V(network_graph_final)$Node_labels <- lab_values[V(network_graph_final)$name]
  
  #Setup the spacing of labels on the graph
  stagger_directions<-sample(seq(1,360), size = length(V(network_graph_final)$Node_labels),replace = TRUE)
  #stagger_directions<-seq(1,360, 360/length(which(!is.na(V(network_graph_final)$Node_labels))))
  
}else{
  V(network_graph_final)$Node_labels<-NA
  stagger_directions<-sample(seq(1,360), size = length(V(network_graph_final)$Node_labels),replace = TRUE)
}

#save pdf
pdf(file = paste0(working_dir, out_dir, "Network_analyses/ERC_network_",fileName, "_", clust_method,"_trimcutoff_", trim_cutoff,".pdf"), width=8, height = 8)

#Plot the graph with all communities
plot(comms_final, network_graph_final,
            vertex.label=V(network_graph_final)$Node_labels,
            vertex.label.cex = 0.5,
            vertex.label.dist=0.5,
            vertex.label.degree = stagger_directions,
            vertex.size=2,
            edge.color=rgb(0,0,0,0.3),
            adjustcolor(V(network_graph_final)$color<-"black", alpha = .5),
            mark.col = comms_plot_col,
            mark.border = comms_plot_border,
            layout=LO_final,
            main=paste0(length(comms_final), " communities clustered using ", algo_name, " algorithm"))

mysubtitle<-paste0("BL method: ",BL_type, "  |  ", "Filter stats: ","R2 " ,RSquared ," P ", PValue , "  |  ")
mtext(side = 3, line = 0, at = 1, adj = 1, mysubtitle)

# legend("topleft",
#        legend = as.factor(1:length(comms_final)),
#        fill = legend_color

dev.off()



### Perform functional categories analysis
#Check to see whether this option was selected by user
if(func_cat_bool == "True"){
  
#Read in func cat files
func_cat_df<-read.table(file = paste0(working_dir, "Functional_categories.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

func_cat_col_df<-read.table(file = paste0(working_dir, "Functional_categories_col_assign.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

if(nrow(func_cat_df)>1){
  print("Functional category file successfully loaded in R")
}else{
  print("Something seems to be wrong with the functional category file")
}

#Merge the HOG info
Func_cats_HOGs<-merge(x = final_network_stats_df_unique_reorder_rows, y = func_cat_df, by = "ID", all.x = TRUE, all.y = FALSE)

#Reorder
Func_cats_HOGs_reorder<-Func_cats_HOGs[order(Func_cats_HOGs$Community_num, -Func_cats_HOGs$Degree_centrality, -Func_cats_HOGs$Eigenvector_centrality), ]

#Write a version of the csv file that includes func data
write.csv(Func_cats_HOGs_reorder, file = paste0(working_dir, out_dir, "Network_analyses/Network_stats_metrics_",fileName, "_", clust_method,"_trimcutoff_", trim_cutoff,"_FUNC_CAT.csv"), row.names = FALSE, quote = FALSE)

# #Assign functional category as a node attribute
#Store as categorical data
category_values <- Func_cats_HOGs_reorder$Functional_category
names(category_values) <- Func_cats_HOGs_reorder$HOG_ID

#Add the attribute
V(network_graph_final)$Functional_category <- category_values[V(network_graph_final)$name]

#Remove NAs
func_cat_col_df <- func_cat_col_df[complete.cases(func_cat_col_df), ]

#Generate color_mapping object
color_mapping <- setNames(func_cat_col_df$Color, func_cat_col_df$Category)

#Plot the network graph
#save pdf
pdf(file = paste0(working_dir, out_dir, "Network_analyses/ERC_network_",fileName, "_", clust_method,"_trimcutoff_", trim_cutoff,"_FUNC_CAT.pdf"), width=8, height = 8)

# Create a color vector based on Functional_category
vertex_colors <- color_mapping[V(network_graph_final)$Functional_category]

# Assign vertex colors
V(network_graph_final)$color <- vertex_colors

# Plot the network
plot(network_graph_final, 
     vertex.label=V(network_graph_final)$Node_labels,
     vertex.label.cex = 0.5,
     vertex.label.dist=0.5,
     vertex.label.degree = stagger_directions,
     vertex.size = 2,
     edge.color = rgb(0, 0, 0, 0.3),
     vertex.color = adjustcolor(V(network_graph_final)$color, alpha = 0.7),  # Set the alpha value for transparency
     layout = LO_final,
     main = paste0(length(comms_final), " communities clustered using ", algo_name, " algorithm"))

     mysubtitle<-paste0("BL method: ",BL_type, "  |  ", "Filter stats: ","R2 " ,RSquared ," P ", PValue , "  |  ")
     mtext(side = 3, line = 0, at = 1, adj = 1, mysubtitle)

dev.off()


#### ASSORTATIVITY ######
#Assign number code for the functional categories
#These codes are needed to run assortativity 
V(network_graph_final)$Functional_cat_code=V(network_graph_final)$Functional_category

#Get list of the different categories
cats<-levels(as.factor(V(network_graph_final)$Functional_cat_code))

#Assign all the codes
for(c in 1:length(cats)){
  V(network_graph_final)$Functional_cat_code[which(V(network_graph_final)$Functional_cat_code==cats[c])]<-c
  }

V(network_graph_final)$Functional_cat_code[which(is.na(V(network_graph_final)$Functional_cat_code))]<-(length(cats)+1)

#assortativity coefficient calculation (observed)
obs_assort<-assortativity.nominal(network_graph_final, types = V(network_graph_final)$Functional_cat_code, directed = FALSE)

#Get a random null distribution
#make function
rando_assort<-function(){
  #Make a copy
  network_graph_rep<-network_graph_final

  #coefficient calculation (randomized)
  assortativity.nominal(network_graph_rep, types = sample(V(network_graph_rep)$Functional_cat_code, replace = FALSE), directed = FALSE)
}

#run the function
assort_reps<-replicate(1000, rando_assort())

#Perform z-test to get p-value
#z-score
z_score<-(obs_assort-mean(assort_reps))/sd(assort_reps)


# One-sided p-value
p_val<-pnorm(q = z_score, lower.tail = FALSE)


#save pdf
pdf(file = paste0(working_dir, out_dir, "Network_analyses/Network_assortativity_",fileName, "_", clust_method,"_trimcutoff_", trim_cutoff,".pdf"), width=6, height = 6)

#plot the curve of the null distribution
plot(density(assort_reps), 
     main = "Random distribution and \nobserved assortativity", 
     xlab = paste0("Observed Assort coefficient: ", 
                   obs_assort, 
                   "\nOne-tailed P-value: ",
                   p_val
                   )
     )
#plot a line where the observed Assortativity is
abline(v=obs_assort)

dev.off()

#Write results to the metadata file
cat(paste0("\n", BL_type,"\t", PValue, "\t", 
           RSquared,"\t",length(unique(c(ERC_hits_df$GeneA_HOG, ERC_hits_df$GeneB_HOG))), "\t",
           nrow(ERC_hits_df),"\t",clust_method, "\t", 
           obs_assort,"\t",z_score,"\t",p_val), 
    file = paste0(working_dir, out_dir, "Network_analyses/Func_categories_rundata.tsv"), append = TRUE)

}#End if statement

### Create files for external use in Cytoscape
#Make a weird c-bound version of ERC results
all_df<-data.frame(HOG=c(ERC_hits_df$GeneA_HOG, ERC_hits_df$GeneB_HOG), ID=c(ERC_hits_df$GeneA_ID, ERC_hits_df$GeneB_ID))

#Get a table that contains gene annotation info.
cyto_df<-data.frame(HOG=unique(c(ERC_hits_df$GeneA_HOG, ERC_hits_df$GeneB_HOG)), ID_comprehensive=NA)

#add gene annotations to HOGs
for(i in 1:nrow(cyto_df)){
  #Find the first row that matches the HOG
  row_num_temp<-which(all_df$HOG == cyto_df$HOG[i])[1]
  #Check if the focal species string is present in ID
  if(length(grep(foc_sp, all_df$ID[row_num_temp]))>0){
    cyto_df$ID_comprehensive[i]<-all_df$ID[row_num_temp]
  }else{
    cyto_df$ID_comprehensive[i]<-all_df$HOG[row_num_temp]
  }
}

#Add the comprehensive ID as an attribute
#Store as categorical data
comp_values <- cyto_df$ID_comprehensive
names(comp_values) <- cyto_df$HOG

#Add the attribute
V(network_graph_final)$Comprehensive_ID <- comp_values[V(network_graph_final)$name]

#Write this graph as a graphML format file
write.graph(network_graph_final, file = paste0(working_dir, out_dir, "Network_analyses/Cytoscape_network_",fileName, "_", clust_method,"_trimcutoff_", trim_cutoff,".graphml"), format = "graphml")
