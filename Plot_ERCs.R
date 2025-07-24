#Script to plot pairwise ERC correlation plots for individual ERC comparisons
#Usage:
#Rscript Plot_ERCs.R <jobname> <BL_type: "r2t" or "bxb"> <GeneA HOG> <GeneB HOG> <Common name for GeneA> <Common name for Gene B>
#Note, the user is responsible for making sure the HOG IDs and common names match up correctly!
#Additional note: the user will need to make sure the ggplot2 library is installed

library("ggplot2")

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

#Get the arguments
jobname<-args[1]
bl_type<-args[2]
hogA<-args[3]
hogB<-args[4]
nameA<-args[5]
nameB<-args[6]
print(jobname)

#Setup file paths
bl_path<-paste0("OUT_", jobname, "/BL_results/", bl_type,"_BLs_normalized.tsv")
out_dir<-paste0("OUT_", jobname, "/Pairwise_ERC_plots/")

#Create the folder if it doesn't already exist
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

#Read in filtered ERC result dataset
bl_df<-read.csv(
  file = bl_path, 
  sep = "\t", header = TRUE)

search_hogs<-c(hogA,hogB)
gene_names<-c(nameA, nameB)

bl_df[which(bl_df$HOG_ID==search_hogs[1]),]
bl_df[which(bl_df$HOG_ID==search_hogs[2]),]

test_df<-data.frame(Branch=names(bl_df),
                    GeneA=paste(bl_df[which(bl_df$HOG_ID==search_hogs[1]),]),
                    GeneB=paste(bl_df[which(bl_df$HOG_ID==search_hogs[2]),]))

test_df<-test_df[-1,]

test_df[,1]<-as.character(test_df[,1])
test_df[,2]<-as.numeric(test_df[,2])
test_df[,3]<-as.numeric(test_df[,3])

all_possible_branches<-test_df$Branch

test_df<-test_df[complete.cases(test_df),]

# Define a custom vector of shapes and repeat it to match the number of unique Branch levels
base_shapes <- c(0, 1, 2)  # Square, Circle, Triangle
shape_values <- rep(base_shapes, length.out = length(all_possible_branches))

# Ensure all_possible_branches is a factor with the correct levels
test_df$Branch <- factor(test_df$Branch, levels = all_possible_branches)

#Fit the linear model
lm_model <- lm(GeneB ~ GeneA, data = test_df)
lm_summary <- summary(lm_model)
r_squared <- lm_summary$r.squared


if(length(shape_values) > 20){
    plot_width=6
}else{
    plot_width=4
}

pdf(file = paste0(out_dir, "ERC_plot_", gene_names[1], "_vs_", gene_names[2], ".pdf"), width = plot_width, height = 4)

# Create the plot
ggplot(test_df, aes(y = GeneB, x = GeneA, color = Branch, shape = Branch)) +
  geom_point(size = 3) +  # Increase the point size
  geom_abline(intercept = coef(lm_model)[1], 
              slope = coef(lm_model)[2], 
              color = "black", linetype = "dashed") +  # Add the regression line
  xlab(gene_names[1]) +
  ylab(gene_names[2]) +
  labs(color = "Branch", shape = "Branch") +
  scale_shape_manual(values = shape_values) +  # Manually specify shapes with repetition
  annotate("text", x = Inf, y = Inf, 
           label = sprintf("R-squared = %.2f", r_squared), 
           hjust = 1.1, vjust = 1.1, size = 3, color = "black") +
  theme(
    plot.title = element_text(size = 8, face = "italic"),
    legend.text = element_text(size = 6, face = "italic"),
    legend.key.size = unit(0.3, "cm")
  )

dev.off()




