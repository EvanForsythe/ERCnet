import argparse
import os
import pandas as pd
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description="Process a job, read a tree, and load a TSV file")
    parser.add_argument("jobname", help="The name of the job")
    parser.add_argument("filename", help="The name of the TSV file")
    parser.add_argument("hog_id1", help="The first HOG ID")
    parser.add_argument("hog_id2", help="The second HOG ID")
    parser.add_argument("sim_tree_path", help="Path to simulations tree")
    return parser.parse_args()

def call_nw_display(tree_path):
    # Call nw_display to print the tree in ASCII format
    try:
        print(f"\nDisplaying phylogenetic tree for: {tree_path}")
        subprocess.run(["nw_display", tree_path], check=True)
    except Exception as e:
        print(f"Error calling nw_display: {e}")

def display_tree(jobname=None, sim_tree_path=None):
    # Determine which tree to display
    if jobname:
        tree_path = os.path.join(f"OUT_{jobname}", "SpeciesTree_rooted_node_labels.txt")
    elif sim_tree_path:
        tree_path = sim_tree_path
    else:
        print("Error: No tree path provided.")
        return
    
    # Call nw_display to display the tree in ASCII format
    call_nw_display(tree_path)

def check_tree_properties(gene_tree_path):
    # Load the tree using Biopython
    from Bio import Phylo
    try:
        tree = Phylo.read(gene_tree_path, "newick")
        # Check if the tree is rooted
        is_rooted = tree.rooted
        print(f"Is the tree rooted? {'Yes' if is_rooted else 'No'}")
        
        # Count the number of multifurcating nodes
        multifurcating_nodes = sum(len(clade.clades) > 2 for clade in tree.get_nonterminals())
        print(f"Number of multifurcating nodes: {multifurcating_nodes}")
        
    except Exception as e:
        print(f"Error loading or checking tree properties: {e}")

def load_tsv(jobname, filename, hog_id1, hog_id2):
    # Construct the TSV file path
    tsv_path = os.path.join(f"OUT_{jobname}", "BL_results", filename)
    
    # Read the TSV file into a dataframe
    try:
        df = pd.read_csv(tsv_path, sep="\t")
        
        # Filter the rows based on HOG_ID1 and HOG_ID2
        df_filtered = df[df['HOG_ID'].isin([hog_id1, hog_id2])]
        
        if df_filtered.empty:
            print(f"No matching rows found for HOG_ID: {hog_id1} or {hog_id2}")
            return None
        
        # Set the HOG_ID as the index for easier plotting
        df_filtered.set_index('HOG_ID', inplace=True)
        
        # Print column values for hog_id1 and hog_id2
        print(f"\nColumn values for HOG_ID {hog_id1}:")
        print(df_filtered.loc[hog_id1].to_string())
        
        print(f"\nColumn values for HOG_ID {hog_id2}:")
        print(df_filtered.loc[hog_id2].to_string())
        
        return df_filtered
    except Exception as e:
        print(f"Error reading TSV file: {e}")
        return None

def correlation_plot(df, hog_id1, hog_id2, output_path):
    if df is not None and len(df.columns) > 1:
        # Transpose the dataframe to plot HOG_IDs as rows
        df_t = df.T
        
        # Remove rows with NaN or Inf values for hog_id1 or hog_id2
        df_t_clean = df_t[[hog_id1, hog_id2]].replace([np.inf, -np.inf], np.nan).dropna()

        # Check if there is enough data left after removing NaNs
        if df_t_clean.empty:
            print(f"Not enough valid data to calculate correlation between {hog_id1} and {hog_id2}")
            return
        
        # Print the table used for calculating the correlation
        print("\nTable used for calculating the correlation:")
        print(df_t_clean)

        # Calculate Pearson and Spearman correlation coefficients
        pearson_corr, pearson_pval = pearsonr(df_t_clean[hog_id1], df_t_clean[hog_id2])
        spearman_corr, spearman_pval = spearmanr(df_t_clean[hog_id1], df_t_clean[hog_id2])

        # Generate a scatter plot to show correlation between two rows
        plt.figure(figsize=(10, 6))
        sns.regplot(x=df_t_clean[hog_id1], y=df_t_clean[hog_id2])

        # Annotate each point with the index (which corresponds to the name of the datapoint)
        for i in range(df_t_clean.shape[0]):
            plt.text(df_t_clean[hog_id1].iloc[i], df_t_clean[hog_id2].iloc[i], 
                     df_t_clean.index[i], fontsize=9, ha='right')

        plt.title(f'Linear Correlation: {hog_id1} vs {hog_id2}')
        plt.xlabel(hog_id1)
        plt.ylabel(hog_id2)

        # Add Pearson and Spearman correlations as text on the plot
        plt.text(0.05, 0.95, f"Pearson r: {pearson_corr:.2f} (p={pearson_pval:.2e})", 
                 transform=plt.gca().transAxes, fontsize=10, verticalalignment='top')
        plt.text(0.05, 0.90, f"Spearman ρ: {spearman_corr:.2f} (p={spearman_pval:.2e})", 
                 transform=plt.gca().transAxes, fontsize=10, verticalalignment='top')

        # Save the plot as a PDF
        plt.savefig(output_path)
        print(f"Correlation plot saved as: {output_path}")
        plt.close()
    else:
        print(f"Not enough data available for correlation plot between {hog_id1} and {hog_id2}")

def main():
    # Parse arguments
    args = parse_args()
    
    # Load and print the phylogenetic tree for jobname using nw_display
    display_tree(jobname=args.jobname)

    # Load and print the simulated tree using nw_display
    display_tree(sim_tree_path=args.sim_tree_path)

    # Check if the trees for HOG_ID1 and HOG_ID2 are rooted and count multifurcating nodes
    gene_tree_path1 = os.path.join(f"OUT_{args.jobname}", "BL_trees", f"RAxML_result.{args.hog_id1}_BL.txt")
    check_tree_properties(gene_tree_path1)

    gene_tree_path2 = os.path.join(f"OUT_{args.jobname}", "BL_trees", f"RAxML_result.{args.hog_id2}_BL.txt")
    check_tree_properties(gene_tree_path2)
    
    # Load the TSV file and filter by HOG_ID1 and HOG_ID2
    df = load_tsv(args.jobname, args.filename, args.hog_id1, args.hog_id2)
    
    # Define the output path for the correlation plot PDF
    output_path = os.path.join(f"OUT_{args.jobname}", "sim_test_plot.pdf")
    
    # Plot the correlation between HOG_ID1 and HOG_ID2 and save as PDF
    correlation_plot(df, args.hog_id1, args.hog_id2, output_path)

if __name__ == "__main__":
    main()
