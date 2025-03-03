import os
import sys
import argparse
import pandas as pd
from collections import defaultdict
from Bio import Phylo
from io import StringIO

#At runtime set working directory to the place where the script lives
working_dir = sys.path[0]+'/'
os.chdir(working_dir)


#Set up an argumanet parser
parser = argparse.ArgumentParser(description='Script for summarizing gene duplication information')

#Add arguments
parser.add_argument('-j', '--JOBname', type=str, metavar='', required=True, help='Unique job name for this run of ERCnet. Avoid including spaces or special characters ("_" is ok)') 
parser.add_argument('-f', '--focal_sp', type=str, metavar='', required=True, help="The species you would like to use for obtaining gene IDs. Provide a string that exactly matches the header of the Filtered_genefam_dataset.csv file") 
parser.add_argument('-r', '--char_range', type=str, metavar='', required=True, help="The range of characters (0-indexed, inclusive) that you'd like to extract from the gene IDs. For example, for 'Atha__AT1G65870', provide '6-14' to extract 'AT1G65870'") 

#Define the parser
args = parser.parse_args()

#Store arguments
JOBname=args.JOBname
focal_sp=args.focal_sp
char_range=args.char_range

#Find the ERCnet dir
out_dir= 'OUT_'+JOBname+'/'

if os.path.isdir(out_dir):
    print("Found the ERCnet output directory!\n")
else:
    print(f"ERROR: ERCnet output, {out_dir}, directory not found! Quitting....\n")
    sys.exit()

#Find the reconciliation dir
rec_dir= out_dir+'DLCpar/'

if os.path.isdir(rec_dir):
    print("Found the reconciliation output directory!\n")
else:
    print(f"ERROR: reconciliation output, {rec_dir}, directory not found! Quitting....\n")
    sys.exit()

#Find the BL results from BXB
BL_table_path=out_dir+'BL_results/bxb_BLs_normalized.tsv'

if os.path.isfile(BL_table_path):
    print("Found the branch lengths file!\n")
else:
    print(f"ERROR: branch lengths file, {BL_table_path}, not found! Quitting....\n")
    sys.exit()


def find_dups(tsv_dir, header_file, output_file, gene_fam_file, focal_sp, char_range):
    print("Looping through HOG files...")

    # Read headers from the first line of "bxb_BLs_normalized.tsv"
    with open(header_file, 'r') as file:
        headers = file.readline().strip().split('\t')

    # Initialize the dataframe with headers and an additional HOG_ID column
    df = pd.DataFrame(columns=headers)

    # Read the gene families file
    gene_fams = pd.read_csv(gene_fam_file, sep=',')

    # Parse the character range
    start_idx, end_idx = map(int, char_range.split('-'))

    # Process each TSV file in the directory
    for filename in os.listdir(tsv_dir):
        if filename.endswith('.dlcpar.locus.recon'):
            filepath = os.path.join(tsv_dir, filename)

            # Read the TSV file into a DataFrame
            file_df = pd.read_csv(filepath, sep='\t', header=None, names=['GT_node', 'ST_node', 'Event'])

            # Initialize a dictionary to count dups for each header
            dup_counts = defaultdict(int)

            # Count "dup" occurrences based on mapping
            for _, row in file_df.iterrows():
                if row['Event'] == 'dup':
                    column_name = next(col for col in df.columns[1:] if str("_" + row['ST_node']) in col)
                    dup_counts[column_name] += 1

            # Prepare a row for the new dataframe
            row_data = {'HOG_ID': str(filename)}
            for header in headers:
                row_data[header] = dup_counts.get(header, 0)

            # Append the row to the dataframe
            df = pd.concat([df, pd.DataFrame([row_data])], ignore_index=True)

            # Extract the specific HOG_ID string
            hog_id = filename.split('_')[0]

            # Explicitly set the HOG_ID for the last row
            df.loc[df.index[-1], 'HOG_ID'] = hog_id

            # Find the corresponding HOG entry in gene_fams
            gene_fam_row = gene_fams[gene_fams['HOG'].str.endswith(hog_id)]
            if not gene_fam_row.empty:
                focal_entry = gene_fam_row[focal_sp].values[0]
                if pd.notna(focal_entry):
                    # Extract the substring based on char_range
                    seq_id = focal_entry[start_idx :end_idx +1]  # 0-indexed, inclusive
                else:
                    seq_id = "NA"
            else:
                seq_id = "NA"

            # Add seqID to the last row
            df.loc[df.index[-1], 'seqID'] = seq_id

    # Write the resulting dataframe to a TSV file
    print(f"Writing file: {output_file}")
    df.to_csv(output_file, sep='\t', index=False)
    return df

# Call the function
dups_counts_df = find_dups(
    tsv_dir=rec_dir,
    header_file=BL_table_path,
    output_file=out_dir + "Duplications_per_gene_with_seqID.tsv",
    gene_fam_file=out_dir + "Filtered_genefam_dataset.csv",
    focal_sp=focal_sp,
    char_range=char_range
)


## Get the sums from the columns
def calculate_column_totals(input_df, output_file):
    """
    Calculates column totals for a DataFrame (excluding a specified column) and writes the results to a TSV file.

    Parameters:
        input_df (pd.DataFrame): The input DataFrame.
        output_file (str): Path to the output TSV file.
    """
    # Exclude the 'HOG_ID' column
    columns_to_sum = [col for col in input_df.columns if col != 'HOG_ID']

    # Calculate totals for each column
    totals = input_df[columns_to_sum].sum()

    # Create a new DataFrame for the results
    result_df = pd.DataFrame({
        'Species_tree_branch': totals.index,
        'Total_duplications_mapped': totals.values
    })

    # Write the result to a TSV file
    print(f"Writing file: {output_file}")
    result_df.to_csv(output_file, sep='\t', index=False)
    return result_df

#Call the function
sums_df = calculate_column_totals(dups_counts_df, out_dir+"Duplications_per_branch.tsv")


# Next read in the species tree and add branch labels

def update_newick_tree(tree_file, df, output_file):
    """
    Updates a Newick tree by removing existing branch lengths and adding new labels based on a DataFrame.

    Parameters:
        tree_file (str): Path to the input Newick tree file.
        df (pd.DataFrame): DataFrame with branch label information.
        output_file (str): Path to the output Newick tree file.
    """
    # Read the Newick tree
    with open(tree_file, 'r') as file:
        tree_str = file.read().strip()

    # Parse the tree using Biopython
    tree = Phylo.read(StringIO(tree_str), 'newick')

    # Remove existing branch lengths
    for clade in tree.find_clades():
        clade.branch_length = None

    # Create a mapping from node names to duplication counts
    branch_label_mapping = {}
    for index, row in df.iterrows():
        column_name = row['Species_tree_branch']
        if '_to_' in column_name:
            _, target_node = column_name.split('_to_')
            branch_label_mapping[target_node] = int(row['Total_duplications_mapped'])

    # Assign new branch labels based on the DataFrame
    for clade in tree.find_clades():
        if clade.name in branch_label_mapping:
            clade.branch_length = branch_label_mapping[clade.name]

    print(f"writing {output_file}")

    # Write the updated tree to a new file
    with open(output_file, 'w') as file:
        Phylo.write(tree, file, 'newick')

#Read in the species tree
SP_tree_path=out_dir+'SpeciesTree_rooted_node_labels.txt'

if os.path.isfile(SP_tree_path):
    print("Found the species tree file!\n")
else:
    print(f"ERROR: species tree file, {SP_tree_path}, not found! Quitting....\n")
    sys.exit()


#Call the function 
update_newick_tree(SP_tree_path, sums_df, out_dir+'SP_tree_with_duplicates.nwk')


