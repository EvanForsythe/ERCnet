import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser(description='Map gene names to ERC results and extract genes of interest.')
print('Starting ERCnet genes of interest extraction...')
parser.add_argument('--job_name', '-j', required=True, help='Jobname of ERCnet run')
parser.add_argument('--goi', '-g', required=True, help='Full path to the gene of interest TSV file with columns Gene_name and Gene_ID. Gene_name should be a common name and Gene_ID should be a string found in the sequence ID of the focal species')
parser.add_argument('--erc_hits', '-r', required=True, help='Full path to the folder containing Filtered_ERC_hits.tsv. The file must be named Filtered_ERC_hits.tsv.')
args = parser.parse_args()

# Validate job_name is in the provided folder path
if args.job_name not in args.erc_hits:
	raise ValueError(f"Job name '{args.job_name}' not found in the ERC hits folder path: {args.erc_hits}")

# Read gene map file
print('Reading gene map file...')
gene_map = pd.read_csv(args.goi, sep='\t')

# Read ERC results file
print('Reading ERC hits file...')
erc_hits_file = os.path.join(args.erc_hits, 'Filtered_ERC_hits.tsv')
erc_df = pd.read_csv(erc_hits_file, sep='\t')

# Example: print first few rows of each
gene_map.columns = [col.strip() for col in gene_map.columns]
erc_df.columns = [col.strip() for col in erc_df.columns]

print('Gene map file:')
print(gene_map.head())
print('\nERC results file:')
print(erc_df.head())


# Filter erc_df to only rows containing any Gene_ID from gene_map in GeneA_ID or GeneB_ID

gene_ids = gene_map['Gene_ID'].astype(str).unique()
print('Filtering ERC hits for genes of interest...')
mask = (
	erc_df['GeneA_ID'].astype(str).apply(lambda x: any(gid in x for gid in gene_ids)) |
	erc_df['GeneB_ID'].astype(str).apply(lambda x: any(gid in x for gid in gene_ids))
)
filtered_erc_df = erc_df[mask].copy()


# Add gene_of_interestA and gene_of_interestB columns
def find_gene_name(gene_id_str, gene_map):
	for _, row in gene_map.iterrows():
		if row['Gene_ID'] in gene_id_str:
			return row['Gene_name']
	return ''

filtered_erc_df['gene_of_interestA'] = filtered_erc_df['GeneA_ID'].astype(str).apply(lambda x: find_gene_name(x, gene_map))
filtered_erc_df['gene_of_interestB'] = filtered_erc_df['GeneB_ID'].astype(str).apply(lambda x: find_gene_name(x, gene_map))

# Move new columns to the front
cols = ['gene_of_interestA', 'gene_of_interestB'] + [col for col in filtered_erc_df.columns if col not in ['gene_of_interestA', 'gene_of_interestB']]
filtered_erc_df = filtered_erc_df[cols]

# Write to TSV
output_file = 'ERC_hits_for_genes_of_interest.tsv'

# Prepare output directory

# Output will be written to the same folder as Filtered_ERC_hits.tsv
output_path = os.path.join(args.erc_hits, output_file)
print(f'Writing filtered results to output file in: {output_path}')

# Write output file in the output directory
filtered_erc_df.to_csv(output_path, sep='\t', index=False)
print(f"\nFiltered ERC results with gene names written to {output_path}")
