

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib as mpl

# Keep text editable in Illustrator
mpl.rcParams['pdf.fonttype'] = 42      # TrueType (keeps text as text)
mpl.rcParams['ps.fonttype']  = 42      # If you ever save PS/EPS
mpl.rcParams['text.usetex'] = False    # Make sure LaTeX isn't turning text into Type 3
mpl.rcParams['font.family']  = 'Arial' # Or another common font installed on the machine opening the PDF


# Will be set after goi_df is loaded
goi_id_to_name = None


def get_pair_names(row):
    # Find which gene IDs from goi_genes are present in GeneA_ID and GeneB_ID
    geneA_names = [goi_id_to_name[gid] for gid in goi_id_to_name if gid in str(row['GeneA_ID'])]
    geneB_names = [goi_id_to_name[gid] for gid in goi_id_to_name if gid in str(row['GeneB_ID'])]
    # Use first match for each, fallback to empty string
    geneA = geneA_names[0] if geneA_names else ''
    geneB = geneB_names[0] if geneB_names else ''
    return f"{geneA}, {geneB}" if geneA or geneB else ''

# File paths

import argparse
parser = argparse.ArgumentParser(description='KDE plot with extracted rows for ERC results.')
parser.add_argument('--pos_control', '-p', required=True, help='Full path to the positive control TSV file, formatted like the gene of interest file, with columns Gene_name and Gene_ID')
parser.add_argument('--job_name', '-j', required=True, help='Jobname of ERCnet run')
parser.add_argument('--goi', '-g', required=True, help='Full path to the gene of interest TSV file with columns Gene_name and Gene_ID. Gene_name should be a common name and Gene_ID should be a string found in the sequence ID of the focal species')
parser.add_argument('--branch_method', '-b', required=True, choices=['BXB', 'R2T'], help='Branch method, must be either BXB or R2T')
parser.add_argument('--test', action='store_true', help='If set, only analyze the first 5 gene of interest and positive control gene pairs')
parser.add_argument('--erc_results_folder', '-e', required=True, help='Full path to the folder containing ERC results files (ERC_results_R2T.tsv or ERC_results_BXB.tsv)')
# New arguments for R and P value cutoffs
parser.add_argument('-pp', '--PearsonP', type=float, required=False, help='P-value cutoff for Pearson correlation. Default is 0.05', default=0.05)
parser.add_argument('-pr', '--PearsonR', type=float, required=False, help='R-value cutoff for Pearson correlation. Default is 0.3', default=0.3)
parser.add_argument('-sp', '--SpearmanP', type=float, required=False, help='P-value cutoff for Spearman correlation. Default is 0.05', default=0.05)
parser.add_argument('-sr', '--SpearmanR', type=float, required=False, help='R-value cutoff for Spearman correlation. Default is 0.3', default=0.3)
parser.add_argument('-kp', '--KendallP', type=float, required=False, help='P-value cutoff for Kendall correlation. Default is 0.05', default=0.05)
parser.add_argument('-kr', '--KendallR', type=float, required=False, help='R-value cutoff for Kendall correlation. Default is 0.3', default=0.3)
args = parser.parse_args()

# Ensure output directory exists
os.makedirs(args.erc_results_folder, exist_ok=True)

# Validate job_name is in the provided ERC results folder path
if args.job_name not in args.erc_results_folder:
    raise ValueError(f"Job name '{args.job_name}' not found in the ERC results folder path: {args.erc_results_folder}")

if args.branch_method == 'R2T':
    ERC_results_path = os.path.join(args.erc_results_folder, 'ERC_results_R2T.tsv')
elif args.branch_method == 'BXB':

    ERC_results_path = os.path.join(args.erc_results_folder, 'ERC_results_BXB.tsv')
else:
    raise ValueError('Invalid branch_method argument. Use R2T or BXB.')

print(f"Reading ERC results from {ERC_results_path} ...")
ERC_results_df = pd.read_csv(ERC_results_path, sep='\t')
print(f"ERC results loaded: {len(ERC_results_df)} rows.")

print(f"Reading positive control gene list from {args.pos_control} ...")
pos_control_df = pd.read_csv(args.pos_control, sep='\t')
pos_control_genes = pos_control_df['Gene_ID'].dropna().unique()
print(f"Loaded {len(pos_control_genes)} unique positive control genes.")

print(f"Reading genes of interest from {args.goi} ...")
goi_df = pd.read_csv(args.goi, sep='\t')
# Keep only one row per unique Gene_name
goi_df_unique = goi_df.drop_duplicates(subset='Gene_name')
goi_genes = goi_df_unique['Gene_ID'].dropna().unique()  # Use unique Gene_IDs for unique Gene_names
goi_id_to_name = dict(zip(goi_df['Gene_ID'].astype(str), goi_df['Gene_name'].astype(str)))

from itertools import combinations
pos_control_pairs = list(combinations(pos_control_genes, 2))
goi_pairs = list(combinations(goi_genes, 2))
print(f"Number of positive control pairs: {len(pos_control_pairs)}")
print(f"Number of gene of interest pairs: {len(goi_pairs)}")

# If test mode, only use first 5 pairs
if args.test:
    print('Test mode: Only analyzing first 5 gene of interest and positive control gene pairs.')
    pos_control_pairs = pos_control_pairs[:5]
    goi_pairs = goi_pairs[:5]


# Helper function for robust substring matching in ERC file columns
def gene_in_field(gene, field):
    # field may be a string with comma-separated IDs
    return any(gene in part.strip() for part in str(field).split(','))


print("Extracting positive control pairs from ERC results ...")
pos_control_rows = pd.DataFrame()
for idx, (gene_a, gene_b) in enumerate(pos_control_pairs):
    if idx % 10 == 0:
        print(f"  Processing positive control gene pair {idx+1}/{len(pos_control_pairs)} ...")
    mask1 = ERC_results_df['GeneA_ID'].apply(lambda x: gene_in_field(gene_a, x)) & ERC_results_df['GeneB_ID'].apply(lambda x: gene_in_field(gene_b, x))
    mask2 = ERC_results_df['GeneA_ID'].apply(lambda x: gene_in_field(gene_b, x)) & ERC_results_df['GeneB_ID'].apply(lambda x: gene_in_field(gene_a, x))
    matched_rows = ERC_results_df[mask1 | mask2]
    pos_control_rows = pd.concat([pos_control_rows, matched_rows], ignore_index=True)
pos_control_rows = pos_control_rows.drop_duplicates()
print(f"Extracted {len(pos_control_rows)} positive control rows.")

print("Extracting gene of interest pairs from ERC results ...")
goi_extracted_rows = pd.DataFrame()
for idx, (gene_a, gene_b) in enumerate(goi_pairs):
    if idx % 10 == 0:
        print(f"  Processing gene of interest pair {idx+1}/{len(goi_pairs)} ...")
    mask1 = ERC_results_df['GeneA_ID'].apply(lambda x: gene_in_field(gene_a, x)) & ERC_results_df['GeneB_ID'].apply(lambda x: gene_in_field(gene_b, x))
    mask2 = ERC_results_df['GeneA_ID'].apply(lambda x: gene_in_field(gene_b, x)) & ERC_results_df['GeneB_ID'].apply(lambda x: gene_in_field(gene_a, x))
    matched_rows = ERC_results_df[mask1 | mask2]
    goi_extracted_rows = pd.concat([goi_extracted_rows, matched_rows], ignore_index=True)
goi_extracted_rows = goi_extracted_rows.drop_duplicates()
print(f"Extracted {len(goi_extracted_rows)} gene of interest rows.")


print("Saving transitional file with all rows of interest (GOI and pos control) ...")

pos_control_rows_copy = pos_control_rows.copy()
pos_control_rows_copy['gene_set'] = 'positive_control'
goi_rows_copy = goi_extracted_rows.copy()
goi_rows_copy['gene_set'] = 'gene_of_interest'
rows_of_interest = pd.concat([pos_control_rows_copy, goi_rows_copy], ignore_index=True)
rows_of_interest = rows_of_interest.drop_duplicates()
rows_filename = f"rows_of_interest_{args.branch_method}.tsv"
rows_filename = os.path.join(args.erc_results_folder, f"Genes_of_interest_and_pos_cont_{args.branch_method}.tsv")
rows_of_interest.to_csv(rows_filename, sep='\t', index=False)
print(f'Saved transitional file: {rows_filename} ({len(rows_of_interest)} rows)')

# Three-panel stacked KDE plot for Pearson_R, Spearman_R, Kendall_Tau, min, and avg
import numpy as np
metrics = [
    ("Pearson_R", "Pearson_R KDE"),
    ("Spearman_R", "Spearman_R KDE"),
    ("Kendall_Tau", "Kendall_Tau KDE")
]

print("Plotting three-panel stacked KDE plot for R values ...")
fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(8, 14), sharex=False)
ERC_results_df_plot = ERC_results_df.sample(n=min(100000, len(ERC_results_df)), random_state=42) if len(ERC_results_df) > 100000 else ERC_results_df
for i, (col, label) in enumerate(metrics):
    ax = axes[i]
    # KDE for full dataset (downsampled)
    sns.kdeplot(ERC_results_df_plot[col], fill=True, color='skyblue', label='Full dataset', ax=ax)
    # KDE for positive control genes
    if not pos_control_rows.empty:
        sns.kdeplot(pos_control_rows[col], fill=True, color='green', label='Positive control genes', ax=ax, alpha=0.5)
    # KDE for gene of interest pairs
    if not goi_extracted_rows.empty:
        sns.kdeplot(goi_extracted_rows[col], fill=True, color='orange', label='Gene of interest pairs', ax=ax, alpha=0.5)
        # Determine R-value threshold for this metric
        if col == 'Pearson_R':
            r_cutoff = args.PearsonR
        elif col == 'Spearman_R':
            r_cutoff = args.SpearmanR
        elif col == 'Kendall_Tau':
            r_cutoff = args.KendallR
        else:
            r_cutoff = 0.3
        # Add vertical lines only for rows passing R threshold
        filtered_rows = goi_extracted_rows[goi_extracted_rows[col] >= r_cutoff]
        y_max = ax.get_ylim()[1]
        for idx, (row_idx, row) in enumerate(filtered_rows.iterrows()):
            value = row[col]
            ax.axvline(value, color='orange', linestyle='-', alpha=0.7, label='Gene of interest pairs' if idx == 0 else "")
            label = get_pair_names(row)
            # Place label at the top of the plot, slightly above y_max
            ax.text(value, y_max*0.5, label, rotation=90, va='bottom', ha='center', fontsize=8, color='orange', clip_on=True)
    ax.set_xlabel(col)
    ax.set_ylabel('Density')
    # No panel title
    # Add more x-axis tick marks
    ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=15))
    handles, labels_ = ax.get_legend_handles_labels()
    by_label = dict(zip(labels_, handles))
    ax.legend(by_label.values(), by_label.keys())
plt.tight_layout()
output_file = os.path.join(args.erc_results_folder, f"Genes_of_interest_and_pos_cont_KDEplot_RVAL_{args.branch_method}.pdf")
plt.savefig(output_file)
print(f"Three-panel KDE plot (Rval) saved as {output_file}")
plt.show()

# Second three-panel figure: Pval values (log10 transformed)
import numpy as np
pval_metrics = [
    ("Pearson_Pval", "Pearson_Pval KDE"),
    ("Spearman_Pval", "Spearman_Pval KDE"),
    ("Kendall_Pval", "Kendall_Pval KDE")
]
# Set minimum value for Pearson_Pval and Spearman_Pval to 1e-10
min_pval = 1e-10
for df in [ERC_results_df, pos_control_rows, goi_extracted_rows]:
    for col in ["Pearson_Pval", "Spearman_Pval"]:
        if col in df.columns:
            df[col] = df[col].clip(lower=min_pval)
# Avoid log(0) by replacing 0 with a very small value for Kendall_Pval
eps = 1e-300
for col in ["Pearson_Pval", "Spearman_Pval", "Kendall_Pval"]:
    ERC_results_df[col + '_log10'] = -np.log10(ERC_results_df[col].replace(0, eps))
    pos_control_rows[col + '_log10'] = -np.log10(pos_control_rows[col].replace(0, eps))
    goi_extracted_rows[col + '_log10'] = -np.log10(goi_extracted_rows[col].replace(0, eps))




pval_metrics_log = [
    ("Pearson_Pval_log10", "-log10(Pearson_Pval) KDE"),
    ("Spearman_Pval_log10", "-log10(Spearman_Pval) KDE"),
    ("Kendall_Pval_log10", "-log10(Kendall_Pval) KDE")
]

print("Plotting three-panel stacked KDE plot for P values (log10 transformed) ...")
fig2, axes2 = plt.subplots(nrows=3, ncols=1, figsize=(8, 14), sharex=False)
# Downsample full dataset for plotting
ERC_results_df_plot = ERC_results_df.sample(n=min(100000, len(ERC_results_df)), random_state=42) if len(ERC_results_df) > 100000 else ERC_results_df
for i, (col, label) in enumerate(pval_metrics_log):
    ax = axes2[i]
    # KDE for full dataset (downsampled)
    sns.kdeplot(ERC_results_df_plot[col], fill=True, color='skyblue', label='Full dataset', ax=ax)
    # KDE for positive control genes
    if not pos_control_rows.empty:
        sns.kdeplot(pos_control_rows[col], fill=True, color='green', label='Positive control genes', ax=ax, alpha=0.5)
    # KDE for gene of interest pairs
    if not goi_extracted_rows.empty:
        sns.kdeplot(goi_extracted_rows[col], fill=True, color='orange', label='Gene of interest pairs', ax=ax, alpha=0.5)
        # Determine P-value threshold for this metric
        if col == 'Pearson_Pval_log10':
            p_cutoff = args.PearsonP
            orig_pval_col = 'Pearson_Pval'
        elif col == 'Spearman_Pval_log10':
            p_cutoff = args.SpearmanP
            orig_pval_col = 'Spearman_Pval'
        elif col == 'Kendall_Pval_log10':
            p_cutoff = args.KendallP
            orig_pval_col = 'Kendall_Pval'
        else:
            p_cutoff = 0.05
            orig_pval_col = None
        # Add vertical lines only for rows passing P threshold
        if orig_pval_col is not None:
            filtered_rows = goi_extracted_rows[goi_extracted_rows[orig_pval_col] <= p_cutoff]
            y_max = ax.get_ylim()[1]
            for idx, (row_idx, row) in enumerate(filtered_rows.iterrows()):
                value = row[col]
                ax.axvline(value, color='orange', linestyle='-', alpha=0.7, label='Gene of interest pairs' if idx == 0 else "")
                label_names = get_pair_names(row)
                ax.text(value, y_max*0.5, label_names, rotation=90, va='bottom', ha='center', fontsize=8, color='orange', clip_on=True)
    ax.set_xlabel(label)
    ax.set_ylabel('Density')
    # No panel title
    # Add integer x-axis tick marks
    ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=15, integer=True))
    handles, labels_ = ax.get_legend_handles_labels()
    by_label = dict(zip(labels_, handles))
    ax.legend(by_label.values(), by_label.keys())
plt.tight_layout()
output_file2 = os.path.join(args.erc_results_folder, f"Genes_of_interest_and_pos_cont_KDEplot_PVAL_{args.branch_method}.pdf")
plt.savefig(output_file2)
print(f"Three-panel KDE plot (Pval, log10) saved as {output_file2}")
plt.show()
