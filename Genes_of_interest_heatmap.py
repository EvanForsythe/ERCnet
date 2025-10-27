# --- imports ---
import argparse
import pandas as pd
import itertools
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib as mpl

# Keep text editable in Illustrator
mpl.rcParams['pdf.fonttype'] = 42      # TrueType (keeps text as text)
mpl.rcParams['ps.fonttype']  = 42      # If you ever save PS/EPS
mpl.rcParams['text.usetex'] = False    # Make sure LaTeX isn't turning text into Type 3



# --- argparse ---

parser = argparse.ArgumentParser(description='Generate ERC heatmap for genes of interest.')
parser.add_argument('--job_name', '-j', required=True, help='Jobname of ERCnet run')
parser.add_argument('--branch_method', '-b', required=True, choices=['BXB', 'R2T'], help='Branch method, must be either BXB or R2T')
parser.add_argument('--HOGs_of_interest', '-g', required=True, help='Full path to a TSV file with the columns "HOG" and "Gene_name" listing HOGs of interest')
parser.add_argument('--test', action='store_true', help='If set, only use the first 5 HOGs_of_interest for testing')
args = parser.parse_args()

# --- output folder setup ---
erc_results_folder = os.path.join("/scratch/forsythe/ERCnet", f"OUT_{args.job_name}", "ERC_results")
output_tsv = os.path.join(erc_results_folder, "Genes_of_interest_heatmap_data.tsv")
output_pdf = os.path.join(erc_results_folder, "Genes_of_interest_heatmap.pdf")
print(f"Saving output files to: {erc_results_folder}")

# --- read HOGs file ---
print(f"Reading HOGs of interest from: {args.HOGs_of_interest}")

df = pd.read_csv(args.HOGs_of_interest, sep="\t")
if args.test:
    print("--test flag is set: using only the first 5 HOGs_of_interest for testing.")
    df = df.iloc[:5]
hog_to_gene = dict(zip(df["HOG"], df["Gene_name"]))
filtered_hogs = df["HOG"].unique().tolist()
print(f"Number of HOGs of interest: {len(filtered_hogs)}")

# --- create HOG pairs ---
hog_pairs = set(frozenset(pair) for pair in itertools.combinations(filtered_hogs, 2))

# --- locate ERC results file ---
if args.branch_method == "BXB":
    input_file = os.path.join(erc_results_folder, "ERC_results_BXB.tsv")
elif args.branch_method == "R2T":
    input_file = os.path.join(erc_results_folder, "ERC_results_R2T.tsv")
else:
    raise ValueError('Invalid branch_method argument. Use BXB or R2T.')
print(f"Reading ERC results from: {input_file}")

chunk_size = 100000  # tune based on available memory
chunks = pd.read_csv(input_file, sep="\t", chunksize=chunk_size)

# --- filter chunks and collect results ---
filtered_results = []
for chunk in chunks:
    chunk = chunk.dropna(subset=["GeneA_HOG", "GeneB_HOG"])
    mask = chunk.apply(
        lambda row: frozenset([row["GeneA_HOG"], row["GeneB_HOG"]]) in hog_pairs,
        axis=1
    )
    filtered_chunk = chunk[mask]
    filtered_results.append(filtered_chunk)
    print(f"Filtering chunk: found {len(filtered_chunk)} filtered rows")

# Concatenate all filtered chunks into a single DataFrame
df = pd.concat(filtered_results, ignore_index=True)
print(f"Total filtered ERC results: {len(df)}")
print("Filtered ERC results shape:", df.shape)
print("First few rows of filtered ERC results:")
print(df.head())
print("Unique GeneA_IDs:", df["GeneA_ID"].unique())
print("Unique GeneB_IDs:", df["GeneB_ID"].unique())

# --- build matrix ---
print("Building heatmap matrix...")
# Use HOGs for matrix construction
# Preserve the order of HOGs from the input `HOGs_of_interest` file (filtered_hogs).
# If the ERC filtered results contain additional HOGs not in the input file,
# append them at the end in sorted order so they are not lost.
data_hogs = list(dict.fromkeys(list(df["GeneA_HOG"]) + list(df["GeneB_HOG"])))
extra_hogs = [h for h in data_hogs if h not in filtered_hogs]
all_hogs = list(filtered_hogs) + sorted(extra_hogs)
print("HOGs in matrix axes (preserving input order, extras appended):", all_hogs)
print("Unique GeneA_HOGs in filtered results:", df["GeneA_HOG"].unique())
print("Unique GeneB_HOGs in filtered results:", df["GeneB_HOG"].unique())
print("Example HOGs from matrix:", all_hogs[:5])
print("Example HOGs from filtered results:", df["GeneA_HOG"].unique()[:5])
heatmap_df = pd.DataFrame(index=all_hogs, columns=all_hogs, dtype=float)


for _, row in df.iterrows():
    a_hog, b_hog = row["GeneA_HOG"], row["GeneB_HOG"]
    Pearson_R = row["Pearson_R"]
    Spearman_R = row["Spearman_R"]
    Kendall_Tau = row["Kendall_Tau"]
    #print(f"P: {Pearson_R}, S: {Spearman_R}, K: {Kendall_Tau}")
    av_R_val = np.mean([Pearson_R, Spearman_R, Kendall_Tau])
    if np.isnan(av_R_val):
        print(f"NaN correlation for pair: {a_hog}, {b_hog}")
    heatmap_df.loc[a_hog, b_hog] = av_R_val
    heatmap_df.loc[b_hog, a_hog] = av_R_val

for hog in all_hogs:
    heatmap_df.loc[hog, hog] = np.nan

non_na_count = np.sum(~np.isnan(heatmap_df.values))
print("Heatmap matrix non-NaN count:", non_na_count)

# Print a sample of the heatmap matrix after filling
print("Heatmap matrix sample after filling:")
print(heatmap_df.iloc[:5, :5])

# Map HOGs to gene names for axis labels in the plot
hog_labels = [hog_to_gene.get(hog, hog) for hog in all_hogs]

mask = np.triu(np.ones_like(heatmap_df, dtype=bool))

# --- save outputs ---

print(f"Saving heatmap data to TSV: {output_tsv}")
heatmap_df.to_csv(output_tsv, sep = "\t")

# --- plot ---
if df.shape[0] == 0:
    print("No filtered ERC results found. PDF will be blank.")
elif heatmap_df.count().sum() == 0:
    print("Heatmap matrix is empty (all NaN). PDF will be blank.")
else:
    print("Creating heatmap plot...")
    plt.figure(figsize=(12, 10))
    sns.heatmap(
        heatmap_df,
        mask=mask,
        cmap="viridis",
        square=True,
        linewidths=0.5,
        linecolor='gray',
        cbar_kws={'label': 'Average_Rvalue'},
        xticklabels=hog_labels,
        yticklabels=hog_labels
    )
    plt.title("Average R value from Peason, Spearman, and Kendall correlations")
    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig(output_pdf, format="pdf")
    plt.show()
    print(f"writing heatmap PDF to: {output_pdf}")
    print("Done! Heatmap figure displayed and saved.")