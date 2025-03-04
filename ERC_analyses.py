#!/usr/bin/env python3
'''
Script for performing BL reconcilation and ERC correlation analyses. This should be run in python3

'''

#Import modules
import os
import re
import sys
import time
import math
import glob
import random
import argparse
import itertools
import subprocess
import numpy as np
import pandas as pd
import ERC_functions as erc
import scipy.stats as stats
from datetime import datetime
from joblib import Parallel, delayed
from scipy.stats import pearsonr, spearmanr, linregress, kendalltau


#At runtime set working directory to the place where the script lives
working_dir = sys.path[0]+'/' 
os.chdir(working_dir)

#Set up an argumanet parser
parser = argparse.ArgumentParser(description='ERC step')

parser.add_argument('-j', '--JOBname', type=str, metavar='', required=True, help='Unique job name for this run of ERCnet. Avoid including spaces or special characters ("_" is ok)') 
parser.add_argument('-m', '--Mult_threads', type=int, metavar='', required=False, default=1, help='Integer: number of threads avilable for parallel computing (default = 1)' )
parser.add_argument('-s', '--FocalSP', type=str, metavar='', required=True, help='The name of the focal species to represent each gene family (should exactly match the tip label of the species tree)') 
parser.add_argument('-b', '--branchMethod', type=str, metavar='', required=False, default="R2T", help='This determins which branch reconcilliation method to use. Enter either "BXB" for branch by branch or "R2T" root to tip.')
parser.add_argument('-t', '--test', type=int, metavar='', required=False, default=0, help='Tests ERC_analyses.py on a small, randomized subset of the data equalling the integer provied.') 

#Define the parser
args = parser.parse_args()

#Store arguments
JOBname=args.JOBname
Mult_threads=args.Mult_threads
FocalSP=args.FocalSP
branchMethod=args.branchMethod
testNum = args.test

#Store output dir as a variable
out_dir= 'OUT_'+JOBname+'/'
fileName = str('ERC_results_' + branchMethod + '.tsv')
bench_fileName = str(JOBname + '_ERC_analyses_benchmark_' + branchMethod + '.tsv')
log_fileName = str(JOBname + '_ERC_analysis_log.tsv')
ERC_command = str('Job: ' + JOBname + ', Threads: ' + str(Mult_threads) + ', FocalSP: ' + str(FocalSP) + ', BranchMethod: ' + str(branchMethod)) 

print('\n')
print(str(branchMethod) + ' chosen for branch method.\n')

#remove previous file (if it exists)
if (erc.CheckFileExists(out_dir+'ERC_results/'+fileName)):
    fileName = erc.GetLatestFileName(out_dir+'ERC_results/', fileName[0:-4])
    fileName = erc.AppendFileName(fileName)
    print('Current ERC_analyses run will be saved to file: ' + str(fileName) + '\n')
    bench_fileName = JOBname + '_ERC_analyses_benchmark_' + fileName[12:]

#Define the time object and folder for optimization testing
timer = time.localtime()
current_time = time.strftime("%a, %d %b %Y %H:%M:%S", timer)

ERC_command = str('Job: ' + JOBname + ', Threads: ' + str(Mult_threads) + ', FocalSP: ' + str(FocalSP) + ', BranchMethod: ' + str(branchMethod) + ', Date: ' + str(current_time))

#Make a directory for Benchmark results
print("Checking or creating folder for benchmark.\n")
erc.CheckAndMakeDir(out_dir, 'benchmark/')

#Adds first timestamp for process start.
with open(out_dir + 'benchmark/' + str(bench_fileName), "a") as bench:
    bench.write("Job Name: " + JOBname + '\t' + "Cores: " + str(Mult_threads) + '\n')
    bench.write("Stage" + '\t' + "Time" + '\n')
    bench.write("Process Start" + '\t' + str(current_time) + '\n')


#Make a directory for ERC results
print("Checking or creating folder for ERC Results.\n")
erc.CheckAndMakeDir(out_dir, 'ERC_results/')

with open(out_dir + 'ERC_results/' + str(log_fileName), "a") as log:
    log.write(ERC_command + '\t' + fileName + '\n')

#Verify files exist in DLCpar/ for BL reconciliation in R. 
if not (erc.CheckFileExists(out_dir + 'DLCpar/*_NODES_BL.txt.dlcpar.locus.recon')):
    print('It appears that something has gone wrong when running GTST_reconciliation.py. Please verify input and output files. ERCnet will now exit.')
    sys.exit()

### Run the BL_reconciliation
print("beginning BL reconciliation in R...\n\n Calling R...\n\n")

#Run the BL reconcilation step.
BL_rec_cmd= 'Rscript Branch_length_reconciliation.R '+JOBname
    
#Run the command (if it contains strings expected in the command, this is a precautin of using shell=True)
if re.search('Branch_length_reconciliation.R', BL_rec_cmd) and re.search('Rscript', BL_rec_cmd):
    print("Running BL reconciliation in R with the folllowing command:")
    print(BL_rec_cmd)
    subprocess.call(BL_rec_cmd, shell=True)

if (erc.CheckFileExists(out_dir+'BL_results/*tsv')):
    print('Finished branch length reconciliation.\n\nResults files written to BL_results/\n\n')

### Run all-by-all correlations
#Read in BL results
if (branchMethod == 'BXB'):
    BLs=pd.read_csv(out_dir+'BL_results/bxb_BLs_normalized.tsv', sep='\t')
else:
    BLs=pd.read_csv(out_dir+'BL_results/r2t_BLs_normalized.tsv', sep='\t')

#Read in gene IDs
gene_fams=pd.read_csv(out_dir+'Filtered_genefam_dataset.csv')

#Clean up the gene IDs dataframe
#Find out what string needs to be removed from HOGs and use replace() to remove it
gene_fams['HOG']=gene_fams['HOG'].str.replace(str(gene_fams['HOG'][0].split(".", 1)[0]+"."), "")

#Verify list of BLs['HOG_ID'] is at least 2. 
if len(BLs['HOG_ID']) < 2:
    print('It appears something has gone wrong and there are no pairs of trees to compare.')
    print('Please verify that the number of gene families from ./Phylogenomics.py is greater than 1 (-t > 1). ERCnet will now exit.')
    sys.exit() 


#Make list of pairwise combinations
pairwise_combos=list(itertools.combinations(BLs['HOG_ID'], 2))


if (testNum > 0):
    print("Testing correlation with " + str(testNum) + " pairwise combinations.")
    pairwise_combos=random.sample(pairwise_combos, testNum)


# Make results file with a header
'''
results_file = os.path.join(out_dir, 'ERC_results', str(fileName))
with open(results_file, "w") as f:
    f.write("GeneA_HOG\tGeneA_ID\tGeneB_HOG\tGeneB_ID\tOverlapping_branches\tSlope\tP_R2\tP_Pval\tS_R2\tS_Pval\tP_FDR_Corrected_Pval\tS_FDR_Corrected_Pval\n")
'''

# Updated function to calculate correlation and extract statistics
def calculate_correlation(geneA, geneB, hog_to_gene_id, BLs, bl_hog_ids):
    """
    Perform correlation calculation for a single pair of genes and return extracted statistics.

    Parameters:
    - geneA, geneB: The HOG IDs for the two genes being compared.
    - hog_to_gene_id: Dictionary mapping HOG IDs to gene IDs for the focal species.
    - BLs: DataFrame containing branch length information.
    - bl_hog_ids: Set containing the HOG_IDs from BLs for faster membership testing.

    Returns:
    - A string with the format:
      "GeneA_HOG\tGeneA_ID\tGeneB_HOG\tGeneB_ID\tNumber of Rows\tSlope\tPearson r-squared\tPearson p-value\tSpearman r-squared\tSpearman p-value"
      If correlation cannot be computed, return None.
    """

    # Use the dictionary to get gene IDs for the focal species
    geneA_ID = hog_to_gene_id.get(geneA)
    geneB_ID = hog_to_gene_id.get(geneB)

    # If either geneA or geneB is not found in the dictionary, skip this pair
    if geneA_ID is None or geneB_ID is None:
        return None

    # Check if geneA and geneB are present in the BLs DataFrame using set membership testing
    if geneA not in bl_hog_ids or geneB not in bl_hog_ids:
        return None

    # Subset the BLs dataframe for the two test HOGs using .isin() for faster filtering
    test_df_temp = BLs[BLs["HOG_ID"].isin([geneA, geneB])]

    # If the resulting subset is empty, skip this pair
    if test_df_temp.empty:
        return None

    # Set 'HOG_ID' as the index and transpose the DataFrame
    test_df_temp = test_df_temp.set_index("HOG_ID")
    test_df_clean = test_df_temp.T

    # Ensure the transposed DataFrame has exactly two columns and assign column names
    if len(test_df_clean.columns) == 2:
        test_df_clean.columns = ['GeneA', 'GeneB']
    else:
        return None

    # Remove rows with NaN values
    test_df_clean = test_df_clean.dropna()

    # Check if we have at least 5 data points for correlation analysis
    if test_df_clean.shape[0] > 4:
        try:
            # Convert to numpy arrays for faster numerical operations
            x = test_df_clean['GeneA'].values
            y = test_df_clean['GeneB'].values

            # Calculate Pearson correlation and p-value
            pearson_corr, pearson_pval = pearsonr(x, y)

            # Calculate Spearman correlation and p-value
            spearman_corr, spearman_pval = spearmanr(x, y)

            # Calculate Kendall's tau correlation and p-value
            kendall_corr, kendall_pval = kendalltau(x, y)
 
            # Calculate slope using linear regression
            slope, _, _, _, _ = linregress(x, y)

            # Store results in the specified order
            results_str = (
                f"{geneA}\t{geneA_ID}\t{geneB}\t{geneB_ID}\t"
                f"{test_df_clean.shape[0]}\t"
                f"{slope}\t"
                f"{pearson_corr**2}\t"
                f"{pearson_pval}\t"
                f"{spearman_corr**2}\t"
                f"{spearman_pval}\t"
                f"{kendall_corr}\t"
                f"{kendall_pval}"
            )

            return results_str
        except ValueError:
            # Skip this pair if there are not enough data points or if input is invalid for correlation
            return None
    else:
        return None


# Function to perform all-by-all correlation analyses with parallelization and incremental result writing
def perform_all_by_all_parallel_incremental(pairwise_combos, gene_fams, BLs, out_dir, fileName, FocalSP, n_jobs, chunk_size):
    """
    Perform all-by-all correlation analyses with parallelization and incremental result writing.

    Parameters:
    - pairwise_combos: List of tuples with all pairwise combinations of HOG_IDs.
    - gene_fams: DataFrame containing gene family information.
    - BLs: DataFrame containing branch length information.
    - out_dir: Output directory where results file will be saved.
    - fileName: Name of the output results file.
    - FocalSP: The column name in gene_fams that stores the gene IDs for the focal species.
    - n_jobs: Number of parallel jobs (default is -1, which uses all available processors).
    - chunk_size: Number of pairwise combinations to process in each chunk (default is 20,000).

    Returns:
    - None, but writes results to a specified file in the output directory.
    """

    # Create a dictionary mapping HOG IDs to gene IDs for the focal species for faster lookup
    hog_to_gene_id = dict(zip(gene_fams['HOG'], gene_fams[FocalSP]))

    # Create a set of HOG IDs in the BLs DataFrame for faster membership testing
    bl_hog_ids = set(BLs['HOG_ID'].values)

    # Prepare results file with header
    results_file = os.path.join(out_dir, 'ERC_results', str(fileName))
    with open(results_file, "w") as f:
        f.write("GeneA_HOG\tGeneA_ID\tGeneB_HOG\tGeneB_ID\tOverlapping_branches\tSlope\tP_R2\tP_Pval\tS_R2\tS_Pval\tK_Tau\tK_Pval\tP_FDR_Corrected_Pval\tS_FDR_Corrected_Pval\tK_FDR_Corrected_Pval\n")
    
    # Divide the pairwise_combos list into chunks of size `chunk_size`
    for chunk_start in range(0, len(pairwise_combos), chunk_size):
        chunk_end = min(chunk_start + chunk_size, len(pairwise_combos))
        chunk = pairwise_combos[chunk_start:chunk_end]

        # Use joblib.Parallel to parallelize the computation across multiple cores for each chunk
        chunk_results = Parallel(n_jobs=n_jobs)(
            delayed(calculate_correlation)(geneA, geneB, hog_to_gene_id, BLs, bl_hog_ids) for geneA, geneB in chunk
        )

        # Filter out None results
        chunk_results = [result for result in chunk_results if result is not None]

        # Write chunk results to the file
        if chunk_results:
            with open(results_file, "a") as f:
                f.write('\n'.join(chunk_results) + '\n')

        # Print progress message
        print(f"Processed chunk {chunk_start} to {chunk_end} of {len(pairwise_combos)}")

    print("All-by-all parallel incremental analysis completed.")

#Timestamp just before function call
erc.benchmarkTime(bench_fileName, out_dir + 'benchmark/', 'start', 'Correlation', timer)

perform_all_by_all_parallel_incremental(pairwise_combos, gene_fams, BLs, out_dir, fileName, FocalSP, n_jobs=int(Mult_threads), chunk_size=100000)

#Timestamp just after function call
erc.benchmarkTime(bench_fileName, out_dir + 'benchmark/', 'end', 'Correlation', timer) 

#Calculate total time of process & lines per minute
num_lines = len(pd.read_csv(out_dir + 'ERC_results/' + str(fileName), sep='\t'))
erc.benchmarkProcess(out_dir + 'benchmark/' + str(bench_fileName), num_lines)

#Progress message
print("Done with all-by-all comparisons...\n")

#Report status
if (erc.CheckFileExists(out_dir+'ERC_results/'+fileName)):
    print('Finished ERC.\n\nResults files written to ERC_results/\n\n')
else:
    print('Something went wrong with ERC analyses...\n')

#Print directions for next steps
print('ERC correlation analyses finished. Exiting....\n\n'\
      'To perform network analyses, run Network_analyses.py. See associated help menu for required arguments.'\
          '\nExample command:\n\n' \
          './Network_analyses.py -j '+JOBname+' -m R2T -c spearman -y fg -s ' +FocalSP+ ' -f <name of .tsv file located in ERC_results/> \n\n' 
      )

