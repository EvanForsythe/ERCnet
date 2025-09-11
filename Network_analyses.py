#!/usr/bin/env python3

#Import modules
import os
import re
import sys
import glob
import time
import subprocess
import argparse
import pandas as pd
import ERC_functions as erc
from scipy import stats

def benchmarkTime(fileName, path, stamp, process, timer):
    timer = time.localtime()
    current_time = time.strftime("%a, %d %b %Y %H:%M:%S", timer)
    with open(path + fileName, "a") as bench:
        bench.write(process + " Call " + stamp  + '\t' + str(current_time) + '\n')
    return 0

#At runtime set working directory to the place where the script lives
working_dir = sys.path[0]+'/' 
os.chdir(working_dir)

#Set up an argumanet parser
parser = argparse.ArgumentParser(description='Network step')

parser.add_argument('-j', '--JOBname', type=str, required=True, help='Unique job name for this run of ERCnet. Avoid including spaces or special characters ("_" is ok)') 
parser.add_argument('-m', '--BLmethod', type=str, required=False, help='Branch length method ERC results to be used in the network. "BXB" for Branch-by-branch. "R2T" for root-to-tip. Default is R2T', default = 'r2t') 
parser.add_argument('-y', '--Clustmeth', type=str, required=True, help='Clustering method to be used to identify communities in network. "fg" for fast-and-greedy (fastest), "eb" for edge-betweenness, "op" for optimal, and "wt" for walktrap. Note that this feature has largely depreciated in favor of using Cytoscape to visualize networks.') 
parser.add_argument('-t', '--Trim_Cutoff', type=int, required=False, help='The user-selected cutoff will be the minimum number of genes necessary for a community to be displayed on the network plot.This is mainly for network visualization and is not recommended for data collection. Must be an integer. 0 (no trimming) is default.', default=0)
parser.add_argument('-s', '--FocalSP', type=str, required=True, help='The name of the focal species to represent each gene family (should exactly match the tip label of the species tree)') 
parser.add_argument('-pp', '--PearsonP', type=float, required=False, help='P-value cutoff for Pearson correlation. Default is 0.05', default=0.05)
parser.add_argument('-pr','--PearsonR', type=float, required=False, help='R-value cutoff for Pearson correlation. Default is 0.5', default=0.5)
parser.add_argument('-sp','--SpearmanP', type=float, required=False, help='P-value cutoff for Spearman correlation. Default is 0.05', default=0.05)
parser.add_argument('-sr','--SpearmanR', type=float, required=False, help='R-value cutoff for Spearman correlation. Default is 0.5', default=0.5)
parser.add_argument('-kp','--KendallP', type=float, required=False, help='P-value cutoff for Kendall correlation. Default is 0.05', default=0.05)
parser.add_argument('-kr','--KendallR', type=float, required=False, help='R-value cutoff for Kendall correlation. Default is 0.5', default=0.5)
parser.add_argument('-F', '--Func_cat', action='store_true', required=False, help='Run a functional clustering analysis with user-provided functional information about genes in the focal species? If selected, youll need to provide two tsv files. See documentation for formatting.') 
parser.add_argument('-L', '--Lab_nodes', action='store_true', required=False, help='Add node labels to the network? If selected, youll need to provide a tsv files of node labels. See documentation for formatting.') 
parser.add_argument('-S', '--Strict', action='store_true', required=False, help='Further filters ERC results based on Benjamini-Hochberg False Discovery Rate')


#Define the parser
args = parser.parse_args()

#Store arguments
JOBname=args.JOBname
BLmethod=args.BLmethod
Clustmeth=args.Clustmeth
Trim_Cutoff=args.Trim_Cutoff
FocalSP=args.FocalSP
pearson_p=args.PearsonP
pearson_r=args.PearsonR
spearman_p=args.SpearmanP
spearman_r=args.SpearmanR
kendall_p=args.KendallP
kendall_r=args.KendallR
func_bool=args.Func_cat
lab_bool=args.Lab_nodes
strict=args.Strict

out_dir = f'OUT_{JOBname}/'

erc.CheckAndMakeDir(out_dir, 'Network_analyses/')

#Define the time object and folder for optimization testing
bench_fileName = JOBname + '_Network_analysis_benchmark.tsv'
timer = time.localtime()
current_time = time.strftime("%H:%M:%S", timer)

benchmarkTime(bench_fileName, os.path.join(out_dir, 'benchmark') + '/', 'start', 'Network_analysis', timer)

#Check to verify the correct methods have been chosen have a relevant filetype
print('Verifying chosen branch and statistical methods have a relevant filetype...')

fileName = "ERC_results_"+BLmethod+".tsv"

erc_results_path = os.path.join(out_dir, 'ERC_results', fileName)
print(f"Checking for ERC results file at: {erc_results_path}")
if os.path.isfile(erc_results_path):
    print('File found. Proceeding to analysis using: ' + str(fileName))
else:
    print(f'No file could be found at the location: {erc_results_path}. Please check -m and -c flags match from ERC_analysis.py step.')
    print('Analysis will now exit...')
    sys.exit()

#Check to see whether Func_cat was specified
if func_bool:
    print("Functional category analysis selected")
    if os.path.isfile("Functional_categories.tsv"):
        print('Functional_categories.tsv file found.')

        
    else:
        print('Cannot find Functional_categories.tsv file. Quitting...\n')
        sys.exit()
    if os.path.isfile("Functional_categories_col_assign.tsv"):
        print('Functional_categories_col_assign.tsv file found.')
    else:
        print('Cannot find Functional_categories_col_assign.tsv file. Quitting...\n')
        sys.exit()
        
if lab_bool:
    if func_bool:
        print("Node labels selected")
        if os.path.isfile("Node_labels.tsv"):
            print('Node_labels.tsv file found.')
        else:
            print('Cannot find Node_labels.tsv file. Quitting...\n')
            sys.exit()

    
#load in the ERCresults file for filtering, prior to sending over to R Network analysis script.
tsvData = erc_results_path
csvData = pd.read_table(tsvData, sep='\t')

#Remove negative Slopes
csvData = csvData[:][csvData['Slope'] > 0]

#Adds FDR data from filtered values
SPs = csvData['Spearman_Pval']
PPs = csvData['Pearson_Pval']

newSPs = erc.false_discovery_control(SPs, axis=0, method='bh')
newPPs = erc.false_discovery_control(PPs, axis=0, method='bh')

csvData['S_FDR_Corrected_Pval'] = newSPs
csvData['P_FDR_Corrected_Pval'] = newPPs

#Check if Kendall's tau is in the dataset
if 'Kendall_Pval' in csvData.columns.tolist():
    KPs = csvData['Kendall_Pval']
    newKPs = erc.false_discovery_control(KPs, axis=0, method='bh')
    csvData['K_FDR_Corrected_Pval'] = newKPs

#Filters based on user selected strictess. 

# Apply filtering for each correlation method with its own cutoffs
if strict:
    csvData = erc.FilterFDR(csvData,
                            pearson_r, pearson_p,
                            spearman_r, spearman_p,
                            kendall_r, kendall_p)
else:
    csvData = erc.FilterSignificance(csvData,
                                     pearson_r, pearson_p,
                                     spearman_r, spearman_p,
                                     kendall_r, kendall_p)

outFileName = fileName.replace('.tsv', '') + f"_PearsonP{pearson_p}_PearsonR{pearson_r}_SpearmanP{spearman_p}_SpearmanR{spearman_r}_KendallP{kendall_p}_KendallR{kendall_r}"

# Create dynamic results subfolder name
if strict:
    results_subfolder = f"Results_{BLmethod}_fdr_pp{pearson_p}_ps{spearman_p}_kp{kendall_p}_pr{pearson_r}_sr{spearman_r}_kr{kendall_r}/"

else:
    results_subfolder = f"Results_{BLmethod}_raw_pp{pearson_p}_ps{spearman_p}_kp{kendall_p}_pr{pearson_r}_sr{spearman_r}_kr{kendall_r}/"


erc.CheckAndMakeDir(os.path.join(out_dir, 'Network_analyses/'), results_subfolder)

results_dir = os.path.join(out_dir+'Network_analyses', results_subfolder)


# Create Func_categories_rundata.tsv in a fixed location outside the dynamic results folder
func_cat_rundata_path = os.path.join(out_dir, 'Network_analyses', 'Func_categories_rundata.tsv')
if func_bool:
    if os.path.isfile(func_cat_rundata_path):
        print(f"metadata about the func category runs will be stored in {func_cat_rundata_path}")
    else:
        print(f"creating file to store metadata about the func category runs: {func_cat_rundata_path}")
        os.makedirs(os.path.dirname(func_cat_rundata_path), exist_ok=True)
        with open(func_cat_rundata_path, "a") as f:
            f.write("Run_name" + "\t" + "N_nodes" + "\t" + "N_connections" + "\t" + "Network_clust_algo" + "\t" +"Obs_assort_coef" + "\t" + "Assort_Zscore" + "\t" + "Assort_Pval") 

# Filtered ERC hits file name
filtered_filename = "Filtered_ERC_hits.tsv"

#Output a filtered version of the ERC_results file
csvData.to_csv(results_dir + filtered_filename, sep='\t', index=False, header=True)

if sum(1 for line in open(results_dir + filtered_filename)) < 2:
    print("It appears that not enough ERC results were retained for further analysis. Consider changing filtering criteria or analysis methods.")
    print("Quitting...")
    sys.exit()

#Run the Network analyses.
Net_cmd = (
    f"Rscript Networks_and_stats.R {JOBname} {BLmethod} {Clustmeth} {Trim_Cutoff} {FocalSP} "
    f"{filtered_filename} {func_bool} {lab_bool}"
)

Net_cmd = (
    f"Rscript Networks_and_stats.R {JOBname} {results_dir} {BLmethod} {Clustmeth} {Trim_Cutoff} {FocalSP} "
    f"{filtered_filename} {func_bool} {lab_bool}"
)
    
#Run the command (if it contains strings expected in the command, this is a precaution of using shell=True)
if re.search('Networks_and_stats.R', Net_cmd) and re.search('Rscript', Net_cmd):

    # Build the Rscript command using the absolute path
    #Net_cmd = (
    #    f"Rscript {r_script_path} {results_dir} {JOBname} {BLmethod} {Clustmeth} {Trim_Cutoff} {FocalSP} "
    #    f"{filtered_filename} {func_bool} {lab_bool}"
    #)
    
    print("Calling R with the following command:")
    print(Net_cmd)
    subprocess.call(Net_cmd, shell=True)

benchmarkTime(bench_fileName, os.path.join(out_dir, 'benchmark') + '/', 'end', 'Network_analysis', timer)


if (erc.CheckFileExists(f'{results_dir}*pdf') and erc.CheckFileExists(f'{results_dir}*csv')):
    print('Finished network analyses.\n\nResults files written to Network_analyses/\n\n')
    print('ERCnet finished! :)') 
else:
    print('Something went wrong with network analyses...\n\n')


