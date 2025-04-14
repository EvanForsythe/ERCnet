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

parser.add_argument('-j', '--JOBname', type=str, metavar='', required=True, help='Unique job name for this run of ERCnet. Avoid including spaces or special characters ("_" is ok)') 
parser.add_argument('-m', '--BLmethod', type=str, metavar='', required=False, help='Branch length method ERC results to be used in the network. "BXB" for Branch-by-branch. "R2T" for root-to-tip. Default is R2T', default = 'r2t') 
parser.add_argument('-p', '--PValue', type=float, metavar='', required=False, help='Cuttoff for P value by which to filter ERC hits in network. Float between 0 and 1. Default value is 0.05', default=0.05) 
parser.add_argument('-r', '--RSquared', type=float, metavar='', required=False, help='Cuttoff R-squared value by which to filter ERC hits in network. Float between 0 and 1. Default value is 0.50', default=0.50) 
parser.add_argument('-y', '--Clustmeth', type=str, metavar='', required=True, help='Clustering method to be used to identify communities in network. "fg" for fast-and-greedy (fastest), "eb" for edge-betweenness, "op" for optimal, and "wt" for walktrap.') 
parser.add_argument('-t', '--Trim_Cutoff', type=int, metavar='', required=False, help='The user-selected cutoff will be the minimum number of genes necessary for a community to be displayed on the network plot.This is mainly for network visualization and is not recommended for data collection. Must be an integer. 0 (no trimming) is default.', default=0)
parser.add_argument('-s', '--FocalSP', type=str, metavar='', required=True, help='The name of the focal species to represent each gene family (should exactly match the tip label of the species tree)') 
parser.add_argument('-c', '--CorrMethod', type=str, metavar='', required=False, help='The type of correlation method you would like to filter P Value and R value by. Should be "pearson", "spearman", "kendall", or "all". Default is pearson.', default='pearson')
parser.add_argument('-f', '--FileName', type=str, metavar='', required=True, help='The filename of ERC_results file you would like to analyze. Should be .tsv file.')
parser.add_argument('-F', '--Func_cat', action='store_true', required=False, help='Run a functional clustering analysis with user-provided functional information about genes in the focal species? If selected, youll need to provide two tsv files. See documentation for formatting.') 
parser.add_argument('-L', '--Lab_nodes', action='store_true', required=False, help='Add node labels to the network? If selected, youll need to provide a tsv files of node labels. See documentation for formatting.') 
parser.add_argument('-S', '--Strict', action='store_true', required=False, help='Further filters ERC results based on Benjamini-Hochberg False Discovery Rate')


#Define the parser
args = parser.parse_args()

#Store arguments
JOBname=args.JOBname
BLmethod=args.BLmethod
PValue=args.PValue
RSquared=args.RSquared
Clustmeth=args.Clustmeth
Trim_Cutoff=args.Trim_Cutoff
FocalSP=args.FocalSP
Corrmethod = args.CorrMethod
fileName=args.FileName
func_bool=args.Func_cat
lab_bool=args.Lab_nodes
strict = args.Strict

#Store output dir as a variable
out_dir= 'OUT_'+JOBname+'/'
#fileName = file.replace('.tsv', '')

# Check if Corrmethod is valid
if Corrmethod not in ['pearson', 'spearman', 'kendall', 'all']:
    print("Invalid option for --CorrMethod/-c")
    exit(1)  # Exit the program with a non-zero status

#Define the time object and folder for optimization testing
bench_fileName = JOBname + '_Network_analysis_benchmark.tsv'
timer = time.localtime()
current_time = time.strftime("%H:%M:%S", timer)

benchmarkTime(bench_fileName, out_dir + 'benchmark/', 'start', 'Network_analysis', timer)



#Make a directory for Network outputs
if not os.path.isdir(out_dir+'Network_analyses/'):
    os.makedirs(out_dir+'Network_analyses/')
    print("created folder: Network_analyses/\n\n")
else: 
    print('ERC results will be stored to Network_analyses/\n\n')

erc.CheckAndMakeDir(out_dir+'ERC_results/', 'Filtered_results')

    
#Make a directory for TSV files from Networks
if not os.path.isdir(out_dir+'Network_analyses/Communities/'):
    os.makedirs(out_dir+'Network_analyses/Communities/')
    print("created folder: Network_analyses/Communities/\n\n")
else: 
    print('ERC results will be stored to Network_analyses/Communities/\n\n')

#Check to verify the correct methods have been chosen have a relevant filetype
print('Verifying chosen branch and statistical methods have a relevant filetype...')
if os.path.isfile(out_dir+'ERC_results/'+fileName):
    print('File found. Proceeding to analysis using: ' + str(fileName))
else:
    print('No file could be found at the location. Please check -m and -c flags match from ERC_analysis.py step.')
    print('Analysis will now exit...')
    sys.exit()

#Check to see whether Func_cat was specified
if func_bool:
    print("Functional category analysis selected")
    if os.path.isfile("Functional_categories.tsv"):
        print('Functional_categories.tsv file found.')
        if os.path.isfile(out_dir+'Network_analyses/Func_categories_rundata.tsv'):
            print("metadata about the func category runs with be stored in Network_analyses/Func_categories_rundata.tsv'")
        else:
            print("creating file to store metadata about the func category runs: Network_analyses/Func_categories_rundata.tsv")
            #Make results file
            with open(out_dir+'Network_analyses/Func_categories_rundata.tsv', "a") as f:
                f.write("Branch_length_method" + "\t" + "Pval_cutoff" + "\t" + "R2_cutoff" + "\t" + "N_nodes" + "\t" + "N_connections" + "\t" + "Network_clust_algo" + "\t" +"Obs_assort_coef" + "\t" + "Assort_Zscore" + "\t" + "Assort_Pval") 
        
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
tsvData = out_dir + 'ERC_results/' + fileName
csvData = pd.read_table(tsvData, sep='\t')

#Remove negative Slopes
csvData = csvData[:][csvData['Slope'] > 0]

#Adds FDR data from filtered values
SPs = csvData['S_Pval']
PPs = csvData['P_Pval']

newSPs = erc.false_discovery_control(SPs, axis=0, method='bh')
newPPs = erc.false_discovery_control(PPs, axis=0, method='bh')

csvData['S_FDR_Corrected_Pval'] = newSPs
csvData['P_FDR_Corrected_Pval'] = newPPs

#Check if Kendall's tau is in the dataset
if 'K_Pval' in csvData.columns.tolist():
    KPs = csvData['K_Pval']
    newKPs = erc.false_discovery_control(KPs, axis=0, method='bh')
    csvData['K_FDR_Corrected_Pval'] = newKPs

#Filters based on user selected strictess. 
if (strict):
    csvData = erc.FilterFDR(csvData, RSquared, PValue, Corrmethod)
else:
    csvData = erc.FilterSignificance(csvData, RSquared, PValue, Corrmethod)

#Changes file output name based on strictness. Allowing for a saved filtered results file in both FDR and regular filtering. 
if (strict):
    outFileName = fileName.replace('.tsv', '') + "_" + str(Corrmethod) + "_" + str(PValue) + "_" + str(RSquared) + "_FDR.tsv"
else:
    outFileName = fileName.replace('.tsv', '') + "_"+ str(Corrmethod) + "_" + str(PValue) + "_" + str(RSquared) + ".tsv"

#Output a filtered version of the ERC_results file
csvData.to_csv(out_dir + "ERC_results/Filtered_results/Filtered_" + outFileName, sep='\t', index=False, header=True) 

if sum(1 for line in open(out_dir+"ERC_results/Filtered_results/Filtered_" + outFileName)) < 2:
    print("It appears that not enough ERC results were retained for further analysis. Consider changing filtering criteria or analysis methods.")
    print("Quitting...")
    sys.exit()

#Run the Network analyses.
#make command.
Net_cmd= 'Rscript Networks_and_stats.R '+JOBname+" "+BLmethod+" "+str(RSquared)+" "+str(PValue)+" "+Clustmeth+" "+str(Trim_Cutoff)+" "+FocalSP+" "+'Filtered_' + outFileName+" "+str(func_bool)+" "+str(lab_bool)
    
#Run the command (if it contains strings expected in the command, this is a precaution of using shell=True)
if re.search('Networks_and_stats.R', Net_cmd) and re.search('Rscript', Net_cmd):
    print("Calling R with the following command:")
    print(Net_cmd)
    subprocess.call(Net_cmd, shell=True)

benchmarkTime(bench_fileName, out_dir + 'benchmark/', 'end', 'Network_analysis', timer)


if (erc.CheckFileExists(out_dir+'Network_analyses/*pdf') and erc.CheckFileExists(out_dir+'Network_analyses/*csv')):
    print('Finished network analyses.\n\nResults files written to Network_analyses/\n\n')
    print('ERCnet finished! :)') 
else:
    print('Something went wrong with network analyses...\n\n')


