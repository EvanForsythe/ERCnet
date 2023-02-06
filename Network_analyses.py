#!/usr/bin/env python3
'''
Script for performing Network analyses

./Network_analyses.py -j TPC_test -m bxb -f pval -c 0.05 -y fg -s Atha

#delete previous runs
#rm -r DLCpar/ BL_results/

'''

#Import modules
import os
import re
import sys
import glob
import subprocess
import argparse
import pandas as pd
from ERC_functions import *

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
parser.add_argument('-c', '--CorrStat', type=str, metavar='', required=False, help='The type of statistical correlation method used from ERC_analyses.py. Enter "spearman" or "pearson".', default='spearman')

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
Corrmethod = args.CorrStat


'''
JOBname ="TEST"
BLmethod="bxb"
Filterstat="pval"
Cutoff=float(0.005)
Clustmeth="fg"
Trim_Cutoff=int(0)
FocalSP="Atha"
'''
#Store output dir as a variable
out_dir= 'OUT_'+JOBname+'/'
fileName = 'ERC_results_' + BLmethod + '_' + Corrmethod + '.tsv'

#Make a directory for Network outputs
if not os.path.isdir(out_dir+'Network_analyses/'):
    os.makedirs(out_dir+'Network_analyses/')
    print("created folder: Network_analyses/\n\n")
else: 
    print('ERC results will be stored to Network_analyses/\n\n')
    
#Make a directory for TSV files from Networks
if not os.path.isdir(out_dir+'Network_analyses/Communities/'):
    os.makedirs(out_dir+'Network_analyses/Communities/')
    print("created folder: Network_analyses/Communities/\n\n")
else: 
    print('ERC results will be stored to Network_analyses/Communities/\n\n')

#Check to verify the correct methods have been chosen have a relevant filetype
print('Verifying chosen branch and statistical methods have a relevant filetype...')
if os.path.isfile(out_dir+'ERC_results/'+fileName):
    print('File found. Proceeding on analysis using: ' + str(fileName))
else:
    print('No file could be found at the location. Please check -m and -c flags match from ERC_analysis.py step.')
    print('Analysis will now exit...')
    sys.exit()

#Run the R script
print("Beginning network analyses using the R package, igraph...\n\n Calling R...\n\n")

#load in the ERCresults file for filtering, prior to sending over to R Network analysis script.
tsvData = out_dir + 'ERC_results/' + fileName
csvData = pd.read_table(tsvData, sep='\t')

#Calls functions from ERC_functions.py to filter the ERC_results file down based on user provided criteria
#FilterBranchType(csvData, branchFilter)
#FilterCorrelationType(csvData, corrFilter)
csvData = FilterSignificance(csvData, RSquared, PValue)

#Output a filtered version of the ERC_results file
csvData.to_csv(out_dir + "ERC_results/Filtered_" + fileName, sep='\t', index=False, header=True) 

if sum(1 for line in open(out_dir+"ERC_results/Filtered_" + fileName)) < 2:
    print("It appears that not enough ERC results were retained for further analysis. Consider changing filtering criteria or analysis methods.")
    print("Quitting...")
    sys.exit()

#Run the Network analyses.
#make command.
Net_cmd= 'Rscript Networks_and_stats.R '+JOBname+" "+BLmethod+" "+str(RSquared)+" "+str(PValue)+" "+Clustmeth+" "+str(Trim_Cutoff)+" "+FocalSP+" "+'Filtered_' + fileName
    
#Run the command (if it contains strings expected in the command, this is a precaution of using shell=True)
if re.search('Networks_and_stats.R', Net_cmd) and re.search('Rscript', Net_cmd):
    print("Calling R with the following command:")
    print(Net_cmd)
    subprocess.call(Net_cmd, shell=True)

if len(glob.glob(out_dir+'Network_analyses/*pdf')) > 0:
    print('Finished network analyses.\n\nResults files written to Network_analyses/\n\n')
else:
    print('Something went wrong with network analyses...\n\n')


