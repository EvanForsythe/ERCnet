#!/usr/bin/env python3
'''
Script for performing Network analyses

./Network_analyses.py -j TPC_test -m bxb -f pval -c 0.05 -y fg -s Atha

#delete previous runs
#rm -r DLCpar/ BL_results/

'''

#Import modules
import os
import sys
import glob
import subprocess
import argparse
import re


#During developent, set working directory:
#working_dir = '/Users/esforsythe/Documents/Work/Bioinformatics/ERC_networks/Analysis/ERCnet_dev/'
#os.chdir('/Users/esforsythe/Documents/Work/Bioinformatics/ERC_networks/Analysis/ERCnet_dev/')

#At runtime set working directory to the place where the script lives
working_dir = sys.path[0]+'/' 
os.chdir(working_dir)

#Set up an argumanet parser
parser = argparse.ArgumentParser(description='Network step')

parser.add_argument('-j', '--JOBname', type=str, metavar='', required=True, help='Unique job name for this run of ERCnet. Avoid including spaces or special characters ("_" is ok)') 
parser.add_argument('-m', '--BLmethod', type=str, metavar='', required=True, help='Branch length method ERC results to be used in the network. "bxb" for Branch-by-branch. "r2t" for root-to-tip.') 
parser.add_argument('-f', '--Filterstat', type=str, metavar='', required=True, help='Correlation statistic to be used to filter ERC hits. "pval" for p-value, "R2" for R-squared') 
parser.add_argument('-c', '--Cutoff', type=float, metavar='', required=True, help='Cuttoff P/R-squared value by which to filter ERC hits in network. float between 0 and 1. Correlations below/above this number will be retained (depending on your choice of pval vs R2).') 
parser.add_argument('-y', '--Clustmeth', type=str, metavar='', required=True, help='Clustering method to be used to identify communities in network. "fg" for fast-and-greedy, "eb" for edge-betweenness (slow).') 
parser.add_argument('-s', '--FocalSP', type=str, metavar='', required=True, help='The name of the focal species to represent each gene family (should exactly match the tip label of the species tree)') 


#Define the parser
args = parser.parse_args()

#Store arguments
JOBname=args.JOBname
BLmethod=args.BLmethod
Filterstat=args.Filterstat
Cutoff=args.Cutoff
Clustmeth=args.Clustmeth
FocalSP=args.FocalSP

'''
JOBname ="TEST"
BLmethod="bxb"
Filterstat="pval"
Cutoff=float(0.005)
Clustmeth="fg"
FocalSP="Atha"
'''
#Store output dir as a variable
out_dir= 'OUT_'+JOBname+'/'

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

#Run the R script
print("Beginning network analyses using the R package, igraph...\n\n Calling R...\n\n")

#Run the Network analyses.
#make command.
Net_cmd= 'Rscript Networks_and_stats.R '+JOBname+" "+BLmethod+" "+Filterstat+" "+str(Cutoff)+" "+Clustmeth+" "+FocalSP
    
#Run the command (if it contains strings expected in the command, this is a precautin of using shell=True)
if re.search('Networks_and_stats.R', Net_cmd) and re.search('Rscript', Net_cmd):
    print("Calling R with the following command:")
    print(Net_cmd)
    subprocess.call(Net_cmd, shell=True)

if len(glob.glob(out_dir+'Network_analyses/*pdf')) > 0:
    print('Finished network analyses.\n\nResults files written to Network_analyses/\n\n')
else:
    print('Something went wrong with network analyses...\n\n')

