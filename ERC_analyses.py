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


#During developent, set working directory:
#working_dir = '/Users/esforsythe/Documents/Work/Bioinformatics/ERC_networks/Analysis/ERCnet_dev/'
#os.chdir('/Users/esforsythe/Documents/Work/Bioinformatics/ERC_networks/Analysis/ERCnet_dev/')

#At runtime set working directory to the place where the script lives
working_dir = sys.path[0]+'/' 
os.chdir(working_dir)

#Set up an argumanet parser
parser = argparse.ArgumentParser(description='ERC step')

parser.add_argument('-j', '--JOBname', type=str, metavar='', required=True, help='Unique job name for this run of ERCnet. Avoid including spaces or special characters ("_" is ok)') 
parser.add_argument('-m', '--Mult_threads', type=int, metavar='', required=False, default=1, help='Integer: number of threads avilable for parallel computing (default = 1)' )
parser.add_argument('-s', '--FocalSP', type=str, metavar='', required=True, help='The name of the focal species to represent each gene family (should exactly match the tip label of the species tree)') 
parser.add_argument('-M', '--Meta_stats',action='store_true', required=False, help='The type of report of metadata from ERC correlations you want (default = False)')
parser.add_argument('-b', '--branchMethod', type=str, metavar='', required=False, default="R2T", help='This determins which branch reconcilliation method to use. Enter either "BXB" for branch by branch or "R2T" root to tip.')
parser.add_argument('-t', '--test', type=int, metavar='', required=False, default=0, help='Tests ERC_analyses.py on a small, randomized subset of the data equalling the integer provied.') 

#Define the parser
args = parser.parse_args()

#Store arguments
JOBname=args.JOBname
Mult_threads=args.Mult_threads
FocalSP=args.FocalSP
Meta_stats=args.Meta_stats
branchMethod=args.branchMethod
testNum = args.test
#JOBname = "Clptest"
#Mult_threads = 1
#FocalSP="A_thaliana_prot"
#Meta_stats="hits"

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

#Make results file
with open(out_dir+'ERC_results/'+str(fileName), "a") as f:
    f.write("GeneA_HOG" + "\t" + "GeneA_ID" + "\t" + "GeneB_HOG" + "\t" + "GeneB_ID" + "\t" + 
            "Overlapping_branches" + "\t" + "Slope" + "\t" + "P_R2" + "\t" + "P_Pval" + '\t' + 
            'S_R2' + '\t' + 'S_Pval' + '\t' + 'P_FDR_Corrected_Pval' + '\t' + 'S_FDR_Corrected_Pval') 

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


#Make functions for parralellizing correlation analyses
def iterate_corr(pairwise_combos):
    for combo in [pairwise_combos]:
        return combo


#Make function for performing actual correlations
def par_corr(i, j):
    #Get the test gene HOG IDs
    geneA=i
    geneB=j

    #geneA="HOG0001953"
    #geneB="HOG0008729"

    #Get the gene ids for the test genes

    # geneA_ID=gene_fams._get_value(list(gene_fams['HOG']).index(geneA), FocalSP)
    # geneB_ID=gene_fams._get_value(list(gene_fams['HOG']).index(geneB), FocalSP)

    ###Could sorting help with this? Should test total time it takes to find both genes as the program completes. 
    rowA = int(np.where(gene_fams['HOG'] == geneA)[0])
    rowB = int(np.where(gene_fams['HOG'] == geneB)[0])

    geneA_ID=gene_fams[FocalSP][rowA]
    geneB_ID=gene_fams[FocalSP][rowB]

    #Branch-by-branch branch lengths
    #subset dataframe to the rows of the two test HOGs
    test_df_temp = BLs[(BLs["HOG_ID"] == geneA) | (BLs["HOG_ID"] == geneB)]


    #transpose the dataframev
    test_df_t=test_df_temp.transpose()

    #Remove the HOG_ID row (first row) and remove rows with NaN
    test_df_clean=test_df_t.iloc[1: , :].dropna()


    #add names to columns
    test_df_clean.columns=['GeneA', 'GeneB']
    xLin = test_df_clean['GeneA']
    yLin = test_df_clean['GeneB']

    #Perform correlation tests if there are at least 5 datapoints
    if test_df_clean.shape[0]>4:

        #Stats order = (0:slope, 1:intercept, 2:pearson_r, 3:pearson_p, 4:stderr, 5:spearman_r, 6:spearman pvalue)
        #This is the old pearon command (used on plastid and methods paper)
        pear_stats=stats.linregress(x=list(test_df_clean['GeneA']), y=list(test_df_clean['GeneB']))

        #This is a new pearson command (used on parastic plant paper)
        pear_stats2=stats.pearsonr(x=test_df_clean['GeneA'], y=test_df_clean['GeneB'])

        #Get spearman stats
        spear_stats=stats.spearmanr(test_df_clean['GeneA'], test_df_clean['GeneB'])
        
        #Full file: "GeneA_HOG", "GeneA_ID",  "GeneB_HOG", "GeneB_ID", "Overlapping_branches", "Slope", "P_R2", "P_Pval", 'S_R2', 'S_Pval', 'P_FDR_Corrected_Pval', 'S_FDR_Corrected_Pval'
        #To be added here: "Overlapping_branches", "Slope", "P_R2", "P_Pval", 'S_R2', 'S_Pval'

        #pear_stats2=stats.pearsonr(x=test_df_clean['GeneA'], y=test_df_clean['GeneB'])
        
        ## TESTING
        '''
        with open(out_dir+'ERC_results/'+str(fileName.replace("ERC_results_", "PEARSON_TEST_")), "a") as f:
                f.write('\n'+ str(geneA) +'\t'+ str(geneA_ID) +'\t'+ str(geneB) +'\t'+ str(geneB_ID) +'\t'+str(pear_stats[2]**2)+'\t'+str(pear_stats2[0]**2)+'\t'+str(abs(pear_stats[2]**2-pear_stats2[0]**2)))
                #print("Appended: " + str(geneA_ID) + " by " + str(geneB_ID))
        '''
        
        #print(f"{pear_stats[2]**2}__{pear_stats2[0]**2}__{abs(pear_stats[2]**2-pear_stats2[0]**2)}")

        #spear_stats2=stats.spearmanr(test_df_clean['GeneA'], test_df_clean['GeneB'])
        
        #print(f"Spearman Pvalue: OG is {spear_stats[1]} and new is {spear_stats2[1]}")


        #For testing
        #results_str= str(test_df_clean.shape[0]) +'\t'+ "NA" +'\t'+ str(pear_stats[0]**2) +'\t'+ str(pear_stats[1]) + '\t' + str(spear_stats[0]**2) + '\t' + str(spear_stats[1])

        #Combine results into a long sting (with tabs seperating values)
        #Old pearson stats indexing
        #results_str= str(test_df_clean.shape[0]) +'\t'+ str(pear_stats[0]) +'\t'+ str(pear_stats[2]**2) +'\t'+ str(pear_stats[3]) + '\t' + str(spear_stats[0]**2) + '\t' + str(spear_stats[1])
        
        #new pearson stats indexing
        results_str= str(test_df_clean.shape[0]) +'\t'+ str(pear_stats[0]) +'\t'+ str(pear_stats2[0]**2) +'\t'+ str(pear_stats2[1]) + '\t' + str(spear_stats[0]**2) + '\t' + str(spear_stats[1])

    else:
        results_str='nan\tnan\tnan\tnan\tnan\tnan'
        #bxb_results_str='NA\tNA\tNA\tNA\tNA\tNA'
    

    #Only write the results if there's some indication of correlation (this keeps the size of the file from inflating)
    if (not results_str.split("\t")[0] == "nan"):
        #write (append) to results file
            with open(out_dir+'ERC_results/'+str(fileName), "a") as f:
                f.write('\n'+ str(geneA) +'\t'+ str(geneA_ID) +'\t'+ str(geneB) +'\t'+ str(geneB_ID) +'\t'+ results_str)
                #print("Appended: " + str(geneA_ID) + " by " + str(geneB_ID))
    
#Timestamp just before Parallel Call
erc.benchmarkTime(bench_fileName, out_dir + 'benchmark/', 'start', 'Correlation', timer)

#Run the correlation analysis (in paralell)
Parallel(n_jobs = int(Mult_threads), max_nbytes=None)(delayed(par_corr)(i,j) for i,j in iterate_corr(pairwise_combos))

#Timestamp just after Parallel Call
erc.benchmarkTime(bench_fileName, out_dir + 'benchmark/', 'end', 'Correlation', timer) 

#Calculate total time of parralel process & lines per minute
num_lines = len(pd.read_csv(out_dir + 'ERC_results/' + str(fileName), sep='\t'))
erc.benchmarkProcess(out_dir + 'benchmark/' + str(bench_fileName), num_lines)

#Progress message
print("Done with all-by-all comparisons...\n")

###Use R to create summary figures of the ERC result

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

