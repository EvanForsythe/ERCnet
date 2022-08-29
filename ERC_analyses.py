#!/usr/bin/env python3
'''
Script for performing BL reconcilation and ERC correlation analyses. This should be run in python3

'''

#Import modules
import os
import re
import sys
import glob
import subprocess
import argparse
import itertools
import pandas as pd
import scipy.stats as stats
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

#Define the parser
args = parser.parse_args()

#Store arguments
JOBname=args.JOBname
Mult_threads=args.Mult_threads
#JOBname = "testBIG"
#Mult_threads = 8

#Store output dir as a variable
out_dir= 'OUT_'+JOBname+'/'

### Run the BL_reconciliation
print("beginning BL reconciliation in R...\n\n Calling R...\n\n")

#Run the BL reconcilation step.
BL_rec_cmd= 'Rscript Branch_length_reconciliation.R '+JOBname
    
#Run the command (if it contains strings expected in the command, this is a precautin of using shell=True)
if re.search('Branch_length_reconciliation.R', BL_rec_cmd) and re.search('Rscript', BL_rec_cmd):
    print("Running BL reconciliation in R with the folllowing command:")
    print(BL_rec_cmd)
    subprocess.call(BL_rec_cmd, shell=True)

if len(glob.glob(out_dir+'BL_results/*tsv')) > 0:
    print('Finished branch length reconciliation.\n\nResults files written to BL_results/\n\n')

### Run all-by-all correlations

#Make a directory for ERC results
if not os.path.isdir(out_dir+'ERC_results/'):
    os.makedirs(out_dir+'ERC_results/')
    print("created folder: ERC_results/\n\n")
else: 
    print('ERC results will be stored to ERC_results/\n\n')
    
#remove previous file (if it exists)
if len(glob.glob(out_dir+'ERC_results/ERC_results.tsv')) > 0:
    os.remove(out_dir+'ERC_results/ERC_results.tsv')

#Make results file
with open(out_dir+'ERC_results/ERC_results.tsv', "a") as f:
    f.write("GeneA_HOG"+"\t"+"GeneB_HOG"+"\t"+"Overlapping_branches_BXB"+"\t"+"Slope_BXB"+"\t"+"Pearson_R2_BXB"+"\t"+"Pearson_P_BXB"+"\t"+"Spearman_R2_BXB"+"\t"+"Spearman_P_BXB"+"\t"+"Overlapping_branches_R2T"+"\t"+"Slope_R2T"+"\t"+"Pearson_R2_R2T"+"\t"+"Pearson_P_R2T"+"\t"+"Spearman_R2_R2T"+"\t"+"Spearman_P_R2T")

#Read in BL results
r2t_BLs=pd.read_csv(out_dir+'BL_results/r2t_BLs_normalized.tsv', sep='\t')
bxb_BLs=pd.read_csv(out_dir+'BL_results/bxb_BLs_normalized.tsv', sep='\t')

#Make list of pairwise combinations
pairwise_combos=list(itertools.combinations(r2t_BLs['HOG_ID'], 2))


#Make functions for parralellizing correlation analyses
def iterate_corr(pairwise_combos):
    for combo in [pairwise_combos]:
        return combo


#Make function for performing actual correlations
def par_corr(i, j):
    
    #Get the test gene HOG IDs
    geneA=i
    geneB=j

    #Branch-by-branch branch lengths
    #subset dataframe to the rows of the two test HOGs
    bxb_test_df_temp = bxb_BLs[(bxb_BLs["HOG_ID"] == geneA) | (bxb_BLs["HOG_ID"] == geneB)]
    
    #transpose the dataframev
    bxb_test_df_t=bxb_test_df_temp.transpose()
    
    #Remove the HOG_ID row (first row) and remove rows with NaN
    bxb_test_df_clean=bxb_test_df_t.iloc[1: , :].dropna()
    
    #add names to columns
    bxb_test_df_clean.columns=['GeneA', 'GeneB']
    
    #Perform correlation tests if there are at least 3 datapoints
    if bxb_test_df_clean.shape[0]>2:
    
        #Stats order = (0:slope, 1:intercept, 2:pearson_r, 3:pearson_p, 4:stderr, 5:spearman_r, 6:spearman pvalue)
        bxb_corr_stats=stats.linregress(x=list(bxb_test_df_clean['GeneA']), y=list(bxb_test_df_clean['GeneB']))+stats.spearmanr(bxb_test_df_clean['GeneA'], bxb_test_df_clean['GeneB'])
        
        #Create string
        bxb_results_str= str(bxb_test_df_clean.shape[0]) +'\t'+ str(bxb_corr_stats[0]) +'\t'+ str(bxb_corr_stats[2]**2) +'\t'+ str(bxb_corr_stats[3]) +'\t'+ str(bxb_corr_stats[5]**2) +'\t'+ str(bxb_corr_stats[6])
        #(note r is squared to get r2)
    else:
        bxb_results_str='NA\tNA\tNA\tNA\tNA\tNA'
    
    #Root to tip branch lengths
    #subset dataframe to the rows of the two test HOGs
    r2t_test_df_temp = r2t_BLs[(r2t_BLs["HOG_ID"] == geneA) | (r2t_BLs["HOG_ID"] == geneB)]
    
    #transpose the dataframe
    r2t_test_df_t=r2t_test_df_temp.transpose()
    
    #Remove the HOG_ID row (first row) and remove rows with NaN
    r2t_test_df_clean=r2t_test_df_t.iloc[1: , :].dropna()
    
    
    #add names to columns
    r2t_test_df_clean.columns=['GeneA', 'GeneB']
    
    #Perform correlation tests if there are at least 3 datapoints
    if r2t_test_df_clean.shape[0]>2:
    
        #Stats order = (0:slope, 1:intercept, 2:pearson_r, 3:pearson_p, 4:stderr, 5:spearman_r, 6:spearman pvalue)
        r2t_corr_stats=stats.linregress(x=list(r2t_test_df_clean['GeneA']), y=list(r2t_test_df_clean['GeneB']))+stats.spearmanr(r2t_test_df_clean['GeneA'], r2t_test_df_clean['GeneB'])
        
        #Create string
        r2t_results_str= str(r2t_test_df_clean.shape[0]) +'\t'+ str(r2t_corr_stats[0]) +'\t'+ str(r2t_corr_stats[2]**2) +'\t'+ str(r2t_corr_stats[3]) +'\t'+ str(r2t_corr_stats[5]**2) +'\t'+ str(r2t_corr_stats[6])
    else:
        r2t_results_str='NA\tNA\tNA\tNA\tNA\tNA'
    
    #write (append) to results file
    with open(out_dir+'ERC_results/ERC_results.tsv', "a") as f:
        f.write('\n'+ geneA +'\t'+ geneB +'\t'+ bxb_results_str +'\t'+ r2t_results_str)


#Run the correlation analysis (in paralell)
Parallel(n_jobs = int(Mult_threads))(delayed(par_corr)(i,j) for i,j in iterate_corr(pairwise_combos))

#Progress message
print("Done with all-by-all comparisons...\n")

###Use R to create summary figures of the ERC results

#Run the R script
ERC_rec_cmd= 'Rscript AllxAll_correlations.R '+JOBname
    
#Run the command (if it contains strings expected in the command, this is a precautin of using shell=True)
if re.search('AllxAll_correlations.R', ERC_rec_cmd) and re.search('Rscript', ERC_rec_cmd):
    print("Creating summary figs in R with the following command:")
    print(ERC_rec_cmd)
    subprocess.call(ERC_rec_cmd, shell=True)

#Report status
if len(glob.glob(out_dir+'ERC_results/*tsv')) > 0:
    print('Finished ERC.\n\nResults files written to ERC_results/\n\n')
else:
    print('Something went wrong with ERC analyses...\n')

#Print wrapup statement
print("ERC correlation analyses finished. Exiting....")

#Print directions for next steps
print('ERC correlation analyses finished. Exiting....\n\n'\
      'To perform network analyses, run Network_analyses.py. See associated help menu for required arguments.'\
          '\nExample command:\n\n' \
          './Network_analyses.py -j '+JOBname+' -m bxb -f pval -c 0.001 -y fg -s <focal_sp>\n\n' 
      )


