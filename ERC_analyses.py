#!/usr/bin/env python3
'''
Script for performing BL reconcilation and ERC correlation analyses. This should be run in python3

'''

#Import modules
import os
import re
import sys
import math
import glob
import argparse
import itertools
import subprocess
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
parser.add_argument('-s', '--FocalSP', type=str, metavar='', required=True, help='The name of the focal species to represent each gene family (should exactly match the tip label of the species tree)') 
parser.add_argument('-M', '--Meta_stats',action='store_true', required=False, help='The type of report of metadata from ERC correlations you want (default = False)')
parser.add_argument('-b', '--branchMethod', type=str, metavar='', required=False, default="R2T", help='This determins which branch reconcilliation method to use. Enter either "BXB" for branch by branch or "R2T" root to tip.')
parser.add_argument('-c', '--corrMethod', type=str, metavar='', required=False, default='spearman', help='This determines which statistical correlation method that will be used. Enter either "spearman" or "pearson".')

#Define the parser
args = parser.parse_args()

#Store arguments
JOBname=args.JOBname
Mult_threads=args.Mult_threads
FocalSP=args.FocalSP
Meta_stats=args.Meta_stats
branchMethod=args.branchMethod
corrMethod=args.corrMethod
#JOBname = "Clptest"
#Mult_threads = 1
#FocalSP="A_thaliana_prot"
#Meta_stats="hits"

#Store output dir as a variable
out_dir= 'OUT_'+JOBname+'/'
fileName = str('ERC_results_' + branchMethod + '_' + corrMethod + '.tsv')

print(str(branchMethod) + ' chosen for branch method.')
print(str(corrMethod) + ' chosen for statistical inference.')

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
if len(glob.glob(out_dir+'ERC_results/'+fileName)) > 0:
    os.remove(out_dir+'ERC_results/'+fileName)

#Make results file
with open(out_dir+'ERC_results/'+str(fileName), "a") as f:
    f.write("GeneA_HOG" + "\t" + "GeneA_ID" + "\t" + "GeneB_HOG" + "\t" + "GeneB_ID" + "\t" + "Overlapping_branches" + "\t" + "Slope" + "\t" + "R2" + "\t" + "Pval") 

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


#Make list of pairwise combinations
pairwise_combos=list(itertools.combinations(BLs['HOG_ID'], 2))


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
    geneA_ID=gene_fams._get_value(list(gene_fams['HOG']).index(geneA), FocalSP)
    geneB_ID=gene_fams._get_value(list(gene_fams['HOG']).index(geneB), FocalSP)
    

    #Branch-by-branch branch lengths
    #subset dataframe to the rows of the two test HOGs
    test_df_temp = BLs[(BLs["HOG_ID"] == geneA) | (BLs["HOG_ID"] == geneB)]
    
    #transpose the dataframev
    test_df_t=test_df_temp.transpose()
    
    #Remove the HOG_ID row (first row) and remove rows with NaN
    test_df_clean=test_df_t.iloc[1: , :].dropna()
    
    #add names to columns
    test_df_clean.columns=['GeneA', 'GeneB']
    
    #Perform correlation tests if there are at least 3 datapoints
    if test_df_clean.shape[0]>2:
            
        if (corrMethod == 'pearson'):
            #Stats order = (0:slope, 1:intercept, 2:pearson_r, 3:pearson_p, 4:stderr, 5:spearman_r, 6:spearman pvalue)
            corr_stats=stats.linregress(x=list(test_df_clean['GeneA']), y=list(test_df_clean['GeneB']))
        else:
            corr_stats=stats.spearmanr(test_df_clean['GeneA'], test_df_clean['GeneB'])

        print(corr_stats)
        #Create string
        results_str= str(test_df_clean.shape[0]) +'\t'+ str(corr_stats[0]) +'\t'+ str(corr_stats[2]**2) +'\t'+ str(corr_stats[3])
        #(note r is squared to get r2)
    else:
        results_str='nan\tnan\tnan\tnan\tnan\tnan'
        #bxb_results_str='NA\tNA\tNA\tNA\tNA\tNA'
    
    #Only write the results if there's some indication of correlation (this keeps the size of the file from inflating)
    if (not results_str.split("\t")[0] == "nan"):
        #write (append) to results file
        with open(out_dir+'ERC_results/'+str(fileName), "a") as f:
            f.write('\n'+ str(geneA) +'\t'+ str(geneA_ID) +'\t'+ str(geneB) +'\t'+ str(geneB_ID) +'\t'+ results_str)

#Run the correlation analysis (in paralell)
Parallel(n_jobs = int(Mult_threads))(delayed(par_corr)(i,j) for i,j in iterate_corr(pairwise_combos))

#Progress message
print("Done with all-by-all comparisons...\n")

###Use R to create summary figures of the ERC results


#Could be redundant now. Come back later!! 
#if (Meta_stats):
    
#    print("Generating plots describing correlation statistic metadata")
    
    #Run the R script
#    ERC_rec_cmd= 'Rscript Corr_meta_stats.R '+JOBname+' '+Meta_stats+' '+fileName
        
    #Run the command (if it contains strings expected in the command, this is a precautin of using shell=True)
#    if re.search('Corr_meta_stats.R', ERC_rec_cmd) and re.search('Rscript', ERC_rec_cmd):
#        print("Creating summary figs in R with the following command:")
#        print(ERC_rec_cmd)
#        subprocess.call(ERC_rec_cmd, shell=True)

#Report status
if len(glob.glob(out_dir+'ERC_results/'+fileName)) > 0:
    print('Finished ERC.\n\nResults files written to ERC_results/\n\n')
else:
    print('Something went wrong with ERC analyses...\n')

#Print directions for next steps
print('ERC correlation analyses finished. Exiting....\n\n'\
      'To perform network analyses, run Network_analyses.py. See associated help menu for required arguments.'\
          '\nExample command:\n\n' \
          './Network_analyses.py -j '+JOBname+' -m r2t -c spearman -y fg -s '+FocalSP+'\n\n' 
      )


