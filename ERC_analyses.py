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
parser.add_argument('-M', '--Meta_stats', type=str, metavar='', required=False, default="none", help='The type of report of metadata from ERC correlations you want. "none" for no reports, "full" for report including all pairwise-stats (very slow/memory-intensive for large dataset), "hits" for trimmed down report of potential hits (default = none)')


#Define the parser
args = parser.parse_args()

#Store arguments
JOBname=args.JOBname
Mult_threads=args.Mult_threads
FocalSP=args.FocalSP
Meta_stats=args.Meta_stats
#JOBname = "Clptest"
#Mult_threads = 1
#FocalSP="A_thaliana_prot"
#Meta_stats="hits"

#Store output dir as a variable
out_dir= 'OUT_'+JOBname+'/'

#Check whether a valid Meta_stats was entered
if not ((Meta_stats == "none") | (Meta_stats == "full") | (Meta_stats == "hits")):
    print("Error: invalid entry for --Meta_stats/-M. Quitting...")
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
    
#remove previous potential hits file (if it exists)
if len(glob.glob(out_dir+'ERC_results/ERC_results_potential_hits.tsv')) > 0:
    os.remove(out_dir+'ERC_results/ERC_results_potential_hits.tsv')

#Make results file
with open(out_dir+'ERC_results/ERC_results.tsv', "a") as f:
    f.write("GeneA_HOG"+"\t"+"GeneA_ID"+"\t"+"GeneB_HOG"+"\t"+"GeneB_ID"+"\t"+"Overlapping_branches_BXB"+"\t"+"Slope_BXB"+"\t"+"Pearson_R2_BXB"+"\t"+"Pearson_P_BXB"+"\t"+"Spearman_R2_BXB"+"\t"+"Spearman_P_BXB"+"\t"+"Overlapping_branches_R2T"+"\t"+"Slope_R2T"+"\t"+"Pearson_R2_R2T"+"\t"+"Pearson_P_R2T"+"\t"+"Spearman_R2_R2T"+"\t"+"Spearman_P_R2T")

#Make results file
with open(out_dir+'ERC_results/ERC_results_potential_hits.tsv', "a") as f:
    f.write("GeneA_HOG"+"\t"+"GeneA_ID"+"\t"+"GeneB_HOG"+"\t"+"GeneB_ID"+"\t"+"Overlapping_branches_BXB"+"\t"+"Slope_BXB"+"\t"+"Pearson_R2_BXB"+"\t"+"Pearson_P_BXB"+"\t"+"Spearman_R2_BXB"+"\t"+"Spearman_P_BXB"+"\t"+"Overlapping_branches_R2T"+"\t"+"Slope_R2T"+"\t"+"Pearson_R2_R2T"+"\t"+"Pearson_P_R2T"+"\t"+"Spearman_R2_R2T"+"\t"+"Spearman_P_R2T")


#Read in BL results
r2t_BLs=pd.read_csv(out_dir+'BL_results/r2t_BLs_normalized.tsv', sep='\t')
bxb_BLs=pd.read_csv(out_dir+'BL_results/bxb_BLs_normalized.tsv', sep='\t')

#Read in gene IDs
gene_fams=pd.read_csv(out_dir+'Filtered_genefam_dataset.csv')

#Clean up the gene IDs dataframe
#Find out what string needs to be removed from HOGs and use replace() to remove it
gene_fams['HOG']=gene_fams['HOG'].str.replace(str(gene_fams['HOG'][0].split(".", 1)[0]+"."), "")


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
    
    #geneA="HOG0001953"
    #geneB="HOG0008729"
    
    #Get the gene ids for the test genes
    geneA_ID=gene_fams._get_value(list(gene_fams['HOG']).index(geneA), FocalSP)
    geneB_ID=gene_fams._get_value(list(gene_fams['HOG']).index(geneB), FocalSP)
    

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
        bxb_results_str='nan\tnan\tnan\tnan\tnan\tnan'
        #bxb_results_str='NA\tNA\tNA\tNA\tNA\tNA'
    
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
        r2t_results_str='nan\tnan\tnan\tnan\tnan\tnan'
        #r2t_results_str='NA\tNA\tNA\tNA\tNA\tNA'
    
    #Only write the results if there's some indication of correlation (this keeps the size of the file from inflating)
    if ((not bxb_results_str.split("\t")[0] == "nan") and (not r2t_results_str.split("\t")[0] == "nan")):
        #write (append) to results file
        with open(out_dir+'ERC_results/ERC_results.tsv', "a") as f:
            f.write('\n'+ str(geneA) +'\t'+ str(geneA_ID) +'\t'+ str(geneB) +'\t'+ str(geneB_ID) +'\t'+ bxb_results_str +'\t'+ r2t_results_str)
        
        if (float(bxb_results_str.split("\t")[1])>0 or float(r2t_results_str.split("\t")[1])>0) and (float(bxb_results_str.split("\t")[3])<0.05 or float(bxb_results_str.split("\t")[5])<0.05 or float(r2t_results_str.split("\t")[3])<0.05 or float(r2t_results_str.split("\t")[5])<0.05):
            #write (append) to results file (for the potential hits)
            with open(out_dir+'ERC_results/ERC_results_potential_hits.tsv', "a") as f:
                f.write('\n'+ str(geneA) +'\t'+ str(geneA_ID) +'\t'+ str(geneB) +'\t'+ str(geneB_ID) +'\t'+ bxb_results_str +'\t'+ r2t_results_str)


#Run the correlation analysis (in paralell)
Parallel(n_jobs = int(Mult_threads))(delayed(par_corr)(i,j) for i,j in iterate_corr(pairwise_combos))

#Progress message
print("Done with all-by-all comparisons...\n")

###Use R to create summary figures of the ERC results

if not (Meta_stats == "none"):
    
    print("Generating plots describing correlation statistic metadata")
    
    #Run the R script
    ERC_rec_cmd= 'Rscript Corr_meta_stats.R '+JOBname+' '+Meta_stats
        
    #Run the command (if it contains strings expected in the command, this is a precautin of using shell=True)
    if re.search('Corr_meta_stats.R', ERC_rec_cmd) and re.search('Rscript', ERC_rec_cmd):
        print("Creating summary figs in R with the following command:")
        print(ERC_rec_cmd)
        subprocess.call(ERC_rec_cmd, shell=True)

#Report status
if len(glob.glob(out_dir+'ERC_results/*tsv')) > 0:
    print('Finished ERC.\n\nResults files written to ERC_results/\n\n')
else:
    print('Something went wrong with ERC analyses...\n')

#Print directions for next steps
print('ERC correlation analyses finished. Exiting....\n\n'\
      'To perform network analyses, run Network_analyses.py. See associated help menu for required arguments.'\
          '\nExample command:\n\n' \
          './Network_analyses.py -j '+JOBname+' -m bxb -f pval -c 0.001 -y fg -s <focal_sp>\n\n' 
      )


