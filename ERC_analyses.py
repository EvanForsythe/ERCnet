'''
Script for performing BL reconcilation and ERC correlation analyses. This should be run in python3

conda activate ERC_networks

python Run_ERC.py -j TEST -o /Users/esforsythe/Documents/Work/Bioinformatics/ERC_networks/Analysis/Orthofinder/Results_Oct15/

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
parser = argparse.ArgumentParser(description='ERC step')

parser.add_argument('-j', '--JOBname', type=str, metavar='', required=True, help='Unique job name for this run of ERCnet. Avoid including spaces or special characters ("_" is ok)') 

#Define the parser
args = parser.parse_args()

#Store arguments
JOBname=args.JOBname
#JOBname = "TEST"

#Store output dir as a variable
out_dir= 'OUT_'+JOBname+'/'

## Run the BL_reconciliation
print("beginning BL reconciliation in R...\n\n Calling R...\n\n")

#Run the BL reconcilation step.
BL_rec_cmd= 'Rscript Branch_length_reconciliation.R '+JOBname
    
#Run the command (if it contains strings expected in the command, this is a precautin of using shell=True)
if re.search('Branch_length_reconciliation.R', BL_rec_cmd) and re.search('Rscript', BL_rec_cmd):
    print("Running BL reconciliation in R with the folllowing command:")
    print(BL_rec_cmd)
    subprocess.call(BL_rec_cmd, shell=True)

if len(glob.glob('BL_results/*tsv')) > 0:
    print('Finished branch length reconciliation.\n\nResults files written to BL_results/\n\n')
    
## Run all-by-all correlations

#Make a directory for ERC results
if not os.path.isdir(out_dir+'ERC_results/'):
    os.makedirs(out_dir+'ERC_results/')
    print("created folder: ERC_results/\n\n")
else: 
    print('ERC results will be stored to ERC_results/\n\n')
    

#Run the R script
#Run the BL reconcilation step.
ERC_rec_cmd= 'Rscript AllxAll_correlations.R '+JOBname
    
#Run the command (if it contains strings expected in the command, this is a precautin of using shell=True)
if re.search('AllxAll_correlations.R', ERC_rec_cmd) and re.search('Rscript', ERC_rec_cmd):
    print("Running ERC comparisons with the following command:")
    print(ERC_rec_cmd)
    subprocess.call(ERC_rec_cmd, shell=True)

#Report status
if len(glob.glob(out_dir+'ERC_results/*tsv')) > 0:
    print('Finished ERC.\n\nResults files written to ERC_results/\n\n')
else:
    print('Something went wrong with ERC analyses...\n')

#Print wrapup statement
print("ERC correlation analyses finished. Exiting....")


print('ERC correlation analyses finished. Exiting....\n\n'\
      'To perform network analyses, run Network_analyses.py. See associated help menu for required arguments.'\
          '\nExample command:\n\n' \
          'python Network_analyses.py -j '+JOBname+' -m bxb -f pval -c 0.001 -y fg -s <focal_sp>\n' 
      )













