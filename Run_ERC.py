'''
Script for performing BL reconcilation. Needs to be run in python 2 env (because dlcpar requires python 2)

conda activate dlcpar_py27 

python Run_ERC.py -o /Users/esforsythe/Documents/Work/Bioinformatics/ERC_networks/Analysis/Orthofinder/Results_Oct15/

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
parser = argparse.ArgumentParser(description='Main ERCnet script')

parser.add_argument('-o', '--OFpath', type=str, metavar='', required=True, help='Full path to the Orthofinder results dir (should contain Species_Tree/, Phylogenetic_Hierarchical_Orthogroups/ etc...)\n Include "/" at the end of the string') 

#Define the parser
args = parser.parse_args()

#Store arguments
OFpath=args.OFpath
#OFpath = "/Users/esforsythe/Documents/Work/Bioinformatics/ERC_networks/Analysis/Orthofinder/Results_Oct15/"

print("beginning BL reconciliation (using R and dlcpar)...\n\n Calling R...\n\n")

#Run the BL reconcilation step.
BL_rec_cmd= 'Rscript BL_reconciliation.R '+OFpath
    
#Run the command (if it contains strings expected in the command, this is a precautin of using shell=True)
if re.search('BL_reconciliation.R', BL_rec_cmd) and re.search('Rscript', BL_rec_cmd):
    subprocess.call(BL_rec_cmd, shell=True)

if len(glob.glob('BL_results/*tsv')) > 0:
    print('Finished branch length reconciliation.\n\nResults files written to BL_results/\n\n')







