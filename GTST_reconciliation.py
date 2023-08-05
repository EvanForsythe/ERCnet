#!/usr/bin/env python2

'''
### Main script for running the GTST reconciliation step in DLCpar ###

conda activate dlcpar_py27

Example command:
    python GTST_reconciliation.py -j TPC_test

'''

#During developent, set working directory:
#import os
#working_dir = '/Users/esforsythe/Documents/Work/Bioinformatics/ERC_networks/Analysis/ERCnet_dev/'
#os.chdir('/Users/esforsythe/Documents/Work/Bioinformatics/ERC_networks/Analysis/ERCnet_dev/')

#Storebought modules
import os
import sys
import glob
import argparse
import subprocess
import time
from datetime import datetime 

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
parser = argparse.ArgumentParser(description='Script to run DLCpar')

parser.add_argument('-j', '--JOBname', type=str, metavar='', required=True, help='Unique job name for this run of ERCnet. Avoid including spaces or special characters ("_" is ok)') 

#Define the parser
args = parser.parse_args()

#Store arguments
JOBname=args.JOBname
#JOBname="TPC_test"

#Store the output dir as a variable
out_dir= 'OUT_'+JOBname+'/'

#Check if dlcpar_search is avialable
proc = subprocess.Popen(['which', 'dlcpar_search'], stdout=subprocess.PIPE)
output = str(proc.stdout.read())

#Define the time object and folder for optimization testing
bench_fileName = JOBname + '_GTST_rec_benchmark.tsv'
timer = time.localtime()
current_time = time.strftime("%H:%M:%S", timer)

benchmarkTime(bench_fileName, out_dir + 'benchmark/', 'start', 'GTST_reconciliation', timer)

if output != 0:
    print('dlcpar appears to be installed.')
else:
    print('could not find dlcpar_search program. Quitting...')
    sys.exit()

#Check if directory exists and contains files
if os.path.isdir(out_dir+'DLCpar/') and len(glob.glob(out_dir+'DLCpar/*NODES_BL.txt'))>0:
    print('Beginning reconciliation using input files located in:')
    print(out_dir+'DLCpar/\n\n')
    #Change directory to the DLCpar folder (this is required to call DLCpar)
    os.chdir(out_dir+'DLCpar/')
else: 
    print('Could not find DLCpar input files. Quitting....')
    sys.exit()

#Get the input files
input_GTs=glob.glob('*NODES_BL.txt')

#Run reconciliation on all gene trees
for file_i, file in enumerate(input_GTs):
    #Run the command
    subprocess.call(["dlcpar_search", "-sSpeciesTree_rooted_node_labels.txt", "-SspeciesIDs.smap", file]) #adding spaces after -s and -S were causing problems
    
    #Report progress
    if file_i % 10 == 0:
        print('%d reconciliations done!' %file_i)

#Change back to original directory
os.chdir(working_dir)

benchmarkTime(bench_fileName, out_dir + 'benchmark/', 'end', 'GTST_reconciliation', timer)



print('IMPORTANT NOTE: for the next step, switch back to the python 3 environement you used in the Phylogenomics steps of the workflow.\n' \
          'For a reminder of the avialable environents:\n' \
              'conda info --envs\n')

print('After switching to python3 environment, run the next step with the following command:\n\n./ERC_analyses.py -j '+JOBname+' -m 1 -s <focal_sp> -b R2T \n\n')
