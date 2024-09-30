#!/usr/bin/env python3

### Main script for running the phylogenomics steps of the ERC workflow###

#Storebought modules
import os
import re
import sys
import glob
import time
import math
import shutil
import argparse
import itertools
import subprocess
import numpy as np
import pandas as pd
from Bio import Phylo
from Bio import SeqIO
from Bio import AlignIO
from datetime import datetime 
from joblib import Parallel, delayed


#Homemade modules
from filterHOGs import make_seq_counts_df, filter_gene_fams, check_alns_for_prune, file_loss_log, bs_reps_check_log
from ERC_functions import benchmarkTime, benchmarkProcess, CheckFileExists, CheckFileNonEmpty, file_counts_log, file_line_count_log, resolve_polytomies


#At runtime set working directory to the place where the script lives
working_dir = sys.path[0]+'/' 
os.chdir(working_dir)

#Set up an argumanet parser
parser = argparse.ArgumentParser(description='Script for running the first step of ERCnet')

parser.add_argument('-j', '--JOBname', type=str, metavar='', required=True, help='Unique job name for this run of ERCnet. Avoid including spaces or special characters ("_" is ok)') 
parser.add_argument('-o', '--OFpath', type=str, metavar='', required=True, help='Full path to the Orthofinder results dir (should contain Species_Tree/, Phylogenetic_Hierarchical_Orthogroups/ etc...)\n Include "/" at the end of the string') 
parser.add_argument('-p', '--MaxP', type=int, metavar='', required=False, default=3, help='Integer: maximum number of paralogs per species allowed in each gene family (default = 3)' )
parser.add_argument('-r', '--MinR', type=int, metavar='', required=False, default=0, help='Integer: minimum number of species represented required in each gene family (default = half of species)' )
parser.add_argument('-t', '--Test_num', type=int, metavar='', required=False, help='Integer: number of gene families to analyze (for testing only)' )
parser.add_argument('-e','--explore_filters', action='store_true', required=False, help='Add this flag to explore filtering options (if selected, program will quit without running downstream steps)')
parser.add_argument('-l', '--Min_len', type=int, metavar='', required=False, default=100, help='Integer: minimum length (amino acid sites) of alignment (after trimming with Gblocks) required to retain gene (default = 100)' )
parser.add_argument('-x', '--Rax_dir', type=str, metavar='', required=True, help='Full path to the location of your raxml install (use which raxmlHPC to locate). Include "/" at the end of the string') 
parser.add_argument('-s','--SPmap', action='store_true', required=False, help='Add this flag to provide a custom species mapping file. This mapping file must be formatted in certian way. See instuctions')
parser.add_argument('-n', '--Node', type=int, metavar='', required=False, help='Integer: node number on orthofinder species tree to be used to obtain HOGs (default = 1)' )
parser.add_argument('-m', '--Mult_threads', type=int, metavar='', required=False, default=1, help='Integer: number of threads avilable for parallel computing (default = 1)' )
parser.add_argument('-a','--Apriori', action='store_true', required=False, help='Add this flag to provide an apriori set of genes to analyze. The input file listing those genes must be formatted in certian way. See instuctions')
parser.add_argument('-c', '--core_distribution', type=int, metavar='', required=False, default=3, help='Integer: Sets the core distribution group which affects number of cores split between front end and back end paralellization of raxml bootstrapping')
parser.add_argument('-P', '--Prune_cutoff', type=float, metavar='', required=False, default=0.9, help='Float: prune seqs from alignments if the proportion of gap sites exceeds this number (default: 0.9)')
parser.add_argument('-T', '--Taper', type=str, metavar='', required=False, default="no", help='Run TAPER trimming of alignments? If selected, the user must include full path to installation of julia (should end in "bin/" (default=no)')


#Define the parser
args = parser.parse_args()

#Store arguments
JOBname=args.JOBname
OFpath=args.OFpath
MaxP_val=args.MaxP
MinR_val=args.MinR
explore_filters=args.explore_filters
Min_len=args.Min_len
Test_num=args.Test_num
Rax_dir=args.Rax_dir
SPmap=args.SPmap
Node=args.Node
Mult_threads=args.Mult_threads
Apriori=args.Apriori
core_dist=args.core_distribution
prune_cutoff=args.Prune_cutoff
taper=args.Taper


#Check to see if the path arguments end with "/"
if not OFpath.endswith('/'):
    print('OFpath does not end with a "/". Quitting...\n')
    sys.exit()

if not Rax_dir.endswith('/'):
    print('Rax_dir does not end with a "/". Quitting...\n')
    sys.exit()

print("Finding the Orthofinder HOG file to be used to designate ingroup subtrees\n")

#Get the status of the variable from the --Node argument
if Node is None:
    print("WARNDING: --Node (-n) not defined. Choosing N1.tsv HOF file by default\nNote that N1.tsv is only appropriate if your species tree contains a single outgroup species.\n")
    Node=1
else:
    if isinstance(Node, int):
        print("--Node (-n) argument used to select N"+str(Node)+".tsv HOG file\n")
    else:
        print("ERROR: Invalid value for --Node (-n). This should be an interger. (e.g. for N2.tsv: '-n 2'')\nQUITTING....\n")
        sys.exit()

#Check availability of taper and julia
if not taper == "no":
    print("-T selected. Checking for TAPER dependencies...")
    if taper.endswith("/"):
        julia_msg= str(subprocess.Popen([taper+'julia', '-v'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate())
        if re.search('version', julia_msg):
            print("Found julia installation!")
            taper_msg=str(subprocess.Popen([taper+'julia', 'correction_multi.jl', '-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate())
            if re.search('fasta', taper_msg):
                print("Found TAPER file (correction_multi.jl)!")
            else:
                print("ERROR: Couldn't find TAPER file (correction_multi.jl). This should be located in the main ERCnet dir. QUITTING...\n")
        else:
            print("ERROR: Couldn't find installation of julia")
            sys.exit()
            
    else:
        print('ERROR: -T argument should be followed by full path to julia installation. Path should end with "/". QUITTING...')
        sys.exit()
else:
    print("Skipping TAPER trimming...")


#Get path to needed Orthofinder files 
sp_tr_path = OFpath+'Species_Tree/SpeciesTree_rooted_node_labels.txt'
HOG_file_path = OFpath+'Phylogenetic_Hierarchical_Orthogroups/N'+str(Node)+'.tsv'
OG_trees_dir = OFpath+'Resolved_Gene_Trees/'
OGseqdir = OFpath+'Orthogroup_Sequences/'

##Check if expected Othofinder files exist
print("Checking if expected input files exist (output file from Orthofinder)...\n")

#Species tree
if os.path.isfile(sp_tr_path):
    print("Found species tree!")
else:
    print("ERROR: species tree not found! Quitting....\n")
    sys.exit()

#HOG file
if os.path.isfile(HOG_file_path):
    print("Found HOG file!")
else:
    print("ERROR: expected HOG file (N"+str(Node)+".tsv) not found! Quitting....\n")
    sys.exit()

#Resolved gene trees
if os.path.isdir(OG_trees_dir):
    print("Found Resolved_Gene_Trees/ directory!")
else:
    print("ERROR: Resolved_Gene_Trees/ directory not found! This should be in the Orthofinder results folder. Check the version of Orthofinder you ran. Quitting....\n")
    sys.exit()

#Orthogroup seuquences
if os.path.isdir(OGseqdir):
    print("Found Orthogroup_Sequences/ directory!\n")
else:
    print("ERROR: Orthogroup_Sequences/ directory not found! This should be in the Orthofinder results folder. Check the version of Orthofinder you ran. Quitting....\n")
    sys.exit()
    
#Create a variable for the output dir
out_dir= 'OUT_'+JOBname+'/'

#Create the output folder
#Make a directory for storing results
if os.path.isdir(out_dir):
    while True:
        user_input = input("This jobname already exists. Would you like to overwrite? (y/n) \n")
        if user_input == 'y':
            print("Clearing contents of " + out_dir + " All output files will be written to this folder\n") 
            shutil.rmtree(out_dir)
            os.makedirs(out_dir)
            break
        if user_input == 'n':
            print('Unique jobname required. Exiting...')
            sys.exit()
            break
        else:
            print("Command not recognized.\n")     

else:
    os.makedirs(out_dir) 
    print('created folder: '+out_dir+'\nAll output files will be written to this folder\n')

if Mult_threads < 4:
    print("Parallel processing is not viable for tree building below 4 cores due to overhead. ERCnet will continue as a linear process for Raxml runs.")

#Created dropped file log
print("Creating Dropped_gene_log.csv file to store information about genes lost during fitering/analysis steps.\n")
dropped_log_handle = open(out_dir+"Dropped_gene_log.csv", "a")
dropped_log_handle.write("HOG,Step,Reason\n")
dropped_log_handle.close()

#Created dropped file log
print("Creating Run_data_counts_log.csv file to store quantitative information from throughout the run.\n")
counts_log_handle = open(out_dir+"Run_data_counts_log.csv", "a")
counts_log_handle.write("Jobname,File_or_folder,Count\n")
counts_log_handle.close()

#Define the time object and folder for optimization testing
bench_fileName = JOBname + '_Phylogenomics_benchmark.tsv'
timer = time.localtime()
current_time = time.strftime("%H:%M:%S", timer)

if not os.path.isdir(out_dir + 'benchmark/'):
    os.makedirs(out_dir + 'benchmark/')
    print("Created benchmarking folder for optimization testing")
else:
    if len(glob.glob(out_dir + 'benchmark/' + str(bench_fileName))) > 0:
        os.remove(out_dir + 'benchmark/' + str(bench_fileName))
    print("ERC_analyses.py benchmark will be stored in ERC_benchmark/ Previous benchmark logs removed.\n\n")

#Adds first timestamp for process start.
with open(out_dir + 'benchmark/' + str(bench_fileName), "a") as bench:
    bench.write("Job Name: " + JOBname + '\t' + "Cores: " + str(Mult_threads) + '\n') 
    bench.write("Stage" + '\t' + "Time" + '\n')
    bench.write("Process Start" + '\t' + str(current_time) + '\n')

#Check HOG file exists 
if os.path.isfile(HOG_file_path):
    print('HOG file for selected node located at:\n'+ HOG_file_path+'\n')
    
else:
    print('OG file for selected node not found. Quitting...\n')
    sys.exit()

#open HOG file
HOG_file=pd.read_csv(HOG_file_path, sep='\t')

#Make a list of species
sp_list_temp = list(HOG_file.columns[3:])

#Generate species mapping file
if SPmap:
    if os.path.isfile("Species_mapping.csv"):
        print('reading user-provided species map. Creating copy in output directory.')
        mapping_table=pd.read_csv("Species_mapping.csv", sep=',')

    else:
        print('Could not find Species_mapping.csv. Quitting...')
        sys.exit()
else:
    print('Attempting to create Species_mapping.csv in output directory.')
    mapping_table = pd.DataFrame(data={'Prefix': sp_list_temp, 'SpeciesID': sp_list_temp})

mapping_table.to_csv(out_dir+'Species_mapping.csv', sep=',' , index=False)

sp_prefix=list(mapping_table['Prefix'])
sp_names=list(mapping_table['SpeciesID'])

print("The following species identifiers will be used: \n"+str(sp_names))
print("These should match the column headers in the HOG file (e.g. N1.tsv) and the tip labels in the species tree")

print("\nThe following species identifiers will be used: \n"+str(sp_prefix))
print("These should be present in the sequence IDs in alignments etc...\n")

#Count number of species
file_line_count_log(JOBname, out_dir+"Species_mapping.csv", out_dir+"Run_data_counts_log.csv")

#Generate a dataframe of the counts data using the my outside module from filterHOGs.py
seq_counts_df=make_seq_counts_df(HOG_file_path, out_dir+'Species_mapping.csv')

#Verify species mapping has been done correctly. 
numSpecies = seq_counts_df.columns[3:]
runningSum = 0

print('Species Mapping Table:')

for spec in numSpecies:
    speciesSum = sum(seq_counts_df[spec])
    runningSum += speciesSum
    print(str(speciesSum) + ' \ttotal genes mapped to \t' + str(spec) + ' genome.')
print('\n')

if (runningSum < 1):
    print('An error has occured in the species mapping table. Please verify species names are consistent. ERCnet will now exit.')
    sys.exit()

#Check if the dataframe was assinged properly
if 'HOG' in list(seq_counts_df.columns):
    print('Sequence counts per species dataframe successfully generated\n')
    print('Writing sequence counts per species table to "Seq_counts_per_species.csv"\n')
    seq_counts_df.to_csv(out_dir+'Seq_counts_per_species.csv', sep=',' , index=False)
else:
    print('ERROR: Sequence counts per sepecies dataframe not properly generated. Quitting...\n')
    sys.exit()


#Count the total number og HOG (before filtering)
file_line_count_log(JOBname, out_dir+"Seq_counts_per_species.csv", out_dir+"Run_data_counts_log.csv")

#If the -e flag was chosen, run the parameter scan to explore the data filtering options
if explore_filters:
    print('-e / --explore_filters flag chosen\nExploring filtering options and quitting\n')
    
    #Set filter ranges
    max_paralogs_vals=list(range(1, 5, 1))
    min_rep_vals=list(range(math.floor(len(sp_names)/2), (len(sp_names)+1), 1))
    
    #Make blank array for parameter scan
    retained_trees=np.zeros((len(max_paralogs_vals), len(min_rep_vals)))
    
    #Start nested loop to fill in the matrix
    for p_count, p_val in enumerate (max_paralogs_vals):
        for s_count, s_val in enumerate (min_rep_vals):
            para_value=p_val
            rep_value=s_val

            passed_only_df = filter_gene_fams(HOG_file, seq_counts_df, sp_names, para_value, rep_value)
            #Number of rows in the dataframe
            retained_trees[p_count,s_count]=passed_only_df.shape[0]
            
    #Write the csv file
    #Convert
    retained_trees_df=pd.DataFrame(retained_trees)
    
    #Set values as integers
    retained_trees_df = retained_trees_df.apply(pd.to_numeric, errors='coerce', downcast='integer')
    
    #Assign new col names
    retained_trees_df.set_axis(max_paralogs_vals, axis=0, inplace=True)
    
    #assign new row names
    retained_trees_df.set_axis(min_rep_vals, axis=1, inplace=True)

    #Print the table
    print('\n\nOutputting table of the total retained trees under different filtering parameter combinations to "Retained_genes_counts.csv"\n')
    retained_trees_df.to_csv(out_dir+'Retained_genes_counts.csv', sep=',' , index=False)
    print('Number of retained genes:')
    print(retained_trees_df)
    print('\nColumns: minimum requred number or species represented\nRows:Maximum allowed number of paralogs/species\n')
    
    print('Parameter exporation complete. To run full analyses, choose -p and -r parameters based on the table above and rerun Phylogenomics.py without the --explore_filters / -e flag.')
    #Exit program (the user will need to rerun after viewing the parameter scan results)
    sys.exit()
    
#Report the filters chosen by user
if (MinR_val == 0):
    print('A value for --MinR was not provided. ERCnet defaults to half the number of species provided in mapping table.')
    MinR_val = max(math.floor(len(sp_names)/2),1)
print('--MaxP_val set to {} and --MinR_val set to {}\nFiltering data....'.format(MaxP_val,MinR_val))

#Use the module from Filter_stats.py to filter the list
Keeper_HOGs_df= filter_gene_fams(out_dir, HOG_file, seq_counts_df, sp_names, MaxP_val, MinR_val)

#Check if any HOGs were retained
if Keeper_HOGs_df.shape[0] > 0:
    print("Dataset filtering successful\n")
else:
    print('ERROR: It appears that zero HOGs were retained. Adjust --MaxP_val and --MinR_val and rerun...\nQuitting....')
    sys.exit()

#Add a filter that can search for an apriori gene set.
if bool(Apriori):
    print('A priori gene dataset option selected (-a/--Apriori). Retrieving genes of interest...\n')
    if isinstance(Test_num, int):
        print('Warning: -t option is also selected. -t could cause problems when used with the -a option. Proceed with caution...\n')
    
    #Read in the file
    apriori_df=pd.read_csv("A_priori_genes.csv", sep=',')
    
    col_header= str(apriori_df.keys()[0])
    
    gene_array=list(apriori_df[col_header])
    
    Keeper_HOGs_df=Keeper_HOGs_df[Keeper_HOGs_df[col_header].str.contains('|'.join(gene_array), case=False, na=False)]
    
    if Keeper_HOGs_df.shape[0] > 0:
        print("Dataset filtering for a priori genes successful\n")
        if len(gene_array) == Keeper_HOGs_df.shape[0]:
            print("It looks like all of the genes you were searching for were present!\n")
        elif len(gene_array) > Keeper_HOGs_df.shape[0]:
            print("It looks like "+str(len(gene_array)-Keeper_HOGs_df.shape[0])+" of your a priori genes were not found")
        else:
            print("ERROR: Something went worng with a priori filter. The dataframe after filtering has more rows than expected. Quitting...\n")
            sys.exit()
        
    else:
        print('ERROR: It appears that zero HOGs were retained. Check to make sure at least some of the exact sequence ID strings are contained in the HOG file. \nQuitting....')
        sys.exit()

#For testing the performance of ERCnet on a small number of genes subset the dataset to a user-defined subset
if isinstance(Test_num, int):
    Keeper_HOGs_df = Keeper_HOGs_df.sample(n=Test_num)

#Write csv file and report the number of gene trees included after filtering
if 'HOG' in list(Keeper_HOGs_df.columns):
    Keeper_HOGs_df.to_csv(out_dir+'Filtered_genefam_dataset.csv', sep=',' , index=False)
    print('Filtered list of gene families successfully generated. \n\nFiltered dataset contains %d genes.\n\nCsv file written to Filtered_genefam_dataset.csv \n' %Keeper_HOGs_df.shape[0])
else:
    print('Error: sequence counts per sepecies dataframe not properly generated. Quitting...\n')
    sys.exit()

#Count the number of HOGs after filtering
file_line_count_log(JOBname, out_dir+"Filtered_genefam_dataset.csv", out_dir+"Run_data_counts_log.csv")

### Test sequence files to look for seq IDs that would cause problems later on.

drop_list=list()

#Loop through all rows
for row_i, row in Keeper_HOGs_df.iterrows():

    #Get the relevant OG and HOG string 
    OGtemp = row['OG']
    HOGtemp = row['HOG']
    
    #Get a list of seqs to retain
    seq_list_temp=[item for item in list(row[sp_names]) if not(pd.isnull(item)) == True]
    
    #Clean up the list by converting to str, removing all the extra stuff then split the str back into a list
    seq_list = str(seq_list_temp).replace(" ", "").replace("'", "").replace("[", "").replace("]", "").split(',')
    
    #open and read fasta file into a dictionary
    OG_handle = open(str(OGseqdir + OGtemp +'.fa'), "r")

    #Loop through the line in the file
    seqID_list_temp=list()
    
    for line in OG_handle:
        if line.startswith(">"):
            id_temp = line.strip() #Removes "\n"
            id_clean = id_temp.replace(">", "") #Removes ">" by replacing with nothing.
            seqID_list_temp.append(id_clean)

    if len(seqID_list_temp) > len(set(seqID_list_temp)):
        print(f"WARNING: found duplicate identical seq IDs in an OG file. The offending OG was {OGtemp} and HOG was {HOGtemp}\n Removing this file from the analysis")
        drop_list.append(OGtemp)
    
    long_seq_counter = 0 
    for i in seqID_list_temp:
        if len(i) > 73:
            long_seq_counter += 1
    
    if long_seq_counter > 0:
        print(f"WARNING: found seq ID(s) in an OG file that are longer than 73 characters. These would create errors in Gblocks later on. The offending OG was {OGtemp} and HOG was {HOGtemp}\n Removing this file from the analysis")
        drop_list.append(OGtemp)

#Remove the any sequences on the drop list
if len(drop_list)>0:
    Keeper_HOGs_df = Keeper_HOGs_df.loc[~Keeper_HOGs_df['OG'].isin(drop_list)]

print(f"Number of HOGs remaining after seq ID formatting problem check: {Keeper_HOGs_df.shape[0]}")

### Extract the subtree sequences to be aligned
#make a dir to store the new fasta files
if not os.path.isdir(out_dir+'HOG_seqs'):
    os.makedirs(out_dir+'HOG_seqs')
    print("created folder : HOG_seqs\n\nfasta files will be written to HOG_seqs/\n")
else: print('Fasta files will be written to HOG_seqs/\n')

#Loop through all rows
for row_i, row in Keeper_HOGs_df.iterrows():

    #Get the relevant OG and HOG string 
    OGtemp = row['OG']
    HOGtemp = row['HOG']
    
    #Get a list of seqs to retain
    seq_list_temp=[item for item in list(row[sp_names]) if not(pd.isnull(item)) == True]
    
    #Clean up the list by converting to str, removing all the extra stuff then split the str back into a list
    seq_list = str(seq_list_temp).replace(" ", "").replace("'", "").replace("[", "").replace("]", "").split(',')
    
    #open and read fasta file into a dictionary
    OG_dict = SeqIO.to_dict(SeqIO.parse(str(OGseqdir + OGtemp +'.fa'), 'fasta'))

    #Subset the dictionary (Turns out this is easy with a one liner)
    HOG_dict= {k: OG_dict[k] for k in OG_dict.keys() & seq_list}
    
    #Write file
    with open(str(out_dir+'HOG_seqs/'+HOGtemp.replace("N"+str(Node)+".", "")+'.fa'), 'w') as handle:
        SeqIO.write(HOG_dict.values(), handle, 'fasta')
    
#Get list of files that were written
seq_file_names = glob.glob(out_dir+'HOG_seqs/HOG*')

#Report the number of files that were written
#Write csv file and report the number of gene trees included after filtering
if len(seq_file_names)>0:
    print('%d HOG sequence files written to HOG_seqs/\n' %int(len(seq_file_names)))
else:
    print('ERROR: writing HOG sequence files failed. Quitting...\n')
    sys.exit()

#Count number of HOG seqs
file_counts_log(JOBname, out_dir, 'HOG_seqs', 'HOG', 'fa', out_dir+'Run_data_counts_log.csv')
    
### Run mafft alignment
mafft_msg= str(subprocess.Popen(['mafft', '--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate())

if re.search('MAFFT', mafft_msg):
    print('Beginning multiple sequence alignment of %d files with MAFFT\n' %len(seq_file_names))
    print('Using MAFFT version %s \n' %str(re.findall(r'MAFFT\s*v\s*\d*.\d*', mafft_msg)))
else:
    print('ERROR: MAFFT not avaialable. Check that mafft is installed and added to $PATH. Quitting... \n')
    sys.exit()
    
#Make a directory for alignments
if not os.path.isdir(out_dir+'Alns'):
    os.makedirs(out_dir+'Alns')
    print("created folder : Alns/\n")
else: 
    print('MAFFT alignments will be written to Alns/\n')

print('Number of alignments finished: ')
##Begin paralellization of alignments

def iterate_maf(seq_file_names):
    for file in [seq_file_names]:
        return file

def par_maf_alns(file):

    temp_out_file_path = str(out_dir+'Alns/ALN_'+file.replace(out_dir+"HOG_seqs/", ""))

    os.system('mafft-linsi --quiet '+file+' > '+temp_out_file_path)
    #os.system('mafft-linsi '+file+' >Alns/ALN_'+file.replace("HOG_seqs/", "")+' 2>&1') #' 2>&1' suppressed stderr from mafft
    #print('mafft-linsi --quiet '+file+' > Alns/ALN_'+file.replace("HOG_seqs/", ""))

    if not CheckFileExists(temp_out_file_path):
        print(f"WARNING: mafft alignment wasn't able to create the file, {temp_out_file_path}. This gene family will be removed from further analysis")
        os.remove(temp_out_file_path)
    
    if not CheckFileNonEmpty(temp_out_file_path):
        print(f"WARNING: mafft alignment yielded an empty file, {temp_out_file_path}. This gene family will be removed from further analysis")
        os.remove(temp_out_file_path)

    ##End paralellization of alignments

#Timestamp just before Parallel Call
benchmarkTime(bench_fileName, out_dir + 'benchmark/', 'start', 'maf_alns', timer)

Parallel(n_jobs = Mult_threads, verbose=0, max_nbytes=None)(delayed(par_maf_alns)(file) for file in iterate_maf(seq_file_names))

#Timestamp just after Parallel Call
benchmarkTime(bench_fileName, out_dir + 'benchmark/', 'end', 'maf_alns', timer)

#Calculate total time of parralel process & lines per minute
num_items = len(seq_file_names)
benchmarkProcess(out_dir + 'benchmark/' + str(bench_fileName), num_items)


#Check alignment status
if len(glob.glob(out_dir+'Alns/ALN*')) == len(seq_file_names):
    print('\n\nDone with multiple sequences alignments\n')
elif len(glob.glob(out_dir+'Alns/ALN*')) < len(seq_file_names):
    print('\nWARNING: some alignments did not finish. Proceeding (with caution)...\n')
elif len(glob.glob(out_dir+'Alns/ALN*')) < len(seq_file_names):
    print('\nWARNING: there are more alignments in the folder than expected. Proceeding (with caution)...\n')

#Check for dropped files
file_loss_log(out_dir, out_dir+'HOG_seqs/', out_dir+'Alns/', '.fa', 'ALN', 'multiple sequence alignment', 'file unsuccesfully aligned for reasons unknown')

#Count number of alns
file_counts_log(JOBname, out_dir, 'Alns', 'ALN', 'fa', out_dir+'Run_data_counts_log.csv')

### Running TAPER Alignment trimming (optionally)
if not taper == "no":
    print("beginning TAPER trimming.")
    
    #Make a directory for TAPER trimmed alignments
    if not os.path.isdir(out_dir+'TAPER_Alns'):
        os.makedirs(out_dir+'TAPER_Alns')
        print("created folder : TAPER_Alns/\n")
    else: 
        print('Trimmed alignments will be written to TAPER_Alns/\n')
    
    #Get list of all alns (that haven't been gblocked/tapered yet)
    aln_file_names = [x for x in glob.glob(out_dir+'Alns/ALN*') if "-gb" not in x]
    
    def iterate_taper(aln_file_names):
        for file in [aln_file_names]:
            return file
        
    def par_taper_trim(file):
        #Build the command used to call taper
        #taper_cmd= taper+'julia correction_multi.jl -m "-" '+file+' > '+file.replace("Alns/", "TAPER_Alns/")
        
        # Build the command without output redirection
        taper_cmd = [taper + 'julia', 'correction_multi.jl', '-m-', file]

        #Run the command
        result = subprocess.run(taper_cmd, input="", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False, text=True)
      
        taper_stdout=result.stdout
        taper_stderr=result.stderr
        
        if not re.search("Version 1.0.0", taper_stderr):
            print(f"Something went wrong with taper on this file: {file}")
            print(f"The following command was responsible: {taper_cmd}")
            print(f"the following error message was generated: {taper_stderr}")
            print("Quitting....\n")
            sys.exit()
        
        else:
            # Open the output file
            with open(file.replace("Alns/", "TAPER_Alns/"), 'w') as output_file:
                # Call the subprocess, redirecting stdout to the output file
                output_file.write(taper_stdout)
            
        #Check to make sure taper worked
        if os.stat(file.replace("Alns/", "TAPER_Alns/")).st_size == 0:
            print("ERROR: something went wrong with the TAPER command. Produced empty file. The following command caused the issue:")
            print(taper_cmd)
            print("Quitting....\n")
            sys.exit()

    Parallel(n_jobs = Mult_threads, verbose=0, max_nbytes=None)(delayed(par_taper_trim)(file) for file in iterate_taper(aln_file_names))
    
    #Check alignment status
    if len(glob.glob(out_dir+'TAPER_Alns/ALN*')) == len(aln_file_names):
        print('\n\nDone with TAPER trimming\n')
    elif len(glob.glob(out_dir+'TAPER_Alns/ALN*')) < len(aln_file_names):
        print('\nWARNING: some TAPER timming did not finish. Proceeding (with caution)...\n')
    elif len(glob.glob(out_dir+'TAPER_Alns/ALN*')) < len(aln_file_names):
        print('\nWARNING: there are more TAPER trimmed alignments in the folder than expected. Proceeding (with caution)...\n')

#Check for dropped files
if not taper == "no":
    file_loss_log(out_dir, out_dir+'Alns/', out_dir+'TAPER_Alns/', 'ALN', 'ALN', 'TAPER trimming', 'file unsuccessfully trimmed with TAPER')
    #Count number of alns
    file_counts_log(JOBname, out_dir, 'TAPER_Alns', 'ALN', 'fa', out_dir+'Run_data_counts_log.csv')


###GBLOCKS cleaning alignments
print('Beginning GBLOCKS\n')

if not taper == "no":
    #Get list of all alns (that haven't been gblocked yet)
    aln_file_names2 = [x for x in glob.glob(out_dir+'TAPER_Alns/ALN*') if "-gb" not in x]

else:    
    #Get list of all alns (that haven't been gblocked yet)
    aln_file_names2 = [x for x in glob.glob(out_dir+'Alns/ALN*') if "-gb" not in x]


#Gblocks
if 'Min_len' in locals():
    print('Beginning gblocks trimming...\n')
else:
    print('Min_len parameter for gblocks filtering not specified. Quitting...\n')

#Make folder for gblocks'd alignments
#Make a directory for gblocks trimmed alignments
if not os.path.isdir(out_dir+'Gb_alns'):
    os.makedirs(out_dir+'Gb_alns')
    print("created folder: Gb_alns/\n")
else: 
    print('Gblocks-trimmed alignments will be stored to Gb_alns/\n')

#Make a directory for gblocks trimmed alignments that are too short 
if not os.path.isdir(out_dir+'Gb_alns/Too_short/'):
    os.makedirs(out_dir+'Gb_alns/Too_short/')
    print("created folder: Gb_alns/Too_short/\n")
else: 
    print('Trimmed alignents that are too short will be stored to Gb_alns/Too_short/\n')

#Loop through the files that need to Gblocks trimmed
##Begin paralellization of gblocks trimming

def iterate_alns(aln_file_names2):
    for aln in [aln_file_names2]:
        return aln

def par_gblocks(aln):

    #Get number of sequences in alignment
    #open connection
    aln_file = open(aln, "r")
    line_count = 0
    for line in aln_file:
        if re.search(r">", line):
            line_count += 1
    #Close the connection
    aln_file.close()
    
    #Create gblocks command
    gblocks_cmd = 'Gblocks '+aln+' -b5=h -b4=5 -p=n -b2='+str(math.floor((line_count/2)+1))
    
    if re.search('Gblocks', gblocks_cmd) and re.search('-b5', gblocks_cmd) and re.search('-p=n', gblocks_cmd):
        #Run cmd and capture std output (in order to know the trimmed alignment length)
        proc = subprocess.Popen(gblocks_cmd, stdout=subprocess.PIPE, shell=True) #apparently shell=True can create 'security issues' (so I put it under an if statement to make sure the cmd is what it should be)
        output = str(proc.stdout.read())
        
        #Extract the length of the trimmed alignment by parsing Gblocks stdout 
        if re.search('Gblocks alignment', output): #putting this in an if statement to make sure it's parsing the right thing

            #Regular expression to pull out the trimmed alignment length    
            match = re.search('\d+|$', re.search(r'(?<=Gblocks alignment:).*|$', output).group()).group() #note that "|$" at the end of both search strings means that match = '' if no hit is found
            
            if match != '':
               trm_aln_ln = int(match)
            else:
                print("ERROR: Problem parsing the Gblocks output for alignment: "+aln+" Something must have gone wrong with Gblocks.")
                print("To trouble shoot Gblocks, try running Gblocks from the command line with:")
                print(gblocks_cmd)
                print("\nQuitting...\n")
                print("Here is the message from Gblocks...")
                print(output)
                sys.exit()
           
            #Move the gblocks files to the appropriate folder
            if trm_aln_ln >= Min_len:
                if not taper == "no":
                    os.replace(str(aln+'-gb'), str(str(aln).replace("ALN_", "GB_ALN_")).replace('TAPER_Alns/', 'Gb_alns/'))
                else:
                    os.replace(str(aln+'-gb'), str(str(aln).replace("ALN_", "GB_ALN_")).replace('Alns/', 'Gb_alns/'))
            else:
                if not taper == "no":
                    print('%s not long enough. Moving to gblocks files to Gb_alns/Too_short/' %aln.replace('TAPER_Alns/', ''))
                    os.replace(str(aln+'-gb'), str(str(aln).replace("ALN_", "GB_ALN_")).replace('TAPER_Alns/', 'Gb_alns/Too_short/'))
                else:
                    print('%s not long enough. Moving to gblocks files to Gb_alns/Too_short/' %aln.replace('Alns/', ''))
                    os.replace(str(aln+'-gb'), str(str(aln).replace("ALN_", "GB_ALN_")).replace('Alns/', 'Gb_alns/Too_short/'))
            
        else:
            print('ERROR: something wrong with Gblocks stdout. Quitting...\n')
            print("Here is the message from Gblocks...")
            print(output)
            sys.exit()

            
    else:
        print('WARNING: Something wrong with Gblocks command...\n')
        print("Here is the Gblocks command that was being attempted\n")
        print(output)
        sys.exit()


#Timestamp just before Parallel Call
benchmarkTime(bench_fileName, out_dir + 'benchmark/', 'start', 'gblocks', timer)

#run the parallel command
Parallel(n_jobs= Mult_threads, verbose=0, max_nbytes=None)(delayed(par_gblocks)(aln) for aln in iterate_alns(aln_file_names2))

#Timestamp just before Parallel Call
benchmarkTime(bench_fileName, out_dir + 'benchmark/', 'start', 'gblocks', timer)

#Calculate total time of parralel process & items per minute
num_items = len(aln_file_names2)
benchmarkProcess(out_dir + 'benchmark/' + str(bench_fileName), num_items)

if not taper == "no":
    file_loss_log(out_dir, out_dir+'TAPER_Alns/', out_dir+'Gb_alns/', 'ALN', 'GB_ALN', 'Gblocks trimming', 'file unsuccessfully Gblocks trimmed or became too short')
else:
    file_loss_log(out_dir, out_dir+'Alns/', out_dir+'Gb_alns/', 'ALN', 'GB_ALN', 'Gblocks trimming', 'file unsuccessfully Gblocks trimmed or became too short')

###Check if there are individual sequences that need to be pruned.
#Make a directory for storing pruning info
if not os.path.isdir(out_dir+'Aln_pruning/'):
    os.makedirs(out_dir+'Aln_pruning/')
    print("created folder: Aln_pruning/\n")
else: 
    print('Alignment pruning info will be stored in Aln_pruning/\n')

#get a list of all the Gblocks output files
alns_to_check=glob.glob(out_dir+'Gb_alns/GB_ALN_*fa')

#Loop through alignments to check and run the pruning function
for a in alns_to_check:
    #run the function (store output as variable os that it'll run)
    store_pruned=check_alns_for_prune(a, prune_cutoff, out_dir)

### Check if alignments still contain at least 4 sequences
#get a list of pruned fasta files
pruned_alns=glob.glob(out_dir+'Aln_pruning/GB_ALN_*fa')

#Make a ticker to keep track of how many files were pruned
prune_replace_ticker=0
prune_remove_ticker=0

for p in pruned_alns:
    #Get the OG
    OG_p=p.replace(out_dir+"Aln_pruning/GB_ALN_", "").replace(".fa", "")
    
    #Delete the original version of the file
    os.remove(out_dir+"Gb_alns/GB_ALN_"+OG_p+".fa")
    
    #Read in alignment
    aln_p = AlignIO.read(p, "fasta")

    #Check how many seqs remain in the alignment after pruning.
    #If at least 4, use the pruned version for downstream analyses, else drop the alignment from downstream analyses
    if len(aln_p) >= 4:
        #move the pruned file to the main GB folder
        os.replace(p, p.replace("Aln_pruning/", "Gb_alns/"))
        
        #Rename the log file to reflect whether the pruned version was retained
        os.replace(out_dir+'Aln_pruning/Prune_IDs_'+OG_p+".txt", out_dir+'Aln_pruning/Prune_IDs_'+OG_p+"_ALN_RETAINED.txt")
    
        prune_replace_ticker += 1
    else:
        #Rename the log file to reflect whether the pruned version was retained
        os.replace(out_dir+'Aln_pruning/Prune_IDs_'+OG_p+".txt", out_dir+'Aln_pruning/Prune_IDs_'+OG_p+"_ALN_DROPPED.txt")
        
        prune_remove_ticker += 1

print("Done pruning super gappy sequences from alignments...\n"+str(prune_replace_ticker)+" files replaced\n"+str(prune_remove_ticker)+" files removed because they were pruned too much\n")

#Count number of g-blocks alignments
file_counts_log(JOBname, out_dir, 'Gb_alns', 'GB_ALN', 'fa', out_dir+'Run_data_counts_log.csv')

#Count number of too-short g-blocks alignments
file_counts_log(JOBname, out_dir, 'Gb_alns/Too_short', 'GB_ALN', 'fa', out_dir+'Run_data_counts_log.csv')

#Count number of files that underwent pruning
file_counts_log(JOBname, out_dir, 'Aln_pruning', 'Prune_IDs_', 'txt', out_dir+'Run_data_counts_log.csv')

### Get the subtrees from orthofinder to use as a constraint tree from bootstrap scoring

#Create folder to store subtrees
if not os.path.isdir(out_dir+'HOG_subtrees/'):
    os.makedirs(out_dir+'HOG_subtrees/')
    print("\nCreated folder: HOG_subtrees/\n")
else: 
    print('HOG subtrees will be stored in HOG_subtrees/\n')

## Extract the subtree corresponding to the HOG from the larger gene tree for each OG

# Loop through all rows in the DataFrame (this replaces previous R code in Get_subtree.R)
for index, row in Keeper_HOGs_df.iterrows():
    # Get the name of the outfile
    HOG_temp = str(row['HOG'])

    # Remove the N1 (or N2, N3, etc.) string
    HOG_temp = re.split(r"\.", HOG_temp)[1]

    # Read in the OG tree
    tree_path = os.path.join(OG_trees_dir, f"{row['OG']}_tree.txt")
    full_tree_temp = Phylo.read(tree_path, "newick")

    # Get list of tips
    raw_strings = row.iloc[3:].tolist()
    raw_strings = [x for x in raw_strings if pd.notnull(x) and x != ""]

    # Split tips into a vector
    keeper_tips = re.split(r"\s*,\s*", ",".join(raw_strings).replace(" ", ""))

    # Check if this was a pruned alignment
    prune_file_pattern = f"{HOG_temp}_ALN_RETAINED.txt"
    prune_files = [f for f in os.listdir(os.path.join(out_dir, "Aln_pruning")) if re.match(prune_file_pattern, f)]
    
    if len(prune_files) == 1:
        prune_seqs = pd.read_table(os.path.join(out_dir, "Aln_pruning", prune_files[0]), header=None)[0].tolist()
        
        # Loop through seqs that need to be pruned from subtree
        keeper_tips = [tip for tip in keeper_tips if tip not in prune_seqs]

    # Get the subtree
    subtree_temp = full_tree_temp.common_ancestor(keeper_tips)

    # Add suffix
    out_file_name = f"{HOG_temp}_tree.txt"

    # Check if the subtree is binary (using `is_bifurcating` method in Biopython)
    if subtree_temp.is_bifurcating():
        # Write the subtree
        Phylo.write(subtree_temp, os.path.join(out_dir, "HOG_subtrees", out_file_name), "newick")
    else:
        # Resolve polytomies and make the subtree bifurcating
        resolved_subtree = resolve_polytomies(subtree_temp)
        
        # Check if the tree has been successfully resolved (optional step)
        if resolved_subtree.is_bifurcating():
            # Write the resolved subtree
            Phylo.write(resolved_subtree, os.path.join(out_dir, "HOG_subtrees", out_file_name), "newick")
        else:
            # If somehow still not bifurcating, log an error
            dropped_log_handle = open(out_dir+"Dropped_gene_log.csv", "a")
            dropped_log_handle.write(f"{HOG_temp},Subtree extraction,Subtree was non bifurcating and could not be resolved\n")

#Check if it worked
if len(glob.glob(out_dir+'HOG_subtrees/*tree.txt')) > 0:
    print('\nFinished running subtree extraction.\n')
else:
    print('ERROR: An error occured during subtree extraction. Quitting...\n')
    sys.exit()

#Check for dropped files
file_loss_log(out_dir, out_dir+'Gb_alns/', out_dir+'HOG_subtrees/', 'GB_ALN', '_tree.txt', 'subtree extraction', 'subtree extraction failed')

#Count number of subtrees
file_counts_log(JOBname, out_dir, 'HOG_subtrees', 'HOG', 'tree.txt', out_dir+'Run_data_counts_log.csv')

#### Perform bootstrap replication with raxml
print("Beginning bootstrapping with raxml.\n")

if Mult_threads == 1:
    print("Number of threads set to 1. Consider adding threads to parallelize if possible.\n")
    print("RAxML requires a minimum of 2 threads to run. ERCnet will not cancel, however this may overtax system resources and cause performance issues. Consider using a minimum of 2 threads.")
elif Mult_threads > 1:
    print("Will attempt to parallelize with "+str(Mult_threads)+" threads. Note that choosing more threads than are avilable can slow raxml performance.\n")
else:
    print("Number of threads argument does not appear to be an interger. Quitting....")
    sys.exit()

if not os.path.isdir(out_dir+'BS_reps/'):
    os.makedirs(out_dir+'BS_reps/')
    print("\nCreated folder: BS_reps/\n")
else: 
    print('Trees with bootstrap values will be stored in BS_reps/\n')

#Get a list of HOGs with retained alignments
HOGs2BS=[x.replace(out_dir+'Gb_alns/GB_ALN_', '').replace('.fa', '') for x in glob.glob(out_dir+'Gb_alns/GB_ALN_*')]

print("Raxml Parallel Group Chosen: " + str(core_dist) + "\n")
 
if (Mult_threads >= 4):
	if(core_dist == 1):
		Rax_front_cores = 2
		Rax_back_cores = int(Mult_threads/2)
	elif(core_dist == 2):
		Rax_front_cores = int(Mult_threads/4)
		Rax_back_cores = 4
	else:
		Rax_front_cores = int(Mult_threads/2)
		Rax_back_cores = 2
else:
	Rax_front_cores = 1
	Rax_back_cores = int(Mult_threads)

if (Rax_back_cores <= 1 and Rax_front_cores >1):
    Rax_front_cores = int(Mult_threads/2)
    
print('Number of trees finished: \n') 

#define rax filepath
file_path = working_dir + out_dir

# Iterate through HOGs for parallel processing
def iterate_HOGS(HOGs2BS):
    for HOG in HOGs2BS:
        yield HOG

# Parallelized RAxML bootstrap inference
def par_raxml_bootstrap(HOG_id, Rax_dir, file_path, cores): 
    # Build the bootstrap command
    raxml_BS_cmd = Rax_dir + 'raxmlHPC-PTHREADS-SSE3 -s ' + file_path + 'Gb_alns/GB_ALN_' + HOG_id + '.fa ' + \
                   '-w ' + file_path + 'BS_reps/ -n ' + HOG_id + '_BS.txt' + \
                   ' -m PROTGAMMALGF -p 12345 -x 12345 -# 100 -T ' + str(cores)

    # Run the command
    if re.search('raxmlHPC', raxml_BS_cmd) and re.search('PROTGAMMALGF', raxml_BS_cmd) and re.search('12345', raxml_BS_cmd):
        subprocess.call(raxml_BS_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

# Parallelized execution of bootstrap inference
def run_parallel_raxml_bootstrap(HOGs2BS, Rax_dir, file_path, Rax_front_cores, Rax_back_cores):
    # Start the parallel RAxML bootstrap inference
    benchmarkTime(bench_fileName, out_dir + 'benchmark/', 'start', 'raxml_bootstrap', timer)
    
    Parallel(n_jobs=Rax_front_cores, verbose=0, max_nbytes=None)(
        delayed(par_raxml_bootstrap)(HOG, Rax_dir, file_path, Rax_back_cores) 
        for HOG in iterate_HOGS(HOGs2BS)
    )
    
    benchmarkTime(bench_fileName, out_dir + 'benchmark/', 'end', 'raxml_bootstrap', timer)

    # Calculate total time and process count
    num_items = len(HOGs2BS)
    benchmarkProcess(out_dir + 'benchmark/' + str(bench_fileName), num_items)
    
    print("Bootstrap tree inference finished.\n")

# Call the function that run raxml bootstrapping
run_parallel_raxml_bootstrap(HOGs2BS, Rax_dir, file_path, Rax_front_cores, Rax_back_cores)

# Check for dropped files and BS replicates that didn't get reach 100 reps
bs_reps_check_log(
    full_out_dir=out_dir, 
    start_dir=out_dir+'Gb_alns/', 
    end_dir=out_dir+'BS_reps/', 
    start_pattern='GB_ALN_', 
    end_pattern='RAxML_bootstrap.', 
    step='Raxml BS replication', 
    error_message='File missing', 
    warning_message='File too short'
)

## Non-parallelized RAxML bootstrap mapping
# Create folder
if not os.path.isdir(out_dir+'BS_trees/'):
    os.makedirs(out_dir+'BS_trees/')
    print("\nCreated folder: BS_trees/\n")
else: 
    print('Trees with bootstrap values will be stored in BS_trees/\n')

def raxml_bootstrap_mapping(HOGs2map, Rax_dir, file_path):
    for HOG_id in HOGs2map:
        # Build the bootstrap mapping command
        raxml_BSmap_cmd = (
            Rax_dir + 'raxmlHPC-PTHREADS-SSE3 -z ' + file_path + 'BS_reps/RAxML_bootstrap.' + HOG_id + '_BS.txt ' +
            '-w ' + file_path + 'BS_trees/ -n ' + HOG_id + '_BS.txt' +
            ' -t ' + file_path + 'HOG_subtrees/' + HOG_id + '_tree.txt -f b -m PROTGAMMALGF -T 1'
        )

        # Run the command
        if re.search('raxmlHPC', raxml_BSmap_cmd):
            subprocess.call(raxml_BSmap_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    print("Bootstrap mapping finished.\n")

# Call the non-parallelized step for bootstrap mapping
raxml_bootstrap_mapping(HOGs2BS, Rax_dir, file_path)

#Check for dropped files
file_loss_log(out_dir, out_dir+'BS_reps/', out_dir+'BS_trees/', 'RAxML_bootstrap.', 'RAxML_bipartitions.', 'Bootstrap score mapping', 'Failed to produce tree file with bootstrap values')

#Count number of subtrees
file_counts_log(JOBname, out_dir, 'BS_trees', 'RAxML_bipartitions.', 'txt', out_dir+'Run_data_counts_log.csv')

## Do the root/rearrange step with Treerecs (GT/ST reconciliation)
#Copy the species tree into outdir
shutil.copy2(sp_tr_path, out_dir)

#Read in the species tree file (as text file)
with open(str(out_dir+'SpeciesTree_rooted_node_labels.txt'), 'r') as file:
  sptr_file = file.read()

#get the names in the speices tree 
#Loop through all rows in the species mapping table
for row_i, row in mapping_table.iterrows():
    
    #Get the relevant OG and HOG string 
    SpeciesID = row['SpeciesID']
    Prefix = row['Prefix']
    
    sptr_file=sptr_file.replace(SpeciesID, Prefix)
    
#Write the new file
with open(str(out_dir+'SpeciesTree_mapped_names.txt'), 'w') as file:
    file.write(sptr_file)
    
#create folder
if not os.path.isdir(out_dir+'Rearranged_trees/'):
    os.makedirs(out_dir+'Rearranged_trees/')
    print("\nCreated folder: Rearranged_trees/\n")
else: 
    print('Rearranged gene trees will be stored in Rearranged_trees/\n')

#Get list of all BS trees
all_bs_trees=glob.glob(out_dir+'BS_trees/RAxML_bipartitions.*')

print('Number of trees rearranged:')

def iterate_trees(all_bs_trees):
    for bs_tree in [all_bs_trees]:
        return bs_tree


##Begin paralellization of rearranement
def par_tree_arrange(bs_tree):
    
    #Build the treerecs command
    treerecs_cmd='treerecs -s '+str(out_dir+'SpeciesTree_mapped_names.txt')+' -g '+bs_tree+' -f -t 80 --output-without-description -r -q -O newick -o '+str(out_dir+'Rearranged_trees/')

    if re.search('treerecs', treerecs_cmd) and re.search(bs_tree, treerecs_cmd):
        #subprocess.call(treerecs_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.call(treerecs_cmd, shell=True)
##End paralellization of rearranement

#Timestamp just before Parallel Call
benchmarkTime(bench_fileName, out_dir + 'benchmark/', 'start', 'tree_arrange', timer)

Parallel(n_jobs = Mult_threads, verbose=0, max_nbytes=None)(delayed(par_tree_arrange)(bs_tree) for bs_tree in iterate_trees(all_bs_trees))

#Timestamp just before Parallel Call
benchmarkTime(bench_fileName, out_dir + 'benchmark/', 'end', 'tree_arrange', timer)

#Calculate total time of parralel process & items per minute
num_items = len(all_bs_trees)
benchmarkProcess(out_dir + 'benchmark/' + str(bench_fileName), num_items)

#Check for dropped files
file_loss_log(out_dir, out_dir+'BS_trees/', out_dir+'Rearranged_trees/', 'RAxML_bipartitions.', 'RAxML_bipartitions.', 'GTSP reconciliation rearrange with treerecs', 'Cant find rearranged file')

#Count number of subtrees
file_counts_log(JOBname, out_dir, 'Rearranged_trees', 'RAxML_bipartitions.', '.nwk', out_dir+'Run_data_counts_log.csv')

### Infer branch-length optimized trees with Raxml (BL trees)
#Get list of gblocks alns (after filters) and subtrees (after filters) and find the overlap

#Make quick function for finding intersect 
def intersection(lst1, lst2):
    temp = set(lst2)
    lst3 = [value for value in lst1 if value in temp]
    return lst3

#Get a list of the IDs that sucessfully created alignments and subtrees
keeperIDs = intersection([x.replace(out_dir+'Gb_alns/GB_ALN_', '').replace('.fa', '') for x in glob.glob(out_dir+'Gb_alns/GB_ALN_*')], [y.replace(out_dir+'Rearranged_trees/RAxML_bipartitions.', '').replace('_BS.txt_recs.nwk', '') for y in glob.glob(out_dir+'Rearranged_trees/*_BS.txt_recs.nwk')])

print('Starting Branch length optimization with RAxML\n')

#Make directory to write output to
if not os.path.isdir(out_dir+'BL_trees/'):
    os.makedirs(out_dir+'BL_trees/')
    print("created folder: BL_trees/\n")
else: 
    print('Branch length optimized trees will be stored in BL_trees/\n')

# Loop through files to process
##Begin paralellization of raxml branch length optimization

def iterate_keeps(keeperIDs):
    for k_ID in [keeperIDs]:
        return k_ID

##Begin paralellization of raxml BL optimization
def par_BL_opt(k_ID, cores): 

    #Build the raxml command
    raxml_bl_cmd=Rax_dir+'raxmlHPC-PTHREADS-SSE3 -s '+str(working_dir+out_dir+'Gb_alns/GB_ALN_'+k_ID+'.fa')+ \
        ' -w '+str(working_dir+out_dir+'BL_trees/')+' -n '+str(k_ID+'_BL.txt')+ \
            ' -t '+str(working_dir+out_dir+'Rearranged_trees/RAxML_bipartitions.'+k_ID+'_BS.txt_recs.nwk')+ \
                ' -m PROTGAMMALGF -p 12345 -f e -T '+str(cores)
                
    
    #Run the command (note, raxml was installed with conda so this wont work in spyder)
    if re.search('raxmlHPC', raxml_bl_cmd) and re.search('PROTGAMMALGF', raxml_bl_cmd) and re.search('12345', raxml_bl_cmd):
        subprocess.call(raxml_bl_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

#Timestamp just before Parallel Call
benchmarkTime(bench_fileName, out_dir + 'benchmark/', 'start', 'keeps', timer)

#Call paralell
Parallel(n_jobs = Rax_front_cores, verbose=0, max_nbytes=None)(delayed(par_BL_opt)(k_ID, Rax_back_cores) for k_ID in iterate_keeps(keeperIDs))

#Timestamp just before Parallel Call
benchmarkTime(bench_fileName, out_dir + 'benchmark/', 'end', 'keeps', timer)

#Calculate total time of parralel process & items per minute
num_items = len(keeperIDs)
benchmarkProcess(out_dir + 'benchmark/' + str(bench_fileName), num_items)

#Check for dropped files
file_loss_log(out_dir, out_dir+'Rearranged_trees/', out_dir+'BL_trees/', 'RAxML_bipartitions.', 'RAxML_result.', 'Branch length optimization', 'Cant find optimized file')

#Count number of subtrees
file_counts_log(JOBname, out_dir, 'BL_trees', 'RAxML_result.', 'txt', out_dir+'Run_data_counts_log.csv')

##End paralellization of raxml branch length optimization
        
print('\n\nDone with RAxML branch length optimization\n')

print("Generating input files for GT/ST reconciliation...\n")

#Make directory to write output to
if not os.path.isdir(out_dir+'DLCpar/'):
    os.makedirs(out_dir+'DLCpar/')
    print("created folder: DLCpar/\n")
else: 
    print('Branch length optimized trees will be stored in DLCpar/\n')
    
#Copy species tree and species mapping file to DLCpar dir
shutil.copy2(sp_tr_path, str(out_dir+'DLCpar/'))

#Create the mapping file
for row_i, row in mapping_table.iterrows():
    #Add the * to the prefix
    row['Prefix'] = str(row['Prefix']).replace(str(row['Prefix']), str(row['Prefix']+'*'))
    
#Write the file for DLCpar
mapping_table.to_csv(out_dir+'DLCpar/speciesIDs.smap', sep='\t' , index=False, header=False)

'''
#Copy trees to DLCpar folder
tree2copy=glob.glob(out_dir+'BL_trees/RAxML_result.*')

for t in tree2copy:
    shutil.copy2(t, out_dir+'DLCpar/')
'''

#Generate the input files for GTST reconciliation.
input_gen_cmd= 'Rscript Generate_rec_inputs.R '+JOBname
#input_gen_cmd= 'Rscript Generate_rec_inputs.R '+JOBname+' '+OFpath    

#Run the command (if it contains strings expected in the command, this is a precautin of using shell=True)
if re.search('Generate_rec_inputs.R', input_gen_cmd) and re.search('Rscript', input_gen_cmd):
    print("Calling R with the folllowing command:")
    print(input_gen_cmd, end="\n\n")
    subprocess.call(input_gen_cmd, shell=True)

if len(glob.glob('DLCpar/*_NODES_BL.txt')) > 0:
    print('Finished generateing input files.\nInput files written to DLCpar/\n')

#Check for dropped files
file_loss_log(out_dir, out_dir+'BL_trees/', out_dir+'DLCpar/', 'RAxML_result.', 'NODES_BL.txt', 'Branch length optimization', 'Cant find optimized file')

#Count number of recs
file_counts_log(JOBname, out_dir, 'DLCpar', 'HOG', 'txt', out_dir+'Run_data_counts_log.csv')

#Count the total number genes lost
file_line_count_log(JOBname, out_dir+"Dropped_gene_log.csv", out_dir+"Run_data_counts_log.csv")

#Make a table of the steps at which genes were lost
# Read the CSV file into a DataFrame
df = pd.read_csv(out_dir+"Dropped_gene_log.csv")

# Count occurrences of each value in the "Step" column
step_counts = df['Step'].value_counts().reset_index()

# Rename columns to match the desired output format
step_counts.columns = ['Step', 'Count']

# Write the result to a new CSV file
step_counts.to_csv(out_dir+"Dropped_gene_table.csv", index=False)

### Finish with messages
print('\nIMPORTANT NOTE: the next step makes use of DLCpar, which requires python 2 (whereas the previous steps are written in python 3).\n' \
      'To run the next step you will need to enter a python 2 anaconda environment and install DLCpar.\n' \
          '\nExample commands:\n' \
            '\nIf this is your first time running, follow these three commands.\n' \
              'conda create --name dlcpar_py27 python=2.7\n' \
                  'conda activate dlcpar_py27\n' \
                      'conda install -c bioconda dlcpar\n')

print("If you've run this step previously, you only need to activate your dlcpar environment.\n" \
        "(eg: conda activate dlcpar_py27\n") 

print('After successfully completing the above steps, run the next step with the following command:\n' \
      './GTST_reconciliation.py -j '+JOBname+'\n\n')

