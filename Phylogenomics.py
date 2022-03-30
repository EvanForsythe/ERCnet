'''
### Main script for running the phylogenomics steps of the ERC workflow###

conda activate ERC_networks

Example command:    
    #Full run: 
    python Phylogenomics.py -j TPC_test -t 50 -p 3 -r 15 -l 100 -s -o /Users/esforsythe/Documents/Work/Bioinformatics/ERC_networks/Analysis/Orthofinder/Plant_cell/Results_Feb15/ -x /opt/anaconda3/envs/ERC_networks/bin/

'''
#During developent, set working directory:
#import os
#working_dir = '/Users/esforsythe/Documents/Work/Bioinformatics/ERC_networks/Analysis/ERCnet_dev/'
#os.chdir('/Users/esforsythe/Documents/Work/Bioinformatics/ERC_networks/Analysis/ERCnet_dev/')

#Storebought modules
import os
import sys
import argparse
import pandas as pd
import numpy as np
import itertools
import glob
import subprocess
import re
import math
import shutil
from Bio import SeqIO

#Homemade modules
from filterHOGs import make_seq_counts_df, filter_gene_fams

#At runtime set working directory to the place where the script lives
working_dir = sys.path[0]+'/' 
os.chdir(working_dir)

#Set up an argumanet parser
parser = argparse.ArgumentParser(description='Main ERCnet script')

parser.add_argument('-j', '--JOBname', type=str, metavar='', required=True, help='Unique job name for this run of ERCnet. Avoid including spaces or special characters ("_" is ok)') 
parser.add_argument('-o', '--OFpath', type=str, metavar='', required=True, help='Full path to the Orthofinder results dir (should contain Species_Tree/, Phylogenetic_Hierarchical_Orthogroups/ etc...)\n Include "/" at the end of the string') 
parser.add_argument('-p', '--MaxP', type=int, metavar='', required=False, help='Integer: maximum number of paralogs per species allowed in each gene family' )
parser.add_argument('-r', '--MinR', type=int, metavar='', required=False, help='Integer: minimum number of species represented required in each gene family' )
parser.add_argument('-t', '--Test_num', type=int, metavar='', required=False, help='Integer: number of gene families to analyze (for testing only)' )
parser.add_argument('-e','--explore_filters', action='store_true', required=False, help='Add this flag to explore filtering options (if selected, program will quit without running downstream steps)')
parser.add_argument('-l', '--Min_len', type=int, metavar='', required=False, help='Integer: minimum length of alignment (after trimming with Gblocks) required to retain gene' )
parser.add_argument('-x', '--Rax_dir', type=str, metavar='', required=True, help='Full path to the location of your raxml install (use which raxmlHPC to locate). Include "/" at the end of the string') 
parser.add_argument('-s','--SPmap', action='store_true', required=False, help='Add this flag to provide a custom species mapping file. This mapping file must be formatted in certian way. See instuctions')


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


'''
#DEV: hardcode arguments
JOBname = "TEST"
OFpath = "/Users/esforsythe/Documents/Work/Bioinformatics/ERC_networks/Analysis/Orthofinder/Results_Oct15/"
MaxP_val=2
MinR_val=8
explore_filters=False
Min_len=100
Test_num=50
Rax_dir= "/opt/anaconda3/envs/ERC_networks/bin/"
#END DEV
'''

'''
#DEV: hardcode arguments
JOBname = "TPC_testing"
OFpath = "/Users/esforsythe/Documents/Work/Bioinformatics/ERC_networks/Analysis/Orthofinder/Plant_cell/Results_Feb15/"
MaxP_val=3
MinR_val=17
explore_filters=False
Min_len=100
Test_num=50
Rax_dir= "/opt/anaconda3/envs/ERC_networks/bin/"
SPmap=True
#END DEV
'''


#Check to see if the path arguments end with "/"
if not OFpath.endswith('/'):
    print('OFpath does not end with a "/". Quitting...\n')
    sys.exit()

if not Rax_dir.endswith('/'):
    print('Rax_dir does not end with a "/". Quitting...\n')
    sys.exit()

#Get path to some needed Orthofinder files 
sp_tr_path = OFpath+'Species_Tree/SpeciesTree_rooted_node_labels.txt'
N1_file_path = OFpath+'Phylogenetic_Hierarchical_Orthogroups/N1.tsv'
OG_trees_dir = OFpath+'Resolved_Gene_Trees/'
OGseqdir = OFpath+'Orthogroup_Sequences/'

#Store output dir as a variable
out_dir= 'OUT_'+JOBname+'/'

#Create the output folder
#Make a directory for storing stats (and log file)
if not os.path.isdir(out_dir):
    os.makedirs(out_dir)
    print('created folder: '+out_dir+'\nAll output files will be written to this folder\n')
else: print('All output files will be written to '+out_dir+'\n')

#Check N1file exists 
if os.path.isfile(N1_file_path):
    print('N1 HOG file located at:\n'+ N1_file_path+'\n')
    
else:
    print('N1 HOG file not found. Quitting...\n')
    sys.exit()


#open N1 file
N1_file=pd.read_csv(N1_file_path, sep='\t')

#Make a list of species
sp_list_temp = list(N1_file.columns[3:])

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
print("These should match the column headers in the N1 file and the tip labels in the species tree")

print("\nThe following species identifiers will be used: \n"+str(sp_prefix))
print("These should be present in the sequence IDs in alignments etc...\n")

#Generate a dataframe of the counts data using the my outside module from filterHOGs.py
seq_counts_df=make_seq_counts_df(N1_file_path, out_dir+'Species_mapping.csv')

#Check if the dataframe was assinged properly
if 'HOG' in list(seq_counts_df.columns):
    print('Sequence counts per sepecies dataframe successfully generated\n')
else:
    print('ERROR: Sequence counts per sepecies dataframe not properly generated. Quitting...\n')
    sys.exit()

#If the -e flag was chosen, run the parameter scan to explore the data filtering options
if explore_filters:
    print('-e / --explore_filters flag chosen\nExploring filtering options and quitting\n')
    
    
    #Set filter ranges
    max_paralogs_vals=list(range(1, 7, 1))
    min_rep_vals=list(range(3, len(sp_names), 1))
    
    #Make blank array for parameter scan
    retained_trees=np.zeros((len(max_paralogs_vals), len(min_rep_vals)))
    
    #Start nested loop to fill in the matrix
    for p_count, p_val in enumerate (max_paralogs_vals):
        for s_count, s_val in enumerate (min_rep_vals):
            para_value=p_val
            rep_value=s_val

            passed_only_df = filter_gene_fams(N1_file, seq_counts_df, sp_names, para_value, rep_value)
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
if 'MaxP_val' in globals() and isinstance(MaxP_val, int) and'MinR_val' in globals() and isinstance(MinR_val, int):
    print('--MaxP_val set to {} and --MinR_val set to {}\nFiltering data....'.format(MaxP_val,MinR_val))
else:
    print('--MaxP_val and --MinR_val required (unless you use the --explore_filters flag)\nexiting...')
    sys.exit()

##Filter results according to filter criteria
#Use the module from Filter_stats.py to filter the list
Keeper_HOGs_df= filter_gene_fams(N1_file, seq_counts_df, sp_names, MaxP_val, MinR_val)

#FOR TESTING: filter the dataset to a user-defined subset
if isinstance(Test_num, int):
    Keeper_HOGs_df = Keeper_HOGs_df.sample(n=Test_num)

#Write csv file and report the number of gene trees included after filtering
if 'HOG' in list(Keeper_HOGs_df.columns):
    Keeper_HOGs_df.to_csv(out_dir+'Filtered_genefam_dataset.csv', sep=',' , index=False)
    print('Filtered list of gene families successfully generated. \n\nFiltered dataset contains %d genes.\n\nCsv file written to Filtered_genefam_dataset.csv \n' %Keeper_HOGs_df.shape[0])
else:
    print('Error: sequence counts per sepecies dataframe not properly generated. Quitting...\n')
    sys.exit()

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
    with open(str(out_dir+'HOG_seqs/'+HOGtemp.replace("N1.", "")+'.fa'), 'w') as handle:
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

print('Alignments finished: ')
for file_i, file in enumerate(seq_file_names):
    if file_i % 10 == 0:
        print(file_i)
    os.system('mafft-linsi --quiet '+file+' > '+out_dir+'Alns/ALN_'+file.replace(out_dir+"HOG_seqs/", ""))
    #os.system('mafft-linsi '+file+' >Alns/ALN_'+file.replace("HOG_seqs/", "")+' 2>&1') #' 2>&1' suppressed stderr from mafft
    #print('mafft-linsi --quiet '+file+' > Alns/ALN_'+file.replace("HOG_seqs/", ""))

#Check alignment status
if len(glob.glob(out_dir+'Alns/ALN*')) == len(seq_file_names):
    print('\n\nDone with multiple sequences alignments\n')
elif len(glob.glob(out_dir+'Alns/ALN*')) < len(seq_file_names):
    print('\nWARNING: some alignments did not finish. Proceeding (with caution)...\n')
elif len(glob.glob(out_dir+'Alns/ALN*')) < len(seq_file_names):
    print('\nWARNING: there are more alignments in the folder than expected. Proceeding (with caution)...\n')

###GBLOCKS cleaning alignments
print('Beginning GBLOCKS\n')

#Get list of all alns (that haven't been gblocked yet)
aln_file_names = [x for x in glob.glob(out_dir+'Alns/ALN*') if "-gb" not in x]


#Gblocks
if 'Min_len' in locals():
    print('Beginning gblocks trimming...\n')
else:
    print('Min_len parameter for gblocks filtering not specified. Quitting....\n')

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
for aln in aln_file_names:
    #aln = aln_file_names[0]
    #Get number of sequences in alignment
    aln_file = open(aln, "r")
    line_count = 0
    for line in aln_file:
        if re.search(r">", line):
            line_count += 1
    
    #Create gblocks command
    gblocks_cmd = 'Gblocks '+aln+' -b5=h -b4=5 -p=n -b2='+str(math.floor((line_count/2)+1))
    
    if re.search('Gblocks', gblocks_cmd) and re.search('-b5', gblocks_cmd) and re.search('-p=n', gblocks_cmd):
        #Run cmd and capture std output (in order to know the trimmed alignment length)
        proc = subprocess.Popen(gblocks_cmd, stdout=subprocess.PIPE, shell=True) #apparently shell=True can create 'security issues' (so I put it under an if statement to make sure the cmd is what it shoudl be)
        output = str(proc.stdout.read())
        
        #Extract the length of the trimmed alignment by parsing Gblocks stdout 
        if re.search('Gblocks alignment', output): #putting this in an if statement to make sure it's parsing the right thing
            trm_aln_ln= int(output.split(r'\n')[6].split(' ')[3]) #This works but I think it could break easily
            
            #Move the gblocks files to the appropriate folder
            if trm_aln_ln >= Min_len:
                os.replace(str(aln+'-gb'), str(aln+'-gb').replace('Alns/', 'Gb_alns/'))
                #os.replace(str(aln+'-gb.htm'), str(aln+'-gb.htm').replace('Alns/', 'Gb_alns/HTML_files/'))
            else:
                print('%s not long enough. Moving to gblocks files to Gb_alns/Too_short/' %aln.replace('Alns/', ''))
                os.replace(str(aln+'-gb'), str(aln+'-gb').replace('Alns/', 'Gb_alns/Too_short/'))
                #os.replace(str(aln+'-gb.htm'), str(aln+'-gb.htm').replace('Alns/', 'Gb_alns/Too_short/'))
            
        else:
            print('ERROR: something wrong with Gblocks stdout. Quitting...\n')
            sys.exit()

            
    else:
        print('WARNING: Something wrong with Gblocks command...\n')
        sys.exit()

### Get the subtree tree files (to use in branch length optimization)
#Create folder to store subtrees
#Make a directory for gblocks trimmed alignments that are too short 
if not os.path.isdir(out_dir+'HOG_subtrees/'):
    os.makedirs(out_dir+'HOG_subtrees/')
    print("\nCreated folder: HOG_subtrees/\n")
else: 
    print('HOG subtrees will be stored in HOG_subtrees/\n')
    
print('Running subtree extraction. Check Non-binary_subtrees.txt for list of trees that are excluded for being non-bifurcating\n')

#Call the R-script used to extract trees
get_st_cmd= 'Rscript Get_subtree.R '+OG_trees_dir+' '+out_dir
    
#Run the command (if it contains strings expected in the command, this is a precautin of using shell=True)
if re.search('Get_subtree.R', get_st_cmd) and re.search(OG_trees_dir, get_st_cmd):
    print("Calling R with the following command:")
    print(get_st_cmd)
    subprocess.call(get_st_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

#Checj if it worked
if len(glob.glob(out_dir+'HOG_subtrees/*tree.txt')) > 0:
    print('\nFinished running subtree extraction. Check "Non-binary_subtrees.txt" for a list trees that were excluded\n')
else:
    print('ERROR: An error occured during subtree extraction. Quitting...\n')
    sys.exit()

### Infer trees with Raxml
#Get list of gblocks alns (after filters) and subtrees (after filters) and find the overlap

#Make quick function for finding intersect 
def intersection(lst1, lst2):
    temp = set(lst2)
    lst3 = [value for value in lst1 if value in temp]
    return lst3

#Get a list of the IDs that sucessfully created alignments and subtrees
keeperIDs = intersection([x.replace(out_dir+'Gb_alns/ALN_', '').replace('.fa-gb', '') for x in glob.glob(out_dir+'Gb_alns/ALN_*')], [y.replace(out_dir+'HOG_subtrees/', '').replace('_tree.txt', '') for y in glob.glob(out_dir+'HOG_subtrees/*_tree.txt')])

print('Starting Branch Length Optimization with RAxML\n')
      
#Make directory to work in
if not os.path.isdir(out_dir+'Trees_working/'):
    os.makedirs(out_dir+'Trees_working/')
    print("created folder: Trees_working/\n")
else: 
    print('HOG subtrees will be stored in Trees_working/\n')
    
#Make directory to write output to
if not os.path.isdir(out_dir+'BL_trees/'):
    os.makedirs(out_dir+'BL_trees/')
    print("created folder: BL_trees/\n")
else: 
    print('HOG subtrees will be stored in BL_trees/\n')

print('Trees finished: ') 
# Loop through files to process
for HOG_i, HOG_id in enumerate(keeperIDs):
    #Track progress
    if HOG_i % 10 == 0:
        print(HOG_i)

    #Copy the needed files to the directory
    #tree
    shutil.copy2(out_dir+'HOG_subtrees/'+HOG_id+'_tree.txt', out_dir+'Trees_working/'+HOG_id+'_tree.txt')
    
    #alignment (remove -gb suffix)
    shutil.copy2(out_dir+'Gb_alns/ALN_'+HOG_id+'.fa-gb', out_dir+'Trees_working/GB_ALN_'+HOG_id+'.fa')
    
    #Build the command
    raxml_cmd= Rax_dir+'raxmlHPC '+'-s '+working_dir+out_dir+'Trees_working/GB_ALN_'+HOG_id+'.fa '+'-w '+working_dir+out_dir+'BL_trees/'+ \
    ' -n '+HOG_id+'_BL.txt'+' -t '+working_dir+out_dir+'Trees_working/'+HOG_id+'_tree.txt'+ \
    ' -m PROTGAMMALGF -p 12345 -f e'
    
    #Run the command (note, raxml was installed with conda so this wont work in spyder)
    if re.search('raxmlHPC', raxml_cmd) and re.search('PROTGAMMALGF', raxml_cmd) and re.search('12345', raxml_cmd):
        subprocess.call(raxml_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        #subprocess.call(raxml_cmd, shell=True)

print('\n\nDone with RAxML branch length optimization\n')

print("Generating input files for GT/ST reconciliation...\n")

#Generate the input files for GTST reconciliation.
input_gen_cmd= 'Rscript Generate_rec_inputs.R '+JOBname+' '+OFpath
    
#Run the command (if it contains strings expected in the command, this is a precautin of using shell=True)
if re.search('Generate_rec_inputs.R', input_gen_cmd) and re.search('Rscript', input_gen_cmd):
    print("Calling R with the folllowing command:")
    print(input_gen_cmd, end="\n\n")
    subprocess.call(input_gen_cmd, shell=True)

if len(glob.glob('DLCpar/*_NODES_BL.txt')) > 0:
    print('Finished generateing input files.\nInput files written to DLCpar/\n')
    

print('\nIMPORTANT NOTE: the next step makes use of DLCpar, which requires python 2 (whereas the previous steps are written in python 3).\n' \
      'To run the next step you will need to enter a python 2 anaconda environment and install DLCpar.\n' \
          '\nExample commands:\n' \
              'conda create --name dlcpar_py27 python=2.7\n' \
                  'conda activate dlcpar_py27\n' \
                      'conda install -c bioconda dlcpar\n')

print('After successfully completing the above steps, run the next step with the following command:\n' \
      'python GTST_reconciliation.py -j '+JOBname
      )

