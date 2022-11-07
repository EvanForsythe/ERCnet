#!/usr/bin/env python3

import pandas as pd
import numpy as np
from Bio import AlignIO, SeqIO


#Hardcoded for development
#HOG_file_path='/Users/esforsythe/Documents/Work/Bioinformatics/ERC_networks/Analysis/Orthofinder/Results_Oct15/Phylogenetic_Hierarchical_Orthogroups/N1.tsv'
#HOG_file_path='/Users/esforsythe/Documents/Work/Bioinformatics/ERC_networks/Analysis/Orthofinder/Plant_cell/Results_Feb15/Phylogenetic_Hierarchical_Orthogroups/N1.tsv'
#mapping_table_path='/Users/esforsythe/Documents/Work/Bioinformatics/ERC_networks/Analysis/ERCnet_dev/OUT_TPC/Species_mapping.csv'

#Function for creating a dataframe with the species counts
def make_seq_counts_df(HOG_file_path, mapping_table_path):

    #open HOG file
    HOG_file=pd.read_csv(HOG_file_path, sep='\t')
    
    #Get a list of the species names (note the first 4 columns are various IDs)
    mapping_table=pd.read_csv(mapping_table_path, sep=',')
    sp_prefix=list(mapping_table['Prefix'])
    sp_names=list(mapping_table['SpeciesID'])
    
    ##Create dataframe that has the counts for each species
    #Create empty df to fill in 
    counts_array=np.empty((HOG_file.shape[0], len(sp_names)))
    counts_array[:] = np.NaN
    
    for row_i, row in HOG_file.iterrows():
        for spec_i, spec in enumerate(sp_prefix):
            if str(row[sp_names[spec_i]]).count(spec)==0:
                count=0
            elif str(row[sp_names[spec_i]]).count(spec)>0:
                count=str(row[sp_names[spec_i]]).count(',')+1
            else:
                count='nan'
            counts_array[row_i, spec_i]=count
    
    #Convert to def
    counts_df=pd.DataFrame(counts_array)
    
    counts_df = counts_df.apply(pd.to_numeric, errors='coerce', downcast='integer')
    
    counts_df.columns=list(sp_names)
    
    counts_df['HOG'] = HOG_file['HOG']
    counts_df['OG'] = HOG_file['OG']
    counts_df['Gene Tree Parent Clade'] = HOG_file['Gene Tree Parent Clade']

    #Rearrnge the columns
    counts_df_rearrange=counts_df.reindex(columns=list(list(counts_df.columns[-3:]) + list(counts_df.columns[0:-3])))

    return(counts_df_rearrange)

if __name__ == "__main__":
    make_seq_counts_df(HOG_file_path, mapping_table_path)
    

#Make function for filtering dataset by criteria
def filter_gene_fams(HOG_file, counts_df_rearrange, sp_names, para_value, rep_value):
    #Max number of paralogs
    def max_paralog_bool(row):
        return row.max() <= para_value
    
    #Apply the function to the rows of the dataframe
    paralog_bool=list(counts_df_rearrange[sp_names].apply(max_paralog_bool, axis=1))
    
    #Minimum number of species represeented
    #rep_value=3
    #row=list(counts_df.loc[50,])
    def min_sp_rep_bool(row):
        def condition(x): 
            return x > 0
        bool_arr = condition(row)
        return len(np.where(bool_arr)[0]) >= rep_value
    
    #Apply the function to the rows of the dataframe
    rep_bool=list(counts_df_rearrange[sp_names].apply(min_sp_rep_bool, axis=1))
    
    #Make blank dataframe to combine the boolian lists above
    conditions_met_df= pd.DataFrame()
    
    #Fill in info
    conditions_met_df['HOG']=HOG_file['HOG']
    conditions_met_df['Paralogs_OK']=paralog_bool
    conditions_met_df['SPrep_OK']=rep_bool
    
    #Get the rows that pass both criteria
    passed_only_df=conditions_met_df[(conditions_met_df['Paralogs_OK']==True) & (conditions_met_df['SPrep_OK']==True)]
    
    #make new df that retains only the keepers 
    Keeper_HOGs_df= HOG_file[HOG_file['HOG'].isin(list(passed_only_df['HOG']))]
    
    #Return the filtered dataframe
    return Keeper_HOGs_df

if __name__ == "__main__":
    filter_gene_fams(HOG_file, counts_df_rearrange, sp_names, para_value, rep_value)

#Make a function that will write the names of ids that should be pruned
def check_alns_for_prune(alns_to_check_temp, prune_cutoff, out_dir):
    #Get the file name
    aln_check_path=str(alns_to_check_temp)
    
    #Get the OG name
    aln_check_OG=aln_check_path.replace(out_dir+"Gb_alns/GB_ALN_", "").replace(".fa", "")
    
    #Read in alignment
    alignment = AlignIO.read(aln_check_path, "fasta")
    
    #Get the number of total sites in the alignment
    n_sites=alignment.get_alignment_length()
    
    #Create blank list
    prune_seqs=[]
    
    #Loop through the seqs in the alignment and count the gap sites (add to list if too many gaps)
    for record in alignment:
        if record.seq.count("-")>(prune_cutoff*n_sites):
            prune_seqs.append(record.id)
    
    #If there are any on the list
    if len(prune_seqs) > 0:
        #Make results file
        with open(out_dir+'Aln_pruning/Prune_IDs_'+aln_check_OG+".txt", "a") as f:
            for item in prune_seqs:
                f.write("%s\n" % item)
        
        #write a new version of the file
        #Open commenction
        FastaDroppedFile = open(out_dir+"Aln_pruning/GB_ALN_"+aln_check_OG+".fa", 'w')
        
        #Loop through the seqs in the alignment and write the keepers to the new file
        for record in alignment:
            if not record.seq.count("-")>(prune_cutoff*n_sites):
                SeqIO.write(record, FastaDroppedFile, 'fasta')
        #Close connection
        FastaDroppedFile.close()
        
    #I need to have a return statement or map() wont run                
    return(aln_check_path)

if __name__ == "__main__":
    check_alns_for_prune(alns_to_check_temp, prune_cutoff, out_dir)
    
    

