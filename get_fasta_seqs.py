'''
Script for extracting individual fasta sequences from a larger fasta file
'''

import sys
import csv
from Bio import SeqIO

#Path to large file that contains all sequences
source_path = "/Users/esforsythe/Documents/Work/Bioinformatics/ERC_networks/Analysis/Orthofinder/Results_Oct15/Orthogroup_Sequences/OG0000002.fa"

#List of sequencxes to extract
seq_list= ("Alyr_AL2G13120","Atha_AT1G62680","Bole_Bol011813","Brap_I01437","Brap_I01189","Bstr_22157s0068","Cgra_1670s0033","Crub_0002s0230","Esal_10023397m","Spar_Sp2g01970")

#open and read fasta file into a dictionary
source_dict = SeqIO.to_dict(SeqIO.parse(source_path, 'fasta'))

#Subset the dictionary (Turns out this is easy with a one liner)
{k: source_dict[k] for k in source_dict.keys() & seq_list}




#Loop through items in list
for i in range(1, len(seq_list)):
    #get sequences
	entry = fastadict[seq_list[i]] #find that particular entry in the dictionary
	newFile.write(">" + row[i] + "\n") #identifier of fasta
    newFile.write(str(entry.seq) + "\n") #sequence