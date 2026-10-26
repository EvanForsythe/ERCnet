# --- IMPORTS ---
import argparse
from pathlib import Path
import shutil
from Bio import SeqIO
import math
import os
import pandas as pd

# --- FUNCTION DECLARATIONS ---
def calculatePercentageGap(window_file):
  totalSequenceChars = 0
  totalGapChars = 0

  records = list(SeqIO.parse(window_file, "fasta"))

  for record in records:
    totalSequenceChars += len(record)

  for record in records:
    totalGapChars += record.seq.count('-')

  return totalGapChars / totalSequenceChars


def sliding_window_fasta(input_fasta, output_dir, num_windows, parent_gene_int):
  """
  Creates smaller alignment files using a non-overlapping sliding window approach.

  Args:
  input_fasta (str): Path to the input FASTA file containing long sequence alignments.
  output_dir (str): Directory to store the output FASTA files.

  Returns:
  None
  """

  logFileColumns = ['Full_Aln_File_Name', 'Window_File_Name', 'Window_Length', 'Window_ Number', 'Start_Index', 'End_Index', 'Percentage_Gap_Across_All_Species']
  logFileRows = []

  # Read the input FASTA file and get all sequences
  records = list(SeqIO.parse(input_fasta, "fasta"))
  
  if len(records) == 0:
    print("No sequences found in the input FASTA file.")
    return

  # Determine the length of the alignment (assuming all sequences have the same length)
  alignment_length = len(records[0].seq)
  
# Calculate window length as 1/1000 of the alignment length
  win_len = math.ceil(alignment_length / num_windows)

  # Iterate over the alignment in windows
  for window_num in range(num_windows):
    logFileData = []
    logFileData.append(input_fasta)

    if parent_gene_int < 10:
      windowName = "GB_ALN_HOG0" + str(parent_gene_int) + "0000" + str(window_num+1) + ".fa"
    else:
      windowName = "GB_ALN_HOG" + str(parent_gene_int) + "0000" + str(window_num+1) + ".fa"

    output_window_file_path = os.path.join(output_dir, windowName)
    logFileData.append(output_window_file_path)
    logFileData.append(win_len)
    logFileData.append(window_num)
    
    start = window_num * win_len
    end = min(start + win_len, alignment_length)  # Ensure the last window doesn't exceed alignment length
    logFileData.append(start)
    logFileData.append(end)

    # Create a new list to store windowed sequences
    window_records = []

    # Extract the windowed sequences for each record
    for record in records:
      window_seq = record.seq[start:end]
      window_record = record[start:end]
      window_records.append(window_record)

    # Write the windowed sequences to a new FASTA file
    SeqIO.write(window_records, output_window_file_path, "fasta")
    logFileData.append(calculatePercentageGap(output_window_file_path))

    logFileRows.append(logFileData)

    # Stop if end has reached the alignment length
    if end >= alignment_length:
      break
    
  print(' > Writing log file...')  
  logFileDF = pd.DataFrame(logFileRows, columns=logFileColumns)
  logWriteDir = os.path.join(output_dir, 'log.tsv')
  logFileDF.to_csv(logWriteDir, sep='\t')

  print("> Done!")

# --- MAIN ---
print('Starting ERCnet windows feature preparation.')
# --- argparse ---
print('Parsing arguments...')
parser = argparse.ArgumentParser()
parser.add_argument('--goi', '-g', required=True, help='Path to tsv file which provides genes of interets for windows analysis.')
parser.add_argument('--job_name', '-j', required=True, help='Job name of standard ERCnet run. This flag string must match job name exactly.')
parser.add_argument('--windows_number', '-n', required=False, type=int, default=4, help='Number of window files to produce for each gene of interest. Defaults to 4 if not specified.')
args = parser.parse_args()

print('Creating new directories...')
original_job_path = 'OUT_' + args.job_name + '/'
if not Path(original_job_path).exists():
  raise SystemExit('No ERCnet run found with provided job name. Please try again. Exiting windows script.')

# Make new windows analysis out directory
windows_job_path = 'OUT_' + args.job_name + '_windows/'
shutil.rmtree(windows_job_path, ignore_errors=True)
Path(windows_job_path).mkdir()

# --- Copy Directories and Files ---
print('Copying directories and files required for Phylogenomics.py from original ERCnet job...')

print('> Dropped_gene_log.csv', flush=True)
shutil.copy(original_job_path + 'Dropped_gene_log.csv', windows_job_path)

print('> Run_data_counts_log.csv', flush=True)
shutil.copy(original_job_path + 'Run_data_counts_log.csv', windows_job_path)

print('> benchmark/', flush=True)
shutil.copytree(original_job_path + 'benchmark/', windows_job_path + 'benchmark/', dirs_exist_ok=True)
#shutil.copy(original_job_path + 'benchmark/' + args.job_name + '_Phylogenomics_benchmark.tsv', windows_job_path)

print('> Species_mapping.csv', flush=True)
shutil.copy(original_job_path + 'Species_mapping.csv', windows_job_path)

print('> Seq_counts_per_species.csv', flush=True)
shutil.copy(original_job_path + 'Seq_counts_per_species.csv', windows_job_path)

print('> Filtered_genefam_dataset.csv', flush=True)
shutil.copy(original_job_path + 'Filtered_genefam_dataset.csv', windows_job_path)

print('> HOG_seqs/', flush=True)
shutil.copytree(original_job_path + 'HOG_seqs/', windows_job_path + 'HOG_seqs/', dirs_exist_ok=True)

print('> Alns/', flush=True)
shutil.copytree(original_job_path + 'Alns/', windows_job_path + 'Alns/', dirs_exist_ok=True)

print('> TAPER_Alns/', flush=True)
shutil.copytree(original_job_path + 'TAPER_Alns/', windows_job_path + 'TAPER_Alns/', dirs_exist_ok=True)

print('> Gb_alns_full_length/', flush=True)
full_Gb_alns_path = windows_job_path + '/Gb_alns_full_length/'
Path(full_Gb_alns_path).mkdir()
shutil.copytree(original_job_path + 'Gb_alns/', full_Gb_alns_path, dirs_exist_ok=True)

print('> Aln_pruning/', flush=True)
shutil.copytree(original_job_path + 'Aln_pruning/', windows_job_path + 'Aln_pruning/', dirs_exist_ok=True)

# --- Create window files and insert into Gb_alns ---
print('Creating Gb_alns directory for windows files...')
windows_Gb_alns_path = windows_job_path + '/Gb_alns/'
Path(windows_Gb_alns_path).mkdir()

print('Reading genes of interest tsv...')
goi_df = pd.read_csv(args.goi, sep='\t')

print('Slicing full gene of interest alignments for windows...')
full_gene_num = 0
for geneHOG in goi_df['HOG']:
  full_gene_num += 1
  print('> Now producing windows for ' + geneHOG)
  sliding_window_fasta(full_Gb_alns_path + 'GB_ALN_' + geneHOG + '.fa', windows_Gb_alns_path, args.windows_number, full_gene_num)
