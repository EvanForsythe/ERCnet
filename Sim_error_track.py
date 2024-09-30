import argparse
import os
import pandas as pd
import re

def parse_args():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Process a jobname and file to extract GeneA_ID and GeneB_ID")
    parser.add_argument("jobname", help="The name of the job")
    parser.add_argument("filename", help="The name of the TSV file")
    return parser.parse_args()

def extract_integer_from_id(gene_id):
    # Extract the integer part from the gene_id using regex
    match = re.search(r"(\d+)$", gene_id)
    if match:
        return int(match.group(1))  # Return the integer value
    else:
        raise ValueError(f"No integer found in {gene_id}")

def process_tsv(jobname, filename):
    # Construct the path to the TSV file
    tsv_path = os.path.join(f"OUT_{jobname}", "ERC_results", "Filtered_results", filename)
    
    # Load the TSV file into a pandas DataFrame
    try:
        df = pd.read_csv(tsv_path, sep="\t")
        
        true_pos_counter=0
        false_pos_counter=0

        # Loop through each row and extract GeneA_ID and GeneB_ID
        for index, row in df.iterrows():
            gene_a_id = row['GeneA_ID']
            gene_b_id = row['GeneB_ID']
            
            # Extract the integer part from GeneA_ID and GeneB_ID
            gene_a_number = extract_integer_from_id(gene_a_id)
            gene_b_number = extract_integer_from_id(gene_b_id)
            
            if gene_a_number <=100 and gene_b_number <=1000:
                true_pos_counter+=1
            else:
                false_pos_counter+=1

    except FileNotFoundError:
        print(f"Error: File {tsv_path} not found.")
    except Exception as e:
        print(f"An error occurred: {e}")
    return true_pos_counter, false_pos_counter



def main():
    # Parse the command line arguments
    args = parse_args()
    
    # Process the TSV file
    true_pos_counter, false_pos_counter = process_tsv(args.jobname, args.filename)

    #Print values
    print(f"Number of true positives: {true_pos_counter}")
    print(f"Number of false positives: {false_pos_counter}")
    
    #Get the expected values
    n_coevolve=100
    n_total=1000
    total_comparisons= (n_total * ( n_total - 1)) // 2
    expected_hits= (n_coevolve * (n_coevolve - 1)) // 2
    expected_misses = total_comparisons - expected_hits
    false_pos_rate=false_pos_counter/expected_misses
    false_neg_rate=1-(true_pos_counter/expected_hits)

    #Print values
    print(f"False positive rate: {false_pos_rate}")
    print(f"False negative rate: {false_neg_rate}")
    
if __name__ == "__main__":
    main()
