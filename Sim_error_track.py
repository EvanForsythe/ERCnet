import argparse
import os
import pandas as pd
import re

def parse_args():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Process a jobname to analyze Network_stats and TSV files")
    parser.add_argument("jobname", help="The name of the job")
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
    tsv_path = os.path.join(f"OUT_{jobname}", "ERC_results", "Filtered_results", f"{filename}.tsv")
    
    # Load the TSV file into a pandas DataFrame
    try:
        df = pd.read_csv(tsv_path, sep="\t")
        
        true_pos_counter = 0
        false_pos_counter = 0

        # Loop through each row and extract GeneA_ID and GeneB_ID
        for index, row in df.iterrows():
            gene_a_id = row['GeneA_ID']
            gene_b_id = row['GeneB_ID']
            
            # Extract the integer part from GeneA_ID and GeneB_ID
            gene_a_number = extract_integer_from_id(gene_a_id)
            gene_b_number = extract_integer_from_id(gene_b_id)
            
            if gene_a_number <= 100 and gene_b_number <= 100:
                true_pos_counter += 1
            else:
                false_pos_counter += 1

    except FileNotFoundError:
        print(f"Error: File {tsv_path} not found.")
    except Exception as e:
        print(f"An error occurred: {e}")
    
    return true_pos_counter, false_pos_counter

def calculate_error_rates(true_pos_counter, false_pos_counter):
    # Calculate expected values and error rates
    n_coevolve = 100
    n_total = 1000
    total_comparisons = (n_total * (n_total - 1)) // 2
    expected_hits = (n_coevolve * (n_coevolve - 1)) // 2
    expected_misses = total_comparisons - expected_hits
    false_pos_rate = false_pos_counter / expected_misses
    false_neg_rate = 1 - (true_pos_counter / expected_hits)

    return false_pos_rate, false_neg_rate

def main():
    # Parse the command line arguments
    args = parse_args()
    
    # Load the Network_stats.csv file into a DataFrame
    network_stats_path = os.path.join(f"OUT_{args.jobname}", "Network_stats.csv")
    df_network_stats = pd.read_csv(network_stats_path)
    
    # Add new columns to the DataFrame
    df_network_stats['true_pos_counter'] = 0
    df_network_stats['false_pos_counter'] = 0
    df_network_stats['false_pos_rate'] = 0.0
    df_network_stats['false_neg_rate'] = 0.0

    # Process each row in the DataFrame
    for index, row in df_network_stats.iterrows():
        filename = row['fileName']  # Get the filename from the 'fileName' column
        
        # Process the TSV file and calculate counters
        true_pos_counter, false_pos_counter = process_tsv(args.jobname, filename)
        
        # Calculate error rates
        false_pos_rate, false_neg_rate = calculate_error_rates(true_pos_counter, false_pos_counter)
        
        # Update the DataFrame with the new values
        df_network_stats.at[index, 'true_pos_counter'] = true_pos_counter
        df_network_stats.at[index, 'false_pos_counter'] = false_pos_counter
        df_network_stats.at[index, 'false_pos_rate'] = false_pos_rate
        df_network_stats.at[index, 'false_neg_rate'] = false_neg_rate

    # Write the updated DataFrame to a new CSV file
    output_csv_path = os.path.join(f"OUT_{args.jobname}", "SIM_error_rates.csv")
    df_network_stats.to_csv(output_csv_path, index=False)
    
    print(f"Results saved to {output_csv_path}")

if __name__ == "__main__":
    main()
