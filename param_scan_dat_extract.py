import os
import sys
import re
from datetime import datetime

def extract_total_time(jobname):
    # Construct the file paths
    base_path = os.path.join("/scratch/forsythe/ERCnet", f"OUT_{jobname}", "benchmark")
    erc_results_path = os.path.join("/scratch/forsythe/ERCnet", f"OUT_{jobname}", "ERC_results")
    data_counts_filepath = os.path.join("/scratch/forsythe/ERCnet", f"OUT_{jobname}", "Run_data_counts_log.csv")

    # Paths for benchmark and ERC results files
    phylo_filepath = os.path.join(base_path, f"{jobname}_Phylogenomics_benchmark.tsv")
    gtst_filepath = os.path.join(base_path, f"{jobname}_GTST_rec_benchmark.tsv")
    r2t_filepath = os.path.join(base_path, f"{jobname}_ERC_analyses_benchmark_R2T.tsv")
    bxb_filepath = os.path.join(base_path, f"{jobname}_ERC_analyses_benchmark_BXB.tsv")
    network_filepath = os.path.join(base_path, f"{jobname}_Network_analysis_benchmark.tsv")
    erc_r2t_filepath = os.path.join(erc_results_path, "ERC_results_R2T.tsv")
    erc_bxb_filepath = os.path.join(erc_results_path, "ERC_results_BXB.tsv")

    # Initialize variables to hold the output values
    values = [jobname]  # Start with jobname as the first value
    total_time_sum = 0.0
    gtst_runtime = None
    r2t_total_time = None
    bxb_total_time = None
    network_runtime = None
    dlcpar_count = None
    erc_r2t_avg_branches = None
    erc_bxb_avg_branches = None

    # Process Phylogenomics benchmark file
    try:
        with open(phylo_filepath, 'r') as file:
            for line in file:
                if "Total time (m)" in line:
                    total_time_sum += float(line.split("\t")[-1].strip())
        values.append(total_time_sum)
    except FileNotFoundError:
        values.append("N/A")
    except Exception as e:
        values.append(f"Error: {e}")

    # Process GTST reconciliation benchmark file
    try:
        with open(gtst_filepath, 'r') as file:
            start_time, end_time = None, None
            for line in file:
                if "GTST_reconciliation Call start" in line:
                    start_time = datetime.strptime(line.split("\t")[-1].strip(), "%a, %d %b %Y %H:%M:%S")
                elif "GTST_reconciliation Call end" in line:
                    end_time = datetime.strptime(line.split("\t")[-1].strip(), "%a, %d %b %Y %H:%M:%S")
            if start_time and end_time:
                gtst_runtime = (end_time - start_time).total_seconds() / 60
        values.append(gtst_runtime if gtst_runtime is not None else "N/A")
    except FileNotFoundError:
        values.append("N/A")
    except Exception as e:
        values.append(f"Error: {e}")

    # Process ERC R2T benchmark file
    try:
        with open(r2t_filepath, 'r') as file:
            for line in file:
                if "Total time (m)" in line:
                    r2t_total_time = float(line.split("\t")[-1].strip())
        values.append(r2t_total_time if r2t_total_time is not None else "N/A")
    except FileNotFoundError:
        values.append("N/A")
    except Exception as e:
        values.append(f"Error: {e}")

    # Process ERC BXB benchmark file
    try:
        with open(bxb_filepath, 'r') as file:
            for line in file:
                if "Total time (m)" in line:
                    bxb_total_time = float(line.split("\t")[-1].strip())
        values.append(bxb_total_time if bxb_total_time is not None else "N/A")
    except FileNotFoundError:
        values.append("N/A")
    except Exception as e:
        values.append(f"Error: {e}")

    # Process Network Analysis benchmark file
    try:
        last_start_time, last_end_time = None, None
        with open(network_filepath, 'r') as file:
            for line in file:
                if "Network_analysis Call start" in line:
                    last_start_time = datetime.strptime(line.split("\t")[-1].strip(), "%a, %d %b %Y %H:%M:%S")
                elif "Network_analysis Call end" in line:
                    last_end_time = datetime.strptime(line.split("\t")[-1].strip(), "%a, %d %b %Y %H:%M:%S")
        if last_start_time and last_end_time:
            network_runtime = (last_end_time - last_start_time).total_seconds() / 60
        values.append(network_runtime if network_runtime is not None else "N/A")
    except FileNotFoundError:
        values.append("N/A")
    except Exception as e:
        values.append(f"Error: {e}")

    # Process Run_data_counts_log.csv to extract DLCpar count
    try:
        with open(data_counts_filepath, 'r') as file:
            for line in file:
                if "DLCpar" in line:
                    dlcpar_count = int(line.split(",")[-1].strip())
                    break
        values.append(dlcpar_count if dlcpar_count is not None else "N/A")
    except FileNotFoundError:
        values.append("N/A")
    except Exception as e:
        values.append(f"Error: {e}")

    # Function to calculate average of Overlapping_branches in large ERC results files
    def calculate_average_overlapping_branches(filepath, filename):
        total = 0
        count = 0
        try:
            with open(filepath, 'r') as file:
                next(file)  # Skip header line
                for line in file:
                    columns = re.split(r'\t(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)', line.strip())
                    if len(columns) >= 5:
                        try:
                            total += int(columns[4])
                            count += 1
                        except ValueError:
                            pass
            return total / count if count > 0 else None
        except FileNotFoundError:
            return "N/A"
        except Exception as e:
            return f"Error: {e}"

    # Calculate averages for ERC_results_R2T.tsv and ERC_results_BXB.tsv
    erc_r2t_avg_branches = calculate_average_overlapping_branches(erc_r2t_filepath, "ERC_results_R2T.tsv")
    erc_bxb_avg_branches = calculate_average_overlapping_branches(erc_bxb_filepath, "ERC_results_BXB.tsv")

    # Append averages to values list
    values.append(erc_r2t_avg_branches if erc_r2t_avg_branches is not None else "N/A")
    values.append(erc_bxb_avg_branches if erc_bxb_avg_branches is not None else "N/A")

    # Print all values in one line, separated by tabs
    print("\t".join(str(value) for value in values))

# Check if jobname argument is provided
if len(sys.argv) != 2:
    print("Usage: python script_name.py <jobname>")
else:
    jobname = sys.argv[1]
    extract_total_time(jobname)
