import pandas as pd
import glob
import re

def extract_numbers(folder_name):
    """Extracts numerical values from folder names like 'OUT_spec10_rep01' for sorting."""
    match = re.search(r"OUT_spec(\d+)_rep(\d+)", folder_name)
    if match:
        spec_num = int(match.group(1))  # Extracts spec number (e.g., 10, 15, 20, 25)
        rep_num = int(match.group(2))   # Extracts rep number (e.g., 01, 02, 03)
        return (spec_num, rep_num)      # Sort by spec first, then rep
    return (float("inf"), float("inf"))  # If no match, send to end

# Get all matching directories and sort them
dirs = sorted(glob.glob("OUT_spec*/"), key=extract_numbers)

# Set a lower per-sheet row limit to avoid exceeding Excel's max
MAX_ROWS = 1_000_000  # Reduced from 1,048,576 to stay safely within limits

# Create a new Excel writer
with pd.ExcelWriter("combined.xlsx", engine="openpyxl") as writer:
    for dir in dirs:
        print(f"Processing {dir}")

        # Get the first matching file for each pattern
        BXB_edges = glob.glob(dir + "Network_analyses/Text_network_edges_Filtered_ERC_results_BXB*")
        R2T_edges = glob.glob(dir + "Network_analyses/Text_network_edges_Filtered_ERC_results_R2T*")
        BXB_nodes = glob.glob(dir + "Network_analyses/Text_network_vertices_Filtered_ERC_results_BXB*")
        R2T_nodes = glob.glob(dir + "Network_analyses/Text_network_vertices_Filtered_ERC_results_R2T*")

        def save_dataframe(df, base_sheet_name):
            """Splits DataFrame into multiple sheets if it exceeds Excel's row limit."""
            num_sheets = (len(df) // MAX_ROWS) + 1  # Calculate number of sheets required
            for i in range(num_sheets):
                start_row = i * MAX_ROWS
                end_row = min((i + 1) * MAX_ROWS, len(df))  # Ensure we don't go beyond df length

                # Add "part1" only if there's more than one sheet
                sheet_suffix = f"_part{i+1}" if num_sheets > 1 else ""
                
                # Sheet name constraint: Excel allows max 31 characters
                sheet_name = f"{base_sheet_name}{sheet_suffix}"[:31]

                # Extract subset, reset index to avoid row number issues
                df_chunk = df.iloc[start_row:end_row].reset_index(drop=True)

                # Write to the Excel file
                df_chunk.to_excel(writer, sheet_name=sheet_name, index=False)

        # Load and save each file while checking row limits
        if BXB_edges:
            df_B_E = pd.read_csv(BXB_edges[0], sep="\t")
            df_B_E.drop(columns=["Node_labels", "color", "Functional_cat_code"], errors="ignore", inplace=True)
            save_dataframe(df_B_E, dir.replace("/", "BXB_edges").replace("OUT_", ""))

        if BXB_nodes:
            df_B_N = pd.read_csv(BXB_nodes[0], sep="\t")
            df_B_N.drop(columns=["Node_labels", "color", "Functional_cat_code"], errors="ignore", inplace=True)
            save_dataframe(df_B_N, dir.replace("/", "BXB_nodes").replace("OUT_", ""))

        if R2T_edges:
            df_R_E = pd.read_csv(R2T_edges[0], sep="\t")
            df_R_E.drop(columns=["Node_labels", "color", "Functional_cat_code"], errors="ignore", inplace=True)
            save_dataframe(df_R_E, dir.replace("/", "R2T_edges").replace("OUT_", ""))

        if R2T_nodes:
            df_R_N = pd.read_csv(R2T_nodes[0], sep="\t")
            df_R_N.drop(columns=["Node_labels", "color", "Functional_cat_code"], errors="ignore", inplace=True)
            save_dataframe(df_R_N, dir.replace("/", "R2T_nodes").replace("OUT_", ""))

print("TSV files combined into 'combined.xlsx'")
