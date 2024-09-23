#!/bin/bash

###Usage:
#sbatch -p forsythe.q -A forsythe --job-name=<jobname> --output=<jobname>.out --error=<jobname>.err auto_ortho_batch.sh sim_ortho_cont <proteomes_path> <ortho_output_path> <root_species>

# Assign command line arguments to variables
job_name=$1
input_path=$2
output_path=$3
outgroup=$4

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=96
#SBATCH --time=96:0:0
#SBATCH --mail-user=evan.forsythe@osucascades.edu
#SBATCH --mail-type=END

submitted_command="sbatch -p forsythe.q -A forsythe --job-name=$job_name --output=${job_name}.out --error=${job_name}.err ortho_batch.sh $job_name $input_path $output_path $outgroup"

echo "Simulated submitted command: $submitted_command"
echo "Job name: $job_name"
echo "Input path: $input_path"
echo "Output path: $output_path"

source activate orthofinder-env

echo "starting orthofinder step1"
echo `date "+%Y-%m-%d %H:%M:%S"`

# Print and run the orthofinder command
echo "Running orthofinder command: orthofinder -f $input_path -y -X -M msa -t 96 -o $output_path -ot"
orthofinder -f $input_path -y -X -M msa -t 96 -o $output_path -ot

echo "finishing orthofinder step 1"
echo `date "+%Y-%m-%d %H:%M:%S"`

# Reroot the species tree
echo "locating the species tree from orthofinder"

# Assume $output_path is defined
results_dir=$(find "$output_path" -maxdepth 1 -type d -name 'Results*' | head -n 1)

# Example use of $results_dir
if [ -d "$results_dir" ]; then
    echo "Processing files in $results_dir"
    # Further commands here that use $results_dir
else
    echo "No Results directory found"
    exit 1
fi


# Reroot the tree
echo "Rerooting the tree with outgroup: $outgroup"
nw_reroot $results_dir/Species_Tree/SpeciesTree_rooted.txt "$outgroup" > $output_path/Newly_rooted_species_tree.txt

if [ $? -eq 0 ]; then
    echo "Tree successfully rerooted and saved"
else
    echo "Error during rerooting"
    exit 1
fi

echo "starting orthofinder step2"
echo `date "+%Y-%m-%d %H:%M:%S"`

# Run the second step
# Print and run the orthofinder command
echo "Running orthofinder command: orthofinder -ft $output_path -s $output_path/Newly_rooted_species_tree.txt -y -X -M msa -t 96"
orthofinder -ft $results_dir -s $output_path/Newly_rooted_species_tree.txt -y -X -M msa -t 96

echo "finishing orthofinder step 2"
echo `date "+%Y-%m-%d %H:%M:%S"`

#Transfer needed files from step 1 orthofinder to step 2 orthofinder

echo "Transferring orthofinder results dir to the new folder"

# Construct the name of the second directory by appending "_1"
results_dir_1="${results_dir}_1"

# Ensure both directories exist
if [ ! -d "$results_dir" ]; then
    echo "Directory $results_dir does not exist."
    exit 1
fi

if [ ! -d "$results_dir_1" ]; then
    echo "Directory $results_dir_1 does not exist."
    exit 1
fi

# Use rsync to copy missing directories from $results_dir to $results_dir_1
rsync -av --ignore-existing "$results_dir/" "$results_dir_1/"

if [ $? -eq 0 ]; then
    echo "Folders from $results_dir copied to $results_dir_1 (if missing)."
else
    echo "Error during the copy process."
    exit 1
fi

