import os
import argparse
import csv
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pyvolve
import random
import copy




# Argument parser to get jobname and optional rand argument from command-line input
def parse_arguments():
    parser = argparse.ArgumentParser(description="Phylogenetic sequence simulation")
    parser.add_argument("-j", "--jobname", required=True, help="Specify the job name for creating output folder")
    parser.add_argument("-r", "--rand", action="store_true", help="Create 100 randomized trees with different sets of accelerated branches")
    parser.add_argument("-i", "--inspect", action="store_true", help="Pause to inspect before proceeding")
    return parser.parse_args()

# Function to create a random bifurcating tree
def create_random_tree(species):
    clades = [Clade(branch_length=random.uniform(0.01, 0.05), name=s) for s in species]

    while len(clades) > 1:
        # Randomly choose two clades to merge
        c1 = clades.pop(random.randint(0, len(clades)-1))
        c2 = clades.pop(random.randint(0, len(clades)-1))
        # Create a new internal clade that joins the two
        new_clade = Clade(branch_length=random.uniform(0.01, 0.05), clades=[c1, c2])
        # Add the new clade back to the list
        clades.append(new_clade)

    # The final remaining clade is the root
    root_clade = clades[0]
    root_clade.branch_length = random.uniform(0.01, 0.05)
    
    # Create the tree
    tree = Tree(root=root_clade)
    
    return tree
'''
# Function to ensure every branch has a branch length (set to 0.0 if missing)
def ensure_branch_lengths(tree):
    for clade in tree.find_clades():
        if clade.branch_length is None or clade.branch_length == 0.0:
            clade.branch_length = random.uniform(0.001, 0.005)  # Add a small default branch length
'''

# Function to set a long branch length for the outgroup (Spec21)
def set_long_branch_for_outgroup(tree, outgroup_name, branch_length):
    # Find Spec21 (the outgroup) and set its own branch length
    for clade in tree.find_clades():
        if clade.name == outgroup_name:
            clade.branch_length = branch_length  # Set long branch length for Spec21 itself
            print(f"Branch length for {outgroup_name} set to {branch_length}")
            return

# Function to save trees with explicit branch lengths, no names
def write_tree_with_explicit_branch_lengths(tree, output_file):
    with open(output_file, 'w') as f:
        for clade in tree.find_clades():
            if clade.branch_length is None:
                clade.branch_length = 0.0  # Ensure every clade has a branch length
        Phylo.write(tree, f, 'newick', branch_length_only=True)

# Function to simulate protein evolution using pyvolve
def simulate_evolution(tree_file, model, seq_length, seqfile):
    # Read the tree in pyvolve format
    pyvolve_tree = pyvolve.read_tree(file=tree_file)
    
    # Define a partition for sequence evolution
    partition = pyvolve.Partition(models=model, size=seq_length)
    
    # Simulate sequence evolution and write to seqfile
    evolver = pyvolve.Evolver(partitions=partition, tree=pyvolve_tree)
    evolver(ratefile=None, seqfile=seqfile)

# Function to split sequences by species and create full proteomes
def split_gene_fam_to_proteomes(gene_fam_dir, full_proteomes_dir, species):
    proteomes = {spec: [] for spec in species}

    for i in range(1, 1001):
        seqfile = os.path.join(gene_fam_dir, f"Gene_fam_{i:04d}.fasta")
        for record in SeqIO.parse(seqfile, "fasta"):
            # Extract species ID from record ID
            spec_id = record.id.split("_")[0]
            # Create a new header SpecXX_Gene_fam_XXXX (no extra ">")
            new_header = f"Spec{spec_id.split('Spec')[1]}_Gene_fam_{i:04d}"
            record.id = new_header
            record.description = ""
            proteomes[spec_id].append(record)

    # Write each species' full proteome to separate files
    for spec, records in proteomes.items():
        output_file = os.path.join(full_proteomes_dir, f"{spec}.fasta")
        SeqIO.write(records, output_file, "fasta")

# Function to create 100 randomized trees with modified branch lengths
def create_randomized_trees(original_tree, rand_dir):
    randomized_trees = []

    for i in range(1, 101):
        # Clone the original tree using deepcopy from the original, unmodified tree
        rand_tree = copy.deepcopy(original_tree)

        # Fetch the clades from the freshly copied tree to ensure independence
        all_clades = list(rand_tree.find_clades())

        random_clades = random.sample(all_clades, 5)  # Randomly select 5 branches

        # Modify branch lengths, starting from the original values
        for clade in random_clades:
            if clade.branch_length:
                original_length = clade.branch_length
                clade.branch_length *= 10  # Multiply branch length by 10
                #print(f"Modified branch {clade.name}: from {original_length} to {clade.branch_length}")

        # Save the randomized tree
        rand_tree_file = os.path.join(rand_dir, f"rand_tree_{i:04d}.nwk")
        write_tree_with_explicit_branch_lengths(rand_tree, rand_tree_file)
        randomized_trees.append(rand_tree_file)

    return randomized_trees


# Function to prompt for user confirmation if --inspect is provided
def inspect_confirmation():
    while True:
        user_input = input("Please inspect the random trees. Do you want to proceed? (y/n): ").strip().lower()
        if user_input == 'y':
            print("Proceeding with the simulation.")
            break
        elif user_input == 'n':
            print("Exiting the script.")
            exit(0)
        else:
            print("Invalid input. Please enter 'y' for yes or 'n' for no.")

# Main function to perform the simulations
def main():
    # Parse the user-defined jobname and optional rand argument
    args = parse_arguments()
    jobname = args.jobname
    rand_mode = args.rand
    output_dir = f"SIM_{jobname}"
    gene_fam_dir = os.path.join(output_dir, "Gene_fam_seqs")
    full_proteomes_dir = os.path.join(output_dir, "Full_proteomes")
    rand_dir = os.path.join(output_dir, "Rand_accelerations")

    # Create output directories
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(gene_fam_dir, exist_ok=True)
    os.makedirs(full_proteomes_dir, exist_ok=True)

    if rand_mode:
        print("Running gene tree accelerations in random mode (evolution across genes not correlated).\n")
        os.makedirs(rand_dir, exist_ok=True)

    # Define the number of species (21)
    species = ["Spec{:02d}".format(i) for i in range(1, 22)]

    # Generate a random tree
    original_tree = create_random_tree(species)

    # Root the tree by Spec21
    original_tree.root_with_outgroup("Spec21")

    # Set the branch leading to the outgroup to a long branch length
    set_long_branch_for_outgroup(original_tree, "Spec21", 0.1)
    
    '''
    # Ensure all branches have branch lengths
    ensure_branch_lengths(original_tree)
    '''
	
    # Save the unmodified tree
    original_tree_file = os.path.join(output_dir, "random_tree_with_branch_lengths.nwk")
    print("Simulated tree with random topology:\n")
    Phylo.draw_ascii(original_tree)
    write_tree_with_explicit_branch_lengths(original_tree, original_tree_file)

    # Initialize tree files for gene family simulation
    modified_tree_files = []

    if rand_mode:
        # Create 100 different randomized trees from the original unmodified tree
        modified_tree_files = create_randomized_trees(original_tree, rand_dir)
    else:
        # Modify 5 random branches on a single tree
        all_clades = list(original_tree.find_clades())
        random_clades = random.sample(all_clades, 5)

        # Modify branch lengths
        for clade in random_clades:
            if clade.branch_length:
                clade.branch_length *= 10

        # Save the modified tree
        modified_tree_file = os.path.join(output_dir, "modified_tree_with_branch_lengths.nwk")
        print("Simulated tree with random topology and modified branch lengths (rate accelerations):\n")
        Phylo.draw_ascii(original_tree)
        write_tree_with_explicit_branch_lengths(original_tree, modified_tree_file)
        modified_tree_files.append(modified_tree_file)
    
    # Check if --inspect is provided
    if args.inspect:
        print("Inspection mode activated.")
        inspect_confirmation()

    # Define model and sequence length for pyvolve
    model = pyvolve.Model("jtt")

    # Simulate protein evolution 1000 times
    for i in range(1, 1001):
        seqfile = os.path.join(gene_fam_dir, f"Gene_fam_{i:04d}.fasta")
        seq_length = random.randint(200, 1000)

        # First 100 simulations from the modified/randomized trees
        if i <= 100:
            if rand_mode:
                # Use randomized trees for the first 100 simulations
                tree_file = modified_tree_files[i-1]
            else:
                # Use the single modified tree for the first 100 simulations
                tree_file = modified_tree_files[0]
            simulate_evolution(tree_file, model, seq_length, seqfile)
        else:
            # Simulate using the unmodified tree for Gene_fam_0101 to Gene_fam_1000
            simulate_evolution(original_tree_file, model, seq_length, seqfile)
        
        #Progress report
        if i % 100 == 0:
        	print(f"Finished simulating {i} genes.\n")

    # Split sequences into full proteomes by species
    print("Creating proteome files from gene families.\n")
    split_gene_fam_to_proteomes(gene_fam_dir, full_proteomes_dir, species)

    # Print summary
    print(f"Simulations completed. All output written to {output_dir}")

if __name__ == "__main__":
    main()
