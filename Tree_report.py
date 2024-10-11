import argparse
import os
from Bio import Phylo

def parse_args():
    parser = argparse.ArgumentParser(description="Report on phylogenetic trees")
    parser.add_argument("jobname", help="The job name")
    parser.add_argument("hog_id", help="The HOG ID")
    return parser.parse_args()

def load_tree(file_path):
    """Load the phylogenetic tree from a file."""
    try:
        tree = Phylo.read(file_path, "newick")
        return tree
    except Exception as e:
        print(f"Error loading tree from {file_path}: {e}")
        return None

def report_on_tree(tree, file_path):
    """Generate a report on the given phylogenetic tree."""
    if tree is None:
        print(f"Cannot generate report for {file_path}. Tree could not be loaded.")
        return

    print(f"\nReport for tree file: {file_path}")
    
    # Whether the tree is rooted
    is_rooted = tree.rooted
    print(f"Is the tree rooted? {'Yes' if is_rooted else 'No'}")
    
    # Count non-bifurcating nodes (multifurcations)
    multifurcating_nodes = sum(len(clade.clades) > 2 for clade in tree.get_nonterminals())
    print(f"Number of non-bifurcating nodes: {multifurcating_nodes}")
    
    # Whether there are node labels
    node_labels_exist = any(clade.name for clade in tree.find_clades())
    print(f"Are there node labels? {'Yes' if node_labels_exist else 'No'}")
    
    # Shortest and longest branch
    branch_lengths = [clade.branch_length for clade in tree.find_clades() if clade.branch_length is not None]
    if branch_lengths:
        shortest_branch = min(branch_lengths)
        longest_branch = max(branch_lengths)
        print(f"Shortest branch length: {shortest_branch}")
        print(f"Longest branch length: {longest_branch}")
    else:
        print("No branch lengths available.")
    
    # Count branches with length = 0
    zero_length_branches = sum(1 for clade in tree.find_clades() if clade.branch_length == 0)
    print(f"Number of branches with length = 0: {zero_length_branches}")

def main():
    # Parse the command-line arguments
    args = parse_args()
    
    # Define the paths for the tree files
    hog_subtree_path = os.path.join(f"OUT_{args.jobname}", "HOG_subtrees", f"{args.hog_id}_tree.txt")
    bs_tree_path = os.path.join(f"OUT_{args.jobname}", "BS_trees", f"{args.hog_id}_BS.treefile")
    re_tree_path=os.path.join(f"OUT_{args.jobname}", "Rearranged_trees", f"{args.hog_id}_BS.treefile_recs.nwk")
    bl_tree_path = os.path.join(f"OUT_{args.jobname}", "BL_trees", f"{args.hog_id}_BL.treefile")


    # Load and report on each tree
    for tree_path in [hog_subtree_path, bs_tree_path, re_tree_path, bl_tree_path]:
        tree = load_tree(tree_path)
        report_on_tree(tree, tree_path)

if __name__ == "__main__":
    main()
