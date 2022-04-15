# ERCnet: a program for running genome-wide ERC analyses in the presence of gene duplication and performing ERC-based network analyses.

## Overview
The full ERCnet workflow consistst of the following steps:
1. *Phylogenomic analyses*
1. *Gene-tree/Species-tree reconciliation*
1. *ERC analyses (including branch-length reconciliation)*
1. *Network analyses (including community classification)*

## Preparing input data for ERCnet

ERCnet directly uses the output from [Orthofinder gene family clustering](https://github.com/davidemms/OrthoFinder)

Orthofinder provides [guidelines](https://davidemms.github.io/orthofinder_tutorials/orthofinder-best-practices.html) to selecting taxa and obtaining and preparing proteome files. In addition, please adhere to the following ERCnet-specific guidelines/recommendations to running Orthofinder upstream of ERCnet:

* For easy installation  of Orthofinder, we recommend creating an [anaconda environment](https://www.anaconda.com/) and using conda to install with the following command `conda install -c bioconda orthofinder=2.5.4`(here we specify the version we used for testing ERCnet data)
    * Even if you haven't used anaconda before, it's worth becoming familiar with anaconda environments because ERCnet relies of conda installations for several dependencies.
    * See below for detailed instructions on setting up conda environments to run ERCnet.

* Include a well-defined outgroup in taxon-sampling.
    * ERCnet uses the outgroup to define and seperate taxon-complete subtrees within larger gene families
    * After running Orthofinder, be sure to inspect the inferred rooted species tree (output by Orthofinder
) to ensure that the expected root was used. 
    * If not, reroot the tree manually and input as a user-defined rooted species tree for another run of Orthofinder. Only proceed to ERCnet with Orthofinder results generated under a believable species tree.

* Run Orthofinder with the -y argument to indicate that subtrees be split into seperate subfamilies (Orthofinder HOGs)

* We would recommend running Orthofinder on a super computer if possible. During testing of moderate-sized datasets (~20 proteomes), Orthofinder finished in <24 hours distributed across 48 cores.

Example Orthofinder run:
```
orthofinder -f <path/to/dir/containing/proteomes/> -y -X -M msa -t <number of threads available on computing system>
```

## Running ERCnet

Different ERCnet steps require different dependencies. Most notably, the *Gene-tree/Species-tree reconciliation* step requires python2 (because DLCpar only supports python2) while the other steps require python3. The need to switch between python 2 and 3 environments is a big part of our recommendation to use anaconda environments.

| ERCnet step  | Python version | Dependencies|
| ------------- |:-------------:|:-------------:|
| *Phylogenomic analyses* | 3 | TBD |
| *Gene-tree/Species-tree reconciliation* | 2 | TBD |
| *ERC analyses* | 3 | TBD |
| *Network analyses* | 3 | TBD |

Below we provide detailed instructions for configuring environments for each step of ERCnet.

### 1. Phylogenomic analyses
#### Installing dependencies

Before running the *Phylogenomic analyses* of ERCnet, we recommend installing all dependencies in an anaconda environment by running the following lines of code:

```
#Create a conda environment with python 3
conda create -n test3 python=3.9.7
conda activate test3

#Install python modules
conda install -c conda-forge pandas
conda install -c conda-forge biopython

#Install phylogenetics-related programs
conda install -c bioconda mafft
conda install -c bioconda gblocks
conda install -c bioconda raxml

#Install R package (specific versions included when important)
conda install -c conda-forge r-base=4.1.2
conda install -c conda-forge r-stringr
conda install -c conda-forge r-ape=5.6
conda install -c conda-forge r-phytools
conda install -c conda-forge r-igraph
```

#### Running Phylogenomic analyses
To see a list of arguments for ERCnet use `./Phylogenomics.py --help`

* Input files:
    * Output from Orthofinder found in the folowing folders:
        * Species_Tree/
        * Phylogenetic_Hierarchical_Orthogroups/
        * Resolved_Gene_Trees/
        * Orthogroup_Sequences/
    * An optional species mapping file which will be used to tell ERCnet when tips on the species tree correspond to sequence identifiers in alignments/gene trees.

Create a comma-seperated file that indicates a prefix text string found in the sequence IDs in alignments/gene trees and the corresponding tip on the species tree (found in the Species_Tree/ directory of Orthofinder output). Please be sure the formatting exactly matches the example provided below and text strings you provide exactly match the strings within your alignments/gene trees/species trees.  

Name the species mapping file "Species_mapping.csv" and save in the main ERC net directory (where all the scripts live).

Example:
```
Prefix,SpeciesID
A_aulacocarpa,A_aulacocarpa_prot
A_thaliana,A_thaliana_prot
A_trichopoda,A_trichopoda_prot
C_sativus,C_sativus_prot
E_grandis,E_grandis_prot
G_maderense,G_maderense_prot
G_raimondii,G_raimondii_prot
H_annuus,H_annuus_prot
L_chinense,L_chinense_prot
L_siphilitica,L_siphilitica_prot
M_acuminata,M_acuminata_prot
O_biennis,O_biennis_prot
O_sativa,O_sativa_prot
P_maritima,P_maritima_prot
P_persica,P_persica_prot
P_trichocarpa,P_trichocarpa_prot
S_lycopersicum,S_lycopersicum_prot
S_noctiflora,S_noctiflora_prot
S_polyrhiza,S_polyrhiza_prot
V_vinifera,V_vinifera_prot
```

To run the *Phylogenomic analyses*, you'll need to set parameters that specify level of stringency with which to filter the Orthofinder gene families to be input into the phylogenomic pipeline. To inform this decision, we recommend running Phylogenomics.py with the --explore_filters (-e) flag, which will perform a parameter scan and return a table indicating the number of gene families retained under different filters.

Example:
```
./Phylogenomics.py -j <jobname> -e -s -o <path/to/orthofinder/results/> -x <path/to/raxml/installation/>
```

All options for Phylogenomics.py:

| Short flag  | Long flag | Description| Required? | Default value |
| ------------- |:-------------:|:-------------:|:-------------:|:-------------:|
| -h | --help | Print help menu | no | NA |
| -o | --OFPath | Full path to the Orthofinder results dir (should containSpecies_Tree/, Phylogenetic_Hierarchical_Orthogroups/ etc...) Include "/" at the end of the string | yes | NA |
| -p | --MaxP | Integer: maximum number of paralogs per species allowed ineach gene family | yes, unless -e is chosen | 3* |
| -r | --MinR | Integer: minimum number of species represented required ineach gene family | yes, unless -e is chosen | 10* |
| -t  | --Test_num | Integer: number of gene families to analyze. This option is intended to help you test whether ERCnet is working on your system by running a small subset of genes before running the full dataset | no | NA |
| -e | --explore_filters | Add this flag to explore filtering options (-p and -r parameters). If selected, program will output parameter scan table and quit without running downstream steps. If -e is chosen it will negate -p and -r | no | NA |
| -l | --Min_len | Integer: minimum length of alignment (after trimming with Gblocks) required to retain gene for downstream analyses | yes | 100 |
| -x | --Rax_dir | Full path to the location of your raxml installation (use which raxmlHPC to locate). Include "/" at the end of the string. | yes | NA |
| -s | --SPmap | Add this flag to provide a custom species mapping file. Not required if the tip labels on the orthofinder species tree exactly match the species prefix in sequence IDs. Mapping file must be formatted in certian way. See instuctions | no | NA |
| -n | --Node |Interger: indicate the node on the species tree that you would like to use to retrieve orthofinder HOGs (subtrees). Assuming your species tree has a single outgroup, you'll probably want N1 (default). However, if you species tree has multiple outgroups (or if you'd just like to perform an ERC analysis on a subset of the species tree), you can indicate which node to use for subtree extracting. E.g. For N2.tsv, "-n 2" or "--Node 2"  | no | 1* |

*Values for these parameters can have a profound impact on analyses so make sure the values make biological sense for your analysis before opting for default values.

Use the table output by running the --explore_filters option (above) to choose reasonable values for -p and -r. 

Use the options described above to run the full *Phylogenomic analyses*

Example:
```
./Phylogenomics.py -j <jobname> -e -s -p 3 -r 10 -l 100 -n 1 -o <path/to/orthofinder/results/> -x <path/to/raxml/installation/>
```





