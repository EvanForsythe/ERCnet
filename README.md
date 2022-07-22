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
| *Phylogenomic analyses* | 3 | pandas, biopython, mafft, gblocks, raxml, stringr, ape, phytools |
| *Gene-tree/Species-tree reconciliation* | 2 | dlcpar |
| *ERC analyses* | 3 | ape, stringr, phytools |
| *Network analyses* | 3 | igraph |

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
conda install -c bioconda treerecs
conda install -c conda-forge joblib

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

Example Species_mapping.csv:
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

You may also optionally provide a list of genes you'd like to study. To do so, use the --A_priori/a flag (see below). To use this option you must create a file named "A_priori_genes.csv", which should have one column containing text strings unique to the sequence ID of your genes of interest. The first line of the file should be the species idenifier (see SpeciesID in the table above) of the species that contains the a priori genes. 

Name the a priori genes file "A_priori_genes.csv" and save in the main ERC net directory (where all the scripts live).

Example A_priori_genes.csv:
```
A_thaliana_prot
AT1G49970
AT1G12410
AT1G09130
AT4G17040
AT3G17000
AT5G17380
AT5G48020
AT1G06550
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
| -j | --JOBname | Unique job name for this run of ERCnet. Avoid including spaces or special characters ("_" is ok) | yes | NA |
| -o | --OFPath | Full path to the Orthofinder results dir (should containSpecies_Tree/, Phylogenetic_Hierarchical_Orthogroups/ etc...) Include "/" at the end of the string | yes | NA |
| -p | --MaxP | Integer: maximum number of paralogs per species allowed ineach gene family | yes, unless -e is chosen | 3* |
| -r | --MinR | Integer: minimum number of species represented required ineach gene family | yes, unless -e is chosen | 10* |
| -t  | --Test_num | Integer: number of gene families to analyze. This option is intended to help you test whether ERCnet is working on your system by running a small subset of genes before running the full dataset | no | NA |
| -e | --explore_filters | Add this flag to explore filtering options (-p and -r parameters). If selected, program will output parameter scan table and quit without running downstream steps. If -e is chosen it will negate -p and -r | no | NA |
| -l | --Min_len | Integer: minimum length of alignment (after trimming with Gblocks) required to retain gene for downstream analyses | yes | 100 |
| -x | --Rax_dir | Full path to the location of your raxml installation (use which raxmlHPC to locate). Include "/" at the end of the string. | yes | NA |
| -s | --SPmap | Add this flag to provide a custom species mapping file. Not required if the tip labels on the orthofinder species tree exactly match the species prefix in sequence IDs. Mapping file must be formatted in certian way. See instuctions | no | NA |
| -n | --Node |Interger: indicate the node on the species tree that you would like to use to retrieve orthofinder HOGs (subtrees). Assuming your species tree has a single outgroup, you'll probably want N1 (default). However, if you species tree has multiple outgroups (or if you'd just like to perform an ERC analysis on a subset of the species tree), you can indicate which node to use for subtree extracting. E.g. For N2.tsv, "-n 2" or "--Node 2"  | no | 1* |
| -m | --Mult_threads |Integer: number of threads avilable for parallel computing (default = 1). Performing a full-genome analyses will likely require supercomputing resources| no | 1 |
| -a | --Apriori | Add this flag to provide an *a priori* list of genes to analyze. The list must be in a file named "A_priori_genes.csv" and formatted in a specific way. See instructions above for more information. When you're using the -a option you probably don't want to use the -t option | no | NA |
| -c | --core_distribution | Sets the group (1/2/3) core_distribution strategy for Raxml processing. This determines how many 'front end' parallel cores are running and how many 'back end' cores are given to each multi-threaded Raxml process. See table below for group definitions | no | 1 |

core_distribution groups for Phylogenomics.py. 

| Group | Front End Cores | Back End Cores | 
| ------------- |:-------------:|:-------------:|
| 1 | num_cores / 2 | 2 |
| 2 | 2 | num_cores / 2 |
| 3 | num_cores / 4 | 4 |

*'Front End' parallelization refers to the paralellization support from Joblib, indicating the number of processes inside any given Phylogenetics computational step are being spawned.

*'Back End' parallelization specifically refers to the number of cores that are passed to RAXML's native multi-threading support. 


*Values for these parameters can have a large impact on analyses so make sure the values make biological sense for your analysis before opting for default values.

Use the table output by running the --explore_filters option (above) to choose reasonable values for -p and -r. 

Use the options described above to run the full *Phylogenomic analyses*

Example:
```
./Phylogenomics.py -j test_job -e -s -p 3 -r 10 -l 100 -n 1 -o <path/to/orthofinder/results/> -x <path/to/raxml/installation/>
```

Brief(ish) walkthrough of what Phylogenomics.py does:
* Filters orthofinder gene families (called HOGs) to remove those that underwent too much gene duplication and/or loss (the user defines what 'too much' means by setting the -r and -p parameters)
   * The user can use the -e argument to make an informed decision on the tradeoffs of -r and -p choice
* Mulitple sequence alignment of gene families (sometimes called 'subtrees' because orthofinder extracts them from larger gene families).
* Gblocks trimming of alignments to remove poorly-aligned sites
* Phylogenetic bootstrap analysis to generate bootstrap support scores for the gene trees generated by Orothfinder
* Rearrange/correct poorly-supported (<80% bootstrap support) branches and root the tree
   * Rearrangement and rooting are both accomlished by gene-tree/species-tree reconciliation in Treerecs
* Optimize branch lengths with raxml using the the Gblocks-trimmed alignments and the new rearranged trees as constraint trees.


## 2. Gene-tree/Species-tree reconciliation
This step uses DLCpar GT/ST reconciliation to 'map' the nodes on the gene trees to appropriate nodes on the species tree. This mapping information is necessary for the subsequence branch-length reconciliation step (which occurs later during the ERC analysis step of ERCnet). Note that GT/ST reconciliation was used in a different context (correcting/rooting trees) in the Phylogenomics step above. 

#### Installing dependencies
This step uses DLCpar, which is written in python2. This means you'll need to setup a python2 environment to run this step of the analysis. Below are instuctions for setting up the environment and installing DLCpar with anaconda
```
#Create python2 env
conda create --name dlcpar_py27 python=2.7
#Activate it
conda activate dlcpar_py27
#Install dlcpar
conda install -c bioconda dlcpar
```
#### Running Gene-tree/Species-tree reconciliation analyses
All options for GTST_reconciliation.py:

| Short flag  | Long flag | Description| Required? | Default value |
| ------------- |:-------------:|:-------------:|:-------------:|:-------------:|
| -h | --help | Print help menu | no | NA |
| -j | --JOBname | Unique job name for this run of ERCnet. This should be the exact same as the jobname used for the Phylogenomics step | yes | NA |

Example command:
```
./GTST_reconciliation.py -j test_job
```

## 3. ERC analyses
Here we extract branch lengths from the gene trees generated above and compare those branch lengths across each pair of genes to ask if those genes show signs of coevolving. Extracting the branch lengths is not trivial, requiring the development of a new approach (described below).

#### Installing dependencies
The rest of ERCnet using python3 so you'll need to switch back to the original anacona env you setup above
```
#List available envs as a reminder
conda info --envs
#Activate the previously setup env (named test3 if you followed along above)
conda activate test3
```
#### Running ERC analyses
All options for ERC_analyses.py:

| Short flag  | Long flag | Description| Required? | Default value |
| ------------- |:-------------:|:-------------:|:-------------:|:-------------:|
| -h | --help | Print help menu | no | NA |
| -j | --JOBname | Unique job name for this run of ERCnet. This should be the exact same as the jobname used in the previous steps | yes | NA |

Example command:
```
./ERC_analyses.py -j test_job
```

What ERC_analyses.py does:
* Branch-length reconciliation
The overall goal of ERCnet is to compare branch lengths (i.e. rates of evolution) between different gene trees. However, when two gene trees have been subject to different histories of duplication/loss (even low levels), they do not have the same branches, so there is not a clear 'apples-to-apples' comparison. In the past, this challenge has been side-stepped by using the length of paths of several branches leading from the root of the tree to a given tip (i.e. root-to-tip approach). This approach is not ideal and can introduce error (both type I and type II). Therefore, we developed a strategy for measuring the rate of evolution along individual branches of the tree (i.e. branch-by-branch). The method is build on similar logic to GT/ST reconciliation, which recognizes that all gene trees have evolved within an over-arching species tree. Therefore, the species tree can be used as a common denominator for making 'apples-to-apples' comparisons between gene trees with vastily different duplication/loss histories. To accomplish this, we use the branches on each gene tree to measure the amount of evolution that occured along branches of the species tree. Sometimes multiple gene tree branches existed within in single species tree branch (e.g. two paralogs both evolving in the same species), meaning we average the branch-lengths. Sometimes a gene tree does not have branches that speak to the evolution in a species tree branch (e.g. the gene was lost in a particular species), meaning the branch legnth is NA for that particular gene. As long and we can successfully map the gene tree branches to the species tree branch, we can extract all available branch length information for use in ERC analyses.
* Root-to-tip branch length extraction
   * Because it's relatively simple, we extract root-to-tip branch lengths in addition to the branch-by-branch approach described above. The user can decide which to use in downstream analyses.
* Normalization of branch lengths
   * branch lengths in each gene are normalized by the genome-wide
* Pair-wise all-by-all ERC analyses
   * Branch length correlation analyses are performed between all combinations of genes

## 4. Network analyses
The all-by-all nature of the ERC results mean that graph theory is a useful framework to explore interactions. Here, we generate networks representing significant 'ERC hits'. The user difines the type of correlation statistic and the cutoff values for what to consider 'significant'.

#### Installing dependencies
This step makes use of the R package igraph. If you followed the instructions above, igraph should already be installed on your python3 environment. 

#### Running Network analyses
All options for Network_analyses.py:

| Short flag  | Long flag | Description| Required? | Default value |
| ------------- |:-------------:|:-------------:|:-------------:|:-------------:|
| -h | --help | Print help menu | no | NA |
| -j | --JOBname | Unique job name for this run of ERCnet. This should be the exact same as the jobname used in the previous steps | yes | NA |
| -m | --BLmethod | Branch length method ERC results to be used in the network. "bxb" for Branch-by-branch. "r2t" for root-to-tip. | yes | NA |
| -f | --Filterstat | Correlation statistic to be used to filter ERC hits. "pval" for p-value, "R2" for R-squared | yes | NA |
| -c | --Cutoff | Cuttoff P/R-squared value by which to filter ERC hits in network. float between 0 and 1. Correlations below/above this number will br retained (depending on your choice of pval vs R2). | yes | NA |
| -y | --Clustmeth | Clustering method to be used to identify communities in network. "fg" for fast-and-greedy (fastest), "eb" for edge-betweenness, "op" for optimal, and "wt" for walktrap. | yes | fg |
| -t | --Trim_Cutoff | Must be an integer. Indicates the minimum number of nodes necessary for a community to be displayed on the network plot. Communitiies smaller than this number will be trimmed from the graph (and associated output tables). This option is mainly for network aesthetics. | no | 0 |
| -s | --FocalSP | The name of the focal species to represent each gene family (should exactly match the tip label of the species tree). See further description below | yes | NA |

Example command:
```
./Network_analyses.py -j test_job -m bxb -f pval -c 0.001 -y fg -t 3 -s A_thaliana_prot
```

What Network_analyses.py does:
* Filter the ERC results to retain only the 'significant' correlations. 
   * The user uses the -m argument to choose which branch lengths should be used for correlation analyses 
   * The user uses the -f argument to choose whether to use the p-value or R-squared value
   * The user uses the -c argument to choose what cutoff value to use. Note that -f=pval means lower values are more stringent whereas -f=R2 means higher values are more stringent. 
* Generate a network diagram from ERC results
   * Nodes represent genes 
   * edges represent significant ERC correlation between genes
* Cluster communities of connected genes
   * Clustering is a complicated task in graph theory. Currently we provide four different algorithms, which the user selects with the -y argument.
* Extract gene names associated with communities
   * A common downstream analysis would be to ask if the genes within a community are enriched for a particular function. To do this, you'll need a gene ID to represnt each gene in the network (technically the nodes represnt gene families/trees). The user can set which species is the best model organism to represent the gene family using the -s argument. We recommend using the species that has the best functional annotations. 
* Extract other global network statistics (TBD) 





