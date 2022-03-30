# ERCnet: a program for running genome-wide ERC analyses in the presence of gene duplication and performing
 network analyses on ERC results

## Overview
The full ERCnet workflow consistst of the following steps:
1. *Phylogenomic analyses*
1. *Gene-tree/Species-tree reconciliation*
1. *ERC analyses (including branch-length reconciliation)*
1. *Network analyses (including community classification)*

## Preparing input data for ERCnet

ERCnet directly uses the output from [Orthofinder gene family clustering](https://github.com/davidemms/Orth
oFinder)

Orthofinder provides [guidelines](https://davidemms.github.io/orthofinder_tutorials/orthofinder-best-practices.html)to selecting taxa and obtaining and preparing proteome files. In addition, please adhere to the f
ollowing ERCnet-specific guidelines/recommendations to running Orthofinder upstream of ERCnet:

* For easy installation  of Orthofinder, we recommend creating an [anaconda environment](https://www.anaconda.com/) and using conda to install with the following command `conda install -c bioconda orthofinder=2.5.4
` (here we specify the version we used for testing ERCnet data)
    * Even if you haven't used anaconda before, it's worth becoming familiar with anaconda environments bec
ause ERCnet relies of conda installations for several dependencies.
    * See below for detailed instructions on setting up conda environments to run ERCnet.

* Include a well-defined outgroup in taxon-sampling.
    * ERCnet uses the outgroup to define and seperate taxon-complete subtrees within larger gene families
    * After running Orthofinder, be sure to inspect the inferred rooted species tree (output by Orthofinder
) to ensure that the expected root was used. 
    * If not, reroot the tree manually and input as a user-defined rooted species tree for another run of O
rthofinder. Only proceed to ERCnet with Orthofinder results generated under a believable species tree.

* Run Orthofinder with the -y argument to indicate that subtrees be split into seperate subfamilies (Orthof
inder HOGs)

* We would recommend running Orthofinder on a super computer if possible. During testing of moderate-sized 
datasets (~20 proteomes), Orthofinder finished in <24 hours distributed across 48 cores.

Example Orthofinder run:
```
orthofinder -f <path/to/dir/containing/proteomes/> -y -
X -M msa -t <number of threads available on computing system>
```

## Running ERCnet

Different ERCnet steps require different dependencies. Most notably, the *Gene-tree/Species-tree reconcilia
tion* step requires python2 (because DLCpar only supports python2) while the other steps require python3. T
he need to switch between python 2 and 3 environments is a big part of our recommendation to use anaconda e
nvironments.

| ERCnet step  | Python version | Dependencies|
| ------------- |:-------------:|:-------------:|
| *Phylogenomic analyses* | 3 | TBD |
| *Gene-tree/Species-tree reconciliation* | 2 | TBD |
| *ERC analyses* | 3 | TBD |
| *Network analyses* | 3 | TBD |

Below we provide detailed instructions for configuring environments for each step of ERCnet.

### 1. Phylogenomic analyses
#### Installing dependencies

Before running the *Phylogenomic analyses* of ERCnet, we recommend installing all dependencies in an anacon
da environment by running the following lines of code:

```
#Create conda environment
conda create -n ERCnet_py3 python=3

#Activate env
conda activate ERCnet_py3

#Install R package (specific versions included when important)
conda install -c conda-forge r-base
conda install -c conda-forge r-stringr
conda install -c conda-forge r-ape=5.6
conda install -c bioconda r-phytools

#Install python modules
conda install -c anaconda pandas
conda install -c conda-forge biopython

#Install additional programs
conda install -c bioconda mafft
conda install -c bioconda gblocks
conda install -c bioconda raxml
```

#### Runnng Phylogenomic analyses
To see a list of arguments for ERCnet use `Phylogenomics.py --help`

* Input files:
    * Output from Orthofinder found in the folowing folders:
        * Species_Tree/
        * Phylogenetic_Hierarchical_Orthogroups/
        * Resolved_Gene_Trees/
        * Orthogroup_Sequences/
    * An optional species mapping file which will be used to tell ERCnet when tips on the species tree corr
espond to sequence identifiers in alignments/gene trees.

Create a comma-seperated file that indicates a prefix text string found in the sequence IDs in alignments/g
ene trees and the corresponding tip on the species tree (found in the Species_Tree/ directory of Orthofinde
r output). Please be sure the formatting exactly matches the example provided below and text strings you pr
ovide exactly match the strings within your alignments/gene trees/species trees.  

Name the species mapping file "Species_mapping.csv" and save in the main ERC net directory (where all the s
cripts live).

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
