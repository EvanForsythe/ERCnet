# ERCnet: a program for running genome-wide ERC analyses in the presence of gene duplication and performing network analyses on ERC results

## Overview
The full ERC net workflow consistst of the fowlling groups of steps:
1. Preprocessing steps
1. Phylogenomic steps
1. Branch length reconciliation steps
1. ERC analysis steps
1. Network analysis steps

## Installing dependencies 

Before running ERCnet, we recommend installing all dependencies in an anaconda environment by running the following lines of code:

```
#Create conda environment
conda create -n ERCnet_env python=3

#Activate env
conda activate ERCnet_env

#Install R package (specific versions included when important)
conda install -c conda-forge r-base
conda install -c conda-forge r-stringr
conda install -c conda-forge r-ape=5.6

#Install python modules
conda install -c anaconda pandas
conda install -c conda-forge biopython

#Install additional programs
conda install -c bioconda mafft
conda install -c bioconda gblocks
conda install -c bioconda raxml
```

## Running ERCnet
To see a list of arguments for ERCnet use `ERCnet_main.py --help`



