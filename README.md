# Supplementary data: evolution of expression noise in gene regulatory networks
This repository contains the code to reproduce the statistical analysis of the following publication (https://doi.org/10.1371/journal.pcbi.1010982):
"Being noisy in a crowd: differential selective pressure on gene expression noise in model gene regulatory networks"

Evolutionary simulations were used to study the evolution of gene-specific expression noise in gene regulatory networks and local and global network properties were found to affect the evolvability of gene-specific expression noise.

## Installation and data download
Clone this repo:

```
$ git clone git@gitlab.gwdg.de:molsysevol/supplementarydata_expressionnoise.git
$ cd supplementarydata_expressionnoise
```

Download the zipped raw data from zenodo (https://zenodo.org/records/6939845) in this folder and unzip.

```
$ unzip results.zip
```

## Reproducing the results
The full code to reproduce the results in the maintext, supplementary information, and additional tests has been compiled into a shell script called pipe_analysis.sh.
Run the analysis script:
```
$ ./pipe_analysis.sh
```
Running the script will knit Rmd files into HTML files and output figures in the plots folder.

## Regenerating the simulation data
The scripts to generate network topologies, perform network analysis, calculate noise metrics and preprocess the data have been compiled into a shell script called pipe_regenerate_data.sh.
Note that to regenerate the simulation data you will have to install the simulator and adapt the scripts to your HPC system.

## Simulator
The code used for running the simulations presented in the paper is not included here, but is available in this repo: https://github.com/NatashaPuzovic/netlings. The version of the program used in this study is tagged as v1.2.0 Release.
