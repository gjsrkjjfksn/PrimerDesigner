
# PrimerDesigner

## Introduction

This repository contains the source code for the PrimerDesigner project. 
The main goal of this project is to create an algorithm for primer design for protein synthesis using microarray oligonucleotide probes. 
The algorithm aims to find the most efficient primer set with complete coverage and no cross hybridization risk.
It uses a primer graph to represent all valid forward and reverse primer combinations and Integer Linear Programming (ILP) with forbidden pair and single path constraints to find an optimal solution

## Requirements

The tool has been tested with the following configuration on a Linux machine:
- Python 3.9.18
- GurobiPy 11.0.2
- NetworkX 3.2.1
- NumPy 1.26.4
- Pandas 2.1.2


## Setting Up the Tool

First, create a file named `gurobi.json` containing the details for the Gurobi license and put it in the main directory.

```json
{
  "WLSACCESSID": "XXXXX",
  "WLSSECRET": "XXXXX",
  "LICENSEID": 12345
}
```

Create a file containing the protein coding-sequences and their names. Each line should contain a protein's name and its DNA coding-sequence, separated by a tab.


## Running the PrimerDesigner Algorithm

To execute the PrimerDesginer algorithm, use the following command:

```bash
python ./tool.py --file_path <file-path> --version <version> 
```
- **file_path**: The file path of the protein coding-sequences
- **version**: Specifies which version of the algorithm to run. Options: Relaxed, Non-relaxed or the Extension (default: Relaxed)
  
The other arguments are optional and include the algorithm parameters:

- **primer_lmin, primer_lmax**: Minimum and maximum primer lengths (default: 18,30).
- **oligo_lmin, oligo_lmax**: Minimum and maximum oligonucleotides lengths (default: 195,205)
- **overlap_lmin, overlap_lmax**: Minimum and maximum overlap length between oligonucleotides  (default: 45,50)
- **allowed_overlap**: Allowed overlap between primer_pairs (default: 6)
- **num_proteins**: Number of variants of the same sequence - used for Relaxed version only (default: 3)
- **apply_threshold**: Boolean flag for applying primer quality threshold (default: False)
- **min_gc, max_gc**: Minimum and maximum threshold on gc content (default: 40,60)
- **min_tm, max_tm**: Minimum and maximum threshold on melting temperature tm (default: 58,65)
- **max_difference**: Maximum difference threshold between the forward and reverse primer in each pair (default: 3)
- **merge_bins**: Boolean flag for merging bins corresponding to identical non-overlapping sequences - used for the Relaxed version (default: False)


Example command:
```bash
python ./tool.py --file_path example_proteins.txt --version Non_relaxed --primer_lmin 20 --primer_lamx 26
```


