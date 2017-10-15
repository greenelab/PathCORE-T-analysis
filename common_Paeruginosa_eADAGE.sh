#!/bin/bash
set -o errexit

# Kathleen Chen 2017

# This file stores the variables used in both the
# ./analysis_Paeruginosa_eADAGE.sh and
# ./web_db_Paeruginosa_eADAGE.sh scripts.

pa_data_dir="./data/pao1_data"
# path to the data compendium file
data_compendium=$pa_data_dir"/all-pseudomonas-gene-normalized.pcl"
# path to the KEGG pathway file
pathway_file=$pa_data_dir"/pseudomonas_KEGG_terms.txt"

# specify the analysis output directory
analysis_dir=$pa_data_dir"/eADAGE_analysis"
mkdir -p $analysis_dir

# network output directory will contain (per model) the generated network file
# and, if applicable, 3 metadata files in a directory `metadata`:
#   - overrepresented pathways in each feature
#   - pathway definitions for each feature after crosstalk removal
#   - positive and negative gene signature sets for each feature
network_output_dir=$analysis_dir"/network_construction"

# specify parameters specific to the permutation test
N_permutations=10000
permutation_output_dir=$analysis_dir"/permutation_test_n="$N_permutations
