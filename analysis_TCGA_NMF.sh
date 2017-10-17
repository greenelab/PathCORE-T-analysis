#!/bin/bash
set -o errexit

# Kathleen Chen 2017

tcga_data_dir="./data/tcga_data"
# path to the PID pathway file
pathway_file=$tcga_data_dir"/PID_pathway_definitions.txt"

nmf_dir=$tcga_data_dir"/NMF_model"

# specify the analysis output directory
analysis_dir=$tcga_data_dir"/NMF_analysis"
mkdir -p $analysis_dir

# network output directory will contain the generated network file
network_output_dir=$analysis_dir"/network_construction"

N_features=300

# number of cores to use in parallel for the Python analysis scripts,
# please modify this to accommodate your need
N_cores=1

# if the gene's weight is above this number of standard deviations from the
# mean, we will include it in the gene signature of the feature
std_cutoff=2.0

# significance level used during both the pathway overrepresentation analysis
# and the edge permutation test
alpha=0.05

###############################################################################
# PathCORE-T analysis
###############################################################################

mkdir -p log

# Generates the significant pathways files and metadata from the
# feature overrepresentation analysis on the model in $network_output_directory.
python run_network_creation.py $nmf_dir $network_output_dir $pathway_file \
                               --n-features=$N_features \
                               --signature=NMF --signature-args=$std_cutoff \
                               --alpha=$alpha --n-cores=$N_cores \
                               --shorten=TCGA_PID \
                               --overlap-correction --all-genes \
> ./log/analysis_TCGA_NMF.log

# Builds a network of pathway co-occurrence relationships for
# each of the models and applies the permutation test to these.
# Results in a single aggregate network where the edges in the
# network were considered significant under the permutation test.
N_permutations=10000

permutation_output_dir=$analysis_dir"/permutation_test_N="$N_permutations

python run_permutation_test.py $network_output_dir $permutation_output_dir \
                               --n-permutations=$N_permutations \
                               --n-features=$N_features \
                               --n-cores=$N_cores \
>> ./log/analysis_TCGA_NMF.log

wait
echo "The TCGA NMF analysis has finished."
exit 0
