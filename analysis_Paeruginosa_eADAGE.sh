#!/bin/bash
set -o errexit

# Kathleen Chen 2017

# These are the filepaths and constants used in both this script
# and the web application setup script (./web_db_Paeruginosa_eADAGE.sh)
source ./common_Paeruginosa_eADAGE.sh

# path to the eADAGE models directory
models_dir=$pa_data_dir"/ensemble_models"

# count number of genes
N_genes=$(expr $(wc -l <$data_compendium) - 1)

# path to the list of genes we have data for in the compendium
pa_genes=$pa_data_dir"/pao1_compendium_genes.txt"
tail -n+2 $data_compendium | awk -F"\t" '{print $1}' > $pa_genes

# the model size that best supports the current P.a expression compendium k=300
N_features=300

# number of cores to use in parallel in run_pathcore_analysis.py,
# please modify this to accommodate your need.
N_cores=1

# if the gene's weight is +/- this number of standard deviations from the
# mean, we will include it in the +/- gene signatures of a feature
std_cutoff=2.5

# significance level used during both the pathway overrepresentation analysis
# and the edge permutation test
alpha=0.05

###############################################################################
# PathCORE-T analysis
###############################################################################

mkdir -p log

# Overrepresentation analysis and network construction:
# Generates the significant pathways files, network files, and metadata from
# running the analysis on all models in $ensemble_directory.
python run_network_creation.py $models_dir $network_output_dir $pathway_file \
                               --n-genes=$N_genes --n-features=$N_features \
                               --signature=eADAGE --signature-args=$std_cutoff \
                               --alpha=$alpha --genes-list=$pa_genes \
                               --n-cores=$N_cores --shorten=PAO1_KEGG \
                               --overlap-correction --all-genes --metadata \
> ./log/analysis_PAO1_eADAGE.log

# Builds a network of pathway co-occurrence relationships for
# each of the models and applies the permutation test to these.
# Results in a single aggregate network where the edges in the
# network were considered significant under the permutation test.
python run_permutation_test.py $network_output_dir $permutation_output_dir \
                               --n-permutations=$N_permutations \
                               --n-features=$N_features \
                               --n-cores=$N_cores \
>> ./log/analysis_PAO1_eADAGE.log

wait
echo "The PAO1 eADAGE analysis has finished."
exit 0
