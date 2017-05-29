#!/bin/bash                                                                                                                                      
set -o errexit

# Kathleen Chen 2017

pa_data_dir="./data/pao1_data"
# path to the data compendium file
data_compendium=$pa_data_dir"/all-pseudomonas-gene-normalized.pcl"
# path to the eADAGE models directory
models_dir=$pa_data_dir"/ensemble_models"
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
#mkdir -p $network_output_dir

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
# PathCORE analysis
###############################################################################

mkdir -p log

# Overrepresentation analysis and network construction:
# Generates the significant pathways files, network files, and metadata from
# running the analysis on all models in $ensemble_directory.
python run_network_creation.py $models_dir $network_output_dir $pathway_file \
--n-genes=$N_genes --n-features=$N_features \
--signature=eADAGE --signature-args=$std_cutoff \
--alpha=$alpha --genes-list=$pa_genes --n-cores=$N_cores --shorten=PAO1_KEGG \
--overlap-correction --all-genes --metadata \
> ./log/analysis_PAO1_eADAGE.log

# Builds a network of pathway co-occurrence relationships for
# each of the models and applies the permutation test to these.
# Results in a single aggregate network where the edges in the
# network were considered significant under the permutation test.
N_permutations=10000

permutation_output_dir=$analysis_dir"/permutation_test_n="$N_permutations

python run_permutation_test.py \
$network_output_dir $permutation_output_dir \
--n-permutations=$N_permutations --n-features=$N_features \
--n-cores=$N_cores \
>> ./log/analysis_PAO1_eADAGE.log

###############################################################################
# PAO1 <> eADAGE <> KEGG web application DB setup
###############################################################################

filtered_network=$permutation_output_dir"/filtered_network.tsv"

# NOTE: this file must be edited after a user has set up an MLab account.
db_credentials_file="./eadage.yml"

python web_initialize_db.py \
$data_compendium $pathway_file $filtered_network $network_output_dir \
$db_credentials_file --models=10 --features=300 --is-pathcore-example

python web_edge_page_data.py $filtered_network $db_credentials_file

wait
echo "The PAO1 eADAGE analysis has finished."
exit 0
