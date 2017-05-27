pa_data_dir="./data/pao1_data"
# path to the normalized dataset file
normalized_expression_dataset=$pa_data_dir"/all-pseudomonas-gene-normalized.pcl"
# path to the KEGG pathway file
pathway_file=$pa_data_dir"/pseudomonas_KEGG_terms.txt"

nmf_dir=$pa_data_dir"/NMF_model"

# specify the analysis output directory
analysis_dir=$pa_data_dir"/NMF_analysis"
mkdir -p $analysis_dir

# network output directory will contain the generated network file
network_output_dir=$analysis_dir"/network_construction"

# count number of genes
N_genes=$(expr $(wc -l <$normalized_expression_dataset) - 1)

# path to the list of genes we have data for in the compendium
pa_genes=$pa_data_dir"/pao1_compendium_genes.txt"
tail -n+2 $normalized_expression_dataset | awk -F"\t" '{print $1}' > $pa_genes

# the model size that best supports the current P.a expression compendium k=300
N_features=300

# number of cores to use in parallel in run_pathcore_analysis.py,
# please modify this to accommodate your need.
N_cores=1

# if the gene's weight is above this number of standard deviations from the
# mean, we will include it in the gene signature of the feature
std_cutoff=2.0

# significance level used during both the pathway overrepresentation analysis
# and the edge permutation test
alpha=0.05

###############################################################################
# PathCORE analysis
###############################################################################

# Generates the significant pathways files and metadata from the
# feature overrepresentation analysis on all models in $ensemble_directory.
python3 build_co_occurrence_network.py \
$nmf_dir $network_output_dir $pathway_file $N_genes \
--n-features=$N_features --signature=NMF --signature-args=$std_cutoff \
--alpha=$alpha --n-cores=$N_cores --shorten=PAO1_KEGG \
--crosstalk-correction --all-genes

# Builds a network of pathway co-occurrence relationships for
# each of the models and applies the permutation test to these.
# Results in a single aggregate network where the edges in the
# network were considered significant under the permutation test.
N_permutations=10000

permutation_output_dir=$analysis_dir"/permutation_test_n="$N_permutations

python3 run_permutation_test.py \
$network_output_dir $permutation_output_dir \
--n-permutations=$N_permutations --n-features=$N_features \
--n-cores=$N_cores
