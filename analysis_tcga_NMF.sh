tcga_data_dir="./data/tcga_data"
# path to the dataset file
expression_dataset=$tcga_data_dir"/HiSeqV2"
# path to the PID pathway file
pathway_file=$tcga_data_dir"/pid_pathway_definitions.txt"

nmf_dir=$tcga_data_dir"/NMF_model"
#mkdir -p $nmf_dir

# specify the analysis output directory
analysis_dir=$tcga_data_dir"/NMF_analysis"
mkdir -p $analysis_dir

# network output directory will contain the generated network file
network_output_dir=$analysis_dir"/network_construction"
#mkdir -p $network_output_dir

# count number of genes
N_genes=$(expr $(wc -l <$normalized_expression_dataset) - 1)

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
# feature overrepresentation analysis on the model in $network_output_directory.
#python3 build_co_occurrence_network.py \
#$nmf_dir $network_output_dir $pathway_file $N_genes \
#--n-features=$N_features --signature=NMF --signature-args=$std_cutoff \
#--alpha=$alpha --n-cores=$N_cores --shorten=TCGA_PID \
#--crosstalk-correction --all-genes

# Builds a network of pathway co-occurrence relationships for
# each of the models and applies the permutation test to these.
# Results in a single aggregate network where the edges in the
# network were considered significant under the permutation test.
N_permutations=10000

permutation_output_dir=$analysis_dir"/permutation_test_N="$N_permutations

python3 run_permutation_test.py \
$network_output_dir $permutation_output_dir \
--n-permutations=$N_permutations --n-features=$N_features \
--n-cores=2
