#!/bin/bash       
set -o errexit

# Kathleen Chen 2017

N_features=300

pa_data_dir="./pao1_data"
compendium_file=$pa_data_dir"/all-pseudomonas-gene-normalized.pcl"

pa_nmf_dir = $pa_data_dir"/NMF_model"
pa_nmf_file = $pa_nmf_dir"/NMF_"$N_features"_features.tsv"

# Generates the NMF k=300 model from the normalized P. aeruginosa
# gene compendium
python generate_nmf_model.py \
$normalized_pancan_file $tcga_nmf_file $N_features

tcga_data_dir="./tcga_data"
pancan_file=$tcga_data_dir"/HiSeqV2"

normalized_pancan_file=$tcga_data_dir"/HiSeqV2_minmaxscale_normalized"
python tcga_dataset_normalization.py \
$pancan_file $normalized_pancan_file

tcga_nmf_dir = $tcga_data_dir"/NMF_model"
tcga_nmf_file=$tcga_nmf_dir"/NMF_"$N_features"_features.tsv"

# Generates the NMF k=300 model from the TCGA Pan-Cancer expression dataset
python generate_nmf_model.py \
$normalized_pancan_file $tcga_nmf_file $N_features



