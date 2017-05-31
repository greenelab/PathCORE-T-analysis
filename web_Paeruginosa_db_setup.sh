#!/bin/bash                                                                                                                                      
set -o errexit

# Kathleen Chen 2017

###############################################################################
# PAO1 <> eADAGE <> KEGG web application DB setup
###############################################################################

pa_data_dir="./data/pao1_data"

filtered_network=$permutation_output_dir"/filtered_network.tsv"

# NOTE: Please create this file after you have set up an mLab account and
# a MongoDB instance with a database user (see README)
db_credentials_file="./db_credentials.yml"

python web_initialize_db.py $data_compendium $pathway_file $filtered_network \
                            $network_output_dir $db_credentials_file \
                            --features=300 --is-pathcore-example

python web_edge_page_data.py $filtered_network $db_credentials_file
