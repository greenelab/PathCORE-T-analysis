"""This script updates an mLab-hosted MongoDB database with the data needed
to compute information for the web application.

Usage:
    web_initialize_db.py <dataset> <pathway-definitions> <network>
        <enrichment-analysis-dir> <db-credentials>
        [--models=<n-models>]
        [--features=<n-features>]
        [--is-pathcore-example]
    web_initialize_db.py -h | --help

Options:
    -h --help                   Show this screen.

    <dataset>                   Path to the gene expression dataset file.

    <pathway-definitions>       Path to the pathway definitions file.
                                formatted as tab-delimited columns:
                                pathway, N (num. genes), geneA;geneB;...geneN

    <network>                   Path to the network file that will be
                                visualized on the web application.

    <enrichment-analysis-dir>   Path to the directory of enrichment analysis
                                files. We will use the `networks` and
                                `metadata` directories created from running
                                `run_network_creation` with the
                                metadata flag (`-m` or `--metadata`).

    <db-credentials>            Path to the .yml file containing the necessary
                                DB credentials to access an mLab-served MongoDB
                                database.

    --features=<n-features>     Number of features constructed in each model.

    --is-pathcore-example       Sample annotations and common names are
                                available for the P. aeruginosa (PAO1) KEGG
                                case study described in the PathCORE paper.
                                This flag runs the functions that incorporate
                                this information into the web interface for the
                                PAO1 case.
                                [default: False]
"""
import os
import sys

from docopt import docopt
import pandas as pd
from pathcore import CoNetwork
from pymongo import ASCENDING  # for indexing the Collection
from pymongo import MongoClient
import yaml

from constants import SHORTEN_PATHWAY_NAMES
from utils import load_expression_dataset, load_pathway_definitions
import utils_setup_PAO1_example as utils_PAO1


def convert_to_gene_object_ids(gene_set, gene_object_id_map):
    """Given a gene set, convert the gene identifiers to their
    corresponding MongoDB ObjectIds.

    Parameters
    -----------
    gene_set : str|set(str)
        Semicolon-separated gene identifiers from a metadata file
        or the gene set directly (as a set of strings). If the input is not
        a str (e.g. it is NaN in a pandas.DataFrame) or set(str),
        this function will return an empty list.
    gene_object_id_map : dict(str -> ObjectId)
        Mapping from the gene identifier to its ObjectId

    Returns
    -----------
    list(ObjectId), the list of gene ObjectIds
    """
    if isinstance(gene_set, str):
        gene_set = set(gene_set.split(";"))
    if not isinstance(gene_set, set):
        return []
    ids_list = []
    for gene in gene_set:
        if gene in gene_object_id_map:
            ids_list.append(gene_object_id_map[gene])
    return ids_list


def _get_metadata_row_records(partial_filename,
                              network_number,
                              metadata_directory):
    """Helper function to read metadata files into a list of dictionaries,
    where each dictionary has column name keys and row values
    """
    filename = "{0}_{1}_pathcore_overrepresentation_analysis.tsv".format(
        network_number, partial_filename)
    path_to_file = os.path.join(metadata_directory, filename)
    row_records = pd.read_table(path_to_file).to_dict("records")
    return row_records


if __name__ == "__main__":
    arguments = docopt(
        __doc__, version="web application DB setup, 1.0")
    expression_data_file = arguments["<dataset>"]
    pathway_definitions_file = arguments["<pathway-definitions>"]
    network_file = arguments["<network>"]
    enrichment_analysis_dir = arguments["<enrichment-analysis-dir>"]
    db_credentials_file = arguments["<db-credentials>"]
    n_features = int(arguments["--features"])

    is_example = arguments["--is-pathcore-example"]

    credentials_info_dict = {}
    with open(db_credentials_file, "r") as credentials_info:
        try:
            credentials_info_dict = yaml.load(credentials_info)
        except yaml.YAMLError as exception:
            sys.exit(exception)

    mongo_client = MongoClient("mongodb://{0}:{1}@{2}/{3}".format(
        credentials_info_dict["user"], credentials_info_dict["password"],
        credentials_info_dict["mLab_URI"], credentials_info_dict["db_name"]))
    db = mongo_client[credentials_info_dict["db_name"]]

    # clearing out all of the collections if any of them already exist
    db.genes.drop()
    db.pathways.drop()
    db.sample_labels.drop()
    db.network_edges.drop()
    db.network_feature_signatures.drop()
    db.network_feature_pathways.drop()

    # after drop, db.<collection> will initialize a new collection
    genes = db.genes
    pathways = db.pathways
    sample_labels = db.sample_labels
    network_edges = db.network_edges
    network_feature_signatures = db.network_feature_signatures
    network_feature_pathways = db.network_feature_pathways

    # get row/colnames & expression values from the gene expression dataset
    expression_data = load_expression_dataset(expression_data_file)
    sample_label_list = expression_data.columns
    gene_list = expression_data.index

    # insert the name of each sample & set the _id as the column index
    index_and_sample_list = []
    for index, sample_label in enumerate(sample_label_list):
        index_and_sample_list.append({"_id": index, "sample": sample_label})
    sample_labels.insert_many(index_and_sample_list, ordered=True)

    # insert the gene's expression value information across those samples
    gene_expression = [list(v) for v in expression_data.values]  # transpose
    gene_and_expression_values_list = []
    for index, gene in enumerate(gene_list):
        gene_expression_values = {"gene": gene,
                                  "expression": gene_expression[index]}
        gene_and_expression_values_list.append(gene_expression_values)
    gene_object_ids = genes.insert_many(
        gene_and_expression_values_list, ordered=True)
    genes.create_index([("gene", ASCENDING)])

    # mapping the gene name to the gene object ID
    gene_to_object_id = {}
    for index in range(len(gene_and_expression_values_list)):
        gene_to_object_id[gene_list[index]] = (
            gene_object_ids.inserted_ids[index])

    # [pa example] update each gene with its common name, if applicable.
    if is_example:
        utils_PAO1.get_gene_common_names(
            genes,
            "./data/pao1_web_info/Pseudomonas_aeruginosa_PAO1.gene_info")

    # insert the pathway annotations
    pathway_definitions_map = load_pathway_definitions(
        pathway_definitions_file)
    pathway_info_list = []
    for pathway, gene_set in pathway_definitions_map.items():
        if is_example:  # [pa example]
            pathway = SHORTEN_PATHWAY_NAMES["PAO1_KEGG"](pathway)
        gene_set_ids = convert_to_gene_object_ids(gene_set, gene_to_object_id)
        pathway_info_list.append(
            {"pathway": pathway, "annotations": gene_set_ids})
    pathway_object_ids = pathways.insert_many(pathway_info_list, ordered=True)

    # mapping the pathway name to the pathway object ID
    pathway_to_object_id = {}
    for index in range(len(pathway_info_list)):
        pathway_to_object_id[pathway_info_list[index]["pathway"]] = (
            pathway_object_ids.inserted_ids[index])

    # include information collected from the overrepresentation analysis
    network_obj = CoNetwork(n_features)
    network_obj.read_network_file(network_file)

    edge_pathways_set = set()
    for (v0, v1), edge in network_obj.edges.items():
        pw0 = network_obj.__getitem__(v0)
        pw1 = network_obj.__getitem__(v1)
        edge_pathways_set.add((pw0, pw1))

    networks_directory = os.path.join(enrichment_analysis_dir, "networks")
    metadata_directory = os.path.join(enrichment_analysis_dir, "metadata")
    for network_num, network_filename in enumerate(
            os.listdir(networks_directory)):
        pathcore_network = os.path.join(networks_directory, network_filename)

        current_network_edges = pd.read_table(
            pathcore_network).to_dict("records")

        metadata_feature_signatures = _get_metadata_row_records(
            "FEATURE_SIGNATURES", network_num, metadata_directory)
        metadata_feature_pathways = _get_metadata_row_records(
            "FEATURE_PATHWAYS", network_num, metadata_directory)

        keep_network_edges = []
        for edge_info in current_network_edges:
            pw0 = edge_info["pw0"]
            pw1 = edge_info["pw1"]
            if (pw0, pw1) in edge_pathways_set:
                edge_info.pop("pw0")
                edge_info.pop("pw1")
                edge_info["edge"] = (
                    pathway_to_object_id[pw0], pathway_to_object_id[pw1])
                edge_info["features"] = [int(float(x)) for x in
                                         edge_info["features"].split(" ")]
                edge_info["network"] = network_num
                keep_network_edges.append(edge_info)

        network_edges.insert_many(keep_network_edges, ordered=True)

        for signature_info in metadata_feature_signatures:
            signature_info["network"] = network_num
            if "positive_signature" in signature_info:
                signature_info["positive_signature"] = \
                    convert_to_gene_object_ids(
                        signature_info["positive_signature"],
                        gene_to_object_id)
            if "negative_signature" in signature_info:
                signature_info["negative_signature"] = \
                    convert_to_gene_object_ids(
                        signature_info["negative_signature"],
                        gene_to_object_id)
            signature_info["feature"] = int(signature_info["feature"])
        network_feature_signatures.insert_many(
            metadata_feature_signatures, ordered=True)

        for pathway_info in metadata_feature_pathways:
            pathway_info["network"] = network_num
            pathway_info["pathway"] = pathway_to_object_id[pathway]
            pathway_info["gene_signature_definition"] = \
                convert_to_gene_object_ids(
                    pathway_info["gene_signature_definition"],
                    gene_to_object_id)
            pathway_info["feature"] = int(pathway_info["feature"])
        network_feature_pathways.insert_many(
            metadata_feature_pathways, ordered=True)

    if is_example:
        db.sample_annotations.drop()
        sample_annotations = db.sample_annotations
        utils_PAO1.get_sample_annotations(
            sample_annotations,
            sample_labels,
            "./data/pao1_web_info/PseudomonasAnnotation.tsv")

    mongo_client.close()
