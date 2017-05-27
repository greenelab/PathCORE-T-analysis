"""
This script iterates through a directory of weight matrices generated
by a feature construction algorithm trained on a transcriptomic dataset.
We allow this analysis to be run on 1 or more weight matrices with
the same shape[*], in the event that multiple models of the data were
produced by the same algorithm initialized with different random seeds.

[*]shape [n,k], where n = all genes for which we have transcriptomic data,
                      k = number of features constructed.

Defaults are set up for the P. aeruginosa eADAGE analysis provided
in the PathCORE paper.

Output (all .tsv files):
    Significant pathways files produced from the pathway overrepresentation
    analysis. For each model, indicates which pathways are enriched in which
    features.

    `networks` directory:
        Networks built using each of the top-level significant pathways files

    If -m | --metadata:
      An additional `metadata` directory is created in the output directory.
      The metadata files contain more information (specific to each
      model) about
        * what the gene signature(s) were for each feature
        * the genes in the gene signature(s) that were annotated to a pathway
          after crosstalk-correction in each feature, if applicable.

Usage:
    pathcore_network_from_model_features.py
        <models-dir> <output-dir> <pathway-definitions> <n-genes>
        [--n-features=<n-features>]
        [--signature=<signature>] [--signature-args=<s-args>]
        [--alpha=<alpha>]
        [--genes-list=<genes-list>]
        [--n-cores=<n-cores>]
        [--replace=<replace>] [--shorten=<pathway-shorten>]
        [-o | --overlap-correction] [-a | --all-genes]
        [-m | --metadata]
    pathcore_network_from_model_features.py -h | --help

Options:
    -h --help                   Show this screen.

    <models-dir>                Path to the directory containing the models
                                generated by an unsupervised feature
                                construction algorithm.

    <output-dir>                Path to the directory that will store the
                                output files. Will be automatically created
                                if the directory does not currently exist.

    <pathway-definitions>       Path to the pathway definitions file.
                                Formatted as tab-delimited columns:
                                pathway, N (num. genes), geneA;geneB;...geneN

    <n-genes>                   Number of genes for which we have recorded
                                expression values in the original dataset

    --n-features=<n-features>   Number of constructed features in each model
                                (currently, they must all contain the same
                                number of features)
                                [default: 300]

    --signature=<signature>     Specify the gene signature definition to use
                                in this analysis. Please review
                                `constants/gene_signature_definitions.py`
                                to select or modify the file to include
                                a partial function that determines the gene
                                signature for a given feature construction
                                algorithm.
                                [default: eADAGE]

    --signature-args=<s-args>   Depends on the gene signature definition
                                selected <signature>. Provide a
                                semicolon-separated list of input arguments
                                needed to use partial function.
                                In the eADAGE signature definition, we only
                                require a user to specify standard deviations
                                from the mean---hence the default.
                                [default: 2.5]

    --alpha=<alpha>             Significance level for pathway enrichment.
                                Overrepresentation is determined by a Fisher's
                                exact test with false-discovery rate
                                correction.
                                [default: 0.05]

    --genes-list=<genes-list>   If the weight matrix files do not include the
                                gene identifiers corresponding to each row,
                                please include the path to a file that contains
                                these identifiers (1 gene id per line).

    --n-cores=<n-cores>         Number of cores used to run the analysis
                                on models in parallel
                                [default: num. available cores - 1]

    --replace=<replace>         The resulting file keeps the naming convention
                                of the input file, with the exception of this
                                substring, which will be replaced by
                                'SigPathways'
                                [default: network]

    --shorten=<pathway-shorten> Specify a function to shorten pathway
                                names for the network. Useful for
                                readability in network visualizations.
                                Please review
                                `constants/shorten_pathway_names.py`
                                and select or modify the file to include
                                a function that takes in a pathway label
                                and outputs a shortened version of it.

    -o --overlap-correction     Donato et al.'s (2013) crosstalk correction
                                algorithm is applied to the input pathway
                                definitions.
                                [default: True]

    -a --all-genes              If `--overlap-correction`, this flag
                                applies the correction to the full
                                gene set annotated to each pathway, as opposed
                                to only genes that are also in a feature's gene
                                signature.
                                [default: True]

    -m --metadata               Collect metadata after overrepresentation
                                analysis of each constructed feature.
                                Used for setting up the PathCORE Flask
                                application.
                                [default: False]
"""
import csv
import multiprocessing
import os
import sys
from time import time

from docopt import docopt
from joblib import Parallel, delayed
import pandas as pd
from pathcore import pathway_enrichment_with_overlap_correction, \
    pathway_enrichment_no_overlap_correction
from pathcore import CoNetwork

from ..constants import GENE_SIGNATURE_DEFINITIONS
from ..constants import SHORTEN_PATHWAY_NAMES
from utils import load_pathway_definitions, load_weight_matrix

# include in the significant pathways output filename
OUTPUT_FILE_NAMING = "SigPathway"
# include in the metadata filenames
METADATA_FILE_ENDS_WITH = "pathcore_overrepresentation_analysis.tsv"


def process(current_weight_matrix, pathway_definitions,
            partial_function_signature,
            alpha=0.05, overlap_correction=True, correct_all_genes=True,
            metadata=False):
    """
    Determines the pathway coverage of a weight matrix

    Parameters
    -----------
    current_weight_matrix : (str, pandas.DataFrame)
        (weight matrix filename, weight matrix dataframe of size [n,k])
    pathway_definitions : dict(str -> set(str))
        Pathway definitions, pre-overlap-correction.
        A pathway (key) is defined by a set of genes (value).
    partial_function_signature : functools.partial
        Accepts a feature weight vector of type pandas.Series(float)
        and returns (set(), set()), the feature's positive and/or negative
        gene signatures. See `gene_signature_definitions.py` in
        `..constants` for examples.
    alpha : float, default=0.05
        Significance level for pathway enrichment.
    overlap_correction : bool (default=True)
        Apply Donato et al.'s (2013) maximum impact estimation algorithm
        on the pathway definitions to remove gene overlap.
    correct_all_genes : bool (default=True)
        For `overlap_correction=True`, the correction is always
        applied to the gene signature(s).
        When `correct_all_genes=True`, the correction is applied to genes
        outside of the signature(s)--termed *remaining* genes.
    metadata : bool (default=False)
        When True, outputs additional files about the per-feature
        overrepresentation analysis, gene signatures, and pathway annotations
        into a `metadata` directory.

    Returns
    -----------
    dict() with the following items:
      - "filename": str filename
      - "significant_pathways": pandas.DataFrame significant
                                pathways dataframe
      - "feature_metadata": dict(int -> dict) if applicable.
                            * keys are the feature numbers
                            * values are the metadata dicts
    """
    filename, weight_matrix = current_weight_matrix
    n_genes, _ = weight_matrix.shape
    print("Processing model {0}".format(filename))
    significant_pathways_df = pd.DataFrame()
    feature_metadata = {}
    for feature in weight_matrix:
        if crosstalk_correction:
            feature_df, additional = \
                pathway_enrichment_with_overlap_correction(
                    weight_matrix[feature], pathway_definitions,
                    partial_function_signature,
                    alpha=alpha, correct_all_genes=correct_all_genes,
                    metadata=metadata)
        else:
            feature_df, additional = \
                pathway_enrichment_no_overlap_correction(
                    weight_matrix[feature], pathway_definitions,
                    partial_function_signature,
                    alpha=alpha, metadata=metadata)
        if feature_df is not None:
            feature_df.loc[:, "feature"] = pd.Series(
                [feature] * len(feature_df.index), index=feature_df.index)
            significant_pathways_df = pd.concat(
                [significant_pathways_df, feature_df], axis=0)
            if additional:
                feature_metadata[feature] = additional
    significant_pathways_df.reset_index(drop=True, inplace=True)
    return {"filename": filename,
            "significant_pathways": significant_pathways_df,
            "feature_metadata": feature_metadata}


def _write_metadata_tsv(partial_filename,
                        rows,
                        model_number,
                        output_directory):
    """A small helper function for writing the metadata information
    to a .tsv file
    """ 
    filepath = os.path.join(
        output_directory,
        "{0}_{1}_{2}".format(
            model_number, partial_filename, METADATA_FILE_ENDS_WITH))
    with open(filepath, "w", newline="") as fp:
        writer = csv.writer(fp, delimiter="\t")
        for row in rows:
            writer.writerow(row)


if __name__ == "__main__":
    arguments = docopt(
        __doc__, version="build co-occurrence network 1.0")
    models_directory = arguments["<models-dir>"]
    output_directory = arguments["<output-dir>"]
    pathway_definitions_file = arguments["<pathway-definitions>"]

    n_genes = int(arguments["<n-genes>"])
    n_features = int(arguments["--n-features"])

    if arguments["--signature"] not in GENE_SIGNATURE_DEFINITIONS:
        raise ValueError("The gene signature definition specified ({0}) "
                         "is not in the `GENE_SIGNATURE_DEFINITIONS` "
                         "dict. Please include your definition in "
                         "`gene_signature_definitions.py` before running "
                         "this script.".format(arguments["--signature"])) 

    which_signature = GENE_SIGNATURE_DEFINITIONS[arguments["--signature"]]
    gene_signature_args = arguments["--signature-args"].split(";")
    # NOTE: no type casting done, all arguments will be strings.
    partial_function_signature = which_signature(*gene_signature_args)

    alpha = float(arguments["--alpha"])
    path_to_genes_file = arguments["--genes-list"]

    substring_to_replace = arguments["--replace"]

    crosstalk_correction = arguments["--crosstalk-correction"]
    correct_all_genes = arguments["--all-genes"]

    network_output_directory = os.path.join(
        output_directory, "networks")
    os.makedirs(output_directory, exist_ok=True)
    os.makedirs(network_output_directory, exist_ok=True)

    metadata = arguments["--metadata"]
    metadata_output_directory = None
    if metadata:
        metadata_output_directory = os.path.join(
            output_directory, "metadata")
        os.makedirs(metadata_output_directory, exist_ok=True)

    if arguments["--n-cores"].isdigit():
        n_cores = int(arguments["--n-cores"])
    else:
        n_cores = multiprocessing.cpu_count() - 1

    t_o = time()

    shorten_function_key = arguments["--shorten"]
    shorten_pathway_names = None
    if shorten_function_key and shorten_function_key in SHORTEN_PATHWAY_NAMES:
        shorten_pathway_names = SHORTEN_PATHWAY_NAMES[shorten_function_key]
    elif shorten_function_key not in SHORTEN_PATHWAY_NAMES:
        raise ValueError("The shorten pathway function specified ({0}) "
                         "is not in the `SHORTEN_PATHWAY_NAMES` "
                         "dict. Please include your definition in "
                         "`shorten_pathway_names.py` before running this"
                         "script.".format(arguments["shorten_function_key"]))
    pathway_definitions = load_pathway_definitions(
        pathway_definitions_file,
        shorten_pathway_names=shorten_pathway_names)

    weight_matrices = []
    for model in os.listdir(models_directory):
        path_to_model_file = os.path.join(models_directory, model)
        weight_matrix = load_weight_matrix(
            path_to_model_file, n_features,
            n_genes=n_genes, path_to_genes_file=path_to_genes_file)
        weight_matrices.append((model, weight_matrix))

    process_model_args = [
        pathway_definitions, partial_function_signature,
        alpha, crosstalk_correction, correct_all_genes, metadata
    ]

    # Please note that we are not accounting for cases where `results`
    # may be too large to fit into memory. The information stored in
    # `results` is not written to any files until after all weight matrices
    # have been processed. (see: line 313)
    with Parallel(n_jobs=n_cores) as parallel:
        results = parallel(
            delayed(process)(weight_matrix_info, *process_model_args)
            for weight_matrix_info in weight_matrices)
    t_f = time() - t_o

    print("{0} models took {1} seconds to run on {2} cores.".format(
        len(results), t_f, n_cores))

    for model_number, model_results in enumerate(results):
        model_filename = model_results["filename"]
        significant_pathways_df = model_results["significant_pathways"]

        # write to significant pathways file
        significant_pathways_out = model_filename.replace(
            substring_to_replace, OUTPUT_FILE_NAMING)
        significant_pathways_out_path = os.path.join(
            output_directory, "{0}_{1}".format(
                model_number, significant_pathways_out))
        significant_pathways_df.to_csv(
            path_or_buf=significant_pathways_out_path, sep="\t", index=False)

        # write to the network file
        model_network = CoNetwork(
            n_features, significant_pathways=significant_pathways_df)
        model_network_df = model_network.to_dataframe()
        network_out = "{0}_pathcore_network.tsv".format(model_number)
        network_out_path = os.path.join(network_output_directory, network_out)
        model_network_df.to_csv(network_out_path, sep="\t", index=False)

        # write to metadata files
        if model_results["feature_metadata"]:
            # both lists will store the rows that are written to their
            # corresponding files.
            # initialize both with the column names as the first row.
            gene_signatures = [
                ("feature", "positive_signature", "negative_signature")]
            signature_pathway_definitions = [
                ("feature", "pathway", "gene_signature_definition")]

            for feature_number, metadata_info in \
                    model_results["feature_metadata"].items():
                gene_signatures.append(
                    (feature_number,
                     ";".join(metadata_info["positive_signature"]),
                     ";".join(metadata_info["negative_signature"])))
                for pathway, definition in \
                        metadata_info["pathway_definitions"].items():
                    definition_genes = ";".join(definition)
                    if definition_genes:
                        signature_pathway_definitions.append(
                             (feature_number, pathway, definition_genes))
            _write_metadata_tsv(
                "FEATURE_SIGNATURES", gene_signatures,
                model_number, metadata_output_directory)
            _write_metadata_tsv(
                "FEATURE_PATHWAYS", signature_pathway_definitions,
                model_number, metadata_output_directory)