"""
Description:
    Generates an NMF model of an input expression dataset. The features
(referred to in scikit-learn as components and in various publications as
patterns) are constructed after taking the transpose of a size [n, m] dataset,
where we assume the rows of the original input are the genes (n) and the
columns are samples (m). This results in a model of size [n, k] that assigns
n genes different weights in each of the k features, k <= m. The PathCORE-T
analysis will use this as the input weight matrix.

Usage:
    data_generate_NMF_model.py
        <data-file> <model-out-file> <n-features>
        [--random-seed=<random-seed>]
    data_generate_NMF_model.py -h --help

Options:
    -h --help                        Show this screen.
    <data-file>                      The path to the expression dataset file
    <model-out-file>                 Specify the output filepath
    <n-features>                     Number of features (NMF components) in
                                     the resulting model
    --random-seed=<random-seed>      Specify the random seed used during
                                     initialization
                                     [default: 42]
"""
from docopt import docopt
import pandas as pd
from sklearn.decomposition import NMF

from utils import load_expression_dataset


def generate_nmf_weight_matrix(dataset,
                               n_features,
                               random_state=42,
                               save_to_file=None):
        """Uses scikit-learn's NMF approximation method
        to produce an NMF model of the input transcriptomic dataset.

        Parameters
        -----------
        dataset : pandas.DataFrame, shape = [n, m]
            Input transcriptomic dataset. Expects a tab-delimited file with
            row and column names, where the rows are the genes and the
            columns are the samples.
        n_features : int, k
            The number of features (NMF components), for dimensionality
            reduction.
        random_state : int
            Set the random state, for random number generator seed control.
        save_to_file : str
            Specify an output filepath to save the resulting NMF components
            as a tab-delimited file.

        Returns
        -----------
        pandas.Dataframe, shape = [n, k]
            where the columns are the constructed features
        """
        model = NMF(n_components=n_features,
                    init="nndsvda",  # for dense datasets
                    solver="cd",
                    random_state=random_state)
        model.fit_transform(dataset.T)

        gene_weights_across_components = {}
        for index, gene_weights in enumerate(model.components_.T):
            gene = dataset.index[index]
            gene_weights_across_components[gene] = list(gene_weights)

        weight_matrix = pd.DataFrame.from_dict(gene_weights_across_components,
                                               orient="index")
        if save_to_file:
            weight_matrix.to_csv(
                save_to_file, sep="\t", header=False, index=True)
        return weight_matrix


if __name__ == "__main__":
    arguments = docopt(
        __doc__, version="generate NMF model 1.0")
    expression_data = arguments["<data-file>"]
    out_filepath = arguments["<model-out-file>"]
    n_features = int(arguments["<n-features>"])

    expression_data_df = load_expression_dataset(expression_data)
    additional_args = {}
    additional_args["save_to_file"] = out_filepath
    if arguments["--random-seed"]:
        additional_args["random_state"] = int(arguments["--random-seed"])
    generate_nmf_weight_matrix(
        expression_data_df, n_features, **additional_args)
