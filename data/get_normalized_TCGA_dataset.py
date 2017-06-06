"""
Description:
    This script is used to zero-one normalize the dataset
    `./tcga_data/HiSeqV2`

Usage:
    get_normalized_TCGA_dataset.py <input-file> <output-file>
    get_normalized_TCGA_dataset.py -h --help

Options:
    -h --help                        Show this screen.
    <input-file>                     Path to the TCGA Pan-Cancer dataset,
                                     filename: HiSeqV2
    <output-file>                    Specify the output filepath
"""
from docopt import docopt
import pandas as pd
from sklearn.preprocessing import MinMaxScaler


def expression_data_minmax_normalization(path_to_file, index_on_col):
    """Used to normalize the TCGA expression dataset.

    Parameters
    -----------
    path_to_file : str
        TThe path to the expression dataset.
        Currently expects a tab-delimited file with row and column names,
        where the rows are the genes and the columns are the samples.
    index_on_col : str
        Name of the column that should be the dataframe's index

    Returns
    -----------
    pandas.DataFrame
    """
    data = pd.read_table(path_to_file)
    data.set_index(index_on_col, inplace=True)
    data = data[-data.index.str.contains('?', regex=False)]
    data = data.sort_index()

    data_normalized = MinMaxScaler().fit_transform(data.T)
    data_normalized = pd.DataFrame(
        data_normalized, index=data.columns, columns=data.index).T
    return data_normalized


if __name__ == "__main__":
    arguments = docopt(__doc__, version="1.0")
    input_file = arguments["<input-file>"]
    output_file = arguments["<output-file>"]
    output_df = expression_data_minmax_normalization(
        input_file, "Sample")
    output_df.to_csv(output_file, sep="\t", header=True, index=True)
