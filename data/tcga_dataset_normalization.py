"""
Usage:
    tcga_data_normalization.py <input-file> <output-file>

Options:
    -h --help                        Show this screen.
    <input-file>                     Path to the TCGA Pan-Cancer dataset,
                                     filename: HiSeqV2
    <output-file>                    Specify the output filepath
"""
from docopt import docopt

from ..utils import expression_data_minmax_normalization

if __name__ == "__main__":
    arguments = docopt(__doc__, version="1.0")
    input_file = arguments["<input-file>"]
    output_file = arguments["<output-file>"]
    output_df = expression_data_minmax_normalization(
        input_file, "Sample")
    output_df.to_csv(output_file, sep="\t", header=True, index=True)
