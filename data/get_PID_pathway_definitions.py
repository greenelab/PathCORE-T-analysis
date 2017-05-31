"""
Description:
    This script parses the file `./tcga_data/c2.cp.v6.0.symbols.gmt`
    and creates the PID pathway definitions file from it
    (`./tcga_data/PID_pathway_definitions.txt`). The output is
    tab-delimited and each row has the following format:
        pathway  N (# of genes in the definition)  gene1;gene2;...geneN

Usage:
    get_PID_pathway_definitions.py <input-file> <output-file>
    get_PID_pathway_definitions.py -h --help

Options:
    -h --help                        Show this screen.
    <input-file>                     Path to the curated canonical pathways
                                     gene sets file (`c2.cp.v6.0.symbols.gmt`)
                                     from MSigDB
    <output-file>                    Specify the output filepath
"""
from docopt import docopt


if __name__ == "__main__":
    arguments = docopt(__doc__, version="1.0")
    input_file = arguments["<input-file>"]
    output_file = arguments["<output-file>"]
    with open(output_file, "w", newline="") as outfile, \
            open(input_file, "r") as infile:
        for row in infile:
            tokens = row.split("\t")
            if "/PID" in tokens[1]:
                pathway = tokens[0]
                genes = tokens[2:]
                n_genes = len(genes)
                if n_genes:
                    genes_str = ";".join(genes)
                    output_row = [pathway, str(n_genes), genes_str]
                    outfile.write("\t".join(output_row))
