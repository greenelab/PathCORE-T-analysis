"""Should be run *after* `web_initialize_db.py`.

For every edge, create a page reporting on its underlying genes.
This script computes and inserts the page information into a MongoDB
collection called `pathcore_edge_data`.

Defaults are the ones used for the PAO1 KEGG eADAGE network demo server.

Usage:
    web_edge_page_data.py <network> <db-credentials>
        [--features=<n-features>]
        [--samples=<keep-n-samples>] [--genes=<keep-n-genes>]
    web_edge_page_data.py -h | --help

Options:
    -h --help                   Show this screen.

    <network>                   Path to the network file that will be
                                visualized on the web application. Expects
                                that this is the network output from the
                                permutation test in the full PathCORE analysis.

    <db-credentials>            Path to the .yml file containing the necessary
                                DB credentials to access an mLab-served MongoDB
                                database.

    --features=<n-features>     Number of features constructed in each model.
                                [default: 300]

    --samples=<keep-n-samples>  The number of most/least expressed samples
                                (each) shown on an edge page.
                                [default: 20]

    --genes=<keep-n-genes>      The number of genes shown on an edge page.
                                [default: 20]

"""
import sys

from docopt import docopt
from pymongo import ASCENDING  # for indexing the Collection
from pymongo import MongoClient
import pandas as pd
import yaml


def index_sort(list_to_sort):
    """
    Parameters
    -----------
    list_to_sort : list(Any)
      List of elements, any type

    Returns
    -----------
    list(tup(int, Any)), a list of tuples, where
      tup[0] is the position of the element in the original list and
      tup[1] is the element itself.
      The tuples in the list are ordered in descending order
      based on tup[1] element comparisons.
    """
    return [index_element_tup for index_element_tup
            in sorted(enumerate(list_to_sort), key=lambda x: -x[1])]


class EdgeInfo:
    """Functions in this class are used to generate the JSON object
    corresponding to each edge page"""

    def __init__(self, db, n_features, keep_n_samples=20, keep_n_genes=20):
        """
        Parameters
        -----------
        db : pymongo.database.Database
          The database retrieved using MongoClient
        n_features : int
          The number of features from which each network was constructed
          (# features constructed in a model)
        keep_n_samples : int (default=20)
          The number of samples to display in each heatmap. Limited so
          the heatmap fits on a web page and can be manually reviewed.
        keep_n_genes : int (default=20)
          The number of genes to display in each heatmap. Limited so
          the heatmap fits on a web page and can be manually reviewed.
        """
        self.db = db
        self.n_features = n_features
        # the final network may be an aggregate of multiple networks
        self.n_networks = len(self.db.network_edges.distinct("network"))
        self.n_total_features = self.n_features * self.n_networks

        self.keep_n_samples = keep_n_samples
        self.keep_n_genes = keep_n_genes

    def get_edge_info(self, pw0_name, pw1_name):
        """
        Parameters
        -----------
        pw0_name : str
          Pathway 0 in the edge
        pw1_name : str
          Pathway 1 in the edge

        Returns
        -----------
        dict(), the dict that will be inserted into the MongoDB collection
        `pathcore_edge_data` as a JSON object. Contains the edge page
        information that will be parsed and displayed on the web application.
        """
        pw0_info = self.db.pathways.find_one({"pathway": pw0_name})
        pw1_info = self.db.pathways.find_one({"pathway": pw1_name})
        pw0_id = pw0_info["_id"]
        pw1_id = pw1_info["_id"]

        pw0_annotations = set(pw0_info["annotations"])
        pw1_annotations = set(pw1_info["annotations"])
        genes_annotated_to_edge = pw0_annotations | pw1_annotations

        # specify the pathway(s) to which a gene is annotated
        genes_which_pathway = {}
        for gene in genes_annotated_to_edge:
            if gene in pw0_annotations and gene not in pw1_annotations:
                genes_which_pathway[gene] = 0
            elif gene in pw1_annotations and gene not in pw0_annotations:
                genes_which_pathway[gene] = 1
            else:  # annotated to both
                genes_which_pathway[gene] = 2

        # track how many features contain both this gene and this edge
        count_gene_edge_common = {}
        for gene in genes_annotated_to_edge:
            count_gene_edge_common[gene] = 0

        edge_in_networks = self.db.network_edges.find(
            {"edge": [pw0_id, pw1_id]})
        n_features_contain_edge = 0
        for edge in edge_in_networks:
            network = edge["network"]
            features = edge["features"]
            n_features_contain_edge += len(features)
            # look at the gene signatures of all features that contain
            # this edge
            all_feature_infos = self.db.network_feature_signatures.find(
                {"network": network, "feature": {"$in": features}})
            for feature_info in all_feature_infos:
                signature = set()
                if "positive_signature" in feature_info:
                    signature |= set(feature_info["positive_signature"])
                if "negative_signature" in feature_info:
                    signature |= set(feature_info["negative_signature"])
                for gene in genes_annotated_to_edge & signature:
                    count_gene_edge_common[gene] += 1
        odds_ratio_information = self._edge_pathway_annotation_odds_ratios(
            genes_which_pathway, n_features_contain_edge,
            count_gene_edge_common)
        return self._edge_page_data(odds_ratio_information)

    def _edge_pathway_annotation_odds_ratios(self,
                                             which_pathway,
                                             n_features_contain_edge,
                                             count_gene_edge_common):
        """Helper function computes the odds ratio of every gene annotated to
        the edge (that is, annotated to one or both of the pathways that make
        up the edge), and returns the odds ratio information for up to
        `self.keep_n_genes` genes with an odds ratio above 1.
        """
        odds_ratios = {}
        odds_ratio_denominator = (n_features_contain_edge /
                                  float(self.n_total_features))
        for gene_id, count in count_gene_edge_common.items():
            # number of times this gene appears in any feature's gene signature
            # (not limited to the features that contain the edge)
            n_features_contain_gene = self.db.network_feature_signatures.find(
                {"$or": [{"positive_signature": gene_id},
                         {"negative_signature": gene_id}]}).count()
            # edge case: divide by 0 occurs if the annotated gene (from the
            # pathway definitions file) was not in the transcriptomic dataset
            if not count and not n_features_contain_gene:
                continue
            # `count` (number of times the edge AND gene appear) is normalized
            # by the number of times the gene appears in any feature
            odds_ratio_numerator = count / float(n_features_contain_gene)
            odds_ratios[gene_id] = (odds_ratio_numerator /
                                    odds_ratio_denominator)

        odds_ratios_above_one = []
        for gene_id, odds_ratio in odds_ratios.items():
            if odds_ratio > 1.:
                gene = self.db.genes.find_one({"_id": gene_id})
                gene_name = gene["gene"]
                if "common_name" in gene:
                    gene_name = gene["common_name"]
                odds_ratios_above_one.append(
                    (str(gene_name), odds_ratio, which_pathway[gene_id]))

        # descending by odds ratio
        odds_ratios_above_one.sort(key=lambda x: -x[1])
        limit_genes_returned = min(
            len(odds_ratios_above_one), self.keep_n_genes)
        return odds_ratios_above_one[:limit_genes_returned]

    def _edge_page_data(self, gene_odds_ratio_list):
        """Helper function creates the dict (MongoDB JSON object) that contains
        all information needed for the PathCORE web application edge page.
        """
        if not gene_odds_ratio_list:  # no genes to display for this edge
            return {"flag": 1}
        norm_sample_matrix = []
        sum_odds_ratios = 0.0

        # `norm_sample_matrix` shape = [n_genes, n_samples]
        for gene_name, odds_ratio, belongs_to in gene_odds_ratio_list:
            norm_sample_matrix.append(
                self.db.genes.find_one(
                    {"$or": [{"gene": gene_name},
                             {"common_name": gene_name}]})["expression"])
            sum_odds_ratios += odds_ratio

        # transpose => shape [n_samples, n_genes]
        norm_sample_matrix = [list(l) for l in zip(*norm_sample_matrix)]

        # for each sample, compute a summary "sample score"
        sample_scores = []
        for sample_data in norm_sample_matrix:
            score = 0.
            for index, gene_expression in enumerate(sample_data):
                _, gene_odds_ratio, _ = gene_odds_ratio_list[index]
                weight = gene_odds_ratio / sum_odds_ratios
                score += weight * gene_expression
            sample_scores.append(score)
        index_sorted_scores = index_sort(sample_scores)

        most_expressed_heatmap, most_expressed_samples = \
            self._get_heatmap_data(
                index_sorted_scores[:self.keep_n_samples], norm_sample_matrix)
        least_expressed_heatmap, least_expressed_samples = \
            self._get_heatmap_data(
                index_sorted_scores[-self.keep_n_samples:], norm_sample_matrix)
        # unzips a list of 3-element tuples
        gene_names, odds_ratios, pw_owner = [list(l) for l in
                                             zip(*gene_odds_ratio_list)]

        return {"most_expressed_heatmap": most_expressed_heatmap,
                "least_expressed_heatmap": least_expressed_heatmap,
                "most_expressed_samples": most_expressed_samples,
                "least_expressed_samples": least_expressed_samples,
                "gene_names": gene_names,
                "odds_ratios": odds_ratios,
                "pathway_owner": pw_owner}

    def _get_heatmap_data(self, index_sample_scores, expression_matrix):
        """Helper function returns 2 lists: one contains the
        positions & expression values of every square in the heatmap and
        the other contains the ordering of the samples (columns) of the heatmap
        """
        heatmap_data = []
        ordered_sample_labels = []
        for rank, (index, score) in enumerate(index_sample_scores):
            sample_label = self.db.sample_labels.find_one(
                {"_id": index})["sample"]
            ordered_sample_labels.append(str(sample_label))
            for gene_index, gene_expression \
                    in enumerate(expression_matrix[index]):
                heatmap_square = {"col_index": rank,
                                  "row_index": gene_index,
                                  "value": gene_expression}
                heatmap_data.append(heatmap_square)
        return (heatmap_data, ordered_sample_labels)


if __name__ == "__main__":
    arguments = docopt(
        __doc__, version="web application setup, 1.0")
    network_file = arguments["<network>"]
    db_credentials_file = arguments["<db-credentials>"]

    n_features = int(arguments["--features"])

    keep_n_samples = int(arguments["--samples"])
    keep_n_genes = int(arguments["--genes"])

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
    db.pathcore_edge_data.drop()

    network = pd.read_table(network_file)
    compute_edge_info = EdgeInfo(db,
                                 n_features,
                                 keep_n_samples=keep_n_samples,
                                 keep_n_genes=keep_n_genes)

    edge_info_list = []
    for index, row in network.iterrows():
        edge_info = compute_edge_info.get_edge_info(
            row["pw0"], row["pw1"])
        edge_info["edge"] = (row["pw0"], row["pw1"])
        edge_info["weight_odds_ratio"] = row["weight"]
        edge_info_list.append(edge_info)
    db.pathcore_edge_data.insert_many(edge_info_list)
    db.pathcore_edge_data.create_index(
         [("edge", ASCENDING),
          ("weight_odds_ratio", ASCENDING)])
    print("Completed the PathCORE demo server insert operations")
