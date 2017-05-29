"""Functions used to process PAO1-specific information we want to
provide in the web interface. A user can modify #TODO to setup a
web application with the information most helpful for analyzing
the edges in the user's PathCORE network.
"""
import pandas as pd


def get_sample_annotations(db_sample_annotations,
                           db_sample_labels,
                           sample_annotations_file):
    """Updates the MongoDB instance with a sample annotations collection.

    Parameters
    -----------
    db_sample_annotations : pymongo.collection.Collection
      The Collection we want to populate with sample annotation information
    db_sample_labels : pymongo.collection.Collection
      The Collection of sample labels (columns in the expression dataset)
      for which we want sample annotations.
    sample_annotations_file : str
      A .tsv file containing annotations for many of the experiments in the
      compendium. The column "CEL file" should match up with the "sample"
      field in `db_sample_labels`--this allows us to only insert annotations
      to existing samples, which we can refer to by their ObjectIds.

    Returns
    -----------
    None, updates the sample annotations Collection.
    """
    sample_annotation_table = pd.read_table(sample_annotations_file)
    sample_annotation_list = []
    for index, row in sample_annotation_table.iterrows():
        retrieve_sample = db_sample_labels.find_one(
            {"sample": row["CEL file"]})
        if retrieve_sample:
            sample_annotation = row.to_dict()
            for key, value in sample_annotation.items():
                if pd.isnull(value):
                    sample_annotation[key] = ""
                if "." in key:
                    new_key = key.replace(".", "")
                    sample_annotation[new_key] = sample_annotation.pop(key)
            # map the annotation to the sample
            sample_annotation["sample_id"] = retrieve_sample["_id"]
            sample_annotation_list.append(sample_annotation)
    db_sample_annotations.insert_many(sample_annotation_list, ordered=True)
    db_sample_annotations.create_index("CEL file")


def get_gene_common_names(db_genes, gene_names_file):
    """Updates the genes collection with common names

    Parameters
    -----------
    db_genes : pymongo.collection.Collection
      The Collection we want to update with a new "common_name" field.
    gene_names_file : str
      A .gene_info file used for matching up the PA gene identifiers to
      common names when they are available.

    Returns
    -----------
    None, updates the genes Collection.
    """
    gene_names_table = pd.read_table(gene_names_file)
    for index, row in gene_names_table.iterrows():
        pa_name = row["LocusTag"]
        pa_or_common_name = row["Symbol"]
        if not pd.isnull(pa_or_common_name) and pa_or_common_name != pa_name:
            db_genes.find_one_and_update(
                {"gene": pa_name},
                {"$set": {"common_name": pa_or_common_name}})
