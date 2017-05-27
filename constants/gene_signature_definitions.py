"""Stores the gene signature definitions that can be used in
`build_co_occurrence_network.py`.

NOTE: any arguments passed in to the top-level
`define_<model>_gene_signature` function must be cast
into the type specified in the function's documentation
because docopt will parse the signature argument input as a string.
(Refers to:
  - lines 31, 68 below
  - lines 266-268 in `build_co_occurrence_network.py`)
"""
from functools import partial


def define_eADAGE_gene_signature(std):
    """Provide a criterion for evaluating whether a gene is part of
    an eADAGE feature's gene signature

    Parameters
    -----------
    std : float
        Only consider "high weight" genes: genes with weights +/- `std`
        standard deviations from the mean.

    Returns
    -----------
    functools.partial that accepts a feature weight vector of type
    pandas.Series(float) and returns (set(), set()), the feature's
    positive and negative gene signatures.
    """
    std = float(std)  # see note on line 4.

    def _gene_signature(feature_weight_vector, std):
        """This definition assumes that the weights in the feature
        are normally distributed and we consider the genes above
        and below some threshold to be more important.
        """
        mean = feature_weight_vector.mean()
        cutoff = std * feature_weight_vector.std()
        positive_gene_signature = set(
            feature_weight_vector[(feature_weight_vector >=
                                   mean + cutoff)].index)
        negative_gene_signature = set(
                feature_weight_vector[(feature_weight_vector <=
                                       mean - cutoff)].index)
        return positive_gene_signature, negative_gene_signature

    return partial(_gene_signature, std=std)


def define_NMF_gene_signature(std):
    """Provide a criterion for evaluating whether a gene is part of
    a NMF feature's gene signature

    Parameters
    -----------
    std : float
        Only consider genes with weights above `std` standard deviations
        from the mean because the weight distribution in an NMF feature
        is nonnegative and right-skewed.

    Returns
    -----------
    functools.partial that accepts a feature weight vector of type
    pandas.Series(float) and returns (set(), set()), the feature's
    positive and negative gene signatures.
    """
    std = float(std)  # see note on line 4.
    def _gene_signature(feature_weight_vector, std):
        mean = feature_weight_vector.mean()
        cutoff = std * feature_weight_vector.std()
        positive_gene_signature = set(
            feature_weight_vector[(feature_weight_vector >=
                                   mean + cutoff)].index)
        return (positive_gene_signature, set())
    return partial(_gene_signature, std=std)

GENE_SIGNATURE_DEFINITIONS = {
	"eADAGE": define_eADAGE_gene_signature,
	"NMF": define_NMF_gene_signature
}