"""Stores functions to shorten pathway names, if desired, during
PathCORE network file generation.
Functions should take in a single input str (pathway name) and
output a shortened version of the pathway name based on known
conventions.
"""


def paeruginosa_shorten_kegg(pathway_name):
    """Based on the naming conventions for PAO1 KEGG pathways.
    """
    REMOVE_SUFFIX = "- Pseudomonas aeruginosa PAO1"
    pathway_short = None
    split_label = pathway_name.split(" ", 1)
    if len(split_label) > 1:
        pathway_short = split_label[1]
    else:
        pathway_short = split_label[0]
    if REMOVE_SUFFIX in pathway_short:
        remove_from_index = pathway_short.index(REMOVE_SUFFIX)
        return "{0}PAO1".format(pathway_short[:remove_from_index])
    return pathway_short.strip()


def tcga_shorten_pid(pathway_name):
    """Based on the naming conventions for Nature-NCI PID pathways.
    """
    REMOVE_SUFFIX = "PATHWAY"
    pathway_short = None
    split_on_underscores = pathway_name.split("_")
    if split_on_underscores[-1] == REMOVE_SUFFIX:
        pathway_short = " ".join(split_on_underscores[:-1])
    else:
        pathway_short = " ".join(split_on_underscores)
    return pathway_short


SHORTEN_PATHWAY_NAMES = {
    "PAO1_KEGG": paeruginosa_shorten_kegg,
    "TCGA_PID": tcga_shorten_pid
}
