# Description

Two [Jupyter](http://jupyter.org/) notebooks in this directory
demonstrate the analyses that support Figure 3 and Supplementary Figure 2
in the PathCORE paper:
- [Figure3_overlap_correction.ipynb](Figure3_overlap_correction.ipynb)
- [Supplemental2_TCGA_NMF_features.ipynb](Supplemental2_TCGA_NMF_features.ipynb)

They must be run after `../ANALYSIS.sh`. Additional Python dependencies
for the notebooks are specified in the `./requirements.txt` file for this
directory.

An additional notebook is included as a supplemental example of the PathCORE
analysis workflow:
- [Supplemental_PAO1_FastICA_example.ipynb](Supplemental_PAO1_FastICA_example.ipynb)

In this notebook, we apply the PathCORE software to a FastICA model
of the normalized _P. aeruginosa_ gene compendium. We apply
[scikit-learn's FastICA implementation](http://scikit-learn.org/stable/modules/generated/sklearn.decomposition.fastica.html). 

Similar to the case studies we describe in our manuscript, we set the number
of features (ICA components) we construct to 300. A user interested in
changing this parameter and examining the differences in the resulting
network can do so in cell [4] of the notebook and re-run the analysis.