#!/bin/bash
set -o errexit

# Kathleen Chen 2017

pip install -r requirements.txt

# TODO: will be available in the next pull request.
# cd data
# ./download_data.sh
# cd ..

wait 

echo "Log information for all analyses is written to *.log files " \
"in the './log' directory"

# Scripts used to generate the results described in the PathCORE paper.
# Ordered based on execution time of each script. They can be modified
# and run independently of each other.

# Supplemental result
./analysis_Paeruginosa_NMF.sh

# Case study 2
./analysis_TCGA_NMF.sh

# Case study 1
./analysis_Paeruginosa_eADAGE.sh

exit 0
