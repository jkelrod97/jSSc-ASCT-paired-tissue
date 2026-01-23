# In this Python script, we run the multimodal sccomposite multiplet detection algorithm for our skin, 
# using the MTX file exported from R (3_skin_preprocessing.Rmd).
#
# Output:
# Cuda is available; Fitting the COMPOSITE model on the RNA modality
# Found 1 GPUs available. Using GPU 0 (NVIDIA GeForce RTX 4090) of compute capability 8.9 with 25.4Gb total memory.

# Warning: too few stable features to provide reliable inference.
# The RNA modality goodness-of-fit score is: 4.372865308179289 
# <3: poor fit 
# 3~5: moderate fit 
# >5: good fit

# Set project path
project_path = "/home/jelrod/scl-CITE-SEQ/preprocessing/multiplet_detection/"
# Paths to RNA and ADT matrices (exported from R)
rna_path = project_path + "skin/RNA.mtx"
# Path where we will save output
save_path = project_path + "skin/multiplet_prediction.csv"

# Load sccomposite package
# See https://github.com/CHPGenetics/COMPOSITE
import sccomposite 
# Run COMPOSITE algorithm for multiplet detection
multiplet_classification, consistency = sccomposite.RNA_modality.composite_rna(rna_path)

# Load pandas
import pandas as pd

# Create dictionary
data = {'multiplet_classification': multiplet_classification}
# Convert to pandas DataFrame
data_file = pd.DataFrame(data)
# Name index explicitly
data_file.index.name = 'index'
# Move the index into a regular column of the DataFrame.
data_file.reset_index(inplace=True)
# Save as CSV
data_file.to_csv(save_path, index=False)
