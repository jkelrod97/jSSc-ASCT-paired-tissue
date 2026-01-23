# In this Python script, we run the multimodal sccomposite multiplet detection algorithm for our PBMCs, 
# using MTX files exported from R (1_pbmc_preprocessing.Rmd).
#
# Output:
# The RNA modality goodness-of-fit score is: 3.9210535234945745 
# <3: poor fit 
# 3~5: moderate fit 
# >5: good fit
# Cuda is available; Fitting the COMPOSITE model on the ADT modality
# Found 1 GPUs available. Using GPU 0 (NVIDIA GeForce RTX 4090) of compute capability 8.9 with 25.4Gb total memory.

# The ADT modality goodness-of-fit score is: 15.113854780774897 
# <3: poor fit 
# 3~5: moderate fit 
# >5: good fit

# Set project path
project_path = "/home/jelrod/scl-CITE-SEQ/preprocessing/multiplet_detection/"
# Paths to RNA and ADT matrices (exported from R)
rna_path = project_path + "PBMC/RNA.mtx"
adt_path = project_path + "PBMC/ADT.mtx"
# Path where we will save output
save_path = project_path + "PBMC/multiplet_prediction.csv"

# Load sccomposite package
# See https://github.com/CHPGenetics/COMPOSITE
import sccomposite 
# Run COMPOSITE algorithm for multiplet detection
multiplet_classification, multiplet_probability = sccomposite.Multiomics.composite_multiomics(RNA = rna_path, ADT =  adt_path)

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
