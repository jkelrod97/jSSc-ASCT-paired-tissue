################################################################################
# RUN_RPCA #####################################################################
################################################################################

# Load dependencies
library(Seurat)
library(stringr)

# Function
run_rpca <- function(seurat_obj, n_dim = 30, k.weight = 100, k.score = 30) {
  # Function to integrate ADT data (already denoised and normalized with DSB) 
  # across batches using reciprocal PCA (RPCA).
  #
  # Arguments:
  # (1) seurat_obj (Seurat): Seurat object with RNA and DSB assays. DSB is denoised and normalized
  #     ADT data from CITE-seq experiment.
  # (2) n_dim (integer): Number of dimensions
  # (3) k.weight (integer): argument for IntegrateData
  # (4) k.score (integer): argument for IntegrateData
  #
  # Output:
  # seurat_obj (Seurat): Updated Seurat object with RNA and DSB_integrated assays. DSB_integrated
  #   is an updated version of DSB, which has been integrated across libraries to correct
  #   for batch effects using RPCA.
  
  # Set default assay
  DefaultAssay(seurat_obj) <- "DSB"
  
  # Split the dataset into a list of seurat objects (one for each library)
  spbmc.list <- SplitObject(seurat_obj, split.by = "LIBRARY")
  # Get number of libraries
  n_lib <- length(spbmc.list)
  # Get names of libraries
  libraries <- names(spbmc.list)
  
  # Remove isotype controls from features to consider
  isotype_ctls <- rownames(seurat_obj[["DSB"]])[str_detect(rownames(seurat_obj[["DSB"]]), fixed("Iso", ignore_case = TRUE))]
  dsb_features <- setdiff(rownames(seurat_obj[["DSB"]]), isotype_ctls)
  
  # For each Seurat object in the list, rescale data and run PCA.
  spbmc.list <- lapply(libraries, FUN = function(SUBJ_NAME) {
    # Extract Seurat object
    x <- spbmc.list[[SUBJ_NAME]]
    # Rescale data
    x <- ScaleData(x, features = dsb_features)
    # Run PCA
    x <- RunPCA(x, features = dsb_features, npcs = n_dim)
    return(x)
  })
  # Name according to library
  names(spbmc.list) = libraries
  
  # Get anchors
  immune_anchors <- FindIntegrationAnchors(object.list = spbmc.list, 
                                           anchor.features = dsb_features, 
                                           reduction = "rpca",
                                           k.score = k.score,
                                           assay = rep("DSB", n_lib),
                                           dims = 1:n_dim)
  
  # Integrate data
  integ_data <- IntegrateData(anchorset = immune_anchors, dims = 1:n_dim, k.weight = k.weight)
  # Add to original seurat object
  seurat_obj[["DSB_integrated"]] <- integ_data[["integrated"]]
  # Set default assay
  DefaultAssay(seurat_obj) <- "DSB_integrated"

  
  # Return seurat object with only two assays--RNA and DSB_integrated
  return(seurat_obj)
}