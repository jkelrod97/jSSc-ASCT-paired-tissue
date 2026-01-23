################################################################################
# RUN_SOUP #####################################################################
################################################################################

# Load dependencies
library(SoupX)
library(Seurat)

# Function
run_soup <- function(raw_counts, filtered_counts, subject) {
  # RUN_SOUP
  # Function that runs SoupX for raw and filtered cellranger counts.
  # Adapted from https://github.com/grabadam-cal/jdm-DECIPHER-2024/blob/main/Preprocessing/process_jdm/run_soupx_well.R
  #
  # Arguments:
  # (1) raw_counts (matrix): raw counts from cellranger
  # (2) filtered_counts (matrix): filtered counts from cellranger
  # (3) subject (string): subject identifier
  #
  # Outputs:
  # soup_outs (list): list containing the following:
  #    - adj_counts: matrix of adjusted counts
  #    - sc: SoupX output
  
  # Extract gene expression
  if (is.list(raw_counts)) {
    raw_rna <- raw_counts$`Gene Expression`
    filtered_rna <- filtered_counts$`Gene Expression`
  } else {
    raw_rna <- raw_counts; filtered_rna <- filtered_counts
  }
  
  # Check row IDs (genes)
  stopifnot(all(rownames(filtered_rna) == rownames(raw_rna)))
  
  # Fix and check the column names
  colnames(raw_rna) <- gsub("-1", "", paste0(subject, "_", colnames(raw_rna)))
  colnames(filtered_rna) <- gsub("-1", "", paste0(subject, "_", colnames(filtered_rna)))
  
  # Store in soup channel--this is just a list with some special properties, 
  # storing all the information associated with a single 10X channel
  sc <- SoupChannel(raw_rna, filtered_rna)
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(filtered_rna)
  
  # Perform normalization and clustering for SoupX
  seurat_obj <- NormalizeData(seurat_obj) %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    RunUMAP(dims = 1:30)
  seurat_obj <- FindNeighbors(seurat_obj)
  seurat_obj <- FindClusters(seurat_obj, res = 0.5)
  
  # Set clusters
  sc <- setClusters(sc, setNames(seurat_obj@meta.data$seurat_clusters,
                                  rownames(seurat_obj@meta.data)))
  
  # Set dimension reduction
  sc <- setDR(sc, seurat_obj@reductions$umap@cell.embeddings[colnames(sc$toc),])
  
  # Estimate contamination fraction
  sc <- autoEstCont(sc)
  
  # Adjust counts
  adj_counts <- adjustCounts(sc)
  
  # Outputs to return
  soup_out <- list(adj_counts = adj_counts,
                   sc = sc)
  
  # Return
  return(soup_out)
}