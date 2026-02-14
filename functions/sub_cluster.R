# Load packages
library(Seurat)
library(dplyr)
library(harmony)

# Source function
setwd(this.path::this.dir())
source("run_rpca.R")

# Custom colors
my_colors <- c("skyblue",
               "#D175A0",
               "#9CCC65",
               "#DE1F6C",
               "violet",
               "steelblue",
               "#99169E",
               "gold",
               "orange",
               "navy",
               "lightpink",
               "tan",
               "red",
               "lightgrey",
               "purple",
               "darkgreen",
               "yellow",
               "slategrey",
               "magenta",
               "coral4",
               "lavender",
               "darkcyan",
               "cornflowerblue",
               "salmon2",
               "bisque1",
               "black",
               "indianred",
               "darkseagreen",
               "gold4",
               "darkmagenta",
               "tomato1",
               "yellow2") 


sub_cluster <- function(obj, 
                        dims = c(30, 18), 
                        res = 1, 
                        cluster.only = FALSE,
                        algorithm = 3,
                        assays = c("RNA", "DSB"), 
                        k.weight = 100,
                        k.score = 30,
                        knn.range = 200) {
  # Store object's original default assay
  orig_assay <- DefaultAssay(obj)
  
  # Remove old clusters
  obj$seurat_clusters <- NULL
  
  # First, do RNA steps if needed
  if (cluster.only == FALSE) {
    # RNA
    DefaultAssay(obj) <- "RNA"
    # Normalize
    obj <- NormalizeData(obj)
    # Variable features
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
    # Scaling
    obj <- ScaleData(obj, features = rownames(obj))
    # PCA
    obj <- RunPCA(obj, features = VariableFeatures(obj), reduction.name = "pca_rna", npcs = dims[1])
    # Harmony (integration)
    nlib <- obj$LIBRARY %>% unique() %>% length()
    if (nlib > 1) {
      obj <- RunHarmony(obj, c("LIBRARY"), verbose = TRUE, reduction.use = "pca_rna", dims.use = 1:dims[1])
    }
  }
  
  
  # Option A: multimodal clustering
  if ("DSB" %in% assays) {
    # Option A1: Everything
    if (cluster.only == FALSE) {
      # DSB ADT
      DefaultAssay(obj) <- "DSB"
      # Remove old integrated DSB
      obj[["DSB_integrated"]] <- NULL
      # Scale DSB output
      obj <- ScaleData(obj, features = rownames(obj))
      # Run RPCA
      if (nlib > 1) {
        obj <- run_rpca(obj, n_dim = dims[2], k.weight = k.weight, k.score = k.score)
      # RNA harmony reduction is ready to be used with WNN.
      # Now, set up integrated DSB values to use in WNN analysis.
      DefaultAssay(obj) <- "DSB_integrated"
      } else {DefaultAssay(obj) <- "DSB"}
      obj <- obj %>%
        ScaleData(features = rownames(obj)) %>%
        RunPCA(reduction.name = 'pca_dsb_int', verbose = FALSE, npcs = dims[2])
      
      # Run WNN
      if (nlib > 1) {reduc_list <- list("harmony", "pca_dsb_int")} else {
        reduc_list <- list("pca_rna", "pca_dsb_int")
      } 
      obj <- FindMultiModalNeighbors(
        obj, reduction.list = reduc_list,
        dims.list = list(1:dims[1], 1:dims[2]),
        modality.weight.name = "RNA.weight",
        verbose = FALSE,
        knn.range = knn.range
      )
      
      # Run UMAP
      obj <- RunUMAP(obj, 
                     nn.name = "weighted.nn", 
                     reduction.name = "wnn.umap", 
                     reduction.key = "wnnUMAP_")
    # Option A2: Clustering only
    }
    # Cluster using SLM algorithm (3)
    obj <- FindClusters(obj, graph.name = "wsnn",
                        algorithm = algorithm,
                        resolution = res,
                        random.seed = 1990)
    
    # Print UMAP
    if (length(my_colors) < length(levels(obj$seurat_clusters))) {
      print(DimPlot(obj, label = T, reduction = "wnn.umap")) # use default colors if there aren't enough in my_colors
    } else {
      print(DimPlot(obj, label = T, cols = my_colors, reduction = "wnn.umap"))
    }
  # Option B: standard RNA clustering
  } else {
    # Option B1: Everything
    if (cluster.only == FALSE) {
      # Run UMAP
      obj <- RunUMAP(obj, reduction = "harmony", dims = 1:dims[1], verbose = TRUE)
      
      # Run WNN
      obj <- FindNeighbors(
        obj, reduction = "harmony",
        dims = 1:dims[1])
    # Option B2: Clustering only
    }
    # Cluster using SLM algorithm (3)
    obj <- FindClusters(obj,
                        algorithm = algorithm,
                        resolution = res,
                        random.seed = 20)
    
    # Print UMAP
    if (length(my_colors) < length(levels(obj$seurat_clusters))) {
      print(DimPlot(obj, label = T, reduction = "umap")) # use default colors if there aren't enough in my_colors
    } else {
      print(DimPlot(obj, label = T, cols = my_colors, reduction = "umap"))
    }
  }
  # Reset to original default assay
  DefaultAssay(obj) <- orig_assay
  return(obj)
}

  