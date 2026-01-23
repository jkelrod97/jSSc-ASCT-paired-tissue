# This script contains the following functions
#  1) run_dsb - runs DSB algorithm
#  2) plot_dsb - plots distribution of cells and empty droplets
#  3) dsb_md - helper for run_dsb that prepares metadata for DSB
#  4) match_seurat - helper for run_dsb that renames matrix from 
#     cell ranger to match cell names in Seurat object


################################################################################
# (1) RUN_DSB ##################################################################
################################################################################

# Load dependencies
library(dsb)

# Function
run_dsb <- function(raw_counts = NULL, 
                    filtered_counts = NULL, 
                    seurat_obj = NULL, 
                    subject = NULL,
                    prot.size.min = 0,
                    prot.size.max = Inf,
                    rna.size.min = 0,
                    rna.size.max = Inf,
                    mt.prop.max = 1,
                    ...) {
  # RUN_DSB
  # Function that runs dsb for raw and filtered cellranger ADT counts.
  # Adapted from https://cran.r-project.org/web/packages/dsb/vignettes/end_to_end_workflow.html
  #
  # Arguments:
  # (1) raw_counts (matrix): raw counts from cellranger
  # (2) filtered_counts (matrix): filtered counts from cellranger
  # (3) seurat_obj (Seurat object): Seurat object
  # (4) subject (string): subject ID
  # (5) prot.size.min (integer - default 0): prot.size is log10(total counts per cell)--set minimum value for background droplet filtering
  # (6) prot.size.max (integer - default infinity): prot.size is log10(total counts per cell)--set maximum value for background droplet filtering
  # (7) rna.size.min (integer - default 0): rna.size is log10(total counts per cell)--set minimum value for background droplet filtering
  # (8) rna.size.max (integer - default infinity): rna.size is log10(total counts per cell)--set maximum value for background droplet filtering
  # (9) mt.prop.max (integer - default 1. Must be between 0 and 1): maximum proportion of mitochondrial RNA allowed for background droplets
  #
  # Outputs:
  # cells.dsb.norm (matrix): denoised and normalized ADT counts matrix
  
  # (1) PREPARE METADATA
  
  # Rename cells in cell ranger output to match Seurat object
  # Also make sure genes & antibodies match those in Seurat object
  raw_counts <- mapply(match_seurat, raw_counts, list("RNA", "ADT"))
  # Apply to filtered counts
  filtered_counts <- mapply(match_seurat, filtered_counts, c("RNA", "ADT"))
  
  # Split the data into separate matrices for RNA and ADT
  prot <- raw_counts$`Antibody Capture`
  rna <- raw_counts$`Gene Expression`
  
  # Now calculate some standard meta data for cells that we will use for quality control 
  # of the background droplets.
  
  # Create metadata of droplet QC stats used in standard scRNAseq processing
  md <- dsb_md(raw_counts, filtered_counts)
  
  # (2) FILTER BACKGROUND DROPLETS
  
  # Filter background droplets based on supplied cutoffs
  
  # Get cell IDs
  background_drops <- md %>% 
    filter(prot.size > prot.size.min &
             prot.size <  prot.size.max &
             rna.size > rna.size.min &
             rna.size < rna.size.max &
             mt.prop < mt.prop.max) %>% 
    rownames()
  
  # Get background ADT matrix
  background_adt_mtx <- prot[, background_drops]
  
  # Make rownames match Seurat object
  rownames(background_adt_mtx) <- gsub("_", "-", rownames(background_adt_mtx))
  
  # (3) FILTER CELLS to only include those in Seurat object (already went through QC)
  
  seurat_adt_counts <- seurat_obj[["ADT"]]$counts

  # Get names of isotype controls
  isotype_ctls <- rownames(seurat_adt_counts)[str_detect(rownames(seurat_adt_counts), fixed("Iso", ignore_case = TRUE))]

  # (5) NORMALIZE AND DENOISE WITH DSB
  
  cells.dsb.norm <- DSBNormalizeProtein(
    cell_protein_matrix = seurat_adt_counts,
    empty_drop_matrix = background_adt_mtx,
    isotype.control.name.vec = isotype_ctls,
    ...
  )
  return(cells.dsb.norm)
}

################################################################################
# (2) PLOT_DSB #################################################################
################################################################################

plot_dsb <- function(raw_counts, filtered_counts) {
  # PLOT_DSB
  # User can view distribution of cells and empty droplets.
  # This is useful for determining appropriate prot.size and rna.size cutoff values
  # for run_dsb.
  #
  # Arguments:
  # (1) raw_counts (matrix): raw counts from cellranger
  # (2) filtered_counts (matrix): filtered counts from cellranger
  #
  # Output:
  # p1 + p2: p1 contains two plots, one for cells and one for empty droplets.
  #   prot.size is plotted on the x axis and rna.size is plotted on the y axis.
  #   Density of droplets with each particular combo is shown with a color scale.
  #   p2 has the same x and y axis, but the color represents the proportion of
  #   mitochondrial RNA. A high proportion of mitochondrial RNA indicates low quality
  #   cells, which should be removed.
  
  # Get meta data
  md <- dsb_md(raw_counts, filtered_counts)
  
  # Plot counts
  p1 <- ggplot(md, aes(x = prot.size, y = rna.size)) +
    geom_bin2d(bins = 300) +
    scale_fill_viridis_c(option = "C") +
    facet_wrap(~drop.class) +
    theme_bw()
  
  # Plot MT proportion
  p2 <- ggplot(md, aes(x = prot.size, y = rna.size, color = mt.prop)) +
    geom_point(size = 0.5, alpha = 0.7) +
    scale_color_viridis_c(option = "C", limits = c(0, 1)) +
    facet_wrap(~drop.class) +
    theme_bw()
  
  # Return plots
  p1 + p2
}

################################################################################
# (3) DSB_MD ###################################################################
################################################################################

dsb_md <- function(raw_counts, filtered_counts) {
  # DSB_MD
  # Helper function for run_dsb
  # Prepares metadata.
  #
  # Arguments:
  # (1) raw_counts (matrix): raw counts from cellranger
  # (2) filtered_counts (matrix): filtered counts from cellranger
  #
  # Output:
  # md (data.frame): data frame with each row representing a cell. 
  #   Contains rna.size, prot.size, n.gene, mt.prop, and
  #   an indicator for cell vs. empty droplet. Droplets with no evidence of capture
  #   are removed. 
  
  # Define cell-containing barcodes
  stained_cells <- colnames(filtered_counts$`Gene Expression`)
  # Define empty droplets
  background <- setdiff(colnames(raw_counts$`Gene Expression`), stained_cells)
  
  # Split the data into separate matrices for RNA and ADT
  prot <- raw_counts$`Antibody Capture`
  rna <- raw_counts$`Gene Expression`
  
  # Now calculate some standard meta data for cells that we will use for quality control 
  # using standard approaches used in scRNAseq analysis pipelines.
  
  # Create metadata of droplet QC stats used in standard scRNAseq processing
  
  # Mitochondrial DNA genes
  mtgene <- grep(pattern = "^MT-", rownames(rna), value = TRUE) # used below
  
  # QC metadata
  md <- data.frame(
    rna.size = log10(Matrix::colSums(rna)), 
    prot.size = log10(Matrix::colSums(prot)), 
    n.gene = Matrix::colSums(rna > 0), 
    mt.prop = Matrix::colSums(rna[mtgene, ]) / Matrix::colSums(rna)
  )
  
  # Add indicator for barcodes Cell Ranger called as cells
  md$drop.class <- ifelse(rownames(md) %in% stained_cells, 'cell', 'background')
  
  # Remove barcodes with no evidence of capture in the experiment
  md <- md[md$rna.size > 0 & md$prot.size > 0, ]
  
  # Return
  return(md)
}

################################################################################
# (4) MATCH_SEURAT #############################################################
################################################################################

# Function for this step
match_seurat <- function(mat, seurat_obj, assay){
  # MATCH_SEURAT:
  # Helper function for run_dsb
  # Rename matrix from cellranger to match cell names in Seurat object
  # Arguments:
  # (1) mat (matrix): counts matrix from cell ranger (can be raw or filtered)
  # (2) seurat_obj (Seurat): Seurat object
  # (3) assay (string): assay that counts matrix corresponds to in Seurat object
  
  # Add subject + "-" to head of colname. Remove trailing 1s.
  colnames(mat) <- paste0(subject, "-", gsub("-1", "", colnames(mat)))
  # Change underscores to dashes. 
  rownames(mat) <- gsub("_", "-", rownames(mat))
  # Make sure order is the same in mat and seurat_obj
  mat <- mat[rownames(seurat_obj[[assay]]), ]
  # Return mat with new colnames/rownames
  mat
  }
