################################################################################
# DISPLAY_TOP ##################################################################
################################################################################

# Load dependencies
library(dplyr)

# Function
display_top <- function(de_df, cluster_id, n = 20) {
  # DISPLAY_TOP: Function that accepts differential expression data frame from
  # Seurat's FindAllMarkers() function and shows the top n (default 20) differentially
  # expressed features for a given cluster.
  #
  # Arguments:
  # (1) de_df (data frame): Output from Seurat's FindAllMarkers().
  # (2) cluster_id (integer): Cluster of interest.
  # (3) n (integer): Number of features to display.
  #
  # Output:
  # de_df_filtered (data frame): Data frame showing results of FindAllMarkers()
  # for the top n features with LFC > 0 for a specified cluster.
  
  # Filter to only include specific luster, LFC > 0, and n features:
  de_df_filtered <- de_df %>% filter(cluster == cluster_id & avg_log2FC > 0) %>% head(n = n)
  
  # Delete row names
  rownames(de_df_filtered) <- NULL
  
  # Create 'feature' column (so that it is applicable to modalities beyond scRNA-seq)
  de_df_filtered$feature <- de_df_filtered$gene
  
  # Delete gene column
  de_df_filtered$gene <- NULL
  
  # Return
  return(de_df_filtered)
}
