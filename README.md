Code accompanying paired tissue (PBMC + skin) longitudinal scRNA-seq jSSc ASCT clinical trial paper (Elrod et al.)

Torok lab, UPMC Children's hospital.

Townes lab, Carnegie Mellon Statistics & Data Science.

Files are executed in this order:
1. **0_A_demultiplex_cellranger_multi_pbmc.sh**: Demultiplex and align PBMC RNA samples with cellranger multi. 
2. **0_B_align_cellranger_multi_pbmc.sh**: Merge demultiplexed PBMC RNA data with other modalities: protein (ADT), T cell receptor (TCR), and B cell receptor (BCR).
3. **0_C_align_cellranger_count_skin.sh**: Align skin scRNA-seq samples (no multiplexing; all samples were sequenced separately).
4. **1_pbmc_preprocessing.Rmd**: Preprocess multimodal PBMC data using Seurat and standard single-cell pipeline.
5. **1a_composite_pbmc.py**: Run scomposite algorithm for multiplet detection. https://github.com/CHPGenetics/COMPOSITE
6. **2_pbmc_cell_type_annotation.Rmd**: Perform multimodal clustering and cell type annotation on PBMCs.
7. **3_skin_preprocessing.Rmd**: Preprocess scRNA-seq skin data using Seurat and standard single-cell pipeline.
8. **3a_composite_skin.py**: Run sccomposite algorithm for multiplet detection. https://github.com/CHPGenetics/COMPOSITE
9. **4_skin_cell_type_annotation.Rmd**: Perform clustering and cell type annotation on skin samples.
10. **5_paired_PBMC_skin_longitudinal_analysis.Rmd**: Implement paired PBMC and skin longitudinal analysis to identify differentially express genes in jSSc patients over time since ASCT. 
