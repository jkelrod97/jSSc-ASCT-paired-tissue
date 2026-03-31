Code accompanying longitudinal scRNA-seq jSSc ASCT monocytes paper (Elrod et al.)

Torok lab, UPMC Children's hospital.

Townes lab, Carnegie Mellon Statistics & Data Science.

Files are executed in this order:
1. **0_A_demultiplex_cellranger_multi_pbmc.sh**: Demultiplex and align PBMC RNA samples with cellranger multi. 
2. **0_B_align_cellranger_multi_pbmc.sh**: Merge demultiplexed PBMC RNA data with other modalities: protein (ADT), T cell receptor (TCR), and B cell receptor (BCR).
3. **0_C_alignment_without_demultiplexing_for_denoising.sh**: Align PBMC RNA data with other modalities without demultiplexing first. Output used for RNA and protein denoising.
4. **1_pbmc_preprocessing.Rmd**: Preprocess multimodal PBMC data using Seurat and standard single-cell pipeline.
5. **1a_composite_pbmc.py**: Run scomposite algorithm for multiplet detection. https://github.com/CHPGenetics/COMPOSITE
6. **2_pbmc_cell_type_annotation.Rmd**: Perform multimodal clustering and cell type annotation on PBMCs.
7. **3_monocyte_subtype_annotation.Rmd**: Sub-cluster monocytes and annotate subtypes.
8. **4_monocyte_longitudinal_analysis.Rmd**: Implement longitudinal analysis to identify differentially express genes in jSSc monocytes over time since ASCT. 
