# In this script, we demultiplex each PBMC library using cellranger multi. 
#
# We run cellranger on a distinct remote server for each library. We use tmux to
# keep the session going while we are logged off. 
#
# Note that the paths and folders referenced here refer to remote paths with a different
# file structure than appears in this repository. Config CSVs referenced here are found in
# cellranger_csvs/demultiplexing. The HTO barcodes used for demultiplexing are defined in
# 5p_hashing_dmux_cmo-set.csv, which is found in the same subfolder. 
#
# This script is followed up by 0_B_align_cellranger_multi_pbmc.sh, where we will merge our
# protein (ADT), T cell receptor (TCR), and B cell receptor (BCR) data with our demultiplexed 
# gene expression (GEX data). 

# Sample 591
ssh jelrod@hydra1.stat.cmu.edu
tmux # session 0
cellranger multi --id=demultiplexed_samples \
--csv=/home/jelrod/scl-CITE-SEQ/cellranger_csvs/591_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/cellranger_output/cellranger_multi_with_bam/591

# Sample 515
ssh jelrod@hydra2.stat.cmu.edu
tmux # session 11
cellranger multi --id=demultiplexed_samples_2 \
--csv=/home/jelrod/scl-CITE-SEQ/cellranger_csvs/515_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/cellranger_output/cellranger_multi_with_bam/515

# Sample 594
ssh jelrod@hydra3.stat.cmu.edu
tmux # session 16
cellranger multi --id=demultiplexed_samples_3 \
--csv=/home/jelrod/scl-CITE-SEQ/cellranger_csvs/594_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/cellranger_output/cellranger_multi_with_bam/594

# Healthy controls
ssh jelrod@hydra4.stat.cmu.edu
tmux # session 2
cellranger multi --id=demultiplexed_samples_4 \
--csv=/home/jelrod/scl-CITE-SEQ/cellranger_csvs/H_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/cellranger_output/cellranger_multi_with_bam/H