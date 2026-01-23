# This script is a follow-up to 0_A_demultiplex_cellranger_multi_pbmc.sh.
# 
# In the alignment step, we match up the protein (ADT), T cell receptor (TCR), and B cell receptor (BCR)
# data to the per-sample gene expression (GEX) data. 
#
# This requires 15 separate runs of cellranger multi, one per sample.
#
# We are starting with sample 515 (recall we are missing the twelve-month benchmark for this sample):
#
# Once again, we use a distinct remote server and tmux session for each cellranger run.
#
# Note that the paths and folders referenced here refer to remote paths with a different
# file structure than appears in this repository. Config CSVs referenced here are found in
# cellranger_csvs/alignment.

#############
## SAMPLE 515
#############

# Baseline
hydra1
# TMUX SESSION 20
tmux
cellranger multi --id=515_baseline --csv=/home/jelrod/scl-CITE-SEQ/cellranger_csvs/final_configs/515_baseline_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/cellranger_output/final_analysis/output/515_baseline

# Six mo.
hydra3
# TMUX SESSION 22
cellranger multi --id=515_six --csv=/home/jelrod/scl-CITE-SEQ/cellranger_csvs/final_configs/515_six_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/cellranger_output/final_analysis/output/515_six

# Twenty-four mo.
hydra4
# TMUX SESSION 5
cellranger multi --id=515_twen4 --csv=/home/jelrod/scl-CITE-SEQ/cellranger_csvs/final_configs/515_twen4_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/cellranger_output/final_analysis/output/515_twen4

#############
## SAMPLE 591
#############

# Baseline
hydra10
# TMUX SESSION 4
tmux
cellranger multi --id=591_baseline --csv=/home/jelrod/scl-CITE-SEQ/cellranger_csvs/final_configs/591_baseline_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/cellranger_output/final_analysis/output/591_baseline

# Six mo.
hydra11
# TMUX SESSION 3
tmux
cellranger multi --id=591_six --csv=/home/jelrod/scl-CITE-SEQ/cellranger_csvs/final_configs/591_six_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/cellranger_output/final_analysis/output/591_six

# Twelve mo.
hydra1
# TMUX SESSION 21
tmux
cellranger multi --id=591_twelve --csv=/home/jelrod/scl-CITE-SEQ/cellranger_csvs/final_configs/591_twelve_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/cellranger_output/final_analysis/output/591_twelve


# Twenty-four mo.
hydra3
# TMUX SESSION 24
tmux
# Delete & recreate output directory (to get rid of failed pipestance)
rm -rf /home/jelrod/scl-CITE-SEQ/cellranger_output/final_analysis/output/591_twen4
mkdir /home/jelrod/scl-CITE-SEQ/cellranger_output/final_analysis/output/591_twen4
cellranger multi --id=591_twen4 --csv=/home/jelrod/scl-CITE-SEQ/cellranger_csvs/final_configs/591_twen4_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/cellranger_output/final_analysis/output/591_twen4


#############
## SAMPLE 594
#############

# Baseline
hydra4
# TMUX SESSION 6
tmux
cellranger multi --id=594_baseline --csv=/home/jelrod/scl-CITE-SEQ/cellranger_csvs/final_configs/594_baseline_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/cellranger_output/final_analysis/output/594_baseline

# Six mo.
hydra7
# TMUX SESSION 5
tmux
cellranger multi --id=594_six --csv=/home/jelrod/scl-CITE-SEQ/cellranger_csvs/final_configs/594_six_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/cellranger_output/final_analysis/output/594_six

# Twelve mo.
hydra12
# TMUX SESSION 2
tmux
cellranger multi --id=594_twelve --csv=/home/jelrod/scl-CITE-SEQ/cellranger_csvs/final_configs/594_twelve_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/cellranger_output/final_analysis/output/594_twelve

# Twenty-four mo.
hydra13
# TMUX SESSION 4
tmux
cellranger multi --id=594_twen4 --csv=/home/jelrod/scl-CITE-SEQ/cellranger_csvs/final_configs/594_twen4_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/cellranger_output/final_analysis/output/594_twen4


#############
## HEALTHY ##
#############

# Baseline
hydra1
# TMUX SESSION 22
tmux
cellranger multi --id=H1 --csv=/home/jelrod/scl-CITE-SEQ/cellranger_csvs/final_configs/H1_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/cellranger_output/final_analysis/output/H1

# Six mo.
hydra4
# TMUX SESSION 7
tmux
cellranger multi --id=H2 --csv=/home/jelrod/scl-CITE-SEQ/cellranger_csvs/final_configs/H2_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/cellranger_output/final_analysis/output/H2

# Twelve mo.
hydra3
# TMUX SESSION 25
tmux
cellranger multi --id=H3 --csv=/home/jelrod/scl-CITE-SEQ/cellranger_csvs/final_configs/H3_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/cellranger_output/final_analysis/output/H3

# Twenty-four mo.
hydra7
# TMUX SESSION 6
tmux
cellranger multi --id=H4 --csv=/home/jelrod/scl-CITE-SEQ/cellranger_csvs/final_configs/H4_config.csv \
--output-dir=/home/jelrod/scl-CITE-SEQ/cellranger_output/final_analysis/output/H4
