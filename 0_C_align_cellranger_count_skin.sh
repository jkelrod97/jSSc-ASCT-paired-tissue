# In this script, we align the skin samples using cellranger count. 
# This is much simpler than for the PBMCs, since we did not use
# multiplexing and only have a single modality (scRNA-seq).
#
# All of this code was executed on a remote server.

cd /home/jelrod/skin/skin-counts-oct-25
cellranger count \
  --id=515_baseline \
  --transcriptome=/home/jelrod/software/refdata-gex-GRCh38-2024-A \
  --fastqs=/home/jelrod/skin/skin-data-oct-25/515_baseline \
  --sample=SC288 \
  --create-bam=false

cd /home/jelrod/skin/skin-counts-oct-25
cellranger count \
  --id=515_six \
  --transcriptome=/home/jelrod/software/refdata-gex-GRCh38-2024-A \
  --fastqs=/home/jelrod/skin/skin-data-oct-25/515_six \
  --sample=SC357 \
  --create-bam=false

cd /home/jelrod/skin/skin-counts-oct-25
cellranger count \
  --id=515_twelve \
  --transcriptome=/home/jelrod/software/refdata-gex-GRCh38-2024-A \
  --fastqs=/home/jelrod/skin/skin-data-oct-25/515_twelve \
  --sample=SC384 \
  --create-bam=false

cd /home/jelrod/skin/skin-counts-oct-25
cellranger count \
  --id=515_twen4 \
  --transcriptome=/home/jelrod/software/refdata-gex-GRCh38-2024-A \
  --fastqs=/home/jelrod/skin/skin-data-oct-25/515_twen4 \
  --sample=SC437GEX_LAF1804A75 \
  --create-bam=false

cd /home/jelrod/skin/skin-counts-oct-25
cellranger count \
  --id=591_baseline \
  --transcriptome=/home/jelrod/software/refdata-gex-GRCh38-2024-A \
  --fastqs=/home/jelrod/skin/skin-data-oct-25/591_baseline \
  --sample=SC367 \
  --create-bam=false

cd /home/jelrod/skin/skin-counts-oct-25
cellranger count \
  --id=591_six \
  --transcriptome=/home/jelrod/software/refdata-gex-GRCh38-2024-A \
  --fastqs=/home/jelrod/skin/skin-data-oct-25/591_six \
  --sample=SC436GEX_LAF1804A74 \
  --create-bam=false

cd /home/jelrod/skin/skin-counts-oct-25
cellranger count \
  --id=591_twelve \
  --transcriptome=/home/jelrod/software/refdata-gex-GRCh38-2024-A \
  --fastqs=/home/jelrod/skin/skin-data-oct-25/591_twelve \
  --sample=SC456Torok_LAF2960A7 \
  --create-bam=false

cd /home/jelrod/skin/skin-counts-oct-25
cellranger count \
  --id=591_twen4 \
  --transcriptome=/home/jelrod/software/refdata-gex-GRCh38-2024-A \
  --fastqs=/home/jelrod/skin/skin-data-oct-25/591_twen4 \
  --sample=SC512_051023_SK_Fresh_TOROK_NRCOS591_3Pv3_GEX_LAF5211A12 \
  --create-bam=false 

cd /home/jelrod/skin/skin-counts-oct-25
cellranger count \
  --id=594_baseline \
  --transcriptome=/home/jelrod/software/refdata-gex-GRCh38-2024-A \
  --fastqs=/home/jelrod/skin/skin-data-oct-25/594_baseline\
  --sample=SC360 \
  --create-bam=false 

cd /home/jelrod/skin/skin-counts-oct-25
cellranger count \
  --id=594_six \
  --transcriptome=/home/jelrod/software/refdata-gex-GRCh38-2024-A \
  --fastqs=/home/jelrod/skin/skin-data-oct-25/594_six\
  --sample=SC400 \
  --create-bam=false 

cd /home/jelrod/skin/skin-counts-oct-25
cellranger count \
  --id=594_twelve \
  --transcriptome=/home/jelrod/software/refdata-gex-GRCh38-2024-A \
  --fastqs=/home/jelrod/skin/skin-data-oct-25/594_twelve\
  --sample=SC429GEX_LAF1804A49 \
  --create-bam=false 

cd /home/jelrod/skin/skin-counts-oct-25
cellranger count \
  --id=594_twen4 \
  --transcriptome=/home/jelrod/software/refdata-gex-GRCh38-2024-A \
  --fastqs=/home/jelrod/skin/skin-data-oct-25/594_twen4\
  --sample=SC513_051023_SK_Frozen_NRCOS594_Nov172022_3Pv3_GEX_LAF5211A13 \
  --create-bam=false 

cd /home/jelrod/skin/skin-counts-oct-25
cellranger count \
  --id=healthy_HSK053 \
  --transcriptome=/home/jelrod/software/refdata-gex-GRCh38-2024-A \
  --fastqs=/home/jelrod/skin/skin-data-oct-25/HC/healthy_HSK053\
  --sample=HSK053 \
  --create-bam=false 

cd /home/jelrod/skin/skin-counts-oct-25
cellranger count \
  --id=healthy_SC296 \
  --transcriptome=/home/jelrod/software/refdata-gex-GRCh38-2024-A \
  --fastqs=/home/jelrod/skin/skin-data-oct-25/HC/healthy_SC296\
  --sample=SC296 \
  --create-bam=false 

cd /home/jelrod/skin/skin-counts-oct-25
cellranger count \
  --id=healthy_SC297 \
  --transcriptome=/home/jelrod/software/refdata-gex-GRCh38-2024-A \
  --fastqs=/home/jelrod/skin/skin-data-oct-25/HC/healthy_SC297\
  --sample=SC297 \
  --create-bam=false 