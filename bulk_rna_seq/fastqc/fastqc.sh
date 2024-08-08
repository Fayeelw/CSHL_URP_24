#!/bin/bash
#PBS -N fastqc
#PBS -l nodes=1:ppn=8
#PBS -l walltime=01:00:00
#PBS -o /grid/beyaz/home/lynchwil/bulk_rna_seq/output/2024_07_27.fastqc.stdout
#PBS -e /grid/beyaz/home/lynchwil/bulk_rna_seq/output/2024_07_27.fastqc.stderr

export LC_ALL=C
module load EBModules-LegacyBNB
module load FastQC/0.11.8-Java-1.8

source /grid/beyaz/home/lynchwil/miniconda3/etc/profile.d/conda.sh
conda activate scrna_seq_env

# Specify the input directories for fastqc
#INPUT_DIR="/grid/beyaz/home/lynchwil/bulk_rna_seq/data/"

# Specify the input SRR ID
SRR_ID="SRR_ID"  # Replace with actual SRR ID

# Specify the output directory for fastqc
OUTPUT_DIR="/grid/beyaz/home/lynchwil/bulk_rna_seq/fastqc/html/"

# Create the output directory
mkdir -p "$OUTPUT_DIR"

# Run fastq-dump on the specified SRR ID
fastq-dump --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip --outdir "$OUTPUT_DIR" "$SRR_ID"

echo "#### FastQC completed `date`"

wait
