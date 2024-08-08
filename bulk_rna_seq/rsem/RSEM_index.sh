#!/bin/bash
#
#SBATCH -p icelake
#SBATCH -A BOEMO-SL3-CPU
#SBATCH --job-name=rsem_index
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16    # 32 cpus per core
#SBATCH --ntasks=1
#SBATCH --output="/home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/output/2024_01_09.rsem_index.GCF_000001405.40_GRCh38.p14_genomic.stdout"
#SBATCH --error="/home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/output/2024_01_09.rsem_index.GCF_000001405.40_GRCh38.p14_genomic.stderr"
#SBATCH --time=08:00:00

cd /home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/software/RSEM

# rsem-prepare-reference
#-gtf: path to hg8 reference genome (downloaded from NIH) in .gtf format
# path to hg8 reference genome (downloaded from NIH) in .fna format
# output path and prefix

rsem-prepare-reference \
--gtf /home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/genomes/GCF_000001405.40_GRCh38.p14_genomic.gtf \
/home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/genomes/GCF_000001405.40_GRCh38.p14_genomic.fna \
/home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/genomes/rsem3/rsem_index

