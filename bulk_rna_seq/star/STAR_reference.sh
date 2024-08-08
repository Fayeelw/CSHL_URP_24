#!/bin/bash
#
#SBATCH -p icelake
#SBATCH -A BOEMO-SL3-CPU
#SBATCH --job-name=star_ref
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --output="/home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/output/2024_01_09.star.GCF_000001405.40_GRCh38.p14_genomic.stdout"
#SBATCH --error="/home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/output/2024_01_09.star.GCF_000001405.40_GRCh38.p14_genomic.stderr"
#SBATCH --time=08:00:00

export PATH=/home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/software/STAR-2.7.11a/source/:$PATH

cd /home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/software/STAR-2.7.11a/bin/Linux_x86_64

#STAR
#-- runThreadN: Number of threads (number of available cores on server)
#--runMode: generates a reference genome for the alignment phase 
#-- genomeDIR: output path
#-- genomeFastaFiles: location of hg38 (downloaed from NIH) in fasta format (.fna)
#-- sjdbGTFfile: location of hg38 (downloaed from NIH) in gtf format (.gtf)
#--sjdbOverhang Reads -1

STAR \
--runThreadN 1 \
--runMode genomeGenerate \
--genomeDir /home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/genomes/star3 \
--genomeFastaFiles /home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/genomes/GCF_000001405.40_GRCh38.p14_genomic.fna \
--sjdbGTFfile /home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/genomes/GCF_000001405.40_GRCh38.p14_genomic.gtf \
--sjdbOverhang 74 # Reads -1