#!/bin/bash
#
#SBATCH -p icelake
#SBATCH -A BOEMO-SL3-CPU
#SBATCH --job-name=6750staralign
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output="/home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/output/2024_01_09.staralignment.SRR8236750.GCF_000001405.40_GRCh38.p14_genomic.stdout"
#SBATCH --error="/home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/output/2024_01_09.staralignment.SRR8236750.GCF_000001405.40_GRCh38.p14_genomic.stderr"
#SBATCH --time=08:00:00

module load samtools/1.14/gcc/v46lwk2d

export PATH=/home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/software/STAR-2.7.11a/source/:$PATH

cd /home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/software/STAR-2.7.11a/bin/Linux_x86_64


#STAR \
#--runThreadN: Number of threads 
#--runMode alignReads \
#--genomeDir: path to where star reference genome is stored 
#--readFilesIn: path to fastq files (downloaded from Hung et al. (2021) and Wu et. al (2022))
#--quantMode TranscriptomeSAM: outputs alignments translated into transcript coordinates (can then be used with RSEM)
#--outSAMtype BAM SortedByCoordinate: sorts the output by coordinate (tidy)
#--outFileNamePrefix: prefix used to name output file


STAR \
--runThreadN 1 \
--runMode alignReads \
--genomeDir /home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/genomes/star3/ \
--readFilesIn /home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/data/RNA_seq/SRR8236750/SRR8236750_1.fastq /home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/data/RNA_seq/SRR8236750/SRR8236750_2.fastq \
--quantMode TranscriptomeSAM \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/data/RNA_seq/SRR8236750/STAR/alignment/SRR8236750_