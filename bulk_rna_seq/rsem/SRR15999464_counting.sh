#!/bin/bash
#
#SBATCH -p icelake
#SBATCH -A BOEMO-SL2-CPU
#SBATCH --job-name=SRR15999464
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --output="/home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/output/2024_01_10.rsem_count.SRR15999464.stdout"
#SBATCH --error="/home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/output/2024_01_10.rsem_count.SRR15999464.stderr"
#SBATCH --time=08:00:00

RSEM_PATH="/home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/software/RSEM"

# go to directory where RSEM executables are
cd $RSEM_PATH

# Run RSEM
# from the help page: rsem-calculate-expression [options] --alignments [--paired-end] input reference_name (= output from rsem-prepare-reference) sample_name (= output path + prefix)
rsem-calculate-expression \
--paired-end \
--alignments /home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/data/RNA_seq/SRR15999464/STAR/alignment/SRR15999464_Aligned.toTranscriptome.out.bam \
/home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/genomes/rsem3/rsem_index \
/home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/data/RNA_seq/SRR15999464/RSEM/SRR15999464_