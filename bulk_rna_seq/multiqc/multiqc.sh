#!/bin/bash
#
#SBATCH -p icelake
#SBATCH -A BOEMO-SL3-CPU
#SBATCH --job-name=multiqc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --output="/home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/output/2024_01_19.multiqc.stdout"
#SBATCH --error="/home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/output/2024_01_19.multiqc.stderr"
#SBATCH --time=08:00:00

# Load the multiqc_env where multiqc has been pip installed
source /home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/output/myenv/bin/activate

STAR=(
	"/home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/data/RNA_seq/SRR15999464/STAR/"
	"/home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/data/RNA_seq/SRR8236750/STAR/"
	"/home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/data/RNA_seq/SRR8236751/STAR/"
	"/home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/data/RNA_seq/SRR5833356/STAR/"	
)
RSEM=(
	"/home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/data/RNA_seq/SRR15999464/RSEM/"
	"/home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/data/RNA_seq/SRR8236750/RSEM/"
	"/home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/data/RNA_seq/SRR8236751/RSEM/"
	"/home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/data/RNA_seq/SRR5833356/STAR/"
)
FASTQC=(
	"/home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/data/RNA_seq/fastqc"

)

output_dir="/home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/data/RNA_seq/multiqc"


# Run multiqc on the fastqc, normalisation and counting steps of the RNAseq pipeline. Store in the output directory.

multiqc "${STAR[@]}" "${RSEM[@]}" "${FASTQC[@]}" -o "$output_dir"
