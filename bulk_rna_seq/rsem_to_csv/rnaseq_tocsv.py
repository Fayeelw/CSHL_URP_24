import pandas as pd

rnaseq = {
    "SRR15999464": "/home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/data/RNA_seq/SRR15999464/RSEM/SRR15999464_.genes.results",
    "SRR8236750": "/home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/data/RNA_seq/SRR8236750/RSEM/SRR8236750_.genes.results",
    "SRR8236751": "/home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/data/RNA_seq/SRR8236751/RSEM/SRR8236751_.genes.results"
}

# Iterate through the dictionary items
for sample, file_path in rnaseq.items():
    rnaseq_final = []

    # Read file into a list of lists 
    with open(file_path, 'r') as rnaseq_file:
        for line in rnaseq_file:
            split_line = line.rstrip().split()
            rnaseq_final.append(split_line)

    columns = ['gene_id', 'transcript_id(s)', 'length', 'effective_length', 'expected_count', 'TPM', 'FPKM']

    rnaseq_csv = pd.DataFrame(rnaseq_final, columns=columns)

    output_file_path = f"/home/fl432/rds/rds-lynch_williams-bQvpQA58m7Y/data/RNA_seq/{sample}/RSEM/{sample}.csv"

    # Output to csv
    rnaseq_csv.to_csv(output_file_path, index=False)
