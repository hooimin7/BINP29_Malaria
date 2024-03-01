import pandas as pd
from Bio import SeqIO

# List of species
species = ['Ht', 'Pb', 'Pc', 'Pf', 'Pk', 'Pv', 'Py', 'Tg']

# Loop over each species
for sp in species:
    # Read the TSV file
    df = pd.read_csv(f'../06_busco/{sp}/run_apicomplexa_odb10/full_table.tsv', sep='\t', comment='#')

    # Filter the dataframe based on the status
    if sp == 'Tg':
        df = df[(df['Status'] == 'Complete') | (df['Status'] == 'Duplicated')]
    else:
        df = df[df['Status'] == 'Complete']

    # Read the FASTA file
    fasta_sequences = SeqIO.parse(open(f'../05_gene_prediction/{sp}_filNoH.faa'),'fasta')

    # Create a dictionary with sequence id as key and sequence as value
    fasta_dict = {fasta.id : str(fasta.seq) for fasta in fasta_sequences}

    # Loop over each row in the dataframe
    for index, row in df.iterrows():
        # Check if the sequence id exists in the FASTA file
        if row['Sequence'] in fasta_dict:
            # Write the sequence to a new FASTA file
            with open(f'f{index+1}.fasta', 'w') as file:
                file.write(f'>{row["Sequence"]}\n{fasta_dict[row["Sequence"]]}')
