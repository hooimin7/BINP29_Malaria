import pandas as pd

# List of species
species = ['Pb', 'Pc', 'Pf', 'Pk', 'Pv', 'Py', 'Tg', 'Ht']

# Read the TSV file
df = pd.read_csv('../05_gene_prediction/myproject.proteinortho.tsv', sep='\t')

# Initialize an array to store sequence IDs
sequence_ids = []
sequence_ids.append(species)

# Loop over each row in the dataframe
for index, row in df.iterrows():
    # Check the conditions and save sequence IDs
    if (row['# Species'] == 8 and row['Genes'] == 8):
        sequence_ids.append(row.iloc[3:11].values.tolist())
    if (row['# Species'] == 8 and row['Genes'] in [9, 10, 11, 12, 13, 14, 15]):
        if any("," in str(value) for value in row.iloc[[9]].values.tolist()) and \
         not(any("," in str(value) for value in row.iloc[[3,4,5,6,7,8,10]].values.tolist())):
            sequence_ids.append(row.iloc[3:11].values.tolist())

outputNumber = 0
# Loop over each species
for sp in sequence_ids[1:]:
    # print(sp)
    outputNumber = outputNumber + 1
    outFileName = 'f{0}.faa'.format(outputNumber)
    with open(outFileName, "w") as out_f:
        # Loop over each sequence ID
        for index, seq_id in enumerate(sp):
            # If the species is 'Tg', split the sequence ID string at each comma and take the first element
            if species[index] == 'Tg':
                seq_id = seq_id.split(",")[0]
            seq_ids = '>' + seq_id
            # Read the FASTA file
            with open(f'../05_gene_prediction/{species[index]}_filNoH.faa', 'r') as fasta_seq_file:
                for line in fasta_seq_file:
                    if line.strip() == seq_ids:
                        out_f.write(line)
                        out_f.write(next(fasta_seq_file))
                    