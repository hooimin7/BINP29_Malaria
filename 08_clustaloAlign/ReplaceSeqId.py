import re
import os

# List of species
species = ['Pb', 'Pc', 'Pf', 'Pk', 'Pv', 'Py', 'Tg', 'Ht']

# Get a list of all *_aligned.faa files in the 08 directory
files = os.listdir('/home/inf-51-2023/malaria/08_clustaloAlign/')

# Loop over each file
for file_path in files:
    # Open the file
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Initialize a counter for the species list
    species_counter = 0

    # Loop over each line in the file
    for i, line in enumerate(lines):
        # If the line starts with '>', replace it with the next species name
        if line.startswith('>'):
            new_header = re.sub(r'\d+_g', species[species_counter], line)
            lines[i] = new_header
            species_counter += 1

    # Write the modified lines back to the file
    with open(file_path, 'w') as file:
        file.writelines(lines)
