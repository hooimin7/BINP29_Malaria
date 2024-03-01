with open("../03_scaffold/scaffold.dat", "r") as f:
    unwanted_contig_ids = set(line.strip() for line in f)

# Read the Ht.genom file and write only the contigs not in scaffold.dat to a new file
with open("../01_remove_scaffold_h_tartakovskyi/Ht.genome", "r") as f, open("Ht.genom_filtered", "w") as out_f:
    write_contig = True
    for line in f:
        if line.startswith(">"):
            contig_id = line.split()[0][1:]  # Remove the ">" prefix
            if contig_id in unwanted_contig_ids:
                write_contig = False
            else:
                write_contig = True
                out_f.write(line)
        elif write_contig:
            out_f.write(line)
