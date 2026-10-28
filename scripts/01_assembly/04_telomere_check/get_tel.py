from Bio import SeqIO

fasta_file = "contigs.fa"
bed_file = "Gifu.tel_unit.bed"
strings = ['CCCTAAA', 'TTTAGGG']
print("Total strings:", len(strings))
x = ['A', 'T', 'G', 'C']

with open(bed_file, 'w') as f:
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq)
        for string in strings:
            for i in range(len(seq) - 6):
                for j in range(7):
                    for m in range(len(x)):
                        modified_string = string[:j] + x[m] + string[j+1:]
                        if seq[i:i + 7] == modified_string:
                            f.write(record.id + "\t" + str(i) + "\t" + str(i + 6) + "\t" + modified_string + "\n")
