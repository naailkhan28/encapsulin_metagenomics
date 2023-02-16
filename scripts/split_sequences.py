from Bio import SeqIO



all_records = list(SeqIO.parse("seqs/filtered_cargo_proteins.fasta", "fasta"))
outfile = SeqIO.write(all_records[:100], "disorder/test_seqs.fasta", "fasta")    